#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NDIM 3
#define PI 3.141592654
#define DR 0.1
#define NTYPES 3 // 2 types + 1, where type 0 represents both.
#define anint(x) ((x >= .5) ? (1.0) : (x <= -.5) ? (-1.0) : (0.0))

void main(int argc, char **argv) {

  int *type, dummy, N, n, i, j, k, t, natoms, nconf, nskip;
  double L, norm, deltar;
  double *xx, **gr, *grsum, *grcross, r[NDIM];
  double **x;
  char temp[100], fname[100], outfile[100];
  FILE *ifp, *ofp;

  /* check that we supplied the necessary information */
  if (argc != 5) {
    printf("usage: gr-xyz <trajectory-file> <box-size> <# configs> <# skiping "
           "configs>\n\n");
    exit(1);
  }

  /* read in number of atoms and box size from configuration file */
  if ((ifp = fopen(argv[1], "r")) == NULL) {
    printf("\nError opening file %s;  Program aborted!\n\n", fname);
    exit(1);
  }

  N = 7600;

  sscanf(argv[2], "%lf", &L);
  sscanf(argv[3], "%d", &nconf);
  sscanf(argv[4], "%d", &nskip);

  printf("%d atoms; Box size L = %lf; # configs = %d; # skiping configs = %d\n",
         N, L, nconf, nskip);

  /* allocate needed memory for positions and g(r) */
  x = (double **)calloc(N, sizeof(double *));
  for (i = 0; i < N; i++) {
    x[i] = (double *)calloc(NDIM, sizeof(double));
  }
    type = (int *)calloc(N, sizeof(int));


  /* allocate extra array space for g(r) in case box grows */
  gr = (double **)calloc(NTYPES, sizeof(double *));
  for (i = 0; i < NTYPES; i++) {
    gr[i] = (double *)calloc((int)(5 * sqrt(3.0) * L / (2.0 * DR)), sizeof(double));
  }
  grcross = (double *)calloc((int)(5 * sqrt(3.0) * L / (2.0 * DR)), sizeof(double));


  printf("\nCalculating g(r) ...");

  /* open file for reading */
  if ((ifp = fopen(argv[1], "r")) == NULL) {
    printf("\nError opening file %s;  Program aborted!\n\n", fname);
    exit(1);
  }

  /* Loop over configurations */
  for (n = 0; n < nconf; n++) {
    printf("\b\b\b%2.0f%%", 100.0 * (double)n / (double)nconf);
    fflush(NULL);

    /* skip first 9 lines of configuration */
    for (i = 0; i < 9; i++)
      fgets(temp, 100, ifp);

    for (i = 0; i < N; i++) {
      fscanf(ifp, "%d%d%d%lf%lf%lf", &dummy, &dummy, &type[i], &x[i][0], &x[i][1], &x[i][2]);
    }
    fgets(temp, 100, ifp);
    if (n >= nskip) {

      /* loop over distinct pairs and calculate separation */
      for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
          for (k = 0; k < NDIM; k++) {
            /* vector separation */
            r[k] = x[i][k] - x[j][k];
            r[k] -= L * anint(r[k] / L); // periodic BC correction
          }
          /* scalar separation */
          deltar = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

          /* add 1 to histogram at this separation for the correct type */
          gr[0][(int)(deltar / DR)]++;

          if (type[i] != type[j]) {
            grcross[(int)(deltar / DR)]++;
          }
          else {
          gr[type[i]][(int)(deltar / DR)]++;
          }
        }
      }
    } /* close if (nskip) */
  }   /* close for (n=0...) */

  /* End loop over configs */

  fclose(ifp);

/* normalization so g(r)->1 at large r */
norm = L * L * L / (2.0 * (double)(N * (N - 1)) * PI * DR) / (double)(nconf - nskip);

ofp = fopen("gr.dat", "w");
fprintf(ofp, "#r\tall\t1-1\t2-2\t1-2\n");
for (i = 1; i < (int)(sqrt(3) * L / (2.0 * DR)); i++) {
  fprintf(ofp, "%lf", (double)i * DR);
  
  for (t = 0; t < NTYPES; t++) {
    fprintf(ofp, "\t%lf", norm * gr[t][i] / ((double)(i * i) * DR * DR));
  }
  
  fprintf(ofp, "\t%lf\n", norm * grcross[i] / ((double)(i * i) * DR * DR));
}
fclose(ofp);
}
