#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NDIM 3
#define PI 3.141592654
#define DR 0.1
#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))

void main(int argc, char **argv) {
  
  int type, idummy, N, n, i, j, k, t, s, natoms, nconf;
  double L[NDIM], invL[NDIM], rij[NDIM], pbc[NDIM], norm, deltar, ddummy;
  double xy, xz, yz;
  double *xx, *gr, *grsum, r[NDIM];
  double **x;
  char temp[100], fname[100], outfile[100];
  FILE *ifp, *ofp;
  
  /* check that we supplied the necessary information */
  if (argc != 10) {
    printf("usage: gr-xyz <trajectory-file> <box-size> <configs-to-skip> <num-configs>\n\n");
    exit(1);
  }

  /* read in number of atoms and box size from configuration file */
  if ((ifp = fopen(argv[1], "r"))==NULL){
    printf("\nError opening file %s;  Program aborted!\n\n", fname);
    exit(1);
  }

  fclose(ifp);
  
  sscanf(argv[2], "%lf", L);
  sscanf(argv[3], "%lf", L+1);
  sscanf(argv[4], "%lf", L+2);
  sscanf(argv[5], "%lf", &xy);
  sscanf(argv[6], "%lf", &xz);
  sscanf(argv[7], "%lf", &yz);
  sscanf(argv[8], "%d", &N); 
  sscanf(argv[9], "%d", &nconf); 

printf("Number of particles = %d\n", N);

  double Lmax=-1;
  for (k=0; k<NDIM; k++) {
	  if (Lmax < L[k]) Lmax = L[k];
	  invL[k] = 1./L[k];
  }

  /* allocate needed memory for positions and g(r) */
  x = (double**) calloc(N, sizeof(double*));
  for (i=0; i<N; i++) {
    x[i] = (double*) calloc(NDIM, sizeof(double));
  }
  //gr = (double*) calloc((int)(sqrt(3.0)*L/(2.0*DR)), sizeof(double));
  /* allocate extra array space for g(r) in case box grows */
  gr = (double*) calloc((int)(5*sqrt(3.0)*Lmax/(2.0*DR)), sizeof(double));


  printf("\nCalculating g(r) for configuration...   ");

  /* open file for reading */
  if ((ifp = fopen(argv[1], "r"))==NULL){
    printf("\nError opening file %s;  Program aborted!\n\n", fname);
    exit(1);
  }

  /* THIS IS WHERE YOU WOULD ADD A LOOP OVER CONFIGURATIONS
     IN THE THE xyz FILE */
  for (s=0; s<nconf; s++){
	  printf("\b\b\b%3d", s);
	  fflush(NULL);
  /* skip first 9 lines of configuration */
  for (i=0; i<9; i++){ 
	  fgets(temp, 1000, ifp);
 // printf("%d", i);
}
  for (i=0; i<N; i++) {
	  fscanf(ifp, "%d %d %lf %lf %lf %lf %lf %lf %lf", &idummy, &idummy, &x[i][0], &x[i][1], &x[i][2], &ddummy, &ddummy, &ddummy, &ddummy);
  }
//    printf(" %lf %lf ", x[0][0], x[N-1][0]);
  fgets(temp, 100, ifp);
  
  /* loop over distinct pairs and calculate separation */
  for (i=0; i<N-1; i++) {
    for (j=i+1; j<N; j++) {
	for (k=0; k<NDIM; k++) {
        rij[k]  = x[i][k] - x[j][k];
        pbc[k]  = anint(rij[k]*invL[k]);
      }
	rij[0] -= (pbc[1]*xy + pbc[2]*xz);
      rij[0] -= L[0]*anint(rij[0]*invL[0]);

      rij[1] -= pbc[2]*yz;
      rij[1] -= L[1]*anint(rij[1]*invL[1]);

	rij[2] -= L[2]*pbc[2];
	    /* scalar separation */
      deltar = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
/* printf("%g",deltar);      
      /* add 1 to histogram at this separation */
      gr[(int)(deltar/DR)]++;
    }
  }
  }  
  /* THIS IS WHERE YOU WOULD END LOOP OVER CONFIGURATIONs */
  
  fclose(ifp);

  /* normalization so g(r)->1 at large r */
  /* will need to be adjusted when you average over many configurations */
  norm = L[0]*L[1]*L[2] / (2.0*(double)(N*(N-1)) * PI * DR) / (double)((nconf)-1);

  sprintf(outfile, "gr.dat");
  ofp = fopen(outfile, "w");
  
  for (i=1; i<(int)(sqrt(3)*Lmax/(2.0*DR)); i++) {
    fprintf(ofp, "%lf\t%lf\n", (double)i * DR, 
	    norm * gr[i]/((double) (i*i) * DR *DR));
  }
  
  
  fclose(ofp);
  
}
