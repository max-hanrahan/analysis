#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NDIM 3
#define PI 3.141592654
#define anint(x) ((x >= .5) ? (1.0) : (x <= -.5) ? (-1.0) : (0.0))
#define QMAX 100
#define NTYPE 3 // 2 types + 1, where type 0 represents both.

typedef struct {
  double r, i;
} complex;

int N, N_per_chain, Npolymer;
// Example of running this code:
// ~/dyn_x/analysis/sq_lin <file of all xyz coords> <boxsize> <lo> <hi> <num particles>
void read_qvectors(int qvector[2 * QMAX][2 * QMAX][NDIM], int nvec[2 * QMAX]);

int main(int argc, char **argv) {
  /* note: all QMAX arrays are over-allocated by 2, becuase some
     wierd corruption was happening that i could not find */
  double **x, rho, L, norm, Lx, sf[NTYPE][2 * QMAX], sfchain[2 * QMAX],
      scross[2 * QMAX];
  complex rhoq[NTYPE], rhoqchain;
  float *xx;
  int i, ii, j, k, t, start, end, skip, vec, nlin, lo, hi, pnum, npow, *type,
      ptype, dummy;
  int q;
  double rx, ry, rz, qx, qy, qz, qmin, qdotr;
  int qvector[2 * QMAX][2 * QMAX][NDIM], nvec[2 * QMAX];
  char fname[100], outfile[100], temp[100], coordfile[100];
  FILE *fp;

  sscanf(argv[1], "%s", coordfile);
  sscanf(argv[2], "%lf", &Lx);
  sscanf(argv[3], "%d", &lo);
  sscanf(argv[4], "%d", &hi);
  sscanf(argv[5], "%d", &pnum);

  printf("Number of Particles: %d\n", pnum);
  printf("Boxsize %lf\n", Lx);
  N_per_chain = 20;
  Npolymer = pnum / N_per_chain;

  rhoq[0].r = 0.0;
  rhoq[0].i = 0.0;
  rhoq[1].r = 0.0;
  rhoq[1].i = 0.0;
  rhoq[2].r = 0.0;
  rhoq[2].i = 0.0;

  /* read qvectors */
  printf("Reading q-vectors to calculate...\n");
  fflush(NULL);
  read_qvectors(qvector, nvec);

  qmin = 2.0 * PI / Lx;

  nlin = hi - lo + 1;

  /* allocate needed memory */
  x = (double **)calloc(pnum, sizeof(double *));
  for (i = 0; i < pnum; i++)
    x[i] = (double *)calloc(NDIM, sizeof(double));
  type = (int *)calloc(pnum, sizeof(int));

  printf("\nCalculating...    ");

  /* read in dump file */

  if ((fp = fopen(argv[1], "r")) == NULL) {
    printf("\nError opening file %s;  Program aborted!\n\n", fname);
    exit(1);
  }

  /* loop over all files */
  for (t = 0; t <= hi; t++) {
    for (i = 0; i < 9; i++) {
      fgets(temp, 100, fp);
    }
    fgets(temp, 100, fp);
    fgets(temp, 100, fp);
    // ii = 0;
    for (i = 0; i < pnum; i++) {
      fscanf(fp, "%d %d %d %lf %lf %lf", &dummy, &dummy, &type[i], &rx, &ry,
             &rz);
      x[i][0] = rx;
      x[i][1] = ry;
      x[i][2] = rz;
    }
    if (t >= lo) {
      for (q = 1; q < QMAX; q++) {
        for (vec = 0; vec < nvec[q]; vec++) {
          for (i = 0; i < NTYPE; i++)
            rhoq[i].r = rhoq[i].i = 0;

          qx = qmin * qvector[q][vec][0];
          qy = qmin * qvector[q][vec][1];
          qz = qmin * qvector[q][vec][2];
          for (i = 0; i < Npolymer; i++) {
            rhoqchain.r = rhoqchain.i = 0.0;
            ii = i * N_per_chain;
            for (j = 0; j < N_per_chain; j++) {
              qdotr = x[ii + j][0] * qx + x[ii + j][1] * qy + x[ii + j][2] * qz;
              rhoqchain.r += cos(qdotr);
              rhoqchain.i += sin(qdotr);
            }
            sfchain[q] += rhoqchain.r * rhoqchain.r + rhoqchain.i * rhoqchain.i;

            rhoq[0].r += rhoqchain.r;
            rhoq[0].i += rhoqchain.i;

            rhoq[type[ii]].r += rhoqchain.r;
            rhoq[type[ii]].i += rhoqchain.i;
          }
          for (i = 0; i < NTYPE; i++)
            sf[i][q] += rhoq[i].r * rhoq[i].r + rhoq[i].i * rhoq[i].i;
          scross[q] += 2 * (rhoq[1].r * rhoq[2].r + rhoq[1].i * rhoq[2].i);
        }
      }
    }
    fgets(temp, 100, fp);
  }

  norm = 1.0 / ((double)nlin);
  norm /= (double)pnum;

  sprintf(outfile, "sq.dat");
  fp = fopen(outfile, "w");
  // todo: modify these headers to be a bit more descriptive.
  // also calculate the number of particles each belongs to
  fprintf(fp, "# Q\tSq_total\tSq_type1\tSq_type_2\tSq_cross\tSq_chain\n");

  for (q = 1; q < QMAX; q++) {
    fprintf(
        fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", (double)(q + 1) * 0.5 * qmin,
        sf[0][q] * norm / (double)nvec[q], sf[1][q] * norm / (double)nvec[q],
        sf[2][q] * norm / (double)nvec[q], scross[q] * norm / (double)nvec[q],
        sfchain[q] * norm / (double)nvec[q]);
  }
  fclose(fp);
  printf("\nStructure factor in file %s\n", outfile);
}

void read_qvectors(int qvector[2 * QMAX][2 * QMAX][NDIM], int nvec[2 * QMAX]) {
  FILE *fp;
  char fname[100];
  int q;

  for (q = 1; q < QMAX; q++) {

    sprintf(fname, "/zfshomes/fstarr/qvector/qvector.%03d", q + 1);
    if ((fp = fopen(fname, "r")) == NULL) {
      printf("\nError opening file %s;  Program aborted!\n\n", fname);
      exit(1);
    }

    nvec[q] = 0;
    fscanf(fp, "%d %d %d", qvector[q][nvec[q]], qvector[q][nvec[q]] + 1,
           qvector[q][nvec[q]] + 2);
    while (!feof(fp)) {
      nvec[q]++;
      fscanf(fp, "%d %d %d", qvector[q][nvec[q]], qvector[q][nvec[q]] + 1,
             qvector[q][nvec[q]] + 2);
    }

    fclose(fp);
  }
}

