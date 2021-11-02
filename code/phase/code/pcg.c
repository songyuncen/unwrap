/*
 *  pcg.c -- function for Preconditioned Conjugate Gradient solution
 *           of weighted least-squares phase-unwrapping problem
 */
#include <stdio.h>
#include <math.h>
#include "solncos.h"
#include "pcg.h"
#include "util.h"
#define SIMIN(x,y) (((x)*(x) < (y)*(y)) ? (x)*(x) : (y)*(y))

/* Main function of PCG algorithm for phase unwrapping.  If dywts */
/* is null, then the dx and dy weights are calculated from dxwts  */
void PCGUnwrap(float *rarray, float *zarray, float *parray,
               float *soln, float *dxwts, float *dywts, int xsize,
               int ysize, int max_iter, double epsi_con)
{
  int      i, j, k, iloop;
  double   sum, alpha, beta, beta_prev, epsi;
  double   *xcos, *ycos;
  for (k=0, sum=0.0; k<xsize*ysize; k++)
    sum += rarray[k]*rarray[k];
  sum = sqrt(sum/(xsize*ysize));
  AllocateDouble(&xcos, xsize, "cosine terms");
  AllocateDouble(&ycos, ysize, "cosine terms");
  for (iloop=0; iloop < max_iter; iloop++) {
    PCGIterate(rarray, zarray, parray, soln, dxwts, dywts, 
               xsize, ysize, xcos, ycos, iloop, sum, &alpha,
               &beta, &beta_prev, &epsi);
    if (epsi < epsi_con) {
      printf("Breaking out of main loop (due to convergence)\n");
      break;
    }
  }
  free(xcos);
  free(ycos);
} 

/* Main function within iterative loop of PCG algorithm. */
void PCGIterate(float *rarray, float *zarray, float *parray,
                float *soln, float *dxwts, float *dywts, int xsize,
                int ysize, double *xcos, double *ycos, int iloop,
                double sum0, double *alpha, double *beta,
                double *beta_prev, double *epsi)
{
  int     i, j, k;
  int     k1, k2, k3, k4;
  double  sum, w1, w2, w3, w4, btemp, delta, avg, scale;
  float   *wts;

  scale = 1.0/(xsize*ysize);
  /* remove constant bias from rarray */
  for (k=0, avg=0.0; k<xsize*ysize; k++)  avg += rarray[k];
  avg *= scale;
  for (k=0; k<xsize*ysize; k++)  rarray[k] -= avg;
  /* compute cosine transform solution of Laplacian in rarray */
  for (k=0; k<xsize*ysize; k++)  {
    zarray[k] = rarray[k];
  }
  DirectSolnByCosineTransform(zarray, xsize, ysize, xcos, ycos);
  /* calculate beta and parray */
  for (k=0, *beta=0.0; k<xsize*ysize; k++) {
    *beta += rarray[k]*zarray[k];
  }
  printf("beta = %lf\n", *beta);
  if (iloop == 0) {
    for (k=0; k<xsize*ysize; k++) {
      parray[k] = zarray[k];
    }
  }
  else {
    btemp = (*beta)/(*beta_prev);
    for (k=0; k<xsize*ysize; k++) {
      parray[k] = zarray[k] + btemp*parray[k];
    }
  }
  /* remove constant bias from parray */
  for (k=0, avg=0.0; k<xsize*ysize; k++)  avg += parray[k];
  avg *= scale;
  for (k=0; k<xsize*ysize; k++)  parray[k] -= avg;
  *beta_prev = *beta;
  /* calculate Qp */
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      k1 = (i<xsize-1) ? k + 1 : k - 1;
      k2 = (i>0) ? k - 1 : k + 1;
      k3 = (j<ysize-1) ? k + xsize : k - xsize;
      k4 = (j>0) ? k - xsize : k + xsize;
      if (dxwts==NULL && dywts==NULL) {  /* unweighted */
        w1 = w2 = w3 = w4 = 1.0;
      }
      else if (dxwts==NULL || dywts==NULL) {  /* one set of wts */
        wts = (dxwts) ? dxwts : dywts;
        w1 = SIMIN(wts[k], wts[k1]);
        w2 = SIMIN(wts[k], wts[k2]);
        w3 = SIMIN(wts[k], wts[k3]);
        w4 = SIMIN(wts[k], wts[k4]);
      }
      else {    /* dxwts and dywts are both supplied */
        w1 = dxwts[k];
        w2 = (i>0) ? dxwts[k-1] : dxwts[k];
        w3 = dywts[k];
        w4 = (j>0) ? dywts[k-xsize] : dywts[k];
      }

      zarray[k] = (w1 + w2 + w3 + w4)*parray[k]
                        - (w1*parray[k1] + w2*parray[k2] 
                                + w3*parray[k3] + w4*parray[k4]);
    }
  }
  /* calculate alpha */
  for (k=0, *alpha=0.0; k<xsize*ysize; k++) {
    *alpha += zarray[k]*parray[k];
  }
  *alpha = *beta/(*alpha);
  printf("alpha = %lf\n", *alpha);
  /* update rarray */
  for (k=0; k<xsize*ysize; k++) {
    rarray[k] -= (*alpha)*zarray[k];
  }
  /* update parray */
  for (k=0; k<xsize*ysize; k++) {
    soln[k] += (*alpha)*parray[k];
  }
  /* remove constant bias from soln */
  for (k=0, avg=0.0; k<xsize*ysize; k++)  avg += soln[k];
  avg *= scale;
  for (k=0; k<xsize*ysize; k++)  soln[k] -= avg;
  /* compute epsi and delta */
  for (k=0, sum=0.0; k<xsize*ysize; k++) {
    sum += rarray[k]*rarray[k];
  }
  *epsi = sqrt(sum/(xsize*ysize))/sum0;
  for (k=0, sum=0.0; k<xsize*ysize; k++) {
    sum += (*alpha)*(*alpha)*parray[k]*parray[k];
  }
  delta = sqrt(sum/(xsize*ysize));
  printf("ITER %d: EPSI %lf DELTA %lf\n", iloop, *epsi, delta);
}
