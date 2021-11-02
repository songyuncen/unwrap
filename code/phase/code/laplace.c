/*
 *  laplace.c -- compute weighted (or unweighted) Laplacian
 */
#include <stdio.h>
#include <math.h>
#include "grad.h"
#include "laplace.h"
#define SIMIN(x,y) (((x)*(x) < (y)*(y)) ? (x)*(x) : (y)*(y))

/* Compute the dx and dy weighted (or unweighted) Laplacian.      */
/* laptype=1 for wrapped phase Laplacian, =0 for usual Laplacian. */
void ComputeLaplacian(float *input, float *laplacian, float *dxwts,
                    float *dywts, int xsize, int ysize, int laptype)
{
  double w1, w2, w3, w4;
  int    i, j, k1, k2, k3, k4, k;
  float  *wts;
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
        w2 = (i>0) ? dxwts[k-1] : dxwts[k];     /* boundary condition */
        w3 = dywts[k];
        w4 = (j>0) ? dywts[k-xsize] : dywts[k];    /* boundary condition */
      }
      if (laptype) {   /* compute wrapped phase Laplacian */
        float *phase = input;
        laplacian[k]  = w1*Gradient(phase[k], phase[k1])
                          + w2*Gradient(phase[k], phase[k2])
                             + w3*Gradient(phase[k], phase[k3])
                                + w4*Gradient(phase[k], phase[k4]);
      }
      else {   /* compute usual (not wrapped) Laplacian */
        float *surface = input;
        laplacian[k]  = w1*(surface[k] - surface[k1])
                          + w2*(surface[k] - surface[k2])
                             + w3*(surface[k] - surface[k3])
                                + w4*(surface[k] - surface[k4]);
      }
    }
  }
}

/* Compute the dx and dy weighted laplacian.  The parameter */
/* e0 is deiscussed in the text.                            */
void ComputeDerivWts(float *phase, float *soln, float *dxwts,
                     float *dywts, float *qual_map,
                     double e0, int xsize, int ysize)
{
  double w, r, eps=1.0e-06;
  int    i, j, k, k1, k2, k3, k4;
  w = 1.0;  /* w must be 1 if qual_map is NULL */
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      k1 = (i < xsize - 1) ? k + 1 : k - 1;
      r = soln[k1] - soln[k] - Gradient(phase[k1], phase[k]);
      dxwts[k] = e0/(r*r + e0);
      if (qual_map) {
        w = qual_map[k];
        if (w > qual_map[k1]) w = qual_map[k1];
      }
      dxwts[k] *= w;  /* must have 0 <= w <= 1 */
    }
  }  
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      k3 = (j < ysize - 1) ? k + xsize : k - xsize;
      r = soln[k3] - soln[k] - Gradient(phase[k3], phase[k]);
      dywts[k] = e0/(r*r + e0);
      if (qual_map) {
        w = qual_map[k];
        if (w > qual_map[k3]) w = qual_map[k3];
      }
      dywts[k] *= w;  /* must have 0 <= w <= 1 */
    }
  }  
}
