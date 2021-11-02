/*
 * gridops.c - prolongation and restriction operators for 
 *             multigrid phase unwrapping.  Also contains
 *             functions for detecting coarsest array size
 *             and initializing an array to zero.
 */
#include <stdio.h>
#include <math.h>
#include "gridops.h"
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))

/* Multigrid prolongation operator: weighted or unweighted. */
/* Resamples coarse grid values and ADDS to fine grid.      */
/* (Note: for unweighted case, pass null pointers for       */
/* coarse_dxwts and coarse_dywts.)                          */  
void ProlongAndAccumulate(float *fine, int wf, int hf,
                          float *coarse, int wc, int hc,
                          float *coarse_dxwts, float *coarse_dywts)
{
  int    i, j, k, k1, k2, k3, k4, a, b, c;
  float  x1, x2, x3, x4, y1, y2, y3, y4, w1, w2, w3, w4, w;
  float  w1x, w2x, w3x, w4x, w1y, w2y, w3y, w4y;

  for (j=0, k=0; j<hc; j++) {
    for (i=0; i<wc; i++, k++) {
      a = 2*i;
      b = 2*j;
      c = b*wf + a;
      k1 = k;
      k2 = (i < wc - 1) ? k1 + 1 : k1;
      k3 = (j < hc - 1) ? k1 + wc : k1;
      k4 = (i < wc - 1) ? ((j < hc - 1) ? k1 + wc + 1 : k2) : k3;
      x1 = coarse[k1];
      x2 = coarse[k2];
      x3 = coarse[k3];
      x4 = coarse[k4];
      if (coarse_dxwts && coarse_dywts) {
        w1x = coarse_dxwts[k1];
        w2x = coarse_dxwts[k2];
        w3x = coarse_dxwts[k3];
        w4x = coarse_dxwts[k4];
        w1y = coarse_dywts[k1];
        w2y = coarse_dywts[k2];
        w3y = coarse_dywts[k3];
        w4y = coarse_dywts[k4];
        w1 = MIN(w1x, w1y);
        w2 = MIN(w2x, w2y);
        w3 = MIN(w3x, w3y);
        w4 = MIN(w4x, w4y);
        y1 = x1;
        w = w1 + w2;
        y2 = (w > 0.0) ? (w1*x1 + w2*x2)/w : 0.5*(x1 + x2);
        w = w1 + w3;
        y3 = (w > 0.0) ? (w1*x1 + w3*x3)/w : 0.5*(x1 + x3);
        w = w1 + w2 + w3 + w4;
        y4 = (w > 0.0) ? (w1*x1 + w2*x2 + w3*x3 + w4*x4)/w
                                  : 0.25*(x1 + x2 + x3 + x4);
      }
      else {
        y1 = x1;                   y2 = 0.5*(x1 + x2);
        y3 = 0.5*(x1 + x3);        y4 = 0.25*(x1 + x2 + x3 + x4);
      }
      fine[c] += y1;
      fine[c + 1] += y2;
      fine[c + wf] += y3;
      fine[c + wf + 1] += y4;
      /* what if wf is odd? fill in extra columns */
      if (i==wc - 1) {
        b = 2*j;
        for (a=2*wc; a<wf; a++) {
          c = b*wf + a;
          fine[c] += y2;
          fine[c + wf] += y4;
        }
      }
      /* what if hf is odd? fill in extra rows */
      if (j==hc - 1) {
        a = 2*i;
        for (b=2*hc; b<hf; b++) {
          c = b*wf + a;
          fine[c] += y3;
          fine[c + 1] += y4;
        }
      }
      /* fill in lower right-hand corner */
      if (i==wc - 1 && j==hc - 1) {
        for (b=2*hc; b<hf; b++) {
          for (a=2*wc; a<wf; a++) {
            c = b*wf + a;
            fine[c] += y4;
          }
        }
      }
    }
  }
}

/* Multigrid restriction operator: weighted or unweighted. */
/* Restricts residuals of derivatives to coarser grid.     */
/* (Note: for unweighted case, pass null pointers for      */
/* dxwts and dywts.)                                       */ 
void Restrict(float *dx2, float *dy2, int wc, int hc,
              float *dx, float *dy, float *dxwts, float *dywts,
              float *soln, int wf, int hf) 
{
  int     a, b, i, j, k;
  int     k1, k2, k3, k4, k5, k6, k7, k8, k9;
  int     k1x, k2x, k3x, k4x, k5x, k6x, k7x, k8x, k9x;
  int     k1y, k2y, k3y, k4y, k5y, k6y, k7y, k8y, k9y;
  float   scale;
  double  f1, f2, f3, f4, f5, f6, f7, f8, f9, wmult;
  double  w1, w2, w3, w4, w5, w6, w7, w8, w9;
  w1 = w2 = w3 = w4 = w5 = w6 = w7 = w8 = w9 = 1.0;

  /*  DX  */
  scale = (wf - 1.0)/(wc - 1.0);
  for (j=0; j<hc; j++) {
    b = 2*j;
    for (i=0; i<wc-1; i++) {   /* note: a < wf - 1 */
      a = 2*i;
      k = b*wf + a;
      /*  Indexes:    */
      /*  k1  k2  k3  */
      /*  k4  k5  k6  */
      /*  k7  k8  k9  */
      k5 = k;
      k4 = (a > 0) ? k5 - 1 : k5 + 1;
      k6 = (a < wf - 1) ? k5 + 1 : k5 - 1;
      k2 = (b > 0) ? k5 - wf : k5 + wf;
      k8 = (b < hf - 1) ? k5 + wf : k5 - wf;
      k1 = (a > 0) ? k2 - 1 : k2 + 1;
      k3 = (a < wf - 1) ? k2 + 1 : k2 - 1;
      k7 = (a > 0) ? k8 - 1 : k8 + 1;
      k9 = (a < wf - 1) ? k8 + 1 : k8 - 1;
      k1x = k2;         k4x = k5;         k7x = k8;
      k2x = k3;         k5x = k6;         k8x = k9;
      k3x = k3 + 1;     k6x = k6 + 1;     k9x = k9 + 1;
      if (dxwts && dywts) {
        /* weights: compute as MIN(k, k-1) */
        if ((w1 = dxwts[k1x]) > dxwts[k1]) w1 = dxwts[k1];
        if ((w2 = dxwts[k2x]) > dxwts[k2]) w2 = dxwts[k2];
        if ((w3 = dxwts[k3x]) > dxwts[k3]) w3 = dxwts[k3];
        if ((w4 = dxwts[k4x]) > dxwts[k4]) w4 = dxwts[k4];
        if ((w5 = dxwts[k5x]) > dxwts[k5]) w5 = dxwts[k5];
        if ((w6 = dxwts[k6x]) > dxwts[k6]) w6 = dxwts[k6];
        if ((w7 = dxwts[k7x]) > dxwts[k7]) w7 = dxwts[k7];
        if ((w8 = dxwts[k8x]) > dxwts[k8]) w8 = dxwts[k8];
        if ((w9 = dxwts[k9x]) > dxwts[k9]) w9 = dxwts[k9];
        wmult = 0.25*w5 + 0.125*(w2 + w8 + w4 + w6) 
                          + 0.0625*(w1 + w3 + w7 + w9);
      }
      /* dx residuals */
      f1 = w1*(dx[k1] - (soln[k1x] - soln[k1])); 
      f2 = w2*(dx[k2] - (soln[k2x] - soln[k2])); 
      f3 = w3*(dx[k3] - (soln[k3x] - soln[k3])); 
      f4 = w4*(dx[k4] - (soln[k4x] - soln[k4])); 
      f5 = w5*(dx[k5] - (soln[k5x] - soln[k5])); 
      f6 = w6*(dx[k6] - (soln[k6x] - soln[k6])); 
      f7 = w7*(dx[k7] - (soln[k7x] - soln[k7])); 
      f8 = w8*(dx[k8] - (soln[k8x] - soln[k8])); 
      f9 = w9*(dx[k9] - (soln[k9x] - soln[k9])); 
      if (dxwts && dywts) {
        if (wmult > 1.0e-6) dx2[j*wc + i] = scale*(0.25*f5 
                 +  0.125*(f4 + f6 + f2 + f8)
                        +  0.0625*(f1 + f3 + f7 + f9))/wmult;
        else dx2[j*wc + i] = dx[k];
      }
      else {
        dx2[j*wc + i] = scale*(0.25*f5 + 0.125*(f4 + f6 + f2 + f8)
                                   +  0.0625*(f1 + f3 + f7 + f9));
      }
    }
  }
  /* correct dx at boundary */
  for (j=0; j<hc; j++) {
    dx2[j*wc + (wc - 1)] = -dx2[j*wc + (wc - 2)];
  }

  /*  DY  */
  scale = (hf - 1.0)/(hc - 1.0);
  for (j=1; j<hc; j++) {   /* note: b > 0 */
    b = 2*j;
    for (i=0; i<wc; i++) {
      a = 2*i;
      k = b*wf + a;
      /*  Indexes:    */
      /*  k1  k2  k3  */
      /*  k4  k5  k6  */
      /*  k7  k8  k9  */
      k5 = k;
      k4 = (a > 0) ? k5 - 1 : k5 + 1;
      k6 = (a < wf - 1) ? k5 + 1 : k5 - 1;
      k2 = (b > 0) ? k5 - wf : k5 + wf;
      k8 = (b < hf - 1) ? k5 + wf : k5 - wf;
      k1 = (a > 0) ? k2 - 1 : k2 + 1;
      k3 = (a < wf - 1) ? k2 + 1 : k2 - 1;
      k7 = (a > 0) ? k8 - 1 : k8 + 1;
      k9 = (a < wf - 1) ? k8 + 1 : k8 - 1;
      k1y = k4;         k2y = k5;         k3y = k6; 
      k4y = k7;         k5y = k8;         k6y = k9;
      k7y = k7 + wf;    k8y = k8 + wf;    k9y = k9 + wf;
      if (dxwts && dywts) {
        /* weights: compute as MIN(k, k-1) */
        if ((w1 = dywts[k1y]) > dywts[k1]) w1 = dywts[k1];
        if ((w2 = dywts[k2y]) > dywts[k2]) w2 = dywts[k2];
        if ((w3 = dywts[k3y]) > dywts[k3]) w3 = dywts[k3];
        if ((w4 = dywts[k4y]) > dywts[k4]) w4 = dywts[k4];
        if ((w5 = dywts[k5y]) > dywts[k5]) w5 = dywts[k5];
        if ((w6 = dywts[k6y]) > dywts[k6]) w6 = dywts[k6];
        if ((w7 = dywts[k7y]) > dywts[k7]) w7 = dywts[k7];
        if ((w8 = dywts[k8y]) > dywts[k8]) w8 = dywts[k8];
        if ((w9 = dywts[k9y]) > dywts[k9]) w9 = dywts[k9];
        wmult = 0.25*w5 + 0.125*(w2 + w8 + w4 + w6) 
                          + 0.0625*(w1 + w3 + w7 + w9);
      }
      /* dy residuals */
      f1 = w1*(dy[k1] - (soln[k1y] - soln[k1])); 
      f2 = w2*(dy[k2] - (soln[k2y] - soln[k2])); 
      f3 = w3*(dy[k3] - (soln[k3y] - soln[k3])); 
      f4 = w4*(dy[k4] - (soln[k4y] - soln[k4])); 
      f5 = w5*(dy[k5] - (soln[k5y] - soln[k5])); 
      f6 = w6*(dy[k6] - (soln[k6y] - soln[k6])); 
      f7 = w7*(dy[k7] - (soln[k7y] - soln[k7])); 
      f8 = w8*(dy[k8] - (soln[k8y] - soln[k8])); 
      f9 = w9*(dy[k9] - (soln[k9y] - soln[k9])); 
      if (dxwts && dywts) {
        if (wmult > 1.0e-6) dy2[j*wc + i] = scale*(0.25*f5 
                 +  0.125*(f4 + f6 + f2 + f8)
                        +  0.0625*(f1 + f3 + f7 + f9))/wmult;
        else dy2[j*wc + i] = scale*dy[k];
      }
      else {
        dy2[j*wc + i] = scale*(0.25*f5 + 0.125*(f4 + f6 + f2 + f8)
                                   +  0.0625*(f1 + f3 + f7 + f9));
      }
    }
  }
  /* correct dy at boundary */
  for (i=0; i<wc; i++) {
    dy2[(hc - 1)*wc + i] = -dy2[(hc - 2)*wc + i];
  }
}

/* returns 1 if w or h is less than 3 */
int Coarsest(int w, int h, int mindim) 
{
  if (mindim==0)
    return (w/2 < 3 || h/2 < 3);
  else
    return (w/2 < mindim || h/2 < mindim);
}

/* initialize array to zeros */
void Zero(float *soln, int w, int h)
{
  int n=w*h, k;
  for (k=0; k<n; k++) {
    soln[k] = 0.0;
  }
}
