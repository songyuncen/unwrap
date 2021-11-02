/* 
 * relax.c - function for solving weighted or unweighted
 *           version of Poisson equation (i.e., weighted
 *           or unweighted least-squares phase unwrapping)
 *           by means of red-black Gauss-Seidel relaxation
 */
#include <stdio.h>
#include <math.h>
#include "relax.h"

/*
 * Black-red Gauss-Seidel relaxation: assumes rho(i,j) 
 *        = 0.25*(dx(i,j) - dx(i+1,j) + dy(i,j) - dy(i,j+1)
 * Thus, f(i,j) 
 *    = 0.25*(f(i-1,j) + f(i+1,j) + f(i,j-1) + f(i,j+1)) + rho(i,j)
 * This version solves the weighted equation.  For the unweighted
 * equation, just pass null pointers for dxwts and dywts.
 */
void Relax(float *soln, float *dx, float *dy, float *dxwts,
           float *dywts, int w, int h, int numit)
{
  int i, j, k, n, ipass, istart;
  float  x1, x2, x3, x4, y, z, r;
  float  w1, w2, w3, w4, norm=0.0;
  double  diff, maxdiff, sum;

  if (w*h < 8*8 && numit < 2) numit = 2*(w + h);
  for (n=0; n<numit; n++) {             /* iteration */
    for (ipass=0; ipass<2; ipass++) {   /* red and black sweeps */
      istart = ipass;
      sum = maxdiff = 0.0;
      for (j=0; j<h; j++) {             /* row */
        for (i=istart; i<w; i+=2) {     /* column */
          k = j*w + i;
          /* weights */
          if (dxwts && dywts) {
            w1 = (i < w - 1) ? dxwts[k+1] : dxwts[k-1];
            if (w1 > dxwts[k]) w1 = dxwts[k];
            w1 = w1*w1;
            w2 = (i > 0) ? dxwts[k-1] : dxwts[k+1];  
            if (w2 > dxwts[k]) w2 = dxwts[k];
            w2 = w2*w2;
            w3 = (j < h - 1) ? dywts[k+w] : dywts[k-w];
            if (w3 > dywts[k]) w3 = dywts[k];
            w3 = w3*w3;
            w4 = (j > 0) ? dywts[k-w] : dywts[k+w];
            if (w4 > dywts[k]) w4 = dywts[k];
            w4 = w4*w4;
            norm = w1 + w2 + w3 + w4;
          } 
          /* surface points */
          x1 = (i < w - 1) ? soln[k+1] : soln[k-1];
          x2 = (i > 0) ? soln[k-1] : soln[k+1];
          x3 = (j < h - 1) ? soln[k+w] : soln[k-w];
          x4 = (j > 0) ? soln[k-w] : soln[k+w];
          /* derivatives */
          y = (i > 0) ? dx[k-1] : -dx[k];
          z = (j > 0) ? dy[k-w] : -dy[k];
          if (norm > 1.0e-6) { 
            r = (w1*(x1 - dx[k]) + w2*(x2 + y) 
                            + w3*(x3 - dy[k]) + w4*(x4 + z))/norm; 
          } 
          else {
            r = 0.25*(x1 - dx[k] + x2 + y + x3 - dy[k] + x4 + z);
          }
          diff = r - soln[k];
          if (diff < 0) diff = -diff;
          maxdiff = (maxdiff < diff) ? diff : maxdiff;
          sum += diff*diff;
          soln[k] = r;
        }
        istart = (istart) ? 0 : 1;
      }
    }
  }
  for (n=1; n<w; n*=2) printf("  ");
  printf("Relax %d %d: rms = %lf\n", w, h, sqrt(sum/(w*h)));
}
