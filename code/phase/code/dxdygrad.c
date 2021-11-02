/*
 * dxdygrad.c -- functions for computing phase derivatives in arrays
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "grad.h"
#include "dxdygrad.h"

/* Compute the wrapped row differences (d/dx partial derivatives */
/* of the phase) and store them in the array "dx".               */
void DxPhaseGradient(float *phase, float *dx, int xsize, int ysize)
{
  int  i, j, k;
  float  grad;
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      if (i < xsize - 1) 
        grad = Gradient(phase[k+1], phase[k]);
      else 
        grad = Gradient(phase[k-1], phase[k]);
      dx[k] = grad;
    }
  }
}

/* Compute the wrapped column differences (d/dy partial  */
/* derivatives) and store them in the array "dy".        */
void DyPhaseGradient(float *phase, float *dy, int xsize, int ysize)
{
  int  i, j, k;
  float  grad;
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      if (j < ysize - 1) 
        grad = Gradient(phase[k+xsize], phase[k]);
      else 
        grad = Gradient(phase[k-xsize], phase[k]);
      dy[k] = grad;
    }
  }
}
