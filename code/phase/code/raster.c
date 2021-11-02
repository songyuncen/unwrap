/*
 *  raster.c -- simple raster function for unwrapping residue-free
 *              phase data by means of Itoh's method.
 */
#include <stdio.h>
#include <math.h>
#include "raster.h"
#include "grad.h"
/* Unwrap phase using simple raster technique based on Itoh's */
/* method.  This method works inly if there are no residues.  */
void RasterUnwrap(float *phase, float *soln, int xsize, int ysize)
{
  int    i, j, k;
  soln[0] = phase[0];
  for (j=1; j<ysize; j++) {  /* unwrap first column */
    k = j*xsize + 0;
    soln[k] = soln[k-xsize] + Gradient(phase[k], phase[k-xsize]);
  }
  for (j=0; j<ysize; j++) {   /* unwrap all the rows */
    for (i=1; i<xsize; i++) {
      k = j*xsize + i;
      soln[k] = soln[k-1] + Gradient(phase[k], phase[k-1]);
    }
  }
}
