/*
 *  solncos.c -- function for solving Poisson's equation by means of
 *               fast cosine transform
 */
#include <stdio.h>
#include <math.h>
#include "dct.h"
#include "pi.h"

/* Main function for direct solution of Poisson's equation for     */
/* unweighted least squares phase unwrapping by means of discrete  */
/* cosine transform (which actually uses FFTs to compute the DCTs) */
void DirectSolnByCosineTransform(float *array, int xsize, int ysize,
                                 double *xcos, double *ycos)
{
  int   i, j, m, n;
  /*  Transform  */
  printf("Forward transform...\n");
  FastCosineTransform(array, xsize, ysize, 1);

  /*  Divide  */
  printf("Scaling...\n");
  for (i=0; i<xsize; i++) {
    xcos[i] = cos(i*PI/(xsize - 1.0));
  }
  for (j=0; j<ysize; j++) {
    ycos[j] = cos(j*PI/(ysize - 1.0));
  }
  array[0] = 0.0;
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      if (i==0 && j==0) {
        array[0] = 0.0;
      }
      else {
        array[j*xsize + i] /= (4.0 - 2.0*xcos[i] - 2.0*ycos[j]);
      }
    }
  }
  /*  Transform  */
  printf("Inverse transform...\n");
  FastCosineTransform(array, xsize, ysize, -1);
}
