/*
 *   residues.c -- function for extracting residues from phase
 */
#include <stdio.h>
#include <math.h>
#include "residues.h"
#include "grad.h"
/* Detect residues in phase data and mark them as positive or  */
/* negative residues in the bitflags array.  Ignore the pixels */
/* marked with avoid_code in the bitflags araay.               */
int Residues(float *phase, unsigned char *bitflags, int posres_code,
             int negres_code, int avoid_code, int xsize, int ysize)
{
  int  i, j, k, NumRes=0;
  double  r;
  for (j=0; j<ysize - 1; j++) {
    for (i=0; i<xsize - 1; i++) {
      k = j*xsize + i;
      if (bitflags && ((bitflags[k] & avoid_code) 
           || (bitflags[k+1] & avoid_code)
             || (bitflags[k+1+xsize] & avoid_code)
                 || (bitflags[k+xsize] & avoid_code))) {
        continue; /* masked region: don't unwrap */
      }
      r = Gradient(phase[k+1], phase[k])
           + Gradient(phase[k+1+xsize], phase[k+1])
              + Gradient(phase[k+xsize], phase[k+1+xsize])
                 + Gradient(phase[k], phase[k+xsize]);
      if (bitflags) {
        if (r > 0.01) bitflags[k] |= posres_code;
        else if (r < -0.01) bitflags[k] |= negres_code;
      }
      if (r*r > 0.01) ++NumRes;
    }
  }
  return NumRes;
}
