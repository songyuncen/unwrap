/*
 *   qualpseu.c -- functions for computing pseudocorrelation 
 *                 quality map
 */
#include <stdio.h>
#include <math.h>
#include "pi.h"
#include "qualpseu.h"
#include "dxdygrad.h"

/* Compute pseudocorrelation quality map in array "result".    */
/* The size of the averaging template is tsize x tsize pizels. */
/* The array temp1 is required for temporary data.             */
void PseudoCorrelation(float *phase, float *result,
                     unsigned char *bitflags, int ignore_code,
                     float *temp1, int tsize, int xsize, int ysize)
{
  int  i, j, k, add_flag;
  printf("Extracting cos's\n");
  for (k=0; k<xsize*ysize; k++) 
    temp1[k] = cos(TWOPI*phase[k]);
  add_flag = 0;
  printf("Filtering cos's\n");
  SqrAvgFilter(temp1, result, xsize, ysize, tsize,
               bitflags, ignore_code, add_flag);
  printf("Extracting sin's\n");
  for (k=0; k<xsize*ysize; k++) 
    temp1[k] = sin(TWOPI*phase[k]);
  add_flag = 1;
  printf("Filtering sin's\n");
  SqrAvgFilter(temp1, result, xsize, ysize, tsize,
               bitflags, ignore_code, add_flag);
  printf("Square root\n");
  for (k=0; k<xsize*ysize; k++) 
    result[k] = 1.0 - sqrt(result[k]);
}

/* Computes average of array 'in' in tsize x tsize windows,  */
/* then squares the values.  If add_code is 1, the results   */ 
/* are added to the values in out, otherwise they overwrite. */
void SqrAvgFilter(float *in, float *out, int xsize, int ysize,
   int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
  int    i, j, k, a, b, aa, bb, cc, n, hs;
  float  r, avg;
  if (tsize < 3 && !add_code) {
    for (i=0; i<xsize*ysize; i++)
      out[i] = 0.0;
  }
  else {
    hs = tsize/2;
    for (j=0; j<ysize; j++) {
      for (i=0; i<xsize; i++) {
        avg = 0.0;
        for (n=0, b=j-hs; b<=j+hs; b++) {
          if ((bb = b) < 0) bb = -bb;
          else if (bb >= ysize) bb = 2*ysize - 2 - bb; 
          for (a=i-hs; a <= i+hs; a++) {  
            if ((aa = a) < 0) aa = -aa;
            else if (aa >= xsize) aa = 2*xsize - 2 - aa; 
            cc = bb*xsize + aa;
            if (aa>=0 && aa<xsize-1 && bb>=0 && bb<ysize-1) {
              r = in[cc];
              avg += r;
              ++n;  
            }
          }
        }
        r = (n>0) ? 1.0/n : 0.0;
        avg *= r;
        if (add_code)
          out[j*xsize + i] += avg*avg;
        else
          out[j*xsize + i] = avg*avg; 
      }
    }   
  }
} 
