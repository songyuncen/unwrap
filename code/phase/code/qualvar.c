/*
 *  qualvar.c -- functions for computing phase derivative variance
 *               quality map
 */
#include <stdio.h>
#include <math.h>
#include "qualvar.h"
#include "dxdygrad.h"

/* Compute variance of phase derivatives in array "result".    */
/* The size of the averaging template is tsize x tsize pizels. */
/* The array temp1 is required for temporary data.             */
void PhaseVariance(float *phase, float *result, 
                   unsigned char *bitflags, int ignore_code,
                   float *temp1, int tsize, int xsize, int ysize)
{
  int  add_flag;
  printf("Extracting dx gradients\n");
  DxPhaseGradient(phase, temp1, xsize, ysize);
  add_flag = 0;
  printf("Extracting dx variances\n");
  DxGradVar(temp1, result, xsize, ysize, tsize, bitflags,
            ignore_code, add_flag);
  printf("Extracting dy gradients\n");
  DyPhaseGradient(phase, temp1, xsize, ysize);
  add_flag = 1;
  printf("Extracting dy variances\n");
  /* add to dx variances */
  DyGradVar(temp1, result, xsize, ysize, tsize, bitflags,
            ignore_code, add_flag);
}

/* Compute variance of dx phase derivs in array "dxvar" in */
/* tsize x tsize windows.  If add_code is 1, the results   */
/* are added to the dxvar array, otherwise they overwrite. */
void DxGradVar(float *dx, float *dxvar, int xsize, int ysize,
  int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
  int    i, j, k, a, b, aa, bb, cc, n, hs;
  float  r, avg, avgsqr;
  if (tsize < 3 && !add_code) {
    for (i=0; i<xsize*ysize; i++)
      dxvar[i] = 0.0;
  }
  else {
    hs = tsize/2;
    for (j=0; j<ysize; j++) {
      for (i=0; i<xsize; i++) {
        avg = avgsqr = 0.0;
        for (n=0, b=j-hs; b<=j+hs; b++) {
          if ((bb = b) < 0) bb = -bb;
          else if (bb >= ysize) bb = 2*ysize - 2 - bb; 
          for (a=i-hs; a < i+hs; a++) {   /* a < to match dx */
            if ((aa = a) < 0) aa = -aa;
            else if (aa >= xsize) aa = 2*xsize - 2 - aa; 
            cc = bb*xsize + aa;
            if (bitflags && (bitflags[cc]&avoid_code)) continue;
            r = dx[cc];
            avg += r;
            avgsqr += r*r;
            ++n;  
          }
        }
        r = (n>0) ? 1.0/n : 0.0;
        avg *= r;
        avgsqr *= r;
        if (add_code)
          dxvar[j*xsize + i] += avgsqr - avg*avg;   /* variance */
        else
          dxvar[j*xsize + i] = avgsqr - avg*avg;   /* variance */
      }
    }   
  }
} 

/* Compute variance of dy phase derivs in array "dyvar" in */
/* tsize x tsize windows.  If add_code is 1, the results   */
/* are added to the dyvar array, otherwise they overwrite. */
void DyGradVar(float *dy, float *dyvar, int xsize, int ysize,
  int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
  int    i, j, k, a, b, aa, bb, cc, n, hs;
  float  r, avg, avgsqr;
  if (tsize < 3 && !add_code) {
    for (i=0; i<xsize*ysize; i++)
      dyvar[i] = 0.0;
  }
  else {
    hs = tsize/2;
    for (j=0; j<ysize; j++) {
      for (i=0; i<xsize; i++) {
        avg = avgsqr = 0.0;
        for (n=0, b=j-hs; b < j+hs; b++) { /* b < to match dy */
          if ((bb = b) < 0) bb = -bb;
          else if (bb >= ysize) bb = 2*ysize - 2 - bb; 
          for (a=i-hs; a<=i+hs; a++) {
            if ((aa = a) < 0) aa = -aa;
            else if (aa >= xsize) aa = 2*xsize - 2 - aa; 
            cc = bb*xsize + aa;
            if (bitflags && (bitflags[cc]&avoid_code)) continue;
            r = dy[cc];
            avg += r;
            avgsqr += r*r;
            ++n;
          }
        }
        r = (n>0) ? 1.0/n : 0.0;
        avg *= r;
        avgsqr *= r;
        if (add_code)
          dyvar[j*xsize + i] += avgsqr - avg*avg;   /* variance */
        else
          dyvar[j*xsize + i] = avgsqr - avg*avg;   /* variance */
      }
    }   
  }
} 
