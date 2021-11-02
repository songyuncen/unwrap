/*
 *  qualgrad.c -- functions for computing max gradient quality map
 */
#include <stdio.h>
#include <math.h>
#include "qualgrad.h"
#include "dxdygrad.h"

/* Compute max gradient quality map in array "result".  Size */
/* of averaging template is tsize x tsize pizels.  The array */
/* temp1 is required for temporary data.                     */
void MaxPhaseGradients(float *phase, float *result,
                    unsigned char *bitflags, int ignore_code,
                    float *temp1, int tsize, int xsize, int ysize)
{
  int  add_flag;
  printf("Extracting dx gradients\n");
  DxPhaseGradient(phase, temp1, xsize, ysize);
  add_flag = 0;
  printf("Extracting dx max's\n");
  DxGradMax(temp1, result, xsize, ysize, tsize, bitflags,
            ignore_code, add_flag);
  printf("Extracting dy gradients\n");
  DyPhaseGradient(phase, temp1, xsize, ysize);
  add_flag = 1;
  printf("Extracting dy max's\n");
  /* add to dx max's */
  DxGradMax(temp1, result, xsize, ysize, tsize, bitflags,
            ignore_code, add_flag);
}

/* Compute max dx gradients in array "dxmax". If add_code */
/* is 1, the results are added to the dxmax array,        */
/* otherwise they overwrite.                              */
void DxGradMax(float *dx, float *dxmax, int xsize, int ysize,
   int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
  int    i, j, a, b, c, hs;
  float  r, dmax;
  if (tsize < 3) {
    for (i=0; i<xsize*ysize; i++)
      dxmax[i] = (add_code) ? dxmax[i] + dx[i] : dx[i];
  }
  else {
    hs = tsize/2;
    for (j=0; j<ysize; j++) {
      for (i=0; i<xsize; i++) {
        for (dmax = 0.0, b=j-hs; b<=j+hs; b++) {
          if (b < 0 || b >= ysize) continue;
          for (a=i-hs; a < i+hs; a++) {  /* a < to match dx */
            if (a < 0 || a >= xsize) continue;
            c = b*xsize + a;
            if (bitflags && (bitflags[c]&avoid_code)) continue;
            r = dx[c];
            if (r < 0) r = -r;
            if (dmax < r) dmax = r; 
          }
        }
        if (add_code)
          dxmax[j*xsize + i] += dmax;
        else
          dxmax[j*xsize + i] = dmax;
      }
    }   
  }
} 

/* Compute max dy gradients in array "dymax". If add_code */
/* is 1, the results are added to the dymax array,        */
/* otherwise they overwrite.                              */
void DyGradMax(float *dy, float *dymax, int xsize, int ysize,
   int tsize, unsigned char *bitflags, int avoid_code, int add_code)
{
  int    i, j, a, b, c, hs;
  float  r, dmax;
  if (tsize < 3) {
    for (i=0; i<xsize*ysize; i++)
      dymax[i] = (add_code) ? dymax[i] + dy[i] : dy[i];
  }
  else {
    hs = tsize/2;
    for (j=0; j<ysize; j++) {
      for (i=0; i<xsize; i++) {
        for (dmax = 0.0, b=j-hs; b < j+hs; b++) {  
          /* b < j+hs to match dy */
          if (b < 0 || b >= ysize) continue;
          for (a=i-hs; a<=i+hs; a++) {
            if (a < 0 || a >= xsize) continue;
            c = b*xsize + a;
            if (bitflags && (bitflags[c] & avoid_code)) continue;
            r = dy[c];
            if (r < 0) r = -r;
            if (dmax < r) dmax = r; 
          }
        }
        if (add_code)
          dymax[j*xsize + i] += dmax;
        else
          dymax[j*xsize + i] = dmax;
      }
    }   
  }
} 
