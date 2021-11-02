/*
 * maskfat.c -- fatten mask
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "util.h"

/* Fattens mask and border by "thickness" pixels */
void FattenMask(unsigned char *array, int mask_code, int thickness,
                int xsize, int ysize)
{
  int i, j, k, ii, jj, kk, t2=thickness;
  unsigned char *temp;
  if (t2 < 1) return;
  AllocateByte(&temp, xsize*ysize, "temp array");
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      temp[k] = 0;
      for (jj=j-t2; jj<=j+t2 && !temp[k]; jj++) {
        for (ii=i-t2; ii<=i+t2 && !temp[k]; ii++) {
          kk = jj*xsize + ii;
          if (ii>=0 && ii<xsize && jj>=0 && jj<ysize 
                                     && (array[kk] & mask_code))
            temp[k] = 1;
        }
      }
    }
  }
  for (k=0; k<xsize*ysize; k++) {
    if (temp[k]) array[k] |= mask_code;
  }
  free(temp);
}

/* Fattens zero values in quality map by "thickness" pixels */
void FattenQual(float *qual_map, int thickness,
                int xsize, int ysize)
{
  int i, j, k, ii, jj, kk, t2=thickness;
  unsigned char *temp;
  if (t2 < 1) return;
  AllocateByte(&temp, xsize*ysize, "temp array");
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      temp[k] = 0;
      for (jj=j-t2; jj<=j+t2 && !temp[k]; jj++) {
        for (ii=i-t2; ii<=i+t2 && !temp[k]; ii++) {
          kk = jj*xsize + ii;
          if (ii>=0 && ii<xsize && jj>=0 && jj<ysize 
                      && qual_map[kk]==0.0)
            temp[k] = 1;
        }
      }
    }
  }
  for (k=0; k<xsize*ysize; k++) {
    if (temp[k]) qual_map[k] = 0.0;
  }
  free(temp);
}
