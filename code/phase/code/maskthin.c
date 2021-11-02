/*
 *   maskthin.c -- function for thinning a mask
 */
#include <stdio.h>
#include <math.h>
#include "list.h"
#include "maskthin.h"

/* Thins mask without disconnecting it. */
void ThinMask(unsigned char *bitflags, int xsize, int ysize,
              int cut_code)
{
  int   i, j, k, n, ii, jj, kk, k1, k2, k3;
  int   numflips, removed_some=1, avoid_code;
  avoid_code = RESIDUE | BORDER;
  while (removed_some) { 
    removed_some = 0;
    for (j=1; j<ysize - 1; j++) {
      for (i=1; i<xsize - 1; i++) {
        k = j*xsize + i;
        /* only consider branch cut pixels for removal */
        if (!(bitflags[k] & cut_code)) continue;
        /* don't remove pixel k if it borders a residue or mask */
        if (bitflags[k-xsize-1] & avoid_code) continue;
        if (bitflags[k-xsize] & avoid_code) continue;
        if (bitflags[k-xsize+1] & avoid_code) continue;
        if (bitflags[k-1] & avoid_code) continue;
        if (bitflags[k] & avoid_code) continue;
        if (bitflags[k+1] & avoid_code) continue;
        if (bitflags[k+xsize-1] & avoid_code) continue;
        if (bitflags[k+xsize] & avoid_code) continue;
        if (bitflags[k+xsize+1] & avoid_code) continue;
        /* walk around border of 3x3 nbhd; make sure removal */
        /* does not leave branch cut pixels disconnected in  */
        /* 3x3 neighborhood  */
        kk = k - xsize;  /* start at upper middle pixel */
        for (n=0, numflips=0; n<4; n++) {
          if (n==0) { k1 = k-xsize; k2 = k1+1; k3 = k2+xsize; }  
          else if (n==1) { k1 = k+1; k2 = k1+xsize; k3 = k2-1; }  
          else if (n==2) { k1 = k+xsize; k2 = k1-1; k3 = k2-xsize; }  
          else if (n==3) { k1 = k-1; k2 = k1-xsize; k3 = k2+1; }  
          if (bitflags[k1] & cut_code) { 
            if (!(bitflags[k3] & cut_code))  
              ++numflips;
          }
          else {
            if (bitflags[k2] & cut_code) {
              ++numflips;
              if (!(bitflags[k3] & cut_code)) ++numflips;
            }
            else if (bitflags[k3] & cut_code) ++numflips;
          }
        }
        if (numflips < 3) { /* remove pixel from branch cut */  
          if (!(bitflags[k-xsize]&cut_code) 
                 || !(bitflags[k+1]&cut_code) 
                      || !(bitflags[k+xsize]&cut_code) 
                           || !(bitflags[k-1]&cut_code)) {
            bitflags[k] &= (~cut_code);
            ++removed_some;
          }
        }
      }
    }
    printf("Trimmed %d mask pixels\n", removed_some);
  }
}
