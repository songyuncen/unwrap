/*
 *   brcut.c -- functions for branch cutting
 */
#include <stdio.h>
#include <math.h>
#include "brcut.h"

/* Place a branch cut in the bitflags array from pixel (a,b) */
/* to pixel (c,d).  The bit for the branch cut pixels is     */
/* given by the value of "code".                             */
void PlaceCut(unsigned char *array, int a, int b, int c, int d,
              int xsize, int ysize, int code)
{
  int  i, j, k, ii, jj, m, n, istep, jstep;
  double  r;

  /* residue location is upper-left corner of 4-square */
  if (c > a && a > 0) a++;
  else if (c < a && c > 0) c++;
  if (d > b && b > 0) b++;
  else if (d < b && d > 0) d++;

  if (a==c && b==d) {
    array[b*xsize + a] |= code;
    return;
  }
  m = (a < c) ? c - a : a - c;
  n = (b < d) ? d - b : b - d;
  if (m > n) {
    istep = (a < c) ? +1 : -1;
    r = ((double)(d - b))/((double)(c - a));
    for (i=a; i!=c+istep; i+=istep) {
      j = b + (i - a)*r + 0.5;
      array[j*xsize + i] |= code;
    }
  }
  else {   /* n < m */
    jstep = (b < d) ? +1 : -1;
    r = ((double)(c - a))/((double)(d - b));
    for (j=b; j!=d+jstep; j+=jstep) {
      i = a + (j - b)*r + 0.5;
      array[j*xsize + i] |= code;
    }
  }
  return;
}

/* Return the squared distance between the pixel (a,b) and the */
/* nearest border pixel.  The border pixels are encoded in the */
/* bitflags array by the value of "border_code".               */
int DistToBorder(unsigned char *bitflags, int border_code,
               int a, int b, int *ra, int *rb, int xsize, int ysize)
{
  int  besta, bestb, found, dist2, best_dist2;
  int  i, j, k, bs;
  *ra = *rb = 0;
  for (bs=0; bs<xsize + ysize; bs++) {
    found = 0;
    best_dist2 = 1000000;  /* initialize to large value */
    /* search boxes of increasing size until border pixel found */
    for (j=b - bs; j<=b + bs; j++) {
      for (i=a - bs; i<=a + bs; i++) {
        k = j*xsize + i;
        if (i<=0 || i>=xsize - 1 || j<=0 || j>=ysize - 1
              || (bitflags[k] & border_code)) {
          found = 1;
          dist2 = (j - b)*(j - b) + (i - a)*(i - a);
          if (dist2 < best_dist2) {
            best_dist2 = dist2;
            besta = i;
            bestb = j;
          }           
        }
      }
    }
    if (found) {
      *ra = besta;
      *rb = bestb;
      break;
    }
  }
  return best_dist2;
} 
