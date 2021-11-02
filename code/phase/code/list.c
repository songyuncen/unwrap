/*
 *   list.c -- functions for managing list of pixels to unwrap
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "list.h"
#include "grad.h"

/* Returns 0 if no pixels left, 1 otherwise */
int GetNextOneToUnwrap(int *a, int *b, int *index_list,
                       int *num_index, int xsize, int ysize)
{
  int index;
  if (*num_index < 1) 
    return 0;   /* return if list empty */
  index = index_list[*num_index - 1];
  *a = index%xsize;
  *b = index/xsize;
  --(*num_index);
  return 1;
}

/* Insert new pixel into the list.  */
/* Note: qual_map can be NULL       */
void InsertList(float *soln, float val, float *qual_map, 
           unsigned char *bitflags, int a, int b, int *index_list,
           int *num_index, int xsize, int ysize, int processed_code,
           int postponed_code, float *min_qual, int max_list_size)
{
  int i, n, index, top, bot, mid;
  double  quality;

  index = b*xsize + a;
  quality = (qual_map) ? qual_map[index] : 0.0;
  if (~(bitflags[index] & postponed_code)) { 
    /* if not postponed, store new unwrapped value */
    if (soln) soln[index] = val;
  }
  else {
    /* if postponed, leave old value */
  }
  /* if quality is too low, postpone it */
  if (qual_map && min_qual && quality < *min_qual) {
    bitflags[index] |= postponed_code;
    return;
  }
  /* otherwise, add to list */
  if (!qual_map) {   /* don't order if qual_map is NULL */
    index_list[*num_index] = index;
    ++(*num_index);
  }
  else {
    /* insert in list in order from lowest to highest quality */
    n = *num_index;
    if (n < 1) {
      (*num_index) = 0;   /* will be incremented below */
      index_list[0] = index;
    }
    else {
      if (quality <= qual_map[index_list[0]]) {
        /* insert at top of list */
        for (i=n; i>0; i--) 
          index_list[i] = index_list[i-1];
        index_list[0] = index;
      }
      else if (quality > qual_map[index_list[n - 1]]) {
        /* insert at bottom */
        index_list[n] = index;
      }
      else {   /* insert in middle */
        top = 0;
        bot = n - 1;
        while (bot - top > 1) {
          mid = (top + bot)/2;
          if (quality <= qual_map[index_list[mid]])  bot = mid;
          else  top = mid;
        }
        for (i=n; i>top+1; i--) 
          index_list[i] = index_list[i-1];
        index_list[top+1] = index;
      }
    }
    ++(*num_index);
  }
  bitflags[index] |= processed_code;
  bitflags[index] &= (~postponed_code);

  /* trim list if it's too big, and increase the quality */
  if (qual_map && min_qual 
          && max_list_size > 0 && *num_index >= max_list_size) {
    n = 0.50*(*num_index);  /* discard 50% */
    for (i=0; i<n; i++) {
      bitflags[index_list[i]] |= postponed_code;
      bitflags[index_list[i]] &= (~processed_code);
    }
    for (i=0; i<*num_index - n; i++)
      index_list[i] = index_list[i + n];
    *num_index -= n;
    *min_qual = qual_map[index_list[0]];
  }

  return;
}

/* Insert the four neighboring pixels of the given pixel */
/* (x,y) into the list.  The quality value of the given  */
/* pixel is "val".                                       */
void UpdateList(float *qual_map, int x, int y, float val,
         float *phase, float *soln, unsigned char *bitflags,
         int xsize, int ysize, int *index_list, int *num_index,
         int ignore_code, int processed_code, int postponed_code, 
         int max_list_size, int dxdy_flag, float *min_qual)
{
  int    i, a, b, k, w;
  float  grad;
  a = x - 1;
  b = y;
  k = b*xsize + a;

  if (a >= 0 
        && !(bitflags[k] & (ignore_code | processed_code))) {
    w = y*xsize + x-1;
    grad = Gradient(phase[w], phase[w+1]);
    if (dxdy_flag && qual_map)
      qual_map[k] = -fabs(grad);
    InsertList(soln, val + grad, qual_map, bitflags, a, b,
               index_list, num_index, xsize, ysize, processed_code,
               postponed_code, min_qual, max_list_size);
  }

  a = x + 1;
  b = y;
  k = b*xsize + a;
  if (a < xsize 
        && !(bitflags[k] & (ignore_code | processed_code))) {
    w = y*xsize + x;
    grad = - Gradient(phase[w], phase[w+1]);
    if (dxdy_flag && qual_map)
      qual_map[k] = -fabs(grad);
    InsertList(soln, val + grad, qual_map, bitflags, a, b,
               index_list, num_index, xsize, ysize, processed_code,
               postponed_code, min_qual, max_list_size);
  }

  a = x;
  b = y - 1;
  k = b*xsize + a;
  if (b >= 0 
        && !(bitflags[k] & (ignore_code | processed_code))) {
    w = (y-1)*xsize + x;
    grad = Gradient(phase[w], phase[w+xsize]);
    if (dxdy_flag && qual_map) 
      qual_map[k] = -fabs(grad);
    InsertList(soln, val + grad, qual_map, bitflags, a, b,
               index_list, num_index, xsize, ysize, processed_code,
               postponed_code, min_qual, max_list_size);
  }

  a = x;
  b = y + 1;
  k = b*xsize + a;
  if (b < ysize 
        && !(bitflags[k] & (ignore_code | processed_code))) {
    w = y*xsize + x;
    grad = - Gradient(phase[w], phase[w+xsize]);
    if (dxdy_flag && qual_map) qual_map[k] = -fabs(grad);
    InsertList(soln, val + grad, qual_map, bitflags, a, b,
               index_list, num_index, xsize, ysize, processed_code, 
               postponed_code, min_qual, max_list_size);
  }
}
