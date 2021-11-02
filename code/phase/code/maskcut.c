/*
 *   maskcut.c -- function for quality map guided mask cut generation
 */
#include <stdio.h>
#include <math.h>
#include "list.h"
#include "util.h"
#include "maskcut.h"

/* Main function for Mask Cut Algorithm for phase unwrapping */
int QualityGuidedMask(float *phase, float *qual_map, 
                      unsigned char *bitflags, int xsize, int ysize, 
                      int mask_code, UnwrapMode unwrap_mode,
                      int debug_flag, char *infile)
{
  int  i, j, k, a, b, c, n, num_pieces=0;
  float  value;
  float  min_qual, small_val = -1.0E+10;
  int    num_index, max_list_size, bench, benchout;
  int    *index_list;
  int    postponed_code=POSTPONED;
  int    avoid_code;
  int    charge, num_pixels;

  bench = xsize*ysize/100;
  if (bench < 1) bench = 1;
  benchout = 10*bench;
  max_list_size = xsize*ysize/2; /* this size can be reduced */ 
  AllocateInt(&index_list, max_list_size + 1,
              "bookkeeping list (index)");

  /* find starting point */
  n = 0;
  num_index = 0;
  avoid_code = BORDER | mask_code;
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      min_qual = small_val;
      num_index = 0;
      k = j*xsize + i;
      if ((bitflags[k] & RESIDUE) && !(bitflags[k] & avoid_code)) {
        charge = (bitflags[k] & POS_RES) ? 1 : -1;
        ++num_pieces;
        num_pixels = 1;
        bitflags[k] |= mask_code;
        if ((bitflags[k] & BORDER) 
            || i==0 || i==xsize-1 || j==0 || j==ysize-1) {
          continue;
        }
        UpdateList(qual_map, i, j, value, phase, NULL, bitflags,
                xsize, ysize, index_list, &num_index,
                postponed_code, VISITED, postponed_code,
                max_list_size, unwrap_mode==dxdygrad, &min_qual);
        while (charge != 0 && num_index > 0) {
          if (n%bench==0) {
            printf("%d ", n/bench);
            fflush(stdout);
          }
          ++n;
          if (!GetNextOneToUnwrap(&a, &b, index_list, &num_index,
                                  xsize, ysize))
            break;   /* no more to unwrap */
          c = b*xsize + a;
          ++num_pixels; 
          if ((bitflags[c] & RESIDUE) && !(bitflags[c] & mask_code)) 
            charge += (bitflags[c] & POS_RES) ? 1 : -1;
          if ((bitflags[c] & BORDER) 
              || a==0 || a==xsize-1 || b==0 || b==ysize-1) {
            charge = 0;
          }
          bitflags[c] |= mask_code;
          UpdateList(qual_map, a, b, value, phase, NULL, bitflags,
                xsize, ysize, index_list, &num_index,
                postponed_code, VISITED, postponed_code,
                max_list_size, unwrap_mode==dxdygrad, &min_qual);
        }
        /* reset pixels on list */
        for (k=0; k<num_index; k++) 
          bitflags[index_list[k]] &= ~(VISITED);
        for (k=0; k<xsize*ysize; k++) 
          bitflags[k] &= ~(VISITED);
      }
    }
  }
  printf("\n");
  free(index_list);
  return num_pieces;
}
