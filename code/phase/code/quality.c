/*
 *  quality.c -- function for quality-guided path-following
 */
#include <stdio.h>
#include <math.h>
#include "quality.h"
#include "list.h"
#include "util.h"

/* Main function for quality-guided algorithm for phase unwrapping */
int QualityGuidedPathFollower(float *phase, float *qual_map,
                            unsigned char *bitflags, float *soln, 
                            int xsize, int ysize, int avoid_code, 
                            UnwrapMode unwrap_mode, int debug_flag,
                            char *infile)
{
  int  i, j, k, a, b, c, n, num_pieces=0;
  float  value;
  float  min_qual, small_val = -1.0E+10;
  int    num_index, max_list_size, bench, benchout;
  int    *index_list;
  int    postponed_code=POSTPONED, unwrapped_code=UNWRAPPED;

  bench = xsize*ysize/100;
  if (bench < 1) bench = 1;
  benchout = 10*bench;
  min_qual = small_val;
  max_list_size = (xsize + ysize);
  AllocateInt(&index_list, max_list_size + 1, 
              "bookkeeping list (index)");
  avoid_code |= unwrapped_code;
  /* find starting point */
  n = 0;
  num_index = 0;
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      if (!(bitflags[k] & avoid_code)) {
        if (bitflags[k] & postponed_code) 
          /* soln[k] already stores the unwrapped value */
          value = soln[k];
        else {
          ++num_pieces;
          value = soln[k] = phase[k];
        } 
        bitflags[k] |= unwrapped_code;
        bitflags[k] &= ~postponed_code;
        UpdateList(qual_map, i, j, value, phase, soln, bitflags,
                   xsize, ysize, index_list, &num_index, avoid_code, 
                   unwrapped_code, postponed_code, max_list_size,
                   unwrap_mode==dxdygrad, &min_qual);
        while (num_index > 0) {
          if (n%bench==0) {
            printf("%d ", n/bench);
            fflush(stdout);
            if (0 && debug_flag && n%benchout==0 && n>0) {
              char filename[300];
              static char file_num = 0;
              sprintf(filename, "%s.path-%d", infile, ++file_num);
              SaveByteToImage(bitflags, "unwrapping path", filename,
                              xsize, ysize, 1, 1, avoid_code);
            }
          }
          ++n;
          if (!GetNextOneToUnwrap(&a, &b, index_list, &num_index,
                xsize, ysize)) break;   /* no more to unwrap */
          c = b*xsize + a;
          bitflags[c] |= unwrapped_code;
          bitflags[c] &= ~postponed_code;
          value = soln[c];
          UpdateList(qual_map, a, b, value, phase, soln, bitflags,
                   xsize, ysize, index_list, &num_index, avoid_code, 
                   unwrapped_code, postponed_code, max_list_size,
                   unwrap_mode==dxdygrad, &min_qual);
          if (num_index <= 0) {
            min_qual = small_val;
            for (c=0; c<xsize*ysize; c++) {
              if ((bitflags[c] & postponed_code) 
                       && !(bitflags[c] & unwrapped_code)) {
                a = c%xsize;    
                b = c/xsize;
                InsertList(soln, soln[c], qual_map, bitflags, a, b,
                           index_list, &num_index, xsize, ysize,
                           unwrapped_code, postponed_code,
                           &min_qual, max_list_size);
              }
            }
            if (num_index > 0) min_qual = qual_map[index_list[0]];
          }
        }  /* while ... */
      }
    }   /* for (b ...  */
  }   /* for (a ... */
  printf("\n");
  free(index_list);
  return num_pieces;
}
