/*
 *   path.c -- function for unwrapping around cuts
 */
#include <stdio.h>
#include <math.h>
#include "list.h"
#include "grad.h"
#include "path.h"

/* Unwrap the phase data (by Itoh's method) without crossing
 * any branch cuts.  Return number of disconnected pieces.
 */
int UnwrapAroundCuts(float *phase, unsigned char *bitflags, 
              float *soln, int xsize, int ysize, int cut_code,
              int debug_flag, char *infile)
{
  int  i, j, k, a, b, c, n, num_pieces=0;
  float  value;
  float  min_qual, small_val = -1.0E+10;
  int    num_index, max_list_size, bench, benchout;
  int    *index_list;
  int    unwrapped_code=UNWRAPPED, postponed_code = POSTPONED;
  int    avoid_code; 

  bench = xsize*ysize/100;
  if (bench < 1) bench = 1;
  benchout = 10*bench;
  min_qual = small_val;
  max_list_size = xsize*ysize; /* this size may be reduced */
  AllocateInt(&index_list, max_list_size + 1, 
              "bookkeeping list (index)");

  avoid_code = cut_code | unwrapped_code | BORDER;

  /* find starting point */
  n = 0;
  num_index = 0;
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      if (!(bitflags[k] & avoid_code)) {
        bitflags[k] |= unwrapped_code;
        if (bitflags[k] & postponed_code)
          /* soln[k] already stores the unwrapped value */
          value = soln[k];
        else {
          ++num_pieces;
          value = soln[k] = phase[k];
        }
        UpdateList(NULL, i, j, value, phase, soln, bitflags, xsize,
                   ysize, index_list, &num_index, avoid_code,
                   unwrapped_code, postponed_code, max_list_size,
                   0, &min_qual);
        while (num_index > 0) {
          if (n%bench==0) {
            printf("%d ", n/bench);
            fflush(stdout);
            if (0 && debug_flag && n%benchout==0 && n>0) {
              char filename[300];
              static char file_num = 0;
              sprintf(filename, "%s.fill-%d", infile, ++file_num);
              SaveByteToImage(bitflags, "fill path", filename,
                              xsize, ysize, 1, 1, avoid_code);
            }
          }
          ++n;
          if (!GetNextOneToUnwrap(&a, &b, index_list,
                                  &num_index, xsize, ysize))
            break;   /* no more to unwrap */
          c = b*xsize + a;
          bitflags[c] |= unwrapped_code;
          value = soln[c];        
          UpdateList(NULL, a, b, value, phase, soln, bitflags,
                     xsize, ysize, index_list, &num_index,
                     avoid_code, unwrapped_code, postponed_code,
                     max_list_size, 0, &min_qual);
        }
      }
    }
  }
  printf("\n");
  free(index_list);
  /* unwrap branch cut pixels */
  for (j=1; j<ysize; j++) {
    for (i=1; i<xsize; i++) {
      k = j*xsize + i;
      if (bitflags[k] & cut_code) {
        if (!(bitflags[k-1] & cut_code)) 
          soln[k] = soln[k-1] + Gradient(phase[k], phase[k-1]);
        else if (!(bitflags[k-xsize] & cut_code))
          soln[k] 
             = soln[k-xsize] + Gradient(phase[k], phase[k-xsize]);
      }
    }
  }
  return num_pieces;
}
