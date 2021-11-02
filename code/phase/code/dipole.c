/*
 *   dipole.c -- function for finding dipoles and eliminating them
 */
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "brcut.h"
#include "list.h"

/* Find the dipoles (adjoining residues of opposite sign) in   */
/* the bitflags array, remove them, and place branch cut.      */
/* The bits for the positive and negative residues are defined */
/* by POS_RES and NEG_RES (defined in brcut.h), while the      */
/* branch cut bits are defined by "branchcut_code".            */
void Dipole(unsigned char *bitflags, int xsize, int ysize,
            int branchcut_code)
{
  int  i, j, k, kk;
  printf("Dipoles...\n");
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      kk = 0;
      if ((bitflags[k] & POS_RES)) {
        if (i<xsize-1 && (bitflags[k+1] & NEG_RES)) kk=k+1;
        else if (j<ysize-1) {
          if ((bitflags[k+xsize] & NEG_RES)) kk=k+xsize;
        }
      }
      else if ((bitflags[k] & NEG_RES)) {
        if (i<xsize-1 && (bitflags[k+1] & POS_RES)) kk=k+1;
        else if (j<ysize-1) {
          if ((bitflags[k+xsize] & POS_RES)) kk=k+xsize;
        }
      }
      if (kk) {
        printf("Connecting dipoles %d,%d - %d,%d\n", i, j, kk%xsize,
               kk/xsize);
        PlaceCut(bitflags, i, j, kk%xsize, kk/xsize,
                 xsize, ysize, branchcut_code);
        bitflags[k] &= (~(RESIDUE)); 
        bitflags[kk] &= (~(RESIDUE)); 
      }
    }
  }
}
