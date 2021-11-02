/*
 * gold.c -- generate branch cuts by Goldstein's algorithm
 */
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "brcut.h"
#include "list.h"

/* Goldstein's phase-unwrapping algorithm.  The bitflags store */
/* the masked pixels (to be ignored) and the residues and      */
/* accumulates other info such as the branch cut pixels.       */
void GoldsteinBranchCuts(unsigned char *bitflags, int MaxCutLen,
               int NumRes, int xsize, int ysize, int branchcut_code)
{
  int            i, j, k, ii, jj, kk, m, n, ri, rj;
  int            charge, boxctr_i, boxctr_j, boxsize, bs2;
  int            dist, min_dist, rim_i, rim_j, near_i, near_j;
  int            ka, num_active, max_active, *active_list;
  int            bench; 
  int            draw_cut_line;
  double         r;
  bench = ysize/100;
  if (bench < 1) bench = 1;
  if (MaxCutLen < 2) MaxCutLen = 2;
  max_active = NumRes + 10;
  AllocateInt(&active_list, max_active + 1, "book keeping data");
  /* branch cuts */
  printf("Computing branch cuts\n");
  for (j=0; j<ysize; j++) {
    if (j%bench==0) {
      printf("%d ", j/bench);
      fflush(stdout);
    }
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      if ((bitflags[k] & (POS_RES | NEG_RES))
               && !(bitflags[k] & VISITED)) {
        bitflags[k] |= VISITED;  /* turn on visited flag */
        bitflags[k] |= ACTIVE;   /* turn on active flag */
        charge = (bitflags[k] & POS_RES) ? 1 : -1;
        num_active = 0;
        active_list[num_active++] = k;
        if (num_active > max_active) num_active = max_active;
        for (boxsize = 3; boxsize<=2*MaxCutLen; boxsize += 2) {
          bs2 = boxsize/2;
          for (ka=0; ka<num_active; ka++) {
            boxctr_i = active_list[ka]%xsize;
            boxctr_j = active_list[ka]/xsize;
            for (jj=boxctr_j - bs2; jj<=boxctr_j + bs2; jj++) {
              for (ii=boxctr_i - bs2; ii<=boxctr_i + bs2; ii++) {
                kk = jj*xsize + ii;
                if (ii<0 || ii>=xsize || jj<0 || jj>=ysize) { 
                  continue;
                }
                else { 
                  if (ii==0 || ii==xsize-1 || jj==0 || jj==ysize-1
                        || (bitflags[kk] & BORDER)) {
                    charge = 0;
                    DistToBorder(bitflags, BORDER, boxctr_i,
                             boxctr_j, &ri, &rj, xsize, ysize);
                    PlaceCut(bitflags, ri, rj, boxctr_i, boxctr_j,
                             xsize, ysize, branchcut_code);
                  }
                  else if ((bitflags[kk] & (POS_RES | NEG_RES))
                                   && !(bitflags[kk] & ACTIVE)) {
                    if (!(bitflags[kk] & VISITED)) {
                      charge += (bitflags[kk] & POS_RES) ? 1 : -1;
                      bitflags[kk] |= VISITED;   /* set flag */
                    }
                    active_list[num_active++] = kk;
                    if (num_active > max_active) 
                           num_active = max_active;
                    bitflags[kk] |= ACTIVE;  /* set active flag */
                    PlaceCut(bitflags, ii, jj, boxctr_i, boxctr_j,
                             xsize, ysize, branchcut_code);
                  }
                  if (charge==0)
                    goto continue_scan;
                }  /* else */
              }   /* for (ii ... */
            }   /* for (jj ... */
          }  /* for (ka ... */  
        }   /* for (boxsize ... */  

        if (charge != 0) {   /* connect branch cuts to rim */
          min_dist = xsize + ysize;  /* large value */
          for (ka=0; ka<num_active; ka++) {
            ii = active_list[ka]%xsize;
            jj = active_list[ka]/xsize;
            if ((dist = DistToBorder(bitflags, BORDER,
                        ii, jj, &ri, &rj, xsize, ysize))<min_dist) {
              min_dist = dist;
              near_i = ii;
              near_j = jj;
              rim_i = ri;
              rim_j = rj;
            }
          } 
          PlaceCut(bitflags, near_i, near_j, rim_i, rim_j,
                   xsize, ysize, branchcut_code);
        } 
        continue_scan :
        /* mark all active pixels inactive */
        for (ka=0; ka<num_active; ka++) 
          bitflags[active_list[ka]] &= ~ACTIVE;  /* turn flag off */
      }  /* if (bitflags ... */
    }  /* for (i ... */
  }  /* for (j ... */
  printf("\n");
  free(active_list);
  return;
} 
