/*
 *  fmg.c -- function for multigrid solution of weighted 
 *           least-squares phase-unwrapping problem
 */
#include <stdio.h>
#include <math.h>
#include "fmg.h"
#include "grid.h"
#include "dxdygrad.h"

/*  Call the functions for performing the multigrid phase-  */
/*  unwrapping algorithm.  If dywts is a null pointer, then */
/*  dxwts are copied into an array for dywts.               */
void MultigridUnwrap(float *soln, float *dx, float *dy, float *dxwts,
  float *dywts, int xsize, int ysize, int num_cycles, int num_iter)
{
  int  n, dywts_was_null=0, coarsest_dim=3;
  if (dywts==NULL) {
    dywts_was_null = 1;
    AllocateDouble(&dywts, xsize*ysize, "dy wts");
    for (n=0; n<xsize*ysize; n++) dywts[n] = dxwts[n];
  }
  for (n=0; n<num_cycles; n++) {
    printf("\nFMG CYCLE %d\n", n+1);
    FullMultigridVcycle(soln, dx, dy, dxwts, dywts,
                        xsize, ysize, num_iter, coarsest_dim);
  }
  if (dywts_was_null) free(dywts);
}
