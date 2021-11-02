/*
 *  unfmg.c -- function for multigrid solution of UNWEIGHTED
 *             least-squares phase-unwrapping problem
 */
#include <stdio.h>
#include <math.h>
#include "unfmg.h"
#include "ungrid.h"
#include "dxdygrad.h"
/* Unweighted multigrid function.                */
/* Initialize soln array to zero before calling. */
void UnweightedMultigridUnwrap(float *soln, float *dx, float *dy,
                int xsize, int ysize, int num_cycles, int num_iter)
{
  int  n, coarsest_dim=3;
  for (n=0; n<num_cycles; n++) {
    printf("\nFMG CYCLE %d\n", n+1);
    Ufmg(soln, dx, dy, xsize, ysize, num_iter, coarsest_dim);
  }
}
