/*
 * ungrid.c - functions for UNWEIGHTED multigrid phase
 *            unwrapping 
 */
#include <stdio.h>
#include <math.h>
#include "ungrid.h"
#include "gridops.h"
#include "gridmem.h"
#include "relax.h"

/* Main function for unweighted multigrid phase unwrapping. */
/* (Called recursively.)                                    */
void Ufmg(float *soln, float *dx, float *dy,
          int w, int h, int numit, int mindim)
{
  float  *dx2=NULL, *dy2=NULL, *soln2=NULL;
  int    w2 = w/2, h2 = h/2;
  if (!Coarsest(w, h, mindim)) {
    dx2 = Allocate(w2, h2, dx_type);
    dy2 = Allocate(w2, h2, dy_type);
    soln2 = Allocate(w2, h2, soln_type);
    Restrict(dx2, dy2, w2, h2, dx, dy, NULL, NULL, soln, w, h);
    Zero(soln2, w2, h2);
    Ufmg(soln2, dx2, dy2, w2, h2, numit, mindim);
    ProlongAndAccumulate(soln, w, h, soln2, w2, h2, NULL, NULL);
  }
  /* perform V-cycle multigrid on fine grid */
  Umv(soln, dx, dy, w, h, numit, mindim);
}

/* V-cycle function for unweighted multigrid phase unwrapping. */
/* (Called recursively.)                                       */
void Umv(float *soln, float *dx, float *dy,
         int w, int h, int numit, int mindim)
{
  float *dx2=NULL, *dy2=NULL, *soln2=NULL;
  int    w2 = w/2, h2 = h/2;
  if (!Coarsest(w, h, mindim)) {
    Relax(soln, dx, dy, NULL, NULL, w, h, numit);
    dx2 = Allocate(w2, h2, dx_type);
    dy2 = Allocate(w2, h2, dy_type);
    soln2 = Allocate(w2, h2, soln_type);
    Restrict(dx2, dy2, w2, h2, dx, dy, NULL, NULL, soln, w, h); 
    Zero(soln2, w2, h2);
    Umv(soln2, dx2, dy2, w2, h2, numit, mindim); 
    ProlongAndAccumulate(soln, w, h, soln2, w2, h2, NULL, NULL);
    Relax(soln, dx, dy, NULL, NULL, w, h, numit);
  }
  else { /* coarsest */
    Relax(soln, dx, dy, NULL, NULL, w, h, 2*w*h); 
  }
}
