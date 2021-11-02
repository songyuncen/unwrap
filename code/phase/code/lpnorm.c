/*
 *  lpnorm.c -- functions for minimum Lp-norm phase-unwrapping
 */
#include <stdio.h>
#include <math.h>
#include "pcg.h"
#include "raster.h"
#include "residues.h"
#include "laplace.h"
#include "lpnorm.h"
#include "util.h"
#include "congruen.h"
#define POS_RES      0x01
#define NEG_RES      0x02
#define BORDER       0x20
#define ZERO_WEIGHT  0x40

/* Main iteration of minimum Lp-norm phase unwrapping algorithm */
void LpNormUnwrap(float *soln, float *phase, float *dxwts,
             float *dywts, unsigned char *bitflags, float *qual_map,
             float *rarray, float *zarray, float *parray, int iter,
             int pcg_iter, double e0, int xsize, int ysize) 
{
  int  i, j, k, n;
  float  *residual;
  residual = rarray;  /* borrow rarray to compute the residual */
  ResidualPhase(residual, phase, soln, xsize, ysize);
  n = Residues(residual, bitflags, POS_RES, NEG_RES, 
               BORDER, xsize, ysize);
  for (k=0; k<iter && n>0; k++) {
    printf("\nIter %d: %d residues\n", k+1, n);
    ComputeDerivWts(phase, soln, dxwts, dywts, qual_map, e0,
                    xsize, ysize);
    ComputeLaplacian(phase, rarray, dxwts, dywts, xsize, ysize, 1);
    /* borrow zarray temporarily to compute Laplacian of soln */
    ComputeLaplacian(soln, zarray, dxwts, dywts, xsize, ysize, 0);
    for (i=0; i<xsize*ysize; i++) rarray[i] -= zarray[i];
    PCGUnwrap(rarray, zarray, parray, soln, dxwts, dywts,
              xsize, ysize, pcg_iter, 0.0);
    ResidualPhase(residual, phase, soln, xsize, ysize);
    for (i=0; i<xsize*ysize; i++) bitflags[i] &= ~(POS_RES | NEG_RES);
    n = Residues(residual, bitflags, POS_RES, NEG_RES, 
                 BORDER, xsize, ysize);
  }
  printf("%d residues after %d iter.  Processing residual...\n", 
         n, k);
  /* see if the entire nonmasked array is residue-free */
  for (i=0; i<xsize*ysize; i++) bitflags[i] &= ~(POS_RES | NEG_RES);
  n = Residues(residual, NULL, POS_RES, NEG_RES, 
               BORDER, xsize, ysize);
  if (n==0) { 
    /* Unwrap residual and add to solution */
    /* (This could be improved by unwrapping only the nonmasked */
    /* pixels, and removing the nonmasked residue test above)   */
    RasterUnwrap(residual, zarray, xsize, ysize); /* borrow zarray */
    for (k=0; k<xsize*ysize; k++)  soln[k] += zarray[k];  
  }
  else {
    /* Make solution congruent to wrapped input phase */
    CongruentSoln(soln, phase, NULL, xsize, ysize);
  }
}
 
/* Compute the residual phase, i.e., the wrapped difference */
/* between the (partially) unwrapped solution "soln" and    */
/* and the wrapped input phase "phase".                     */
void ResidualPhase(float *resid, float *phase, float *soln,
                   int xsize, int ysize)
{
  int    n, k;
  double r;
  for (k=0; k<xsize*ysize; k++) {
    r = phase[k] - soln[k];
    if (r < 0) r += (2 - ((int)r));  /* make r positive */
    r -= (int)r;  /* get fractional part of r */
    if (r > 0.5) r -= 1.0;  /* wrap to range (-0.5, 0.5) */
    resid[k] = r;
  }
}
