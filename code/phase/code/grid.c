/*
 * grid.c - functions for (weighted) multigrid phase unwrapping 
 */
#include <stdio.h>
#include <math.h>
#include "grid.h"
#include "gridops.h"
#include "gridmem.h"
#include "relax.h"

/* Main function for weighted multigrid (called recursively) */
void FullMultigridVcycle(float *soln, float *dx, float *dy,
   float *dxwts, float *dywts, int w, int h, int numit, int mindim)
{
  float  *dx2=NULL, *dy2=NULL, *soln2=NULL;
  float  *dxwts2=NULL, *dywts2=NULL;
  int    w2 = w/2, h2 = h/2;
  if (!Coarsest(w, h, mindim)) {
    dxwts2 = Allocate(w2, h2, dxwts_type);
    dywts2 = Allocate(w2, h2, dywts_type);
    dx2 = Allocate(w2, h2, dx_type);
    dy2 = Allocate(w2, h2, dy_type);
    soln2 = Allocate(w2, h2, soln_type);
    RestrictDxwts(dxwts2, w2, h2, dxwts, w, h);
    RestrictDywts(dywts2, w2, h2, dywts, w, h);
    Restrict(dx2, dy2, w2, h2, dx, dy, dxwts, dywts, soln, w, h); 
    Zero(soln2, w2, h2);
    FullMultigridVcycle(soln2, dx2, dy2, dxwts2, dywts2, w2, h2,
                        numit, mindim);
    ProlongAndAccumulate(soln, w, h, soln2, w2, h2, dxwts2, dywts2);
  }
  /* perform V-cycle multigrid on fine grid */
  Vcycle(soln, dx, dy, dxwts, dywts, w, h, numit, mindim);
}

/* V-cycle multigrid algorithm (called recursively) */
void Vcycle(float *soln, float *dx, float *dy, float *dxwts,
            float *dywts, int w, int h, int numit, int mindim)
{
  float *dx2=NULL, *dy2=NULL, *soln2=NULL;
  float *dxwts2=NULL, *dywts2=NULL;
  int    w2 = w/2, h2 = h/2;
  if (!Coarsest(w, h, mindim)) {
    Relax(soln, dx, dy, dxwts, dywts, w, h, numit);
    dxwts2 = Allocate(w2, h2, dxwts_type);
    dywts2 = Allocate(w2, h2, dywts_type);
    dx2 = Allocate(w2, h2, dx_type);
    dy2 = Allocate(w2, h2, dy_type);
    soln2 = Allocate(w2, h2, soln_type);
    RestrictDxwts(dxwts2, w2, h2, dxwts, w, h);
    RestrictDywts(dywts2, w2, h2, dywts, w, h);
    Restrict(dx2, dy2, w2, h2, dx, dy, dxwts, dywts, soln, w, h); 
    Zero(soln2, w2, h2);
    Vcycle(soln2, dx2, dy2, dxwts2, dywts2, w2, h2, numit, mindim); 
    ProlongAndAccumulate(soln, w, h, soln2, w2, h2, dxwts2, dywts2);
    Relax(soln, dx, dy, dxwts, dywts, w, h, numit);
  }
  else { /* coarsest */
    Relax(soln, dx, dy, dxwts, dywts, w, h, 2*w*h); 
  }
}

/* multigrid restriction operation for dx weights */
void RestrictDxwts(float *coarse, int wc, int hc, 
                   float *fine, int wf, int hf) 
{
  int     a, b, i, j, k, m, n;
  int     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k20;
  float   *f, *c;
  f = fine;         c = coarse;
  /*  Indexes:       */
  /*  k1   k2    k3  */
  /*  k4   k5=k  k6  */
  /*  k7   k8    k9  */
  for (j=0; j<hc; j++) {
    b = 2*j;
    for (i=0; i<wc; i++) {
      a = 2*i;
      k = b*wf + a;
      k5 = k;
      k6 = (a < wf - 1) ? k5 + 1 : k5 - 1;
      k4 = (a > 0) ? k5 - 1 : k5 + 1;
      if (b < hf - 1) {
        k7 = k4 + wf;   k8 = k5 + wf;    k9 = k6 + wf;
      }
      else {
        k7 = k4 - wf;   k8 = k5 - wf;    k9 = k6 - wf;
      }
      if (b > 0) {
        k1 = k4 - wf;   k2 = k5 - wf;   k3 = k6 - wf;
      }
      else {
        k1 = k4 + wf;   k2 = k5 + wf;   k3 = k6 + wf;
      }
      if (f[k5]==0.0 && f[k9]==0.0) c[j*wc + i] = 0.0;
      else if (f[k6]==0.0 && f[k8]==0.0) c[j*wc + i] = 0.0;
      else if (f[k6]==0.0 && f[k9]==0.0) c[j*wc + i] = 0.0;
      else if (f[k5]==0.0 && f[k8]==0.0) c[j*wc + i] = 0.0;
      else if ((i==0 || c[j*wc + i - 1] != 0.0) &&
        ((f[k4]==0.0 && f[k8]==0.0) 
          ||  (f[k5]==0.0 && f[k7]==0.0)))
                                     c[j*wc + i] = 0.0;
      else c[j*wc + i] = 0.25*(f[k5] + f[k6] + f[k8] + f[k9]);
    }
  }
}

/* multigrid restriction operation for dy weights */
void RestrictDywts(float *coarse, int wc, int hc, 
                   float *fine, int wf, int hf) 
{
  int     a, b, i, j, k, m, n;
  int     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k20;
  float  *f, *c;
  f = fine;         c = coarse;
  /*  Indexes:       */
  /*  k1   k2    k3  */
  /*  k4   k5=k  k6  */
  /*  k7   k8    k9  */
  for (j=0; j<hc; j++) {
    b = 2*j;
    for (i=0; i<wc; i++) {
      a = 2*i;
      k = b*wf + a;
      k5 = k;
      k6 = (a < wf - 1) ? k5 + 1 : k5 - 1;
      k4 = (a > 0) ? k5 - 1 : k5 + 1;
      if (b < hf - 1) {
        k7 = k4 + wf;    k8 = k5 + wf;    k9 = k6 + wf;
      }
      else {
        k7 = k4 - wf;    k8 = k5 - wf;    k9 = k6 - wf;
      }
      if (b > 0) {
        k1 = k4 - wf;    k2 = k5 - wf;    k3 = k6 - wf;
      }
      else {
        k1 = k4 + wf;    k2 = k5 + wf;    k3 = k6 + wf;
      }
      if (f[k5]==0.0 && f[k9]==0.0) c[j*wc + i] = 0.0;
      else if (f[k6]==0.0 && f[k8]==0.0) c[j*wc + i] = 0.0;
      else if (f[k8]==0.0 && f[k9]==0.0) c[j*wc + i] = 0.0;
      else if (f[k5]==0.0 && f[k6]==0.0) c[j*wc + i] = 0.0;
      else if ((j==0 || c[(j - 1)*wc + i] != 0.0) &&
          ((f[k2]==0.0 && f[k6]==0.0)
            ||  (f[k3]==0.0 && f[k5]==0.0)))
                                     c[j*wc + i] = 0.0;
      else c[j*wc + i] = 0.25*(f[k5] + f[k6] + f[k8] + f[k9]);
    }
  }
}     
