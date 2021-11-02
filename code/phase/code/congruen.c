/*
 *  congruen.c -- make surface congruent to phase 
 */
#include <stdio.h>
#include <math.h>
#include "grad.h"
#include "util.h"

/* Make the surface in the array "surf" congruent to the wrapped */
/* phase (scaled to 0-1) in the array "phase".  Ignore the       */
/* zero-weights in the quality map array qual_map.               */
void CongruentSoln(float *surf, float *phase, float *qual_map,
                   int xsize, int ysize) 
{
  int     i, j, k, n, npos;
  int     kk, num, rnum, rnum_min;
  double  rmin, rr, rr_min, r;
  float   *ss;
  AllocateFloat(&ss, xsize*ysize, "temp array");
  /* ensure surface is positive */
  rmin = surf[0];
  for (k=0; k<xsize*ysize; k++) {
    if (rmin > surf[k]) rmin = surf[k];
  }
  npos = (rmin < 0) ? -rmin + 2.0 : 0.0;
  for (k=0; k<xsize*ysize; k++) {
    surf[k] += npos;
  }
  /* find best additive offset (to minimize discontinuities) */
  for (kk=0; kk<=10; kk++) {
    rr = 0.1*kk;
    /* create temporary surface by adding a constant rr */
    for (k=0; k<xsize*ysize; k++) {
      if (qual_map && qual_map[k]==0.0) continue; /* ignore 0-wts */ 
      n = (int)(surf[k] + rr);
      r = (surf[k] + rr) - n;
      r -= phase[k];
      if (r < -0.5) ss[k] = n - 1 + phase[k];
      else if (r > 0.5) ss[k] = n + 1 + phase[k];
      else  ss[k] = n + phase[k];
    }
    /* measure the discontinuities in the temporary surface */
    for (k=0, rnum=0; k<xsize*ysize; k++) {
      if (qual_map && qual_map[k]==0.0) continue; /* ignore 0-wts */
      i = k%xsize;
      j = k/xsize;
      if (i>0) {
        if (qual_map && qual_map[k-1]==0.0) continue; /* ignore */
        r = ss[k] - ss[k-1];
        if (r > 0.5 || r < -0.5) rnum++;
      }
      if (j>0) {
        if (qual_map && qual_map[k-xsize]==0.0) continue; /*ignore*/
        r = ss[k] - ss[k-xsize];
        if (r > 0.5 || r < -0.5) rnum++;
      }
    }
    /* save the rr that gives the minimum error (rnum) */
    if (kk==0 || rnum < rnum_min) {
      rnum_min = rnum;
      rr_min = rr;
    }
    printf("Offset %lf yields %d discontinuities\n", rr, rnum);
  }
  free(ss);
  printf("Adding offset %lf to surface\n", rr_min);
  for (k=0; k<xsize*ysize; k++) {
    surf[k] += rr_min; 
  }
  for (k=0; k<xsize*ysize; k++) {
    n = (int)surf[k];
    r = surf[k] - n;
    r -= phase[k];
    if (r < -0.5) surf[k] = n - 1 + phase[k];
    else if (r > 0.5) surf[k] = n + 1 + phase[k];
    else  surf[k] = n + phase[k];
    surf[k] -= npos;  /* also un-do addition of npos */
  }
}
