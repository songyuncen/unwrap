/*
 * view.c -- utility functions for viewing surface data
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pi.h"
#include "view.h"

/* Generate shaded bird's eye view of surface.              */
/* elev_angle = elevation angle of light vector in radians  */
/* (elev_angle=0.0 results in a default angle of 45 deg.)   */
void Shade(float *surf, char *image, int xsize, int ysize,
           double elev_angle)
{
  int    i, j, k;
  double angle, ht, r, ca, sa, ch, sh, lightness;
  ca = cos(elev_angle);
  sa = sin(elev_angle);
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      k = j*xsize + i;
      if (i==xsize - 1) {
        image[k] = 0;
      }
      else {
        ht = surf[k] - surf[k+1];
        if (ht < 0.0) {
          image[k] = 0;
        } 
        else {
          /* Dot product of surface normal with light vector. */
          /* Assumes light vector perpendicular to array columns. */
          r = sqrt(1.0 + ht*ht);
          ch = 1.0/r;
          sh = ht/r;
          lightness = ch*ca + sh*sa;
          image[k] = 255.0*lightness; 
        }
      }
    }
  }
}

/* Generate image of surface elevations.                  */
void SurfImage(float *surf, char *image, int xsize, int ysize)
{
  int     i, j, k;
  double  rmin, rmax, rscale;
  for (k=0, rmax = -(rmin = 1.0e+20); k<xsize*ysize; k++) {
    if (rmin > surf[k]) rmin = surf[k];
    if (rmax < surf[k]) rmax = surf[k];
  }
  rscale = (rmin < rmax) ? 1.0/(rmax - rmin) : 0.0;
  for (k=0, rmax = -(rmin = 1.0e+20); k<xsize*ysize; k++) {
    image[k] = 255.0*rscale*(surf[k] - rmin);
  } 
}

/* Generate image of rewrapped phase.              */
/* Assumes surface values are scaled so that the   */
/* interval 0-TWOPI is scaled to the interval 0-1. */
void Rewrap(float *surf, char *image, int xsize, int ysize)
{
  int     i, j, k;
  double  r;
  for (k=0; k<xsize*ysize; k++) {
    r = surf[k];
    if (r < 0) r += (int)(-r) + 2;  /* ensure r > 0 */
    r -= (int) r;
    image[k] = 256.0*r;
  }
}

