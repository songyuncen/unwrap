/*
 * getqual.c -- generate or process quality map and set quality mode
 */
#include <stdio.h>
#include <math.h>
#include "getqual.h"
#include "file.h"
#include "util.h"
#include "pi.h"
#include "qualvar.h"
#include "qualpseu.h"
#include "qualgrad.h"

/* Generate a new quality map from the phase data, or process */
/* the quality map that was input.                            */
void GetQualityMap(int mode, float *qual_map, float *phase, 
                   unsigned char *bitflags, int border_code,
                   int tsize, int xsize, int ysize)
{
  float  *temp;
  double rmin, rmax, rscale;
  int    i, j, k;
  /* process phase gradients */
  if (mode==variance) {
    AllocateFloat(&temp, xsize*ysize, "temp data");
    PhaseVariance(phase, qual_map, bitflags, border_code, temp,
                  tsize, xsize, ysize);
    free(temp);
    /* convert from cost to quality, and scale to interval (0,1) */
    for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k++) {
      if (rmin > qual_map[k]) rmin = qual_map[k];
      if (rmax < qual_map[k]) rmax = qual_map[k];
    }
    printf("Min & max phase derivative variance = %lf, %lf\n",
           rmin, rmax);
    rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0;
    for (k=0; k<xsize*ysize; k++) {
      qual_map[k] = (rmax - qual_map[k])*rscale;
      if (bitflags && (bitflags[k]&border_code)) qual_map[k] = 0.0;
    }
  }
  else if (mode==gradient) {
    AllocateFloat(&temp, xsize*ysize, "temp data");
    MaxPhaseGradients(phase, qual_map, bitflags, border_code, temp,
                      tsize, xsize, ysize);
    free(temp);
    /* convert from cost to quality, and scale to interval (0,1) */
    for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k++) {
      if (rmin > qual_map[k]) rmin = qual_map[k];
      if (rmax < qual_map[k]) rmax = qual_map[k];
    }
    printf("Min&max 'max phase gradient' = %lf %lf\n", rmin,rmax);
    rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0;
    for (k=0; k<xsize*ysize; k++) {
      qual_map[k] = (rmax - qual_map[k])*rscale;
      if (bitflags && (bitflags[k]&border_code)) qual_map[k] = 0.0;
    }
  }
  else if (mode==pseudocorr) {
    AllocateFloat(&temp, xsize*ysize, "temp data");
    PseudoCorrelation(phase, qual_map, bitflags, border_code, temp,
                      tsize, xsize, ysize);
    for (k=0; k<xsize*ysize; k++) {
      if (bitflags && (bitflags[k]&border_code)) qual_map[k] = 0.0;
    }
    free(temp);
  }
  else if (mode==none) {
    for (k=0; k<xsize*ysize; k++) {
      qual_map[k] = 1.0;
      if (bitflags && (bitflags[k]&border_code)) qual_map[k] = 0.0;
    }
  }
  else { /* quality map was input */
    /* scale to interval (0,1) */
    for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k++) {
      if (rmin > qual_map[k]) rmin = qual_map[k];
      if (rmax < qual_map[k]) rmax = qual_map[k];
    }
    printf("Min & max corr. coeff. = %lf, %lf\n", rmin, rmax);
    rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0;
    for (k=0; k<xsize*ysize; k++) {
      qual_map[k] = (qual_map[k] - rmin)*rscale;
      if (bitflags && (bitflags[k]&border_code)) qual_map[k] = 0.0;
    }
  }
}

/* Determine the quality mode based on a keyword, and return a */
/* corresponding integer.  If keyword unrecognized, return -1.  */
/* Allow the keyword "none" only if allow_none = 1.            */
int SetQualityMode(char *modekey, char *qualfile, int allow_none)
{
  int mode;
  if (Keyword(modekey, "none")) {
    mode = none;
    printf("No weights supplied\n");
    if (!allow_none) {
      fprintf(stderr, "Invalid mode: %s\n", modekey);
      mode = -1;  /* "none" not allowed */
    }
  }
  else if (Keyword(modekey, "max_corr")) {
    mode = corr_coeffs;
    if (Keyword(qualfile, "none")) {
      fprintf(stderr, "Must supply quality image for this mode\n");
      mode = -1;
    }
  }
  else if (Keyword(modekey, "min_grad")) {
    if (!Keyword(qualfile, "none"))
      printf("Correlation image file will be ignored\n");
    mode = gradient;
  }
  else if (Keyword(modekey, "min_var")) {
    if (!Keyword(qualfile, "none"))
      printf("Correlation image file will be ignored\n");
    mode = variance;
  }
  else if (Keyword(modekey, "max_pseu")) {
    if (!Keyword(qualfile, "none"))
      printf("Correlation image file will be ignored\n");
    mode = pseudocorr;
  }
  else {
    fprintf(stderr, "Invalid mode = '%s'.  Must be\n", modekey); 
    fprintf(stderr, "max_corr, min_var, max_grad, max_pseu "
                    "or none.\n");
    mode = -1;  /* invalid */
  }
  return mode;
}
