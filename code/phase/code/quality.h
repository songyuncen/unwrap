#ifndef __QUALITY
#define __QUALITY
#include "getqual.h"
int QualityGuidedPathFollower(float *phase, float *qual_map,
      unsigned char *bitflags, float *soln, int xsize, int ysize,
      int avoid_code, UnwrapMode unwrap_mode, int debug_mode,
      char *infile);
#endif
