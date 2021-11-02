#ifndef __MASKCUT
#define __MASKCUT
#include "getqual.h"
int QualityGuidedMask(float *phase, float *qual_map,
                      unsigned char *bitflags, int xsize, int ysize,
                      int mask_code, UnwrapMode unwrap_mode,
                      int debug_flag, char *infile);
#endif
