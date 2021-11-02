#ifndef __QUALGRAD
#define __QUALGRAD
void MaxPhaseGradients(float *phase, float *result,
                      unsigned char *bitflags, int ignore_code,
                      float *soln, int tsize, int xsize, int ysize);
void DxGradMax(float *dx, float *dxvar, int xsize, int ysize,
   int size, unsigned char *bitflags, int avoid_code, int add_flag);
void DyGradMax(float *dy, float *dyvar, int xsize, int ysize,
   int size, unsigned char *bitflags, int avoid_code, int add_flag);
#endif
