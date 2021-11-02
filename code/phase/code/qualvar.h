#ifndef __QUALVAR
#define __QUALVAR
void PhaseVariance(float *phase, float *result, 
                   unsigned char *bitflags, int ignore_code, 
                   float *soln, int tsize, int xsize, int ysize);
void DxGradVar(float *dx, float *dxvar, int xsize, int ysize,
  int size, unsigned char *bitflags, int avoid_code, int add_flag);
void DyGradVar(float *dy, float *dyvar, int xsize, int ysize,
  int size, unsigned char *bitflags, int avoid_code, int add_flag);
#endif
