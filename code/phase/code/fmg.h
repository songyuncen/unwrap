#ifndef __FMG
#define __FMG
void MultigridUnwrap(float *soln, float *dx, float *dy,
                     float *dxwts, float *dywts, int xsize,
                     int ysize, int num_cycles, int num_iter);
#endif
