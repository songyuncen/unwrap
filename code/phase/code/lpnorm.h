#ifndef __LPNORM
#define __LPNORM
void LpNormUnwrap(float *soln, float *phase, float *dxwts, 
          float *dywts, unsigned char *bitflags, float *qual_map, 
          float *rarray, float *zarray, float *parray, int iter,
          int pcg_iter, double e0, int xsize, int ysize);
void RasterUnwrap(float *phase, float *soln, int xsize, int ysize);
void ResidualPhase(float *resid, float *phase, float *soln,
                   int xsize, int ysize);
#endif
