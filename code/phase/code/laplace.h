#ifndef __LAPLACE
#define __LAPLACE
void ComputeLaplacian(float *phase, float *laplacian, float *dxwts,
                   float *dywts, int xsize, int ysize, int laptype);
void ComputeDerivWts(float *phase, float *soln, float *dxwts,
                     float *dywts, float *qual_map, double e0,
                     int xsize, int ysize);
#endif
