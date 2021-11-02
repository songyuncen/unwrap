#ifndef __PCG
#define __PCG
void PCGUnwrap(float *rarray, float *zarray, float *parray, 
               float *soln, float *dxwts, float *dywts, int xsize,
               int ysize, int max_iter, double epsi_con);
void PCGIterate(float *rarray, float *zarray, float *parray,
                float *soln, float *dxwts, float *dywts, int xsize,
                int ysize, double *xcos, double *ycos, int iloop,
                double sum0, double *alpha, double *beta,
                double *beta_prev, double *epsi);
#endif
