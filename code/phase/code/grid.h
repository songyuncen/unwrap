#ifndef __GRID
#define __GRID
void FullMultigridVcycle(float *soln, float *dx, float *dy,
                         float *dxwts, float *dywts, int w, int h,
                         int numit, int coarsest_dim);
void Vcycle(float *soln, float *dx, float *dy, float *dxwts,
            float *dywts, int w, int h, int numit, int coarsest_dim);
void RestrictDxwts(float *coarse, int wc, int hc, float *fine,
                   int wf, int hf);
void RestrictDywts(float *coarse, int wc, int hc, float *fine,
                   int wf, int hf);
#endif
