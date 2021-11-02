#ifndef __GRIDOPS
#define __GRIDOPS
void Restrict(float *dx2, float *dy2, int wc, int hc,
              float *dx, float *dy, float *dxwts,
              float *dywts, float *soln, int wf, int hf);
void ProlongAndAccumulate(float *fine, int wf, int hf,
                          float *coarse, int wc, int hc,
                          float *coarse_dxwts, float *coarse_dywts);
int Coarsest(int w, int h, int coarsest_dim);
void Zero(float *soln, int w, int h);
#endif
