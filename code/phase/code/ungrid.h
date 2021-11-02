#ifndef __UNGRID
#define __UNGRID
void Ufmg(float *soln, float *dx, float *dy, int w, int h,
          int numit, int coarsest_dim);
void Umv(float *soln, float *dx, float *dy, int w, int h,
         int numit, int coarsest_dim);
#endif
