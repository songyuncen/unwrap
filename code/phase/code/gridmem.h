#ifndef __GRIDMEM
#define __GRIDMEM
#define MAX_ARRAYS   100
#define NUM_TYPES      5
typedef enum {
  dx_type, dy_type, soln_type, dxwts_type, dywts_type
} ArrayType;
float *Allocate(int w, int h, ArrayType type);
void FreeAll();
#endif
