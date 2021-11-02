#ifndef __FLYNN
#define __FLYNN
#define isign(A,B) ((B) >= 0 ? (A) : (-(A)))
#define nint(A) ((A)>0. ? (int)((A)+0.5) : (int)((A)-0.5))
void FlynnMinDisc(float *phase, int *value, unsigned char *bitflags,
                  float *qual_map, short *vjump, short *hjump,
                  int xsize, int ysize);
#endif
