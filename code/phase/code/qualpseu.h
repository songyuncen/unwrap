#ifndef __QUALPSEU
#define __QUALPSEU
void PseudoCorrelation(float *phase, float *result,
                     unsigned char *bitflags, int ignore_code,
                     float *temp1, int tsize, int xsize, int ysize);
void SqrAvgFilter(float *in, float *out, int xsize, int ysize,
                  int size, unsigned char *bitflags, int avoid_code,
                  int add_flag);
#endif
