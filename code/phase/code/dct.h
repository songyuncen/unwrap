#ifndef __DCT
#define __DCT
#include "pi.h"
#define BLOCKSIZE     50
#define MAX_NUM_ROWS  2050
void FastCosineTransform(float *array, int xsize, int ysize,
                         int sign);
void RowDcts(float *array, int xsize, int ysize, int sign);
void ColDcts(float *array, int xsize, int ysize, int sign,
             float *vector, float *buffer);
void FastDct(float *array, int size, int sign);
void RealFft(float *array, int size, int sign);
void Fft(float *array, int size, int sign);
#endif
