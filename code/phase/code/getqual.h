#ifndef __GETQUAL
#define __GETQUAL
typedef enum {none, gradient, variance, pseudocorr,
              corr_coeffs, dxdygrad
} UnwrapMode;
void GetQualityMap(int mode, float *qual_map, float *phase,
                   unsigned char *bitflags, int border_code,
                   int tsize, int xsize, int ysize);
int SetQualityMode(char *mode_key, char *qual_file, int allow_none);
#endif
