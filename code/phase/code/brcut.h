#ifndef __BRCUT
#define __BRCUT
void PlaceCut(unsigned char *array, int a, int b, int c, int d,
              int xsize, int ysize, int code);
int DistToBorder(unsigned char *bitflags, int border_code,
            int a, int b, int *ra, int *rb, int xsize, int ysize);
#endif
