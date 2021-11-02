#ifndef __VIEW
#define __VIEW
void Shade(float *surf, char *image, int xsize, int ysize,
           double elev_angle);
void SurfImage(float *surf, char *image, int xsize, int ysize);
void Rewrap(float *surf, char *image, int xsize, int ysize);
#endif
