/*
 * extract.c -- extract phase from input data
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "extract.h"
#include "pi.h"
/*
 *  Get the phase data from the input data.  Allocate the temporary
 *  memory required and read the input file.
 *
 *  in_format = 0 for 8-byte complex data,
 *              1 for 4-byte complex data,
 *              2 for 1-byte phase data,
 *              3 for 4-byte float phase data
 *
 *  Output: normalized (0 - 1) phase values in array "phase" 
 */
void GetPhase(int in_format, FILE *ifp, char *infile, float *phase,
              int xsize, int ysize)
{
  void  *in_data;
  if (in_format==0)  /* floating-pt real & imaginary parts */
    AllocateFloat((float**)(&in_data), 2*xsize*ysize, "input data");
  else if (in_format==1)  /* short integer real & imaginary parts */
    AllocateShort((short**)(&in_data), 2*xsize*ysize, "input data");
  else if (in_format==2)  /* 1-byte (quantized) phase values */
    AllocateByte((unsigned char **)(&in_data), xsize*ysize,
                 "input data");
  else                    /* floating-pt phase values */
    AllocateFloat((float **)(&in_data), xsize*ysize, "input data");

  printf("Reading input data...\n");
  if (in_format==1)      /* short int real and imaginary parts */
    ReadShort(ifp, (short *)in_data, 2*xsize*ysize, infile);
  else if (in_format==0) /* floating-pt real and imag parts */
    ReadFloat(ifp, (float *)in_data, 2*xsize*ysize, infile);
  else if (in_format==3) /* floating-pt phase */
    ReadFloat(ifp, (float *)in_data, xsize*ysize, infile);
  else                   /* 1-byte quantized phase */
    ReadByte(ifp, (unsigned char *)in_data, xsize*ysize, infile);

  ExtractPhase(in_format, in_data, phase, xsize, ysize, 1);
  free(in_data);
}

/* Extract the phase from the input data.  If status_flag is 0 */
/* then do not print any status.                               */
void ExtractPhase(int in_format, void *in_data, float *phase,
                  int xsize, int ysize, int status_flag)
{
  int            i, j;
  double         x, y, r, angle;
  static double  one_over_twopi = 1.0/TWOPI;
  static double  scale;
  float          *in8_data = in_data;
  short          *in4_data = in_data;
  unsigned char  *quantized_phase = in_data;
  float          *float_phase = in_data;
  if (in_format==0 || in_format==1) {
    for (j=0; j<ysize; j++) {
      if (status_flag && j%100==0) {
        printf("%d ", j);
        fflush(stdout);
        if (j+100 >= ysize) printf("\n");
      }
      for (i=0; i<xsize; i++) {
        if (in_format==1) {
          x = in4_data[2*(j*xsize + i)];
          y = in4_data[2*(j*xsize + i) + 1];
        }
        else {
          x = in8_data[2*(j*xsize + i)];
          y = in8_data[2*(j*xsize + i) + 1];
        }
        r = sqrt((double)(x*x + y*y));
        r = (r==0.0) ? 0.0 : 1.0/r;
        x *= r;
        y *= r;
        angle = atan2((double)y, (double)x);
        if (angle < 0) angle += TWOPI;
        if (angle >= TWOPI) angle -= TWOPI;
        angle *= one_over_twopi;
        phase[j*xsize + i] = angle;
      }
    }
  }
  else if (in_format==2) {     /* quantized phase */
    scale = 1.0/256.0;
    for (j=0; j<xsize*ysize; j++) {
      phase[j] = quantized_phase[j]*scale;
    }
  }
  else {    /* 4-byte float phase */
    scale = one_over_twopi;
    for (j=0; j<xsize*ysize; j++) { 
      /* re-scale phase to interval (0,1) */
      phase[j] = float_phase[j]*scale;
    }
  }
}
