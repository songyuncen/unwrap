#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "dct.h"

/*
 * Two-dimensional cosine tranform of array.
 * Array must be dimensioned from 0 to xsize * ysize - 1.
 * (xsize & ysize each must be one more than a power of two)
 * For forward transform, sign = 1.  For inverse transform,
 * sign = -1.
 */
void FastCosineTransform(float *array, int xsize, int ysize,
                         int sign)
{
  int          i, j;
  double       tmpr;
  float        *vector, *buffer;
  vector = (float*) malloc(ysize*sizeof(float));
  if (!vector) {
    fprintf(stderr, "Cannot allocate memory for cosine transf.\n");
    exit(MEMORY_ALLOCATION_ERROR);
  }
  buffer = (float *) malloc(BLOCKSIZE*ysize*sizeof(float));
  if (!buffer) {
    fprintf(stderr, "Cannot allocate buffer for cosine transf.\n");
    exit(MEMORY_ALLOCATION_ERROR);
  }
  RowDcts(array, xsize, ysize, sign);
  ColDcts(array, xsize, ysize, sign, vector, buffer);
  free(vector);
  free(buffer);
}

/* Perform 1-dim fast DCT's on rows of array */
void RowDcts(float *array, int xsize, int ysize, int sign)
{
  int  i, j;
  for (j=0; j<ysize; j++) {
    FastDct(&array[j*xsize], xsize, sign);
  }
}

/* Perform 1-dim fast DCT's on columns of array */
void ColDcts(float *array, int xsize, int ysize, int sign,
             float *vector, float *buffer)
{
  int  i, j;
  int  ii, istop, blocksize=BLOCKSIZE;
  if (blocksize<=0) {
    for (i=0; i<xsize; i++) {
      for (j=0; j<ysize; j++) {
        vector[j] = array[j*xsize + i];
      }
      FastDct(vector, ysize, sign);
      for (j=0; j<ysize; j++) {
        array[j*xsize + i] = vector[j];
      } 
    }
  }
  else {
    for (ii=0; ii<xsize; ii+=blocksize) {
      istop = (ii + blocksize <= xsize) ? ii + blocksize : xsize;
      for (j=0; j<ysize; j++) {
        for (i=ii; i<istop; i++) {
          buffer[(i - ii)*ysize + j] = array[j*xsize + i];
        }
      }
      for (i=ii; i<istop; i++) {        
        FastDct(&buffer[(i - ii)*ysize], ysize, sign);
      }
      for (j=0; j<ysize; j++) {
        for (i=ii; i<istop; i++) {
          array[j*xsize + i] = buffer[(i - ii)*ysize + j];
        }
      }
    }
  }
}

/* Fast Discrete Cosine Transform 
 *   The array must be dimensioned from 0 to size-1 
 *   (size must be one more than a power of two).
 *   If sign = -1, the inverse transform is computed.
 */
void FastDct(float *array, int size, int sign)
{
  int     j;
  double  rmult, sum, p1, p2;
  double  angle, zr, zi, zpi, zpr, r;
  --array;            --size;
  angle = PI/size;    r = sin(0.5*angle);
  zpr = -2.0*r*r;     zpi = sin(angle);
  zr = 1.0;           zi = 0.0;
  sum =      0.5*(array[1] - array[size + 1]);
  array[1] = 0.5*(array[1] + array[size + 1]);
  for (j=2; j<=size/2; j++) {
    r = zr;
    zr =  r*zpr - zi*zpi + zr;
    zi = zi*zpr +  r*zpi + zi;
    p1 = 0.5*(array[j] + array[size + 2 - j]);
    p2 =     (array[j] - array[size + 2 - j]);
    array[j] = p1 - zi*p2;
    array[size + 2 - j] = p1 + zi*p2;
    sum += zr*p2;
  }
  RealFft(array, size, 1);
  array[size + 1] = array[2];
  array[2] = sum;
  for (j=4; j<=size; j+=2) {
    sum += array[j];
    array[j] = sum;
  }
  if (sign < 0) rmult = 1.0/size;
  else rmult = 2.0;
  for (j=1; j<=size+1; j++) array[j] *= rmult;
}

/* FFT of real function (size must be a power of two) */
void RealFft(float *array, int size, int sign)
{
  int     i, i1, i2, i3, i4;
  double  s1, s2, qr, qi, rr, ri;
  double  zr, zi, zpr, zpi, r, angle;

  angle = PI/(size/2); 
  s1 = 0.5;
  if (sign > 0) {
    s2 = -0.5;        
    Fft(array, size/2, 1);
  }
  else {
    s2 = 0.5;         
    angle = -angle;
  }
  r = sin(0.5*angle);
  zpr = -2.0*r*r;     zpi = sin(angle);
  zr = 1.0 + zpr;     zi = zpi;
  for (i=2; i<= size/4; i++) {
    i1 = 2*i - 1;
    i2 = i1 + 1;
    i3 = size + 3 - i2;
    i4 = i3 + 1;
    qr = s1*(array[i1] + array[i3]);
    qi = s1*(array[i2] - array[i4]);
    rr = -s2*(array[i2] + array[i4]);
    ri =  s2*(array[i1] - array[i3]);
    array[i1]  =  qr + zr*rr - zi*ri;
    array[i2]  =  qi + zr*ri + zi*rr;
    array[i3]  =  qr - zr*rr + zi*ri;
    array[i4]  = -qi + zr*ri + zi*rr;
    r = zr;
    zr = r*zpr - zi*zpi + zr;
    zi = zi*zpr + r*zpi + zi;
  }
  qr = array[1];
  if (sign > 0) {
    array[1] = qr + array[2];
    array[2] = qr - array[2];
  }
  else {
    array[1] = s1*(qr + array[2]);
    array[2] = s1*(qr - array[2]);
    Fft(array, size/2, -1);
  }
}

/* FFT (size must be a power of two) */
void Fft(float *array, int size, int sign)
{
  int           n, nmax, m, j, istep, i;
  double        r, zr, zpr, zpi, zi, angle;
  double        tmp, tmpr, tmpi;
  j = 1;   n = 2*size;
  for (i=1; i<n; i+=2) {
    if (j > i) {
      /* perform swaps */
      tmp = array[j]; array[j] = array[i]; array[i] = tmp;
      tmp = array[j+1]; array[j+1] = array[i+1]; array[i+1] = tmp;
    }
    m = n/2;
    while (j > m && m > 1) {
      j -= m;    m /= 2;
    }
    j += m;
  }
  nmax = 2;
  while (n > nmax) {
    istep = 2*nmax;
    angle = TWOPI/nmax;
    if (sign < 0) angle = -angle;
    r = sin(0.5*angle);
    zpr = -2.0*r*r;    zpi = sin(angle);
    zr = 1.0;          zi = 0.0;
    for (m=1; m<nmax; m+=2) {
      for (i=m; i<=n; i+=istep) {
        j = i + nmax;
        tmpr = zr*array[j] - zi*array[j+1];
        tmpi = zr*array[j+1] + zi*array[j];
        array[j]   = array[i]   - tmpr;
        array[j+1] = array[i+1] - tmpi;
        array[i]   += tmpr;
        array[i+1] += tmpi;
      }
      r = zr;
      zr = r*zpr - zi*zpi + zr;
      zi = zi*zpr + r*zpi + zi;
    }
    nmax = istep;
  }
}
