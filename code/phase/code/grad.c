/*
 *   grad.c -- function for computing phase derivative 
 *             (wrapped phase difference)
 */
#include <stdio.h>
#include <math.h>
#include "grad.h"
/* return wrapped phase difference */
float Gradient(float p1, float p2)
{
  float  r;
  r = p1 - p2;
  if (r > 0.5) r -= 1.0;
  if (r < -0.5) r += 1.0;
  return r;
}
