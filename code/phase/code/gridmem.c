/*
 * gridmem.c - dynamic memory management functions
 *             for multigrid phase unwrapping 
 */
#include <stdio.h>
#include <math.h>
#include "gridmem.h"

static int sizes[NUM_TYPES][MAX_ARRAYS];
static int num[NUM_TYPES];
static float *arrays[NUM_TYPES][MAX_ARRAYS];

/* allocate array: if already allocated, return pointer to it */
float *Allocate(int w, int h, ArrayType type)
{
  int size = w*h, i, j, k;
  k = (int) type;
  for (i=0; i<num[k]; i++) {
    if (sizes[k][i]==size) 
      return arrays[k][i];
  }
  sizes[k][num[k]] = size;
  arrays[k][num[k]] = (float *) malloc(size*sizeof(float));
  if (!arrays[k][num[k]]) {
    printf("Error: cannot allocate memory (%d bytes)\n", size*4);
    return 0;  /* error */
  }
  return arrays[k][num[k]++];
}

/* free all allocated arrays */
void FreeAll()
{ 
  int i, j;
  for (j=0; j<NUM_TYPES; j++) {
    for (i=0; i<num[j]; i++) {
      free(arrays[j][i]);
    }
  }
}
