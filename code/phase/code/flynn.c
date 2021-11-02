/*
 *  flynn.c -- functions for Flynn's min. discontinuity algorithm
 */
#include <stdio.h>
#include <math.h>
#include "flynn.h"
#include "trees.h"
/* define large value-increase */
#define BIG  25500

/* Main iteration (and initialization steps) for Flynn's minimum
 * discontinuity phase-unwrapping algorithm.  The quality values
 * (weights) are stored in the qual_map array.  This array and the
 * phase array are xsize x ysize pixels in size.  The bitflags array,
 * which stores the change and direction info, and the other arrays
 * (value, vjump and hjump) are working arrays for managing the data.
 * The dimensions of these arrays are (xsize + 1) x (ysize + 1).
 * If there are masked (border) pixels, they should be coded as
 * zero weights in qual_map.
 */
void FlynnMinDisc(float *phase, int *value, unsigned char *bitflags,
                  float *qual_map, short *vjump, short *hjump,
                  int xsize, int ysize)
{
  double          eps = 1.0e-06, phd;
  short           jmpold;
  int             j, i, jnext, inext;
  int             ilast, jlast;
  int             new_loops, new_edges;
  int             val, value_incr;
  int             loop_found;
  int             value_change;

  /* Compute phase jump counts.
   * If dx(i,j) and dy(i,j) represent the wrapped derivatives
   * phase(i+1,j) - phase(i,j) and phase(i,j+1) - phase(i,j),
   * respectively, then the jump count hjump(i,j) corresponds to
   * dy(i-1,j-1), and vjump(i,j) corresponds to dx(i-1,j-1).
   */  
  for (j=1; j <= ysize-1; j++) {
    for (i=1; i <= xsize; i++) {
      phd = phase[j*xsize + i-1] - phase[(j-1)*xsize + i-1] + eps;
      hjump[j*(xsize+1) + i] = (short) nint(phd);
    }
  }
  for (j=1; j <= ysize; j++) {
    for (i=1; i <= xsize-1; i++) {
      phd = phase[(j-1)*xsize + i] - phase[(j-1)*xsize + i-1] + eps;
      vjump[j*(xsize+1) + i] = (short) nint(phd);
    }
  }
  /* Make add nodes bitflags initially */
  for (j=0; j <= ysize; j++) {
    for (i=0; i <= xsize; i++) {
      bitflags[j*(xsize+1) + i] |= THIS_TIME;
    }
  }
  /* Main iteration */
  do {
    new_loops = 0;
    new_edges = 0;
    /* Add left-to-right edges to tree */
    for (j=0; j <= ysize-1; j++) {
      for (i=0; i <= xsize; i++) {
        jnext = j+1;
        inext = i;
        if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
              || (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)) {
          if (jnext <= 1 || jnext >= ysize 
               || inext <= 0 || inext >= xsize) {
            value_incr = -isign(1, -vjump[jnext*(xsize+1) + inext]);
          }
          else {
            /* quality map values must be between 0 and 1 */
            val = (int)(1.0 + BIG*qual_map[jnext*xsize + inext]);
            value_incr  
               = -isign(val, -vjump[jnext*(xsize+1) + inext]);
          }
          value_change = value[j*(xsize+1) + i] + (int)value_incr
                                   - value[jnext*(xsize+1) + inext];
          if (value_change > 0) {
            new_edges++;
            /* revise values in subtree of [jnext][inext] */
            /* and check for loop */
            ilast = i;  jlast = j;  loop_found = 0;
            ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                        value_change, value, bitflags, xsize, ysize);
            if (loop_found) {
              RemoveLoop(inext, jnext, i, j, value, bitflags,
                         vjump, hjump, xsize, ysize);
              ChangeOrphan(inext, jnext, -value_change, value,
                           bitflags, xsize, ysize);
              new_loops++;
            }
            else {
              /* add edge and remove edges to [jnext][inext] */
              bitflags[j*(xsize+1) + i] |= RIGHT;
              if (jnext<ysize) 
                bitflags[(jnext+1)*(xsize+1) + inext] &= ~LEFT;
              if (inext<xsize) 
                bitflags[jnext*(xsize+1) + inext+1] &= ~UP;
              if (inext > 0)
                bitflags[jnext*(xsize+1) + inext-1] &= ~DOWN;
            }
          }
        }
      }
    }
    /* Add top-to-bottom edges to tree */
    for (j=0; j <= ysize; j++) {
      for (i=0; i <= xsize-1; i++) {
        jnext = j;
        inext = i+1;
        if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
             || (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)){
          if (jnext <= 0 || jnext >= ysize 
               || inext <= 1 || inext >= xsize) {
            value_incr = -isign(1, hjump[jnext*(xsize+1) + inext]);
          }
          else {
            /* quality map values must be between 0 and 1 */
            val = (int)(1.0 + BIG*qual_map[jnext*xsize + inext]);
            value_incr 
                = -isign(val, hjump[jnext*(xsize+1) + inext]);
          }
          value_change = value[j*(xsize+1) + i] + (int)value_incr
                                 - value[jnext*(xsize+1) + inext];
          if (value_change > 0) {
            new_edges++;
            ilast = i; jlast = j; loop_found = 0;
            ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                        value_change, value, bitflags, xsize, ysize);
            if (loop_found){
              RemoveLoop(inext, jnext, i, j, value, bitflags,
                         vjump, hjump, xsize, ysize);
              ChangeOrphan(inext, jnext, -value_change, value,
                           bitflags, xsize, ysize);
              new_loops++;
            }
            else {
              bitflags[j*(xsize+1) + i] |= DOWN;
              if (jnext<ysize) 
                bitflags[(jnext+1)*(xsize+1) + inext] &= ~LEFT;
              if (jnext > 0) 
                bitflags[(jnext-1)*(xsize+1) + inext] &= ~RIGHT;
              if (inext<xsize) 
                bitflags[jnext*(xsize+1) + inext+1] &= ~UP;
            }
          }
        }
      }
    }
    /* Add right-to-left edges to tree */
    for (j=ysize; j >=1 ; j--) {
      for (i=xsize; i >=0; i--) {
        jnext = j-1;
        inext = i;
        if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
               || (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)){
          if (j <= 1 || j >= ysize || i <= 0 || i >= xsize) {
            value_incr = -isign(1, vjump[j*(xsize+1) + i]);
          }
          else {
            /* quality map values must be between 0 and 1 */
            val = (int)(1.0 + BIG*qual_map[j*xsize + i]);
            value_incr = -isign(val, vjump[j*(xsize+1) + i]);
          }
          value_change = value[j*(xsize+1) + i] + (int)value_incr
                                 - value[jnext*(xsize+1) + inext];
          if (value_change > 0) {
            new_edges++;
            ilast = i; jlast = j; loop_found = 0;
            ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                        value_change, value, bitflags, xsize, ysize);
            if (loop_found){
              RemoveLoop(inext, jnext, i, j, value, bitflags, 
                         vjump, hjump, xsize, ysize);
              ChangeOrphan(inext, jnext, -value_change, value,
                           bitflags, xsize, ysize);
              new_loops++;
            }
            else {
              bitflags[j*(xsize+1) + i] |= LEFT;
              if (jnext > 0)
                bitflags[(jnext-1)*(xsize+1) + inext] &= ~RIGHT;
              if (inext < xsize)
                bitflags[jnext*(xsize+1) + inext+1] &= ~UP;
              if (inext > 0) 
                bitflags[jnext*(xsize+1) + inext-1] &= ~DOWN;
            }
          }
        }
      }
    }
    /* Add bottom-to-top edges to tree */
    for (j=ysize; j >=0; j--) {
      for (i=xsize; i >= 1; i--) {
        jnext = j;
        inext = i-1;
        if ((bitflags[j*(xsize+1) + i] & THIS_TIME)
               || (bitflags[jnext*(xsize+1) + inext] & THIS_TIME)){
          if (j <= 0 || j >= ysize || i <= 1 || i >= xsize) {
            value_incr = -isign(1, -hjump[j*(xsize+1) + i]);
          }
          else {
            /* quality map values must be between 0 and 1 */
            val = (int)(1.0 + BIG*qual_map[j*xsize + i]);
            value_incr = -isign(val, -hjump[j*(xsize+1) + i]);
          }
          value_change = value[j*(xsize+1) + i] + (int)value_incr
                                 - value[jnext*(xsize+1) + inext];
          if (value_change > 0) {
            new_edges++;
            ilast = i; jlast = j; loop_found = 0;
            ChangeExten(inext, jnext, ilast, jlast, &loop_found,
                        value_change, value, bitflags, xsize, ysize);
            if (loop_found){
              RemoveLoop(inext, jnext, i, j, value, bitflags,
                         vjump, hjump, xsize, ysize);
              ChangeOrphan(inext, jnext, -value_change, value,
                           bitflags, xsize, ysize);
              new_loops++;
            }
            else {
              bitflags[j*(xsize+1) + i] |= UP;
              if (jnext < ysize) 
                bitflags[(jnext+1)*(xsize+1) + inext] &= ~LEFT;
              if (jnext > 0) 
                bitflags[(jnext-1)*(xsize+1) + inext] &= ~RIGHT;
              if (inext > 0) 
                bitflags[jnext*(xsize+1) + inext-1] &= ~DOWN;
            }
          }
        }
      }
    }
    printf("New edges: %d  New loops: %d\n",new_edges,new_loops);
    /* Update testability matrix */
    for (j=0; j <= ysize; j++) {
      for (i=0; i <= xsize; i++) {
        if (bitflags[j*(xsize+1) + i] & NEXT_TIME) {
          bitflags[j*(xsize+1) + i] |= THIS_TIME;
        }
        else {
          bitflags[j*(xsize+1) + i] &= ~NEXT_TIME;
          bitflags[j*(xsize+1) + i] &= ~THIS_TIME;
        }
      }
    }
  } while (new_edges > 0);

  /* Compute new unwrapped phase and exit */
  printf("Computing revised unwrapped phase...\n");
  for (i = 0; i < xsize - 1; i++) {
    phd = phase[0*xsize + (i + 1)] - phase[0*xsize + i] + eps;
    jmpold = (short) nint(phd);
    phase[0*xsize + (i + 1)] 
        += (float)(vjump[1*(xsize + 1) + (i + 1)] - jmpold);
  }
  for (j = 0; j < ysize - 1; j++) {
    for (i = 0; i < xsize; i++) {
      phd = phase[(j + 1)*xsize + i] 
                    - phase[j*xsize + i] + eps;
      jmpold = (short) nint(phd);
      phase[(j + 1)*xsize + i] 
        += (float)(hjump[(j + 1)*(xsize + 1) + (i + 1)] - jmpold);
    }
  }
}
