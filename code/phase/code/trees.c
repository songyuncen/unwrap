/*
 *  trees.c -- functions for managing trees, loops, edges, etc.
 *             in Flynn's min. discontinuity algorithm
 */
#include <stdio.h>
#include <math.h>
#include "trees.h"

/* Remove the loop and update all of the edges. */  
void RemoveLoop(int ibase, int jbase, int ilast, int jlast,
                int *value, unsigned char *bitflags, short *vjump,
                short *hjump, int xsize, int ysize)
{
  int   jtip,itip;
  int   value_shift;
  /* Remove edge from last node to base of loop */
  if (jbase>jlast)       --vjump[jbase*(xsize+1) + ibase];
  else if (jbase<jlast)  ++vjump[jlast*(xsize+1) + ilast];
  else if (ibase>ilast)  ++hjump[jbase*(xsize+1) + ibase];
  else if (ibase<ilast)  --hjump[jlast*(xsize+1) + ilast];
  /* Back up over loop: remove edges and change subtree values */
  jtip = jlast;
  itip = ilast;
  do {
    value_shift = -value[jtip*(xsize+1) + itip];
    ChangeOrphan(itip, jtip, value_shift, value, bitflags,
                 xsize, ysize);
    if (jtip > 0 
           && (bitflags[(jtip-1)*(xsize+1) + itip] & RIGHT)) {
      vjump[jtip*(xsize+1) + itip]--;
      bitflags[(jtip-1)*(xsize+1) + itip] &= ~RIGHT;
      jtip--;
    }
    else if ((jtip<ysize) &&
        (bitflags[(jtip+1)*(xsize+1) + itip] & LEFT)) {
      vjump[(jtip+1)*(xsize+1) + itip]++;
      bitflags[(jtip+1)*(xsize+1) + itip] &= ~LEFT;
      jtip++;
    }
    else if (itip > 0 
                && (bitflags[jtip*(xsize+1) + itip-1] & DOWN)) {
      hjump[jtip*(xsize+1) + itip]++;
      bitflags[jtip*(xsize+1) + itip-1] &= ~DOWN;
      itip--;
    }
    else if ((itip<xsize) 
                && (bitflags[jtip*(xsize+1) + itip+1] & UP)) {
      hjump[jtip*(xsize+1) + itip+1]--;
      bitflags[jtip*(xsize+1) + itip+1] &= ~UP;
      itip++;
    }
    else {
      printf("Error: broken loop at col %d, row %d\n",itip,jtip);
      exit(-1);
    }
  } while ((jtip != jbase) || (itip != ibase));
}

/*
 * Add 'value_change' (which is > 0) to the values of all nodes
 * of the tree rooted at the node (i,j).  Set a flag if the tree
 * includes (ilast,jlast). Make (i,j) child in the next iteration.
 */
void ChangeExten(int i, int j, int ilast, int jlast,
                 int *loop_found, int value_change, int *value,
                 unsigned char *bitflags, int xsize, int ysize)
{
  unsigned char kids;
  *loop_found |= ((j == jlast) && (i == ilast));
  kids = bitflags[j*(xsize+1) + i];
  if (kids & LEFT)  
    ChangeExten(i, j-1, ilast, jlast, loop_found, value_change,
                value, bitflags, xsize, ysize);
  if (kids & RIGHT)  
    ChangeExten(i, j+1, ilast, jlast, loop_found, value_change,
                value, bitflags, xsize, ysize);
  if (kids & UP)    
    ChangeExten(i-1, j, ilast, jlast, loop_found, value_change,
                value, bitflags, xsize, ysize);
  if (kids & DOWN) 
    ChangeExten(i+1, j, ilast, jlast, loop_found, value_change,
                value, bitflags, xsize, ysize);
  value[j*(xsize+1) + i] += value_change;
  bitflags[j*(xsize+1) + i] |= (NEXT_TIME | THIS_TIME);
}

/*
 * Add 'value_shift' (which is < 0) to the values of all nodes
 * of the tree rooted at the node (i,j).  If the change would make
 * 'value' < 0, make 'value' 0 instead  and split the node from
 * its parent.  Make (i,j) child in the next iteration.
 */
void ChangeOrphan(int i, int j, int value_shift, int *value,
                  unsigned char *bitflags, int xsize, int ysize)
{
  unsigned char kids;
  if (value[j*(xsize+1) + i] + (value_shift) < 0) {
    value_shift = -value[j*(xsize+1) + i];
    if (j>0)     bitflags[(j-1)*(xsize+1) + i] &=  ~RIGHT;
    if (j<ysize) bitflags[(j+1)*(xsize+1) + i] &=  ~LEFT;
    if (i>0)     bitflags[j*(xsize+1) + i-1]   &=  ~DOWN;
    if (i<xsize) bitflags[j*(xsize+1) + i+1]   &=  ~UP;
  }
  kids = bitflags[j*(xsize+1) + i];
  if (kids & LEFT) ChangeOrphan(i, j-1, value_shift, value,
                                bitflags, xsize, ysize);
  if (kids & RIGHT) ChangeOrphan(i, j+1, value_shift, value,
                                 bitflags, xsize, ysize);
  if (kids & UP) ChangeOrphan(i-1, j, value_shift, value,
                              bitflags, xsize, ysize);
  if (kids & DOWN) ChangeOrphan(i+1, j, value_shift, value,
                                bitflags, xsize, ysize);
  bitflags[j*(xsize+1) + i] |= (NEXT_TIME | THIS_TIME);
  value[j*(xsize+1) + i] += value_shift;
}
