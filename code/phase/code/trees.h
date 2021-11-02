#ifndef __TREES
#define __TREES
#define LEFT    (0x08)
#define RIGHT   (0x04)
#define UP      (0x02)
#define DOWN    (0x01)
#define THIS_TIME  (0x40)
#define NEXT_TIME  (0x80)
void RemoveLoop(int ibase, int jbase, int ilast, int jlast,
                int *value, unsigned char *bitflags, short *vjump,
                short *hjump, int xsize, int ysize);
void ChangeExten(int i, int j, int ilast, int jlast,
                 int *loop_found, int value_change, int *value,
                 unsigned char *bitflags, int xsize, int ysize);
void ChangeOrphan(int i, int j, int value_shift, int *value,
                  unsigned char *bitflags, int xsize, int ysize);
#endif
