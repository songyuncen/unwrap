/*
 * mainjump.c - Generate image of discontinuities 
 *
 * Source code files required:
 *      maskfat.c    mainjump.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "util.h"
#include "maskfat.h"

main (int argc, char *argv[])
{
  int            i, j, k, m, n, a, b;
  int            ii, jj, kk;
  FILE           *ifp, *ofp, *mfp=NULL;
  float          *surf;     /* array */ 
  unsigned char  *mask=NULL, *out;
  char           infile[200], outfile[200];
  char           maskfile[200];
  int            xsize, ysize;   /* dimensions of arrays */ 
  double         r, r1, r2, sum, xsum;

  printf("Generate discontinuity map of surface data\n");
  if (argc!=6) {
    fprintf(stderr,"Usage: %s surf image xsize ysize mask\n",
                  argv[0]);
    fprintf(stderr,"\nwhere surf = input surface (float) file\n");
    fprintf(stderr,"        image = output image (byte) file\n");
    fprintf(stderr,"        xsize, ysize = array dimensions\n");
    fprintf(stderr,"        mask = mask file (0=mask).\n");
    fprintf(stderr,"        (If no mask, enter 'none'.)\n");
    exit(BAD_USAGE);
  }
  sscanf(argv[1], "%s", infile);
  sscanf(argv[2], "%s", outfile);
  sscanf(argv[3], "%d", &xsize);
  sscanf(argv[4], "%d", &ysize);
  sscanf(argv[5], "%s", maskfile);

  printf("Input file = %s\n", infile);
  printf("Output file = %s\n", outfile);
  printf("Mask file = %s\n", maskfile);
  printf("Array size = %d x %d\n", xsize, ysize);
  if (Keyword(maskfile, "none")) printf("No mask file\n");

  /*  OPEN FILES, ALLOCATE MEMORY      */
  OpenFile(&ifp, infile, "r"); 
  OpenFile(&ofp, outfile, "w"); 
  if (!Keyword(maskfile, "none")) OpenFile(&mfp, maskfile, "r"); 
  AllocateFloat(&surf, xsize*ysize, "surface data");
  AllocateByte(&out, xsize*ysize, "binary image");
  if (mfp) AllocateByte(&mask, xsize*ysize, "mask data");
 
  /*  READ DATA  */
  printf("Reading data...\n");
  ReadFloat(ifp, surf, xsize*ysize, infile);
  if (mfp) ReadByte(mfp, mask, xsize*ysize, maskfile);

  /* FATTEN MASK BY 1 PIXEL (OPTIONAL) */
  if (mask) {
    for (k=0; k<xsize*ysize; k++) mask[k] = !mask[k];
    FattenMask(mask, 1, 1, xsize, ysize);
    for (k=0; k<xsize*ysize; k++) mask[k] = !mask[k];
  }

  /* Generate image of discontinuities (where derivs > PI) */
  for (j=0; j<ysize-1; j++) {
    for (i=0; i<xsize-1; i++) {
      k = j*xsize + i;
      out[k] = 255; 
      if (mask && (!mask[k] || !mask[k+1] || !mask[k+xsize])) {
        out[k] = 0;
      }
      else {
        r = surf[k] - surf[k+1];
        if (r > PI || r < -PI) out[k] = 0;
        r = surf[k] - surf[k+xsize];
        if (r > PI || r < -PI) out[k] = 0;
      }  
    }
  } 
  /* Compute measure (sum) of discontinuities */
  for (j=0, sum=xsum=0.0; j<ysize-1; j++) {
    for (i=0; i<xsize-1; i++) {
      k = j*xsize + i;
      if (mask && (!mask[k] || !mask[k+1] || !mask[k+xsize]))
        continue;
      if (!out[k]) ++sum;
      ++xsum;
    }
  }
  printf("Saving image to %s\n", outfile);
  printf("L0-measure (normalized to 0-1) of discontinuities %lf\n",
         sum/xsum);
  printf("(%lf%% of %d nonmasked pixels)\n", 
         100.0*sum/xsum, (int)xsum); 
  WriteByte(ofp, out, xsize*ysize, outfile);
  free(surf);
  fclose(ifp);
  free(out);
  fclose(ofp);
  if (mfp) free(mask);
  if (mfp) fclose(mfp);
}

