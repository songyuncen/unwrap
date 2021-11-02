/*
 * maindiff.c - difference between two surfaces & 
 *              RMS and absolute difference measures
 *
 * Source code files required:
 *   maindiff.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "util.h"

main (int argc, char *argv[])
{
  int            k, num;
  FILE           *ifp1, *ifp2, *ofp, *mfp;
  float          *surf1, *surf2;     /* array */ 
  float          *outsurf;      /* array */
  unsigned char  *mask;
  char           ifile1[200], ifile2[200], outfile[200];
  char           maskfile[200];
  int            xsize, ysize;   /* dimensions of arrays */ 
  double         sum, avg, sumsqr, avgsqr, rms, sumabs, avgabs;
  double         r, r1, r2;

  printf("Surface difference and RMS measure\n");
  if (argc!=7) {
    fprintf(stderr, "Usage: %s surf1 surf2 outfile xsize ysize"
                    " maskfile\n", argv[0]);
    fprintf(stderr, "\nwhere surf1 = first surface file\n");
    fprintf(stderr, "        surf2 = second surface file\n");
    fprintf(stderr, "        outfile = output surface file\n");
    fprintf(stderr, "        xsize, ysize = array dimensions\n");
    fprintf(stderr, "        maskfile = mask file (0=masked).\n");
    fprintf(stderr, "        (If no mask, enter 'none')\n");
    exit(BAD_USAGE);
  }
  sscanf(argv[1], "%s", ifile1);
  sscanf(argv[2], "%s", ifile2);
  sscanf(argv[3], "%s", outfile);
  sscanf(argv[4], "%d", &xsize);
  sscanf(argv[5], "%d", &ysize);
  sscanf(argv[6], "%s", maskfile);
  printf("Input file 1 = %s\n", ifile1);
  printf("Input file 2 = %s\n", ifile2);
  printf("Output file = %s\n", outfile);
  printf("Mask file = %s\n", maskfile);
  printf("Array size = %d x %d\n", xsize, ysize);
  if (Keyword(maskfile, "none")) printf("No mask file\n");

  /*  OPEN FILES, ALLOCATE MEMORY      */
  mfp = ifp1 = ifp2 = ofp = NULL;
  OpenFile(&ifp1, ifile1, "r"); 
  OpenFile(&ifp2, ifile2, "r"); 
  OpenFile(&ofp, outfile, "w"); 
  if (!Keyword(maskfile, "none")) OpenFile(&mfp, maskfile, "r"); 

  surf1 = surf2 = outsurf = NULL;
  mask = NULL;
  AllocateFloat(&surf1, xsize*ysize, "surface 1 data");
  AllocateFloat(&surf2, xsize*ysize, "surface 2 data");
  AllocateFloat(&outsurf, xsize*ysize, "output data");
  if (mfp) AllocateByte(&mask, xsize*ysize, "mask data");

  /*  READ DATA  */
  printf("Reading data...\n");
  ReadFloat(ifp1, surf1, xsize*ysize, ifile1);
  ReadFloat(ifp2, surf2, xsize*ysize, ifile2);
  if (mfp) ReadByte(mfp, mask, xsize*ysize, maskfile);

  /*  PROCESS  */
  for (k=0, num=0, sum=sumsqr=0.0; k<xsize*ysize; k++) {
    r1 = (surf1) ? surf1[k] : 0.0;
    r2 = (surf2) ? surf2[k] : 0.0;
    r = r1 - r2;
    if (mask && !mask[k]) {
      r = 0;
    }
    else {
      sum += r;
      sumsqr += r*r;
      ++num;
    }
  }  
  avg = (num > 0) ? sum/num : 0.0;
  avgsqr = (num > 0) ? sumsqr/num : 0.0;
  rms = sqrt(avgsqr - avg*avg);
  for (k=0, sumabs=0.0; k<xsize*ysize; k++) {
    r1 = (surf1) ? surf1[k] : 0.0;
    r2 = (surf2) ? surf2[k] : 0.0;
    r = (mask && !mask[k]) ? 0.0 : r1 - r2;
    outsurf[k] = r - avg;  /* subtract average */
    if (!mask || mask[k])  sumabs += (r > avg) ? r - avg : avg - r; 
  }
  avgabs = (num > 0) ? sumabs/num : 0.0;
  printf("Subtracting average value %lf from differences.\n", avg);
  printf("RMS value of differences = %lf.\n", rms);
  printf("Absolute difference = %lf.\n", avgabs);
  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, outsurf, xsize*ysize, outfile);
  free(surf1);
  free(surf2);
  free(outsurf);
  if (mask) free(mask);
  fclose(ifp1);
  fclose(ifp2);
  fclose(ofp);
  if (mfp) fclose(mfp);
}
