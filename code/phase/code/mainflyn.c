/*
 * mainflyn.c -- phase unwrapping by Flynn's min discontinuity alg.
 *
 * Source code files required:
 *     dxdygrad.c     extract.c       flynn.c     getqual.c
 *         grad.c       histo.c    mainflyn.c     maskfat.c
 *     qualgrad.c    qualpseu.c     qualvar.c       trees.c
 *         util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "histo.h"
#include "maskfat.h"
#include "pi.h"
#include "util.h"
#include "extract.h"
#include "getqual.h"
#include "flynn.h"
#define BORDER  (0x20)

main (int argc, char *argv[])
{
  int            i, j, k, tsize;
  FILE           *ifp, *ofp, *mfp=0, *qfp=0;
  float          *phase;     /* array */ 
  float          *qual_map;  /* array */
  unsigned char  *bitflags;  /* array */
  short          *hjump;     /* array */
  short          *vjump;     /* array */
  int            *value;     /* array */
  char           buffer[200], tempstr[200];
  char           infile[200], outfile[200];
  char           bmaskfile[200], qualfile[200];
  char           format[200], modekey[200];
  int            in_format, debug_flag;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            avoid_code, thresh_flag, fatten, guess_mode;
  UnwrapMode     mode;
  double         rmin, rmax, rscale, one_over_twopi = 1.0/TWOPI;
  /* define "use" statement */
  char           use[] =     /* define usage statement */
    "Usage: program-name -input file -format fkey -output file\n"
    "  -xsize x -ysize y [ -mode mkey -bmask file -corr file\n"
    "  -tsize size -debug yes/no -thresh yes/no -fat n\n"
    "  -guess yes/no]\n"
    "where 'fkey' is a keyword designating the input file type\n"
    "(key = complex8, complex4, float or byte), 'x' and 'y' are\n"
    "the dimensions of the file, bmask is an optional byte-file\n"
    "of masks for masking out undefined phase values, corr\n"
    "is an optional byte-file of cross-correlation values,\n"
    "tsize is the size of the square template for averaging the\n"
    "corr file or quality values (default = 1), and 'mkey' is\n"
    "a keyword designating the source of the quality values\n"
    "that guide the unwrapping path.  The value of mkey may be\n"
    "'min_grad' for Minimum Gradient unwrapping, 'min_var' for\n"
    "Minimum Variance unwrapping, 'max_corr' for Maximum\n"
    "Correlation unwrapping, 'max_pseu' for Maximum Pseudo-\n"
    "correlation unwrapping, or 'none'.  All files are simple\n"
    "raster files, and the output file consists of floating\n"
    "point numbers that define the heights of the unwrapped\n"
    "surface.  If the 'debug' parm is 'yes', then the\n"
    "intermediate byte-files are saved (quality map, etc.)\n"
    "To apply an automatic threshold to the quality map to make\n"
    "a quality mask, the 'thresh' parm should be yes.  To\n"
    "thicken the quality mask, the 'fat' parm should be the\n"
    "number of pixels by which to thicken.  If the 'guess'\n"
    "parameter is yes, then the input file must be a floating\n"
    "point array that represents an intermediate solution.\n"
    "This solution must be congruent to the wrapped phase.\n";

  printf("Phase Unwrapping by Flynn's Min. Discontinuity Method\n");

  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use);
  CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use);
  if (!CommandLineParm(argc, argv, "-mode", StringParm,
        modekey, 0, use))   strcpy(modekey, "none");
  if (!CommandLineParm(argc, argv, "-bmask", StringParm, 
        bmaskfile, 0, use)) strcpy(bmaskfile, "none");
  if (!CommandLineParm(argc, argv, "-corr", StringParm,
        qualfile, 0, use)) strcpy(qualfile, "none");
  if (!CommandLineParm(argc, argv, "-tsize", IntegerParm, 
        &tsize, 0, use)) tsize = 1;
  if (!CommandLineParm(argc, argv, "-debug", StringParm,
        tempstr, 0, use)) debug_flag = 0;
  else debug_flag = Keyword(tempstr, "yes");
  if (!CommandLineParm(argc, argv, "-thresh", StringParm,
        tempstr, 0, use)) thresh_flag = 0;
  else thresh_flag = Keyword(tempstr, "yes");
  if (!CommandLineParm(argc, argv, "-fat", IntegerParm,
        &fatten, 0, use)) fatten = 0;
  if (!CommandLineParm(argc, argv, "-guess", StringParm,
        tempstr, 0, use)) guess_mode = 0;
  else guess_mode = Keyword(tempstr, "yes");

  if (Keyword(format, "complex8"))  in_format = 0;
  else if (Keyword(format, "complex4"))  in_format = 1;
  else if (Keyword(format, "byte"))  in_format = 2;
  else if (Keyword(format, "float"))  in_format = 3;
  else {
    fprintf(stderr, "Unrecognized format: %s\n", format);
    exit(BAD_PARAMETER);
  }

  printf("Input file =  %s\n", infile);
  printf("Input file type = %s\n", format);
  printf("Output file =  %s\n", outfile);
  printf("File dimensions = %dx%d (cols x rows).\n", xsize, ysize);

  if (Keyword(bmaskfile, "none")) printf("No border mask file.\n");
  else printf("Border mask file = %s\n", bmaskfile);
  printf("Quality mode = %s\n", modekey);
  if (modekey==none) {
    printf("No quality map.\n");
    strcpy(qualfile, "none");
  }
  else {
    if (Keyword(qualfile, "none")) printf("No correlation file.\n");
    else printf("Correlation image file = %s\n", qualfile);
    printf("Averaging template size = %d\n", tsize);
    if (tsize < 0 || tsize > 30) {
      fprintf(stderr, "Illegal size: must be between 0 and 30\n");
      exit(BAD_PARAMETER);
    }
  }
  mode = SetQualityMode(modekey, qualfile, 1);
  if (mode < 0) exit(BAD_PARAMETER);  /* error msg already printed */




  /*  OPEN FILES, ALLOCATE MEMORY   */
  OpenFile(&ifp, infile, "r");
  OpenFile(&ofp, outfile, "w");
  if (!Keyword(bmaskfile, "none")) OpenFile(&mfp, bmaskfile, "r");
  if (mode==corr_coeffs) OpenFile(&qfp, qualfile, "r");
  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateByte(&bitflags, (xsize+1)*(ysize+1), "bitflags array");
  AllocateFloat(&qual_map, xsize*ysize, "quality map");
  AllocateShort(&hjump, (xsize+1)*(ysize+1), "hjump array");
  AllocateShort(&vjump, (xsize+1)*(ysize+1), "vjump array");
  AllocateInt(&value, (xsize+1)*(ysize+1), "value array");

  /*  READ AND PROCESS DATA  */
  if (guess_mode) {
    printf("Guess Mode: Inputting initial guess...\n");
    ReadFloat(ifp, phase, xsize*ysize, infile);
    for (k=0; k<xsize*ysize; k++) phase[k] *= one_over_twopi;
  }
  else {
    printf("Reading phase data...\n");
    GetPhase(in_format, ifp, infile, phase, xsize, ysize);
  }
  
  if (qfp) {
    printf("Reading quality data...\n");
    /* borrow the bitflags array temporarily */
    ReadByte(qfp, bitflags, xsize*ysize, qualfile);
    /* process data and store in quality map array */
    AverageByteToFloat(bitflags, qual_map, tsize, xsize, ysize);
  }

  /* border mask data */
  printf("Processing border mask data...\n");
  if (mfp) {
    ReadByte(mfp, bitflags, xsize*ysize, bmaskfile);
  }
  else {
    for (k=0; k<xsize*ysize; k++)
      bitflags[k] = 255;
  }
  for (k=0; k<xsize*ysize; k++) {
    bitflags[k] = (!bitflags[k]) ? BORDER : 0;
  }
  if (mfp) FattenMask(bitflags, BORDER, 1, xsize, ysize);
  GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                tsize, xsize, ysize);
  if (thresh_flag) {
    HistoAndThresh(qual_map, xsize, ysize, 0.0, 0, 0.0,
                   bitflags, BORDER);
    if (fatten > 0)
      FattenQual(qual_map, fatten, xsize, ysize);
  }
  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.qual", outfile);
    SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                     0, 0, 0);
  }

  /* embed bitflags array in larger array for use by Flynn's routines */
  for (j=ysize - 1; j>=0; j--) {
    for (i=xsize - 1; i>=0; i--) {
      bitflags[j*(xsize+1) + i] = bitflags[j*xsize + i];
    }
  }
  for (j=0; j<=ysize; j++) bitflags[j*(xsize + 1) + xsize] = 0;
  for (i=0; i<=xsize; i++) bitflags[ysize*(xsize + 1) + i] = 0;

  /*  UNWRAP  */
  printf("Unwrapping...\n");
  FlynnMinDisc(phase, value, bitflags, qual_map, vjump, hjump,
               xsize, ysize);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, phase, "solution");

  /*  SAVE RESULT  */
  for (k=0; k<xsize*ysize; k++)
    phase[k] *= TWOPI;
  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, phase, xsize*ysize, outfile);
  free(phase);
  free(bitflags);
  if (qual_map) free(qual_map);
  free(hjump);
  free(vjump);
  free(value);
}
