/*
 * mainpcg.c -- phase unwrapping by means of PCG algorithm
 *
 * Source code files required:
 *     congruen.c         dct.c    dxdygrad.c     extract.c
 *      getqual.c        grad.c       histo.c     laplace.c
 *      mainpcg.c     maskfat.c         pcg.c    qualgrad.c
 *     qualpseu.c     qualvar.c     solncos.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "histo.h"
#include "pi.h"
#include "getqual.h"
#include "congruen.h"
#include "maskfat.h"
#include "util.h"
#include "extract.h"
#include "pcg.h"
#include "laplace.h"
#define BORDER   0x20
#define DEFAULT_NUM_ITER  20

main (int argc, char *argv[])
{
  int            i, j, k, n;
  FILE           *ifp, *ofp, *mfp=0, *qfp=0;
  float          *phase;     /* array */ 
  float          *soln;      /* array */ 
  float          *qual_map;  /* array */
  float          *rarray;  /* array */
  float          *parray;  /* array */
  float          *zarray;  /* array */
  unsigned char  *bitflags;
  double         *xcos, *ycos;
  char           buffer[200], tempstr[200];
  char           infile[200], outfile[200];
  char           bmaskfile[200], qualfile[200];
  char           format[200], modekey[200];
  int            in_format, debug_flag;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            xsize_actual, ysize_actual;
  int            xsize_dct, ysize_dct;
  int            avoid_code, thresh_flag, fatten;
  int            tsize, num_iter=DEFAULT_NUM_ITER;
  double         rmin, rmax, rscale, epsi_con;
  double         one_over_twopi = 1.0/TWOPI;
  UnwrapMode     mode;
  char           use[] =       /* define usage statement */
   "Usage: prog-name -input file -format fkey -output file\n"
   "  -xsize x -ysize y [ -mode mkey -bmask file -corr file\n"
   "  -tsize size -debug yes/no -iter num -converge epsi\n"
   "  -thresh yes/no -fat n ]\n"
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
   "Correlation unwrapping, or 'max_pseu' for Maximum\n"
   "Pseudocorrelation unwrapping.  All files are simple\n"
   "raster files, and the output file consists of floating\n"
   "point numbers that define the heights of the unwrapped\n"
   "surface.  If the 'debug' parm is 'yes', then intermediate\n"
   "byte-files are saved (quality map, etc.)  The maximum number\n"
   "of iterations is given by 'num' (default 20), and the con-\n"
   "vergence tolerance is given by 'epsi' (default = 0).  To\n"
   "apply an automatic threshold to the quality map to make\n"
   "a quality mask, the 'thresh' parm should be yes.  To\n"
   "thicken the quality mask, the 'fat' parm should be the\n"
   "number of pixels by which to thicken.\n";

  printf("Phase Unwrapping by PCG algorithm\n");

  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use);
  CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use);
  if (!CommandLineParm(argc, argv, "-mode", StringParm,
        modekey, 0, use)) strcpy(modekey, "none");
  if (!CommandLineParm(argc, argv, "-bmask", StringParm,    
        bmaskfile, 0, use)) strcpy(bmaskfile, "none");
  if (!CommandLineParm(argc, argv, "-corr", StringParm,
        qualfile, 0, use)) strcpy(qualfile, "none");
  if (!CommandLineParm(argc, argv, "-tsize", IntegerParm,
        &tsize, 0, use)) tsize = 1;
  if (!CommandLineParm(argc, argv, "-debug", StringParm,
        tempstr, 0, use)) debug_flag = 0;
  else debug_flag = Keyword(tempstr, "yes");
  CommandLineParm(argc, argv, "-iter", IntegerParm,
                  &num_iter, 0, use);
  if (num_iter < 0) num_iter = 1;
  if (!CommandLineParm(argc, argv, "-toler", DoubleParm,
        &epsi_con, 0, use)) epsi_con = 0.0;
  if (!CommandLineParm(argc, argv, "-thresh", StringParm,
        tempstr, 0, use)) thresh_flag = 0;
  else thresh_flag = Keyword(tempstr, "yes");
  if (!CommandLineParm(argc, argv, "-fat", IntegerParm,
        &fatten, 0, use)) fatten = 0;

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
  if (Keyword(qualfile, "none")) printf("No quality file.\n");
  else printf("Correlation image file = %s\n", qualfile);

  printf("Averaging template size = %d\n", tsize);
  if (tsize < 0 || tsize > 30) {
    fprintf(stderr, "Illegal size: must be between 0 and 30\n");
    exit(BAD_PARAMETER);
  }

  printf("Quality mode = %s\n", modekey);
  mode = SetQualityMode(modekey, qualfile, 0);
  if (mode < 0) exit(BAD_PARAMETER);  /* error msg already printed */

  /* Increase dimensions to power of two (plus one) */
  xsize_actual = xsize;
  ysize_actual = ysize;
  for (xsize_dct = 1; xsize_dct + 1 < xsize_actual; xsize_dct*=2)
    ;
  xsize_dct = xsize_dct + 1;
  for (ysize_dct = 1; ysize_dct + 1 < ysize_actual; ysize_dct*=2)
    ;
  ysize_dct = ysize_dct + 1;
  if (xsize_dct != xsize_actual || ysize_dct != ysize_actual) {
    printf("Dim's increased from %dx%d to %dx%d for FFT/DCT's\n",
           xsize_actual, ysize_actual, xsize_dct, ysize_dct);
  }

  /*  OPEN FILES, ALLOCATE MEMORY   */
  /* note: xsize_dct >= xsize;  ysize_dct >= ysize */
  xsize = xsize_dct;
  ysize = ysize_dct;
  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateFloat(&soln, xsize*ysize, "scratch data");
  AllocateFloat(&qual_map, xsize*ysize, "quality map");
  AllocateByte(&bitflags, xsize*ysize, "bitflags array");
  OpenFile(&ifp, infile, "r");
  OpenFile(&ofp, outfile, "w");
  if (!Keyword(bmaskfile, "none")) OpenFile(&mfp, bmaskfile, "r");
  if (mode==corr_coeffs) OpenFile(&qfp, qualfile, "r");

  /*  READ AND PROCESS DATA  */
  xsize = xsize_actual;
  ysize = ysize_actual;
  printf("Reading input data...\n");
  GetPhase(in_format, ifp, infile, phase, xsize, ysize);

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
  /* Set border (masked) pixels to 0-wts and free bitflags */
  for (k=0; k<xsize*ysize; k++) {
    if (bitflags[k] & BORDER) qual_map[k] = 0.0;
  }
  free(bitflags);   /* bitflags array no longer needed */

  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.qual", outfile);
    SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                     0, 0, 0);
  }

  /* embed arrays in possibly larger FFT/DCT arrays */
  for (j=ysize_dct-1; j>=0; j--) {
    for (i=xsize_dct-1; i>=0; i--) {
      if (i<xsize_actual && j<ysize_actual)
        phase[j*xsize_dct + i] = phase[j*xsize_actual + i];
      else phase[j*xsize_dct + i] = 0.0;
    }
  }
  if (qual_map) {
    for (j=ysize_dct-1; j>=0; j--) {
      for (i=xsize_dct-1; i>=0; i--) {
        if (i<xsize_actual && j<ysize_actual)
          qual_map[j*xsize_dct + i] = qual_map[j*xsize_actual + i];
        else
          qual_map[j*xsize_dct + i] = 0.0;
      }
    }
  }

  /* Set dimensions to DCT dimensions */
  xsize = xsize_dct;
  ysize = ysize_dct;
  /* Allocate more memory */
  AllocateFloat(&rarray, xsize*ysize, "r array data");
  AllocateFloat(&zarray, xsize*ysize, "z array data");
  AllocateFloat(&parray, xsize*ysize, "p array data");

  /*  UNWRAP  */
  printf("Unwrapping...\n");
  for (k=0; k<xsize*ysize; k++)  soln[k] = 0.0;
  ComputeLaplacian(phase, rarray, qual_map, NULL, xsize, ysize, 1);
  PCGUnwrap(rarray, zarray, parray, soln, qual_map, NULL,
            xsize, ysize, num_iter, epsi_con);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /* restore dimensions to input sizes */
  xsize = xsize_actual;
  ysize = ysize_actual;
  /* extract results from enlarged arrays */
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      soln[j*xsize + i] = soln[j*xsize_dct + i];
    }
  }
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      qual_map[j*xsize + i] = qual_map[j*xsize_dct + i];
    }
  }
  for (j=0; j<ysize; j++) {
    for (i=0; i<xsize; i++) {
      phase[j*xsize + i] = phase[j*xsize_dct + i];
    }
  }
  /* make result congruent to input phase */
/***  CongruentSoln(soln, phase, qual_map, xsize, ysize); ***/
  /* scale and save result */
  for (k=0; k<xsize*ysize; k++) soln[k] *= TWOPI;
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  printf("Wrote congruent surface to file %s\n", outfile);
  free(qual_map);
  free(phase);
  free(soln);
  free(rarray);
  free(parray);
  free(zarray);
}
