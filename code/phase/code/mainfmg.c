/*
 * mainfmg.c -- phase unwrapping by means of multigrid algorithm
 *
 * Source code files required:
 *     congruen.c    dxdygrad.c     extract.c         fmg.c
 *      getqual.c        grad.c        grid.c     gridmem.c
 *      gridops.c       histo.c     mainfmg.c     maskfat.c
 *     qualgrad.c    qualpseu.c     qualvar.c       relax.c
 *         util.c
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
#include "fmg.h"
#define BORDER           0x20
#define DEFAULT_NUM_ITER    2
#define DEFAULT_NUM_CYCLES  2

main (int argc, char *argv[])
{
  int            i, j, k, n;
  FILE           *ifp, *ofp, *mfp=0, *qfp=0;
  float          *phase;     /* array */ 
  float          *soln;      /* array */ 
  float          *qual_map;  /* array */
  float          *dx;        /* array */
  float          *dy;        /* array */
  unsigned char  *bitflags;
  char           buffer[200], tempstr[200];
  char           infile[200], outfile[200];
  char           bmaskfile[200], qualfile[200];
  char           format[200], modekey[200];
  int            in_format, debug_flag;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            avoid_code, thresh_flag, fatten, tsize;
  int            num_iter=DEFAULT_NUM_ITER;
  int            num_cycles=DEFAULT_NUM_CYCLES;
  double         rmin, rmax, rscale, one_over_twopi = 1.0/TWOPI;
  UnwrapMode     mode;
  char           use[] =   /* define usage statement */
    "Usage: prog-name -input file -format fkey -output file\n"
    "  -xsize x -ysize y [ -mode mkey -bmask file -corr file\n"
    "  -tsize size -debug yes/no -cycles numc -iter num\n"
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
    "byte-files are saved (quality map, etc.)  The number of\n"
    "multigrid cycles is given by 'numc'.  The number of Gauss-\n"
    "Seidel iterations is given by 'num'.  To apply an auto-\n"
    "matic threshold to the quality map in order to produce\n"
    "a quality mask, the 'thresh' parm should be yes.  To\n"
    "thicken the quality mask, the 'fat' parm should be the\n"
    "number of pixels by which to thicken.\n";

  printf("Phase Unwrapping by Weighted Multigrid Algorithm\n");
      
  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use);
  CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use);
  if (!CommandLineParm(argc, argv, "-mode", StringParm, modekey,
        0, use)) strcpy(modekey, "none");
  if (!CommandLineParm(argc, argv, "-bmask", StringParm, bmaskfile,
        0, use)) strcpy(bmaskfile, "none");
  if (!CommandLineParm(argc, argv, "-corr", StringParm, qualfile,
        0, use)) strcpy(qualfile, "none");
  if (!CommandLineParm(argc, argv, "-tsize", IntegerParm, &tsize,
        0, use)) tsize = 1;
  if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
        0, use)) debug_flag = 0;
  else debug_flag = Keyword(tempstr, "yes");
  CommandLineParm(argc, argv, "-iter", IntegerParm, &num_iter, 
                  0, use);
  if (num_iter < 0) num_iter = 1;
  CommandLineParm(argc, argv, "-cycles", IntegerParm, &num_cycles,
                  0, use);
  if (num_cycles < 0) num_cycles = 1;
  if (!CommandLineParm(argc, argv, "-thresh", StringParm, tempstr,
        0, use)) thresh_flag = 0;
  else thresh_flag = Keyword(tempstr, "yes");
  if (!CommandLineParm(argc, argv, "-fat", IntegerParm, &fatten,
        0, use)) fatten = 0;

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

  /*  OPEN FILES, ALLOCATE MEMORY   */
  OpenFile(&ifp, infile, "r");
  OpenFile(&ofp, outfile, "w");
  if (!Keyword(bmaskfile, "none"))
    OpenFile(&mfp, bmaskfile, "r");
  if (mode==corr_coeffs) 
    OpenFile(&qfp, qualfile, "r");
  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateFloat(&soln, xsize*ysize, "scratch data");
  AllocateByte(&bitflags, xsize*ysize, "bitflags array");
  AllocateFloat(&qual_map, xsize*ysize, "quality map");

  /*  READ AND PROCESS DATA  */
  printf("Reading phase data...\n");
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

  for (k=0; k<xsize*ysize; k++) {
    if (bitflags[k] & BORDER) qual_map[k] = 0.0;
  }
  free(bitflags);

  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.qual", outfile);
    SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                     0, 0, 0);
  }

  /*  UNWRAP  */
  printf("Unwrapping...\n");
  AllocateFloat(&dx, xsize*ysize, "dx derivs");
  AllocateFloat(&dy, ysize*ysize, "dy derivs");
  DxPhaseGradient(phase, dx, xsize, ysize);
  DyPhaseGradient(phase, dy, xsize, ysize);
  for (k=0; k<xsize*ysize; k++)  soln[k] = 0.0;
  MultigridUnwrap(soln, dx, dy, qual_map, NULL, xsize, ysize,
                  num_cycles, num_iter);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /* make result congruent to wrapped input phase */
  CongruentSoln(soln, phase, qual_map, xsize, ysize);
  /*  scale and save result  */
  for (k=0; k<xsize*ysize; k++) soln[k] *= TWOPI;
  printf("Saving unwrapped congruent surface to %s\n", outfile);
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  free(dx);
  free(dy);
  free(phase);
  free(qual_map);
  free(soln);
}
