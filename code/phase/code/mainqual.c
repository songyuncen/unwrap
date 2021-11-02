/*
 * mainqual.c -- phase unwrapping by quality-guided path following
 *
 * Source code files required:
 *     dxdygrad.c     extract.c     getqual.c        grad.c
 *         list.c    mainqual.c     maskfat.c     quality.c
 *     qualgrad.c    qualpseu.c     qualvar.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "quality.h"
#include "util.h"
#include "extract.h"
#include "getqual.h"

#define BORDER (0x20)
main (int argc, char *argv[])
{
  int            i, j, k;
  FILE           *ifp, *ofp, *mfp=0, *qfp=0;
  float          *phase;     /* array */ 
  float          *soln;      /* array */ 
  float          *qual_map;  /* array */
  unsigned char  *bitflags;  /* array */
  char           buffer[200], string[200];
  char           infile[200], outfile[200];
  char           bmaskfile[200], qualfile[200];
  char           tempstr[200], format[200], modekey[200];
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            tsize, in_format, debug_flag, avoid_code;
  double         rmin, rmax, rscale;
  UnwrapMode     mode;
  char           use[] =      /* define usage statement */
    "Usage: prog-name -input file -format fkey -output file\n"
    "  -xsize x -ysize y -mode mkey [ -bmask file -corr file\n"
    "  -tsize size -debug yes/no]\n"
    "where 'fkey' is a keyword designating the input file type\n"
    "(key = complex8, complex4, float or byte), 'x' and 'y' are\n"
    "the dimensions of the file, bmask is an optional byte-file\n"
    "of masks for masking out undefined phase values, corr\n"
    "is an optional byte-file of cross-correlation values,\n"
    "tsize is the size of the square template for averaging the\n"
    "corr file or quality values (default = 1), and 'mkey' is\n"
    "a keyword designating the source of the quality values\n"
    "that guide the unwrapping path.  The value of mkey must be\n"
    "'min_grad' for Minimum Gradient unwrapping, 'min_var' for\n"
    "Minimum Variance unwrapping, 'max_corr' for Maximum\n"
    "Correlation unwrapping, or 'max_pseu' for Maximum\n"
    "Pseudocorrelation unwrapping.  All files are simple\n"
    "raster files, and the output file consists of floating\n"
    "point numbers that define the heights of the unwrapped\n"
    "surface.  If the 'debug' parm is 'yes', then intermediate\n"
    "byte-files are saved (quality map, unwrapping paths, etc.)\n";

  printf("Phase Unwrapping by Quality Map Guided Path Following\n");

  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use);
  CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use);
  CommandLineParm(argc,argv, "-mode", StringParm, modekey, 1, use);
  if (!CommandLineParm(argc, argv, "-bmask", StringParm,
        bmaskfile, 0, use)) strcpy(bmaskfile, "none");
  if (!CommandLineParm(argc, argv, "-corr", StringParm,
        qualfile, 0, use)) strcpy(qualfile, "none");
  if (!CommandLineParm(argc, argv, "-tsize", IntegerParm,
        &tsize, 0, use)) tsize = 1;
  if (!CommandLineParm(argc, argv, "-debug", StringParm,
        tempstr, 0, use)) debug_flag = 0;
  else debug_flag = Keyword(tempstr, "yes");

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
  AllocateFloat(&soln, xsize*ysize, "unwrapped data");
  AllocateByte(&bitflags, xsize*ysize, "bitflag data");
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
  for (k=0; k<xsize*ysize; k++)
    bitflags[k] = (!bitflags[k]) ? BORDER : 0;
  FattenMask(bitflags, BORDER, 1, xsize, ysize);

  GetQualityMap(mode, qual_map, phase, bitflags, BORDER,
                tsize, xsize, ysize);
  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.qual", outfile);
    SaveFloatToImage(qual_map, "quality", filename, xsize, ysize,
                     0, 0, 0);
  }

  /* initialize soln array */
  for (k=0; k<xsize*ysize; k++)
    soln[k] = 0;

  /*  UNWRAP  */
  avoid_code = BORDER;
  if (mode==gradient && tsize==1) mode = dxdygrad;
  k = QualityGuidedPathFollower(phase, qual_map, bitflags, soln,
              xsize, ysize, avoid_code, mode, debug_flag, outfile);
  if (k > 1) printf("\n%d pieces\n", k);

  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /*  SAVE RESULT  */
  /* scale output */
  for (k=0; k<xsize*ysize; k++) 
    soln[k] *= TWOPI;

  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  free(soln);
  free(phase);
  free(bitflags);
  free(qual_map);
}
