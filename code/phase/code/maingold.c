/*
 * maingold.c - phase unwrapping by means of residues & branch cuts
 *
 * Source code files required:
 *        brcut.c      dipole.c     extract.c        gold.c
 *         grad.c        list.c    maingold.c     maskfat.c
 *         path.c    residues.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "residues.h"
#include "util.h"
#include "extract.h"
#include "list.h"
#include "gold.h"
#include "dipole.h"

main (int argc, char *argv[])
{
  int            i, j, k, m, n, a, b;
  FILE           *ifp, *ofp, *mfp;
  float          *phase;     /* array */ 
  float          *soln;      /* array */
  unsigned char  *bitflags;
  char           buffer[200], string[200];
  char           infile[200], outfile[200];
  char           maskfile[200], tempstr[200], format[200];
  int            in_format, debug_flag, dipole_flag;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            NumRes, MaxCutLen;
  char           use[] =    /* define usage statement */
    "Usage: prog-name -input file -format fkey -output file\n"
    "  -xsize x -ysize y [ -mask file -cutlen c -debug yes/no\n"
    "  -dipole yes/no ]\n"
    "where 'fkey' is a keyword designating the input file type\n"
    "(key = complex8, complex4, float or byte), 'x' and 'y' are\n"
    "the dimensions of the file, mask is an optional byte-file\n"
    "of masks for masking out undefined phase values and cutlen\n"
    "is the max branchcut length allowed.  All files are simple\n"
    "raster files, and the output file consists of floating\n"
    "point numbers that define the heights of the unwrapped\n"
    "surface.  If the 'debug' parm is 'yes', then intermediate\n"
    "byte-files are saved (residues, branch cuts, etc.)  If\n"
    "the 'dipole' parm is 'yes', then dipole-residues are\n"
    "eliminated before unwrapping.\n";

  printf("Phase Unwrapping By Means of Goldstein's Algorithm\n");

  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc, argv, "-xsize", IntegerParm, &xsize, 1,use);
  CommandLineParm(argc, argv, "-ysize", IntegerParm, &ysize, 1,use);
  if (!CommandLineParm(argc, argv, "-mask", StringParm, maskfile,
        0, use)) strcpy(maskfile, "none");
  if (!CommandLineParm(argc, argv, "-cutlen", IntegerParm,
     &MaxCutLen, 0, use)) MaxCutLen = 0; /* default defined below */
  if (!CommandLineParm(argc, argv, "-debug", StringParm, tempstr,
        0, use)) debug_flag = 0;
  else debug_flag = Keyword(tempstr, "yes");
  if (!CommandLineParm(argc, argv, "-dipole", StringParm, tempstr,
        0, use)) dipole_flag=0;
  else dipole_flag = Keyword(tempstr, "yes");

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

  if (Keyword(maskfile, "none")) printf("No mask file.\n");
  else printf("Mask file = %s\n", maskfile);

  if (dipole_flag) printf("Dipole-residues will be eliminated.\n");
  if (debug_flag) 
    printf("Intermediate files will be saved (i.e. debug on).\n");

  /*  OPEN FILES, ALLOCATE MEMORY      */
  OpenFile(&ifp, infile, "r");
  OpenFile(&ofp, outfile, "w");
  mfp = NULL;
  if (!Keyword(maskfile, "none"))
    OpenFile(&mfp, maskfile, "r");

  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateFloat(&soln, xsize*ysize, "unwrapped data");
  AllocateByte(&bitflags, xsize*ysize, "bitflag data");

  /*  READ AND PROCESS DATA  */
  printf("Reading phase data...\n");
  GetPhase(in_format, ifp, infile, phase, xsize, ysize);

  /* mask data */
  printf("Processing mask data...\n");
  if (mfp) {
    ReadByte(mfp, bitflags, xsize*ysize, maskfile);
  }
  else {
    for (k=0; k<xsize*ysize; k++)
      bitflags[k] = 255;
  }
  for (k=0; k<xsize*ysize; k++)
    bitflags[k] = (!bitflags[k]) ? BORDER : 0;

  FattenMask(bitflags, BORDER, 1, xsize, ysize);

  /*  LOCATE AND PROCESS RESIDUES  */
  /* compute residues and store in bitflags array */
  NumRes = Residues(phase, bitflags, POS_RES, NEG_RES,
                    BORDER, xsize, ysize);
  printf("%d Residues\n", NumRes);

  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.res", outfile);
    SaveByteToImage(bitflags, "residues", filename, xsize, ysize,
                    1, 1, 0);
  }

  /*  GENERATE BRANCH CUTS  */
  if (dipole_flag) {  /* elimate dipole-residues first */
    Dipole(bitflags, xsize, ysize, BRANCH_CUT);
    i = Residues(phase, bitflags, 0, 0, BORDER | BRANCH_CUT,
                 xsize, ysize);
    printf("%d Residues are left\n", i);
  }

  if (MaxCutLen==0) MaxCutLen = (xsize + ysize)/2;
  GoldsteinBranchCuts(bitflags, MaxCutLen, NumRes, xsize, ysize,
                      BRANCH_CUT);
  if (debug_flag) {
    char filename[300];
    sprintf(filename, "%s.brc", outfile);
    SaveByteToImage(bitflags, "branch cuts", filename, 
                    xsize, ysize, 1, 1, BRANCH_CUT | BORDER);
  }

  /*  UNWRAP AROUND CUTS */
  for (k=0, i=0; i<xsize*ysize; i++) 
    if (bitflags[i] & BRANCH_CUT) k++;
  printf("%d BRANCH CUT PIXELS\n", k);
  printf("Unwrapping around branch cuts\n");
  k = UnwrapAroundCuts(phase, bitflags, soln, xsize, ysize, AVOID,
                       0, NULL);
  if (k > 1) printf("%d disconnected pieces.\n", k);
  else printf("%d piece\n", k);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /*  SAVE RESULT  */
  for (k=0; k<xsize*ysize; k++)
    soln[k] *= TWOPI; /* scale output */
  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  free(soln);
  free(phase);
  free(bitflags);
}
