/*
 * mainunmg.c -- unweighted least-squares phase unwrapping
 *               by means of unweighted multigrid algorithm
 *
 * Source code files required:
 *     dxdygrad.c     extract.c        grad.c     gridmem.c
 *      gridops.c    mainunmg.c       relax.c       unfmg.c
 *       ungrid.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "file.h"
#include "pi.h"
#include "util.h"
#include "extract.h"
#include "unfmg.h"
#define DEFAULT_NUM_ITER    2
#define DEFAULT_NUM_CYCLES  2

main (int argc, char *argv[])
{
  int            i, j, k, n;
  FILE           *ifp, *ofp;
  float          *phase;     /* array */ 
  float          *soln;      /* array */ 
  float          *dx;        /* array */
  float          *dy;        /* array */
  char           infile[200], outfile[200];
  char           format[200];
  int            in_format;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            num_iter=DEFAULT_NUM_ITER;
  int            num_cycles=DEFAULT_NUM_CYCLES;
  char           use[] =   /* define usage statement */
    "Usage: prog-name -input file -format fkey -output file\n"
    "  -xsize x -ysize y [ -cycles numc -iter num]\n"
    "where 'fkey' is a keyword designating the input file type\n"
    "(key = complex8, complex4, float or byte) and 'x' and 'y'\n"
    "are the dimensions of the file.  All files are simple\n"
    "raster files, and the output file consists of floating\n"
    "point numbers that define the heights of the unwrapped\n"
    "surface.  The number of cycles is given by 'numc'.  The\n"
    "number of Gauss-Seidel iterations is given by 'num'.\n";

  printf("Unweighted Phase Unwrapping by Multigrid Algorithm\n");
      
  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use);
  CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use);
  CommandLineParm(argc, argv, "-iter", IntegerParm, &num_iter, 
                  0, use);
  if (num_iter < 0) num_iter = 1;
  CommandLineParm(argc, argv, "-cycles", IntegerParm, &num_cycles,
                  0, use);
  if (num_cycles < 0) num_cycles = 1;

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

  /*  OPEN FILES, ALLOCATE MEMORY   */
  OpenFile(&ifp, infile, "r");
  OpenFile(&ofp, outfile, "w");
  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateFloat(&soln, xsize*ysize, "solution array");

  /*  READ AND PROCESS DATA  */
  printf("Reading phase data...\n");
  GetPhase(in_format, ifp, infile, phase, xsize, ysize);

  /*  UNWRAP  */
  printf("Unwrapping...\n");
  AllocateFloat(&dx, xsize*ysize, "dx derivs");
  AllocateFloat(&dy, ysize*ysize, "dy derivs");
  DxPhaseGradient(phase, dx, xsize, ysize);
  DyPhaseGradient(phase, dy, xsize, ysize);
  for (k=0; k<xsize*ysize; k++)  soln[k] = 0.0;
  UnweightedMultigridUnwrap(soln, dx, dy, xsize, ysize,
                            num_cycles, num_iter);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /*  SAVE RESULT  */
  for (k=0; k<xsize*ysize; k++) 
    soln[k] *= TWOPI;
  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  free(dx);
  free(dy);
  free(soln);
}
