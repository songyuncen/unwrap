/*
 * mainunwt.c -- unweighted least-squares phase unwrapping
 *               by means of direct transforms
 *
 * Source code files required:
 *          dct.c     extract.c        grad.c     laplace.c
 *     mainunwt.c     solncos.c        util.c
 */
#include <stdio.h>
#include <math.h>
#include "solncos.h"
#include "laplace.h"
#include "file.h"
#include "pi.h"
#include "util.h"
#include "extract.h"

main (int argc, char *argv[])
{
  int            i, j, k, n;
  FILE           *ifp, *ofp;
  float          *phase;     /* array */ 
  float          *soln;      /* array */ 
  double         *xcos, *ycos;
  char           infile[200], outfile[200];
  char           format[200];
  int            in_format;
  int            xsize, ysize;   /* dimensions of arrays */ 
  int            xsize_dct, ysize_dct;
  double         rmin, rmax, rscale;
  double         one_over_twopi = 1.0/TWOPI;
  char           use[] =       /* define usage statement */
   "Usage: prog-name -input file -format fkey -output file\n"
   "  -xsize x -ysize y\n"
   "where 'fkey' is a keyword designating the input file type\n"
   "(key = complex8, complex4, float or byte).  The dimensions\n"
   "x and y must be powers of two plus 1 (e.g., 257, 513, etc.)\n";

  printf("Unweighted Phase Unwrapping by DCT/FFT\n");

  /* GET COMMAND LINE PARAMETERS AND CHECK */
  CommandLineParm(argc,argv, "-input", StringParm, infile, 1, use);
  CommandLineParm(argc,argv, "-format", StringParm, format, 1, use);
  CommandLineParm(argc,argv, "-output", StringParm, outfile, 1,use);
  CommandLineParm(argc,argv, "-xsize", IntegerParm, &xsize, 1, use);
  CommandLineParm(argc,argv, "-ysize", IntegerParm, &ysize, 1, use);

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

  /* Verify dimensions are of form 2**n + 1 */
  for (xsize_dct = 1; xsize_dct + 1 < xsize; xsize_dct*=2)
    ;
  xsize_dct = xsize_dct + 1;
  for (ysize_dct = 1; ysize_dct + 1 < ysize; ysize_dct*=2)
    ;
  ysize_dct = ysize_dct + 1;
  if (xsize_dct != xsize || ysize_dct != ysize) {
    fprintf(stderr, "Error: dims %dx%d must be 2**n + 1\n",
            xsize, ysize);
    exit(BAD_PARAMETER);
  }

  /*  OPEN FILES, ALLOCATE MEMORY   */
  OpenFile(&ifp, infile, "r");
  OpenFile(&ofp, outfile, "w");
  AllocateFloat(&phase, xsize*ysize, "phase data");
  AllocateFloat(&soln, xsize*ysize, "scratch data");

  /*  READ AND PROCESS DATA  */
  printf("Reading phase data...\n");
  GetPhase(in_format, ifp, infile, phase, xsize, ysize);

  /*  UNWRAP  */
  printf("Computing Laplacian...\n");
  ComputeLaplacian(phase, soln, NULL, NULL, xsize, ysize, 1);
  free(phase);  
  AllocateDouble(&xcos, xsize, "cosine terms");
  AllocateDouble(&ycos, ysize, "cosine terms");
  printf("Performing direct transform...\n");
  DirectSolnByCosineTransform(soln, xsize, ysize, xcos, ycos);
  free(xcos);
  free(ycos);
  printf("\nFinished\n");
  PrintMinAndMax(xsize, ysize, soln, "solution");

  /*  SAVE RESULT  */
  for (k=0; k<xsize*ysize; k++)
    soln[k] *= TWOPI;
  printf("Saving unwrapped surface to file '%s'\n", outfile);
  WriteFloat(ofp, soln, xsize*ysize, outfile);
  free(soln);
}
