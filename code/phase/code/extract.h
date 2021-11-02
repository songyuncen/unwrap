#ifndef __EXTRACT
#define __EXTRACT
void GetPhase(int in_format, FILE *ifp, char *in_file, float *phase,
              int xsize, int ysize);
void  ExtractPhase(int in_format, void *in_data, float *phase,
                   int xsize, int ysize, int status_flag);
#endif
