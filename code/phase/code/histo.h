#ifndef __HISTO
#define __HISTO
void HistoAndThresh(float *qual, int xdim, int ydim,
               double thresh, int percent_flag, double percentage, 
               unsigned char *ignore, int ignore_code);
void Histogram(float *qual, int xdim, int ydim, int *histo, 
               int num_bins, unsigned char *ignore, int ignore_code);
#endif
