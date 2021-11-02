/*
 *  histo.c -- functions for histogramming and thresholding
 */
#include <stdio.h>
#include <math.h>
#include "histo.h"
#define NUM_BINS           10000
#define NUM_BINS_TO_PRINT  10

/*
 * Compute histogram, apply histogram equalization and
 * apply threshold.  If thresh=0, the threshold is determined
 * automatically.  If percent_flag is 1, the threshold is
 * selected to threshold the percentage of the pixels
 * defined by "percentage".  The bitflags array is "ignore",
 * and ignore_code defines those pixels to ignore in the
 * array "qual" of values to be histogrammed and thresholded.
 * The results overwrite the values in "qual".
 */
void HistoAndThresh(float *qual, int xsize, int ysize,
               double thresh, int percent_flag, double percentage, 
               unsigned char *ignore, int ignore_code)
{
  int     i, j, k;
  int     num_zero=0, sum, num;
  int     histo[NUM_BINS + 1];
  double  first_val, last_val, rmult;
  printf("Computing histogram, equalizing & applying threshold\n");
  /* compute and print histogram */
  Histogram(qual, xsize, ysize, histo, NUM_BINS,
            ignore, ignore_code);
  /* stretch histogram so at least 5% in first & last bins */
  for (num=0.0, i=0; i<NUM_BINS; i++) { /* find number of pixels */
    num += histo[i];
  }
  for (sum=0.0, i=0; i<NUM_BINS; i++) {
    sum += histo[i];
    if (sum >= 0.05*num) break;
  }
  first_val = (1.0*i)/NUM_BINS;
  for (sum=0.0, i=0; i<NUM_BINS; i++) {
    sum += histo[i];
    if (sum >= 0.95*num) break;
  }
  last_val = (1.0*i)/NUM_BINS;
  if (last_val <= first_val) last_val = first_val + 1.0/NUM_BINS;
  rmult = 1.0/(last_val - first_val);
  /* stretch weights */
  for (j=0; j<xsize*ysize; j++) {
    if (ignore[j] & ignore_code) continue;   /* ignore these */
    if (qual[j] < first_val)
      qual[j] = 0.0;
    else if (qual[j] > last_val)
      qual[j] = 1.0;
    else
      qual[j] = (qual[j] - first_val)*rmult;
  }
  /* compute histogram again */
  printf("Stretching histogram...\n");
  Histogram(qual, xsize, ysize, histo,
            NUM_BINS, ignore, ignore_code);
  for (num=0.0, i=0; i<NUM_BINS; i++) {  /* find number of pixels */
    num += histo[i];
  }
  if (thresh==0.0) {
    /* Find optimal threshold by finding the trough of the */
    /* histogram.  Failing that, choose percentage = 50.   */
    int  downhill=0, uphill=0, num_coarse_bins = 10;
    double  coarse_histo[NUM_BINS];
    printf("Looking for trough of histogram to set threshold...\n");
    while (num_coarse_bins <= 100) {
      for (k=0; k<num_coarse_bins; k++)
        coarse_histo[k] = 0.0;
      for (i=0, k=0, sum=0.0; i<NUM_BINS; i++) {
        if (i==NUM_BINS - 1) sum += histo[i];
        if ((int)((num_coarse_bins*1.0*i)/NUM_BINS) > k 
                                        || i==NUM_BINS - 1) {
          /* coarse_histo holds the percentages */
          coarse_histo[k] = (100.0*sum)/(num);
          sum = 0.0;
          ++k;
        }
        sum += histo[i];
      }
      /* find first minima */
      for (i=1, downhill=0; i<num_coarse_bins; i++) {
        if (downhill) {
          if (coarse_histo[i] > coarse_histo[i-1] + 0.5) {
            printf("Uphill at %d\n", i);
            break;
          }
        }
        else {
          if (coarse_histo[i] < coarse_histo[i-1] - 1.0) {
            downhill = 1;
            printf("Downhill at %d\n", i);
          }
        }
      }
      if (i < num_coarse_bins) {
        /* thresh = (i - 1.0)/num_coarse_bins; */
        for (j=0, percentage=0.0; j<i; j++) {
          percentage += coarse_histo[j];
        }
        for (sum=0.0, i=0; i<NUM_BINS; i++) {
          sum += histo[i];
          if (sum >= percentage*num*0.01) break;
        }
        if (i<NUM_BINS - 1) thresh = (i + 1.0)/NUM_BINS;
        else thresh = (NUM_BINS - 1.0)/NUM_BINS;
        printf("Found trough.  Setting threshold to %g (%g%%).\n",
               thresh, percentage);
        break;
      }
      else {
        /* repeat with more bins in histogram */
        num_coarse_bins += 10;
      }
    }
    if (thresh==0.0) {
      percentage = 50.0;
      for (sum=0.0, i=0; i<NUM_BINS; i++) {
        sum += histo[i];
        if (sum >= percentage*num*0.01) break;
      }
      if (i<NUM_BINS - 1) thresh = (i + 1.0)/NUM_BINS;
      else thresh = (NUM_BINS - 1.0)/NUM_BINS;
      printf("Could not find min of histogram.\n");
      printf("Setting thresh at 50%%.\n");
    }
  }

  /* compute threshold from desired percentage */
  if (percent_flag) {
    for (sum=0.0, i=0; i<NUM_BINS; i++) {
      sum += histo[i];
      if (sum >= percentage*num*0.01) break;
    }
    if (i<NUM_BINS - 1) thresh = (i + 1.0)/NUM_BINS;
    else thresh = (NUM_BINS - 1.0)/NUM_BINS;
  }
  /* apply threshold */
  printf("Threshold = %lf\n", thresh);
  for (j=0, num_zero=0, num=0; j<xsize*ysize; j++) {
    if (ignore[j] & ignore_code) continue;   /* ignore these */
    ++num;
    if (qual[j]<thresh) ++num_zero;
    qual[j] = (qual[j]<thresh) ? 0.0 : 1.0;
  }
  printf("%lg percent of the weights are now zero-weights\n",
         100.0*((double)num_zero)/((double)(num)));
}

/* Generate histogram (in array histo) of the values in array */
/* "qual" using num_bins bins.  The bigflags array is array   */
/* "ignore", and the pixels marked with the ignore_code must  */
/* be ignored (these are mask pixels, for example).           */
void Histogram(float *qual, int xsize, int ysize, int *histo,
               int num_bins, unsigned char *ignore, int ignore_code)
{
  double min_qual, max_qual;
  int   sum, num, i, j, k, off=0;
  /* zero histogram bins */
  min_qual = 1.0e+10;
  max_qual = -1.0e+10;
  for (i=0; i<=num_bins; i++) histo[i] = 0;
  /* find min & max weight, and histogram of weights */
  for (j=0, num=0; j<xsize*ysize; j++) {
    if (ignore[j] & ignore_code) continue;   /* ignore these */
    ++num;
    if (qual[j] < min_qual) min_qual = qual[j];
    if (qual[j] > max_qual) max_qual = qual[j];
    if ((qual[j] < -1.0e-5 || qual[j] > 1.0 + 1.0e-5) && !off) {
      off = 1;
      printf("WARNING: Quality values must be between 0 and 1!\n");
      printf("(Found a magnitude of %f.)\n", qual[j]);
      printf("Results of this run may be invalid!\n\n");
    }
    if (qual[j] < 0.0) qual[j] = 0.0;
    if (qual[j] > 1.0) qual[j] = 1.0;
    ++histo[(int)(qual[j]*num_bins)];
  }
  histo[NUM_BINS - 1] += histo[NUM_BINS];  /* ignore last bin */
  /* print histogram */
  printf("Min & max qual = %lf, %lf\nHistogram (percentages) = \n", 
         min_qual, max_qual);
  /* quantize histogram to fewer bins and print */
  for (i=0, k=0, sum=0.0; i<NUM_BINS; i++) {
    if (i==NUM_BINS - 1) sum += histo[i];
    if ((int)((NUM_BINS_TO_PRINT*1.0*i)/NUM_BINS) > k 
                                         || i==NUM_BINS - 1) {
      printf("Bins %.3f-%.3f:  %d\n", (k + 0.0)/NUM_BINS_TO_PRINT,
         (k + 1.0)/NUM_BINS_TO_PRINT, (int)((100.0*sum)/(num)+0.5));
      sum = 0.0;
      ++k;
    }
    sum += histo[i];
  }
}
