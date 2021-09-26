#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "OTFParameters.h"
#include "Stats.h"

// 30-dec-2020:    I slapped these together from NEMO's median.c and moment.c
//                 which has a a few options if the current one doesn't work out.
//



// helper function for qsort()
int compar_float(const void *va, const void *vb)
{
  float *a = (float *) va;
  float *b = (float *) vb;
  return *a < *b ? -1 : *a > *b ? 1 : 0;
}

/*
 *   rms_stats:   a dynamic rms_cutoff evaluator
 *                rms_cutoff > 0  : classic absolute cutoff
 *                           < 0  : robust  mean + |rms_cutoff|*std
 */

void rms_stats(int n, float *RMS, int *Pixel,
	       int npix, int *use_pixel, float *RMS_cut,
	       float rms_cutoff)
{
  int i, j, k, l;
  float *dat;
  float mean, median, std;
  float s0, s1, s2;
  float m1, m2, m3, iqr, dlo, dhi;
  float frob = 1.5;   // hardcoded for a robust mean/std
  
  if (n<0 || npix<0) return;
  printf("rms_stats: n=%d rms_cutoff=%g\n",n, rms_cutoff);
  dat = (float *) malloc(n * sizeof(float));  

  for (i=0; i<MAXPIXEL; i++) {
    // ignore pixels we don't use
    if (use_pixel[i] == 0) {
      RMS_cut[i] = -1.0;
      continue;
    }
    // accumate data for this pixel
    for (j=0, k=0; j<n; j++) {
      if (Pixel[j] == i)
	if (rms_cutoff > 0 && RMS[j] < rms_cutoff)
	  dat[k++] = RMS[j];
        else if (rms_cutoff < 0)
	  dat[k++] = RMS[j];
	  
    }
    if (k==0) {
      printf("Warning: no pixels found for pixel %d\n",i);
      RMS_cut[i] = -1.0;
      continue;
    }
      
    // sort the accumulated data since we need to use the quartiles
    // first compute a straight mean/median/std for all points
    qsort(dat,k,sizeof(float),compar_float);
    s0=s1=s2=0.0;
    for(j=0; j<k; j++) {
      s1 += dat[j];
      s2 += dat[j] * dat[j];
    }
    mean = s1/k;
    std = sqrt(s2/k - mean*mean);
    m1 = dat[k/4];
    m2 = dat[k/2];      // median
    m3 = dat[(k*3)/4];
    // printf("RMS count %d %d  %g %g %g \n",i, k, mean, m2, std);

    // 
    iqr = m3-m1;
    dlo = m1 - frob*iqr;
    dhi = m3 + frob*iqr;
    // now only take all data between dlo and dhi
    // and recompute the mean,median/std again
    s0=s1=s2=0.0;
    for(j=0, l=0; j<k; j++) {
      if (dat[j]<dlo || dat[j]>dhi) continue;
      l++;
      s1 += dat[j];
      s2 += dat[j] * dat[j];
    }
    mean = s1/l;
    std = sqrt(s2/l - mean*mean);
    if (rms_cutoff > 0)
      RMS_cut[i] = rms_cutoff;                       // classic hardcoded
    else
      RMS_cut[i] = mean - rms_cutoff * std;          // robust cutoff
      
    //printf("RMS count %d %d  %g %g %g   %g\n",i,l, mean, mean, std, RMS_cut[i]);    

    //printf("\n");
    
  }
  free(dat);
}

#if 0
// MAD = Median Absolute Deviation
//      (sigma = 1.4826 * mad for a normal distribution)
real mad_moment(Moment *m)
{
  real median, x;
  int i, n; 
  Moment tmp;

  if (m->ndat==0)
    error("mad_moment cannot be computed with ndat=%d",m->ndat);
  median = median_moment(m);
  n = MIN(m->n, m->ndat);
  ini_moment(&tmp,1,n);
  for (i=0; i<n; i++) {
    x = m->dat[i] - median;
    if (x > 0)
      accum_moment(&tmp,x,1.0);
    else
      accum_moment(&tmp,-x,1.0);
  }
  median = median_moment(&tmp);
  free_moment(&tmp);
  return median;
}
// MARD = Mean Absolute Relative Difference

real mard_moment(Moment *m)
{
  real mean, x;
  int i, n;
  Moment tmp;

  if (m->ndat==0)
    error("mard_moment cannot be computed with ndat=%d",m->ndat);
  mean = sum1/sum0;
  n = MIN(m->n, m->ndat);
  ini_moment(&tmp,1,n);
  for (i=0; i<n; i++) {
    x = m->dat[i] - mean;
    if (x > 0)
      accum_moment(&tmp,x,1.0);
    else
      accum_moment(&tmp,-x,1.0);
  }
  mean = mean_moment(&tmp);
  free_moment(&tmp);
  return mean;
}
#endif
