#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

void rms_stats(int n, float *RMS, int *Pixel,
	       int npix, int *use_pixel, float *RMS_cut,
	       float rms_cutoff);
