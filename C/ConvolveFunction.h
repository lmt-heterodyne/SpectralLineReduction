#ifndef _CONVOLVE_H_
#define _CONVOLVE_H_

#define CONVOLVE_UNDEFINED -1
#define CONVOLVE_BOX        0
#define CONVOLVE_JINC       1
#define CONVOLVE_GAUSS      2
#define CONVOLVE_TRIANGLE   3

typedef struct
{
  int type;
  int npts;
  int n_cells;
  float rmax, delta, resolution_size, cell_size;
  float jinc_a, jinc_b, jinc_c, gauss_b;
  float *array;
} ConvolveFunction;

void initialize_convolve_function(ConvolveFunction *CF, float resolution_size, float cell_size, float rmax, int npts);
float get_weight(ConvolveFunction *CF, float r);
void initialize_box_filter(ConvolveFunction *CF, float a);
void initialize_gauss_filter(ConvolveFunction *CF, float hpw);
void initialize_jinc_filter(ConvolveFunction *CF, float a, float b, float c);
void initialize_triangle_filter(ConvolveFunction *CF, float a);

#endif
