#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ConvolveFunction.h"

double jinc(double r);

void initialize_convolve_function(ConvolveFunction *CF, float resolution_size, float cell_size, float rmax, int npts)
{
  CF->type = CONVOLVE_UNDEFINED;
  CF->resolution_size = resolution_size;     // this is lambda/D  (PJT: or should it be 1.15 lambda/D ?)
  CF->cell_size = cell_size;                 // actual cell size
  CF->rmax = rmax;                           // this is number of lambda/D's to cut off
  CF->npts = npts;                           // OTF.nsamples (usually 256)
  CF->delta = rmax*resolution_size/(npts-1); // in units of arcsec
  CF->n_cells = (int) floor(rmax*CF->resolution_size/CF->cell_size) + 1;
  CF->array = (float*)malloc(CF->npts*sizeof(float));
}

float get_weight(ConvolveFunction *CF, float r)
{
  int index;
  float result;
  index = (int) floor(r/CF->delta);
  if(index>=CF->npts)
    result = 0.0;
  else
    result = CF->array[index];
  return(result);
}

void initialize_box_filter(ConvolveFunction *CF, float a)
{
  int i;
  float x;
  CF->type = CONVOLVE_BOX;
  for(i=0;i<CF->npts;i++)
    {
      x = (i * CF->delta);
      if(x<=a)
	CF->array[i] = 1.0;
      else
	CF->array[i] = 0.0;
    }
}

void initialize_triangle_filter(ConvolveFunction *CF, float a)
{
  int i;
  float x;
  CF->type = CONVOLVE_BOX;
  for(i=0;i<CF->npts;i++)
    {
      x = (i * CF->delta);
      if(x<=a)
	CF->array[i] = 1.0 - x/a;
      else
	CF->array[i] = 0.0;
    }
}

void initialize_gauss_filter(ConvolveFunction *CF, float a)
{
  int i;
  float x;

  CF->type = CONVOLVE_GAUSS;
  CF->gauss_b = a;  // at gauss_b=1 this follows exactly the FWHM 
  for(i=0;i<CF->npts;i++)
    {
      // i*delta is distance in arcsec, so we must normalize by lambda/D
      x = (i * CF->delta)/CF->resolution_size / CF->gauss_b;
      CF->array[i] = exp(-2.77258872*x*x);   //  2.77.. = 4*ln(2)
    }
}

void initialize_jinc_filter(ConvolveFunction *CF, float a, float b, float c)
{
  int i;
  float x;
  float w1, w2, w3;

  CF->type = CONVOLVE_JINC;
  CF->jinc_a = a;
  CF->jinc_b = b;
  CF->jinc_c = c;

  for(i=0;i<CF->npts;i++)
    {
      // i*delta is distance in arcsec so we must normalize by lambda/D
      x = (i*CF->delta)/CF->resolution_size;
      w1 = 2.0*jinc(6.28318531*x/a);
      w2 = exp(-pow(2.0*x/b,c));
      w3 = 2.0*jinc(3.831706*x/CF->rmax);
      CF->array[i] = w1*w2*w3;
    }
}

// This is the jinc (J1(x)/x) function from the FCRAO program
double jinc(double r)
{
  double arg,f1,theta,j;

  if (r <=3.0) 
    {
      arg=pow(r/3.0,2.0);
      j=0.5+arg*(-0.56249985+arg*(0.21093573+arg*(-0.03954289+arg*(0.00443319+arg*(-0.00031761+arg*0.00001109)))));
    }
  else 
    {
      arg=3.0/r;
      f1=0.79788456+arg*(0.00000156+arg*(0.01659667+arg*(0.00017105+arg*(-.00249511+arg*(0.00113653-arg*0.00020033)))));
      theta=r-2.35619449+arg*(0.12499612+arg*(0.00005650+arg*(-.00637879+arg*(0.00074348+arg*(0.00079824-arg*0.00029166))))) ;
   
      j=f1*cos(theta)/(r*sqrt(r));
    }
  return(j);
}

