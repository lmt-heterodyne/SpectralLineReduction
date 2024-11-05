#ifndef _CUBE_H_
#define _CUBE_H_

#define X_AXIS 1
#define Y_AXIS 0
#define Z_AXIS 2

#include "Version.h"

typedef struct 
{
  int obsnum[MAXOBS];
  int nobsnum;
  char source[40];
  char receiver[40];   // FITS 'INSTRUME'
  int map_coord;       // 0=azel  1=radec  2=latlon
  float x_position, y_position;
  float *cube;
  float *caxis[3];
  float crval[3], crpix[3], cdelt[3];
  char ctype[3][16], units[3][16];
  int n[3],ncube,nplane;
  int nchan0, chan0;        // provenance from original cube
  // new for provenance    ; or use a OTF pointer
  // still todo: HISTORY
  double restfreq;
  float vlsr;
  float zsource;
  char date_obs[40];        // eg 2021-09-14T07:25:23.370   = 23 chars
  float resolution_size;
  char history1[MAXHIST];
  char history2[MAXHIST];
} Cube;

void initialize_cube(Cube* C, int *n);

void initialize_cube_axis(Cube* C, int axis, float crval, float crpix, float cdelt, char *ctype, char *units);

int cube_axis_index(Cube* C, int axis, float value);

int cube_z_index(Cube* C, float x, float y);

void write_netcdf_cube(Cube* C, char *filename); 

void write_fits_cube(Cube *C, char *filename, char *comment);

void print_fits_error(int);

#endif
