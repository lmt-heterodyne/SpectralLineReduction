#ifndef _PLANE_H_
#define _PLANE_H_

#define PLANE_X_AXIS 1
#define PLANE_Y_AXIS 0

typedef struct
{
  float *plane;
  float x_position, y_position;
  int   map_coord;
  float *caxis[2];
  float crval[2], crpix[2], cdelt[2];
  char ctype[2][16];
  char cunit[2][16];
  int n[2],nplane;
} Plane;

/** Plane methods
 */
  
void initialize_plane(Plane* P, int *n);

void initialize_plane_axis(Plane *P, int axis, float crval, float crpix, float cdelt, char *ctype, char *cunit);

int plane_axis_index(Plane *P, int axis, float value);

int plane_index(Plane *P, float x, float y);

void write_fits_plane(Plane *P, char *filename);
void read_fits_plane(Plane *P, char *filename);
float get_value(Plane *P, float x, float y);
void print_fits_error(int);
#endif
