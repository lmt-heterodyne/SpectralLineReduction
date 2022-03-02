/** plane.c - methods for handling data plane
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Plane.h"
#include "fitsio.h"

/*************************  PLANE METHODS ***********************************/

/** initialize_plane - allocated data for the plane
 */
void initialize_plane(Plane* P, int *n)
{
  int i;
  for(i=0;i<2;i++)
    P->n[i] = n[i];
  P->nplane = n[0]*n[1];
  P->plane = (float*)malloc(P->nplane*sizeof(float));
  if(P->plane == NULL)
    fprintf(stderr,"Plane: Failed to Allocate Plane\n");
  else
    for(i=0;i<P->nplane;i++)
      P->plane[i] = 0.0;
}

/** initialize_axis - initializes the axis data 
 */
void initialize_plane_axis(Plane *P, int axis, float crval, float crpix, float cdelt, char *ctype, char *cunit)
{
  int i;

  P->crval[axis] = crval;
  P->crpix[axis] = crpix;
  P->cdelt[axis] = cdelt;
  strcpy(P->ctype[axis],ctype);
  strcpy(P->cunit[axis],cunit);

  P->caxis[axis] = (float*)malloc(P->n[axis]*sizeof(float));
  if(P->caxis[axis] == NULL)
    fprintf(stderr,"Plane: Failed to allocate axis %d\n",axis);

  for(i=0;i<P->n[axis];i++)
      P->caxis[axis][i] = (i+1-P->crpix[axis])*P->cdelt[axis]+P->crval[axis];

  printf("P-AXIS %d starts at %g\n",axis,  P->caxis[axis][0]);  
}

/** axis_index - finds the index in an axis array given a value
 */
int plane_axis_index(Plane *P, int axis, float value)
{
  float result_f;
  int result_i;
  
  result_f = (value - P->crval[axis])/P->cdelt[axis] + P->crpix[axis] - 0.5;
  result_i = (int)floor(result_f);
  
  if((result_i<0) || (result_i>=P->n[axis])) return -1;

  return result_i;
}

/** plane_index - finds the index of an element in the plane according to 
    given x and y positions
*/
int plane_index(Plane *P, float x, float y)
{
  int ix, iy;
  int result;

  ix = plane_axis_index(P, PLANE_X_AXIS, x);
  iy = plane_axis_index(P, PLANE_Y_AXIS, y);

  result = iy*P->n[PLANE_X_AXIS] + ix;
  return result;
}


void write_fits_plane(Plane *P, char *filename)
{
  int i,j,k,ii,ic;
  int retval, status;
  int naxis;
  long naxes[3], obsnum;
  float equinox;
  char radesys[20];
  float *buffer;
  fitsfile *fptr;

  char ctype[20], cunit[20], comment[60];
  float crval, cdelt, crpix, bmaj, bpa;
  
  naxis = 2;
  naxes[0] = P->n[PLANE_X_AXIS];
  naxes[1] = P->n[PLANE_Y_AXIS];

  // create the buffer to reorder the cube to FITS standard
  buffer = (float*)malloc(P->nplane*sizeof(float));


  ic = 0;
  i=0;
  for(j=0;j<P->n[PLANE_Y_AXIS];j++)
    for(k=0;k<P->n[PLANE_X_AXIS];k++)
      {
	ii = i + j*P->n[PLANE_X_AXIS] + (P->n[PLANE_X_AXIS]-k-1);
	buffer[ic] = P->plane[ii];
	ic++;
      }


  
  // you MUST initialize status
  status = 0;
  if((retval=fits_create_file(&fptr, filename, &status)) != 0)
    print_fits_error(status);

  if((retval=fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) != 0)
    print_fits_error(status);

  // scale axes to standards
  strcpy(ctype,"RA---SFL");          // nominal projection Sanson-Flamsteed
  crval = P->x_position;             // degrees
  cdelt = -P->cdelt[PLANE_X_AXIS] / 3600.;  // degrees - we flipped the RA axis
  crpix = P->crpix[PLANE_X_AXIS];
  strcpy(cunit,"deg     ");
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE1  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE1 %s\n",ctype);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL1  ", &crval, cunit, &status)) != 0)
    {
      printf("CRVAL1 %f\n",crval);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT1  ", &cdelt, cunit, &status)) != 0)
    {
      printf("CDELT1 %f\n",cdelt);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRPIX1  ", &crpix, " ", &status)) != 0)
    {
      printf("CRPIX1 %f\n",crpix);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "CUNIT1  ", cunit, " ", &status)) != 0)
    {
      printf("CUNIT1 %s\n",cunit);
      print_fits_error(status);
    }

  strcpy(ctype,"DEC--SFL");          // nominal projection Sanson-Flamsteed
  crval = P->y_position;             // degrees  
  cdelt = P->cdelt[PLANE_Y_AXIS] / 3600.;  // degrees 
  crpix = P->crpix[PLANE_Y_AXIS];
  strcpy(cunit,"deg     ");
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE2  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE2 %s\n",ctype);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL2  ", &crval, cunit, &status)) != 0)
    {
      printf("CRVAL2 %f\n",crval);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT2  ", &cdelt, cunit, &status)) != 0)
    {
      printf("CDELT2 %f\n",cdelt);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRPIX2  ", &crpix, " ", &status)) != 0)
    {
      printf("CRPIX2 %f\n",crpix);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "CUNIT2  ", cunit, " ", &status)) != 0)
    {
      printf("CUNIT2 %s\n",cunit);
      print_fits_error(status);
    }


  // write the data cube
  if((retval=fits_write_img(fptr, TFLOAT, 1, P->nplane, buffer, &status)) != 0)
    print_fits_error(status);

  // close the file
  if((retval=fits_close_file(fptr, &status)) != 0)
    print_fits_error(status);
}

//  The routines below are support to read a model, and for a given SpecFile
//  simulate how LMT would observe it. The -a flag (or --model) is used for this.
//
//  Important notes and deficiencies
//  - input map needs to be a jy/pixel convolved with the LMT beam.
//    and CRPIX is ignored (see next item)
//  - WCS is that of the original specfile
//  - X and Y confusion for non-sq. maps, this needs a real fix
//  - get_value uses "nearest" pixel, should do an interpolation

float get_value(Plane *P, float x, float y)
{
  int ix = plane_axis_index(P, PLANE_X_AXIS, x);
  int iy = plane_axis_index(P, PLANE_Y_AXIS, y);
  if (ix<0 || iy<0) return 0.0;
  int izp = plane_index(P, x, y);
  // should never hapen
  if (izp < 0) printf("PIX: BAD %g %g -> %d %d  %d\n",x,y,ix,iy,izp);
  return  P->plane[izp];
}

void read_fits_plane(Plane *P, char *filename)
{
  int i,j,k,ii,ic;
  int retval, status;
  int naxis;
  long naxes[3], obsnum;
  float equinox;
  char radesys[20];
  float *buffer;
  float nulval = 0.0;
  int anynul;
  fitsfile *fptr;

  char ctype[20], cunit[20], comment[60];
  float crval, cdelt, crpix, bmaj, bpa;

  // code not done yet
  printf("Warning: reading a model in channel-0 is an experimental feature\n");
  printf("There are still several limitations\n");

  status = 0;
  if((retval=fits_open_file(&fptr, filename, 0, &status)) != 0)
    print_fits_error(status);      

  naxis = 2;
  naxes[0] = naxes[1] = 0;
  if((retval=fits_get_img_size(fptr, naxis, naxes, &status)) != 0)
    print_fits_error(status);
  printf("NAXIS: %ld %ld\n",naxes[0],naxes[1]);
  
  P->n[PLANE_X_AXIS] = naxes[0];
  P->n[PLANE_Y_AXIS] = naxes[1];
  P->nplane = naxes[0] * naxes[1];

  // create the buffer to reorder the cube to FITS standard
  buffer = (float*)malloc(P->nplane*sizeof(float));
  P->plane = (float*)malloc(P->nplane*sizeof(float));

  //int fits_read_img / ffgpv
  //(fitsfile *fptr, int  datatype, long firstelem, long nelements,
  // DTYPE *nulval, > DTYPE *array, int *anynul, int *status)

  // read the data cube
  if((retval=fits_read_img(fptr, TFLOAT, 1, P->nplane, &nulval, buffer, &anynul, &status)) != 0)
    print_fits_error(status);


  if((retval=fits_read_key(fptr, TFLOAT, "CRPIX1", &P->crpix[0], comment, &status)) != 0)  
    print_fits_error(status);
  if((retval=fits_read_key(fptr, TFLOAT, "CRPIX2", &P->crpix[1], comment, &status)) != 0)  
    print_fits_error(status);
  printf("CRPIX: %g %g\n",P->crpix[0], P->crpix[1]);
  
  if((retval=fits_read_key(fptr, TFLOAT, "CDELT1", &P->cdelt[0], comment, &status)) != 0)  
    print_fits_error(status);
  if((retval=fits_read_key(fptr, TFLOAT, "CDELT2", &P->cdelt[1], comment, &status)) != 0)  
    print_fits_error(status);
  P->cdelt[0] *= 3600;
  P->cdelt[1] *= 3600;
  printf("CDELT: %g %g\n",P->cdelt[0], P->cdelt[1]);
  
  P->crval[0] = 0.0;
  P->crval[1] = 0.0;
  printf("CRVAL: %g %g\n",P->crval[0], P->crval[1]);


  ic = 0;
  i=0;
  for(j=0;j<P->n[PLANE_Y_AXIS];j++)
    for(k=0;k<P->n[PLANE_X_AXIS];k++)
      {
	ii = i + j*P->n[PLANE_X_AXIS] + (P->n[PLANE_X_AXIS]-k-1);
	P->plane[ii] = buffer[ic] ;
	ic++;
      }

  // close the file
  if((retval=fits_close_file(fptr, &status)) != 0)
    print_fits_error(status);

} 
