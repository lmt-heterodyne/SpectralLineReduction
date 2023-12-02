/** cube.c - methods for handling data cube
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <netcdf.h>
#include "Cube.h"
#include "SpecFile.h"
#include "Version.h"
#include "fitsio.h"

/*************************  CUBE METHODS ***********************************/

/** initialize_cube - allocated data for the cube
 */
void initialize_cube(Cube* C, int *n)
{
  int i;
  for(i=0;i<3;i++)
    C->n[i] = n[i];
  C->ncube = n[0]*n[1]*n[2];
  C->nplane = n[0]*n[1];
  C->cube = (float*)malloc(C->ncube*sizeof(float));
  if(C->cube == NULL)
    {
      fprintf(stderr,"Cube: Failed to Allocate Cube\n");
      exit(1);
    }
  else
    for(i=0;i<C->ncube;i++)
      C->cube[i] = 0.0;
}

/** initialize_axis - initializes the axis data 
 */
void initialize_cube_axis(Cube *C, int axis, float crval, float crpix, float cdelt, char *ctype, char *units)
{
  int i, nn;

  C->crval[axis] = crval;
  C->crpix[axis] = crpix;
  C->cdelt[axis] = cdelt;
  for(i=0;i<16;i++)
    {
      C->ctype[axis][i] = '\0';
      C->units[axis][i] = '\0';
    }
  strcpy(C->ctype[axis],ctype);
  strcpy(C->units[axis],units);

  C->caxis[axis] = (float*)malloc(C->n[axis]*sizeof(float));
  if(C->caxis[axis] == NULL)
    {
      fprintf(stderr,"Cube: Failed to allocate axis %d\n",axis);
      exit(1);
    }

  for(i=0;i<C->n[axis];i++)
      C->caxis[axis][i] = (i+1-C->crpix[axis])*C->cdelt[axis]+C->crval[axis];
  
  printf("C-AXIS %d starts at %g\n",axis,  C->caxis[axis][0]);
}

/** axis_index - finds the index in an axis array given a value
 */
int cube_axis_index(Cube *C, int axis, float value)
{
  float result_f;
  int result_i;

  result_f = (value - C->crval[axis])/C->cdelt[axis] + C->crpix[axis] - 0.5;
  result_i = (int)floor(result_f);
  
  if((result_i<0) || (result_i>=C->n[axis])) return -1;

  return result_i;
}

/** cube_z_index - finds the index of the first element of z array in cube
    given x and y positions
*/
int cube_z_index(Cube *C, float x, float y)
{
  int ix, iy;
  int result;

  ix = cube_axis_index(C, X_AXIS, x);
  iy = cube_axis_index(C, Y_AXIS, y);

  result = iy*C->n[X_AXIS]*C->n[Z_AXIS] + ix*C->n[Z_AXIS];
  return result;
}

/** write netcdf data cube
 *  not used anymore, we use fits cubes
 */
void write_netcdf_cube(Cube *C, char *filename)
{
  int ncid;
  int retval;
  int dimid_x, dimid_y, dimid_z, dimid_label, dimid_T;
  int T_dimids[3], l_dimids[3];
  int varid_obsnum, varid_source;
  int varid_NAXIS1, varid_NAXIS2, varid_NAXIS3;
  int varid_CRVAL1, varid_CRPIX1, varid_CDELT1, varid_CTYPE1, varid_CAXIS1; 
  int varid_CRVAL2, varid_CRPIX2, varid_CDELT2, varid_CTYPE2, varid_CAXIS2; 
  int varid_CRVAL3, varid_CRPIX3, varid_CDELT3, varid_CTYPE3, varid_CAXIS3; 

  if ((retval = nc_create(filename, NC_NOCLOBBER, &ncid)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_dim(ncid, "x", (size_t)C->n[X_AXIS], &dimid_x)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_dim(ncid, "y", (size_t)C->n[Y_AXIS], &dimid_y)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_dim(ncid, "z", (size_t)C->n[Z_AXIS], &dimid_z)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_dim(ncid, "label", (size_t)16, &dimid_label)) != NC_NOERR)
    ERR(retval);
  l_dimids[0] = dimid_label;
  l_dimids[1] = dimid_label;
  l_dimids[2] = dimid_label;


  if((retval = nc_def_var(ncid, "Header.Obs.ObsNum", NC_INT, 0, &dimid_y, &varid_obsnum)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "Header.Obs.SourceName", NC_CHAR, 1, &dimid_label, &varid_source)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_var(ncid, "NAXIS1", NC_INT, 0, &dimid_z, &varid_NAXIS1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRVAL1", NC_FLOAT, 0, &dimid_z, &varid_CRVAL1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRPIX1", NC_FLOAT, 0, &dimid_z, &varid_CRPIX1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CDELT1", NC_FLOAT, 0, &dimid_z, &varid_CDELT1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CTYPE1", NC_CHAR, 1, &dimid_label, &varid_CTYPE1)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CAXIS1", NC_FLOAT, 1, &dimid_z, &varid_CAXIS1)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_var(ncid, "NAXIS2", NC_INT, 0, &dimid_x, &varid_NAXIS2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRVAL2", NC_FLOAT, 0, &dimid_x, &varid_CRVAL2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRPIX2", NC_FLOAT, 0, &dimid_x, &varid_CRPIX2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CDELT2", NC_FLOAT, 0, &dimid_x, &varid_CDELT2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CTYPE2", NC_CHAR, 1, &dimid_label, &varid_CTYPE2)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CAXIS2", NC_FLOAT, 1, &dimid_x, &varid_CAXIS2)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_def_var(ncid, "NAXIS3", NC_INT, 0, &dimid_y, &varid_NAXIS3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRVAL3", NC_FLOAT, 0, &dimid_y, &varid_CRVAL3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CRPIX3", NC_FLOAT, 0, &dimid_y, &varid_CRPIX3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CDELT3", NC_FLOAT, 0, &dimid_y, &varid_CDELT3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CTYPE3", NC_CHAR, 1, &dimid_label, &varid_CTYPE3)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_def_var(ncid, "CAXIS3", NC_FLOAT, 1, &dimid_y, &varid_CAXIS3)) != NC_NOERR)
    ERR(retval);


  T_dimids[0] = dimid_y;
  T_dimids[1] = dimid_x;
  T_dimids[2] = dimid_z;
  if((retval = nc_def_var(ncid, "T", NC_FLOAT, 3, T_dimids, &dimid_T)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_enddef(ncid)) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_int(ncid, varid_obsnum, &C->obsnum[0])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_source, C->source)) != NC_NOERR)
    ERR(retval);



  if((retval = nc_put_var_int(ncid, varid_NAXIS1, &C->n[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRVAL1, &C->crval[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRPIX1, &C->crpix[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CDELT1, &C->cdelt[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_CTYPE1, C->ctype[Z_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CAXIS1, C->caxis[Z_AXIS])) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_int(ncid, varid_NAXIS2, &C->n[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRVAL2, &C->crval[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRPIX2, &C->crpix[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CDELT2, &C->cdelt[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_CTYPE2, C->ctype[X_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CAXIS2, C->caxis[X_AXIS])) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_int(ncid, varid_NAXIS3, &C->n[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRVAL3, &C->crval[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CRPIX3, &C->crpix[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CDELT3, &C->cdelt[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var(ncid, varid_CTYPE3, C->ctype[Y_AXIS])) != NC_NOERR)
    ERR(retval);
  if((retval = nc_put_var_float(ncid, varid_CAXIS3, C->caxis[Y_AXIS])) != NC_NOERR)
    ERR(retval);

  if((retval = nc_put_var_float(ncid, dimid_T, C->cube)) != NC_NOERR)
    ERR(retval);
	   
  nc_close(ncid);  
}

void write_fits_cube(Cube *C, char *filename)
{
  int i,j,k,ii,ic;
  int retval, status;
  int naxis;
  long naxes[3]; // obsnum;
  float equinox;
  char radesys[20];
  float *buffer;
  fitsfile *fptr;

  char ctype[20], cunit[20], comment[60];
  float crval, cdelt, crpix, bmaj, bpa;
  
  naxis = 3;
  naxes[0] = C->n[X_AXIS];
  naxes[1] = C->n[Y_AXIS];
  naxes[2] = C->n[Z_AXIS];

  // create the buffer to reorder the cube to FITS standard
  buffer = (float*)malloc(C->ncube*sizeof(float));
  if(buffer == NULL)
    {
      fprintf(stderr,"write_fits_cube: Failed to allocate buffer\n");
      exit(1);
    }
  
  ic = 0;
  for(i=0;i<C->n[Z_AXIS];i++)
      for(j=0;j<C->n[Y_AXIS];j++)
	  for(k=0;k<C->n[X_AXIS];k++)
	    {
	      ii = i + j*C->n[X_AXIS]*C->n[Z_AXIS] + (C->n[X_AXIS]-k-1)*C->n[Z_AXIS];
	      buffer[ic] = C->cube[ii];
	      ic++;
	    }

  // you MUST initialize status
  status = 0;
  if((retval=fits_create_file(&fptr, filename, &status)) != 0)
    print_fits_error(status);

  if((retval=fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) != 0)
    print_fits_error(status);

  // Write the header variables
  if((retval=fits_update_key(fptr, TSTRING, "TELESCOP", "LMT", " ", &status)) != 0)
    {
      printf("TELESCOP\n");
      print_fits_error(status);
    }
  // INSTRUME is not in https://heasarc.gsfc.nasa.gov/docs/fcg/common_dict.html
  strcpy(comment,"SEQUOIA, 1MM, or OMAYA");  
  if((retval=fits_update_key(fptr, TSTRING, "INSTRUME", C->receiver, comment, &status)) != 0)  
    {
      printf("INSTRUME\n");
      print_fits_error(status);
    }
  
  // @todo OBSERVER (the PI ?)
  
  if((retval=fits_update_key(fptr, TSTRING, "OBJECT  ", C->source, " ", &status)) != 0)
    {
      printf("OBJECT\n");
      print_fits_error(status);
    }
  strcpy(comment,"Data.TelescopeBackend.TelTime");
  if((retval=fits_update_key(fptr, TSTRING, "DATE-OBS", C->date_obs, comment, &status)) != 0)
    {
      printf("DATE-OBS\n");
      print_fits_error(status);
    }
  strcpy(comment,"LMTSLR Software version");
  if((retval=fits_update_key(fptr, TSTRING, "ORIGIN  ", LMTSLR_VERSION, comment, &status)) != 0)
    {
      printf("ORIGIN\n");
      print_fits_error(status);
    }

  // BUNIT
  strcpy(cunit,"K       ");
  if((retval=fits_update_key(fptr, TSTRING, "BUNIT   ", cunit, " ", &status)) != 0)
    {
      printf("BUNIT %s\n",cunit);
      print_fits_error(status);
    }

  // scale axes to standards
  if (C->map_coord == 2)
    strcpy(ctype,"GLON-SFL");          // nominal projection Sanson-Flamsteed
  else
    strcpy(ctype,"RA---SFL");          // nominal projection Sanson-Flamsteed
  
  crval = C->x_position;               // degrees
  cdelt = -C->cdelt[X_AXIS] / 3600.;   // degrees - we flipped the RA axis
  crpix = C->crpix[X_AXIS];
  strcpy(cunit,"deg     ");
  sprintf(comment,"map_coord=%d", C->map_coord);
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE1  ", ctype, comment, &status)) != 0)
    {
      printf("CTYPE1 %s\n",ctype);
      print_fits_error(status);
    }
  strcpy(comment,"deg Header.Obs.XPosition");
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL1  ", &crval, comment, &status)) != 0)
    {
      printf("CRVAL1 %f\n",crval);
      print_fits_error(status);
    }
  sprintf(comment,"deg (%g arcsec)", C->cdelt[X_AXIS]);
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT1  ", &cdelt, comment, &status)) != 0)
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

  if (C->map_coord == 2)
    strcpy(ctype,"GLAT-SFL");          // nominal projection Sanson-Flamsteed
  else
    strcpy(ctype,"DEC--SFL");          // nominal projection Sanson-Flamsteed
  crval = C->y_position;             // degrees
  cdelt = C->cdelt[Y_AXIS] / 3600.;  // degrees 
  crpix = C->crpix[Y_AXIS];
  strcpy(cunit,"deg     ");
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE2  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE2 %s\n",ctype);
      print_fits_error(status);
    }
  strcpy(comment,"deg Header.Obs.YPosition");  
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL2  ", &crval, comment, &status)) != 0)
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

  // MIRIAD says the VELO is an AIPS convention
  // CASA wants VRAD
  // strcpy(ctype,"VELO-LSR");        // this becomes VOPT in CASA, unless you add VELREF=257 ,then it's ok
  // VELO-F2V (but VELO-V2F fails) should be same, but is not 
  strcpy(ctype,"VRAD");       
  crval = C->crval[Z_AXIS]*1000.;    // m/s
  cdelt = C->cdelt[Z_AXIS]*1000.;    // m/s
  crpix = C->crpix[Z_AXIS];
  strcpy(cunit,"m/s     ");
  
  if((retval=fits_update_key(fptr, TSTRING, "CTYPE3  ", ctype, " ", &status)) != 0)
    {
      printf("CTYPE3 %s\n",ctype);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRVAL3  ", &crval, cunit, &status)) != 0)
    {
      printf("CRVAL3 %f\n",crval);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CDELT3  ", &cdelt, cunit, &status)) != 0)
    {
      printf("CDELT3 %f\n",cdelt);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "CRPIX3  ", &crpix, " ", &status)) != 0)
    {
      printf("CRPIX3 %f\n",crpix);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "CUNIT3  ", cunit, " ", &status)) != 0)
    {
      printf("CUNIT3 %s\n",cunit);
      print_fits_error(status);
    }

  // Alternate CTYPE3F = 'FREQ' could be made here
  // SPECSYSF = 'LSRK'
  // SSYSOBSF = 'TOPOCENT'
  // CDELT3F,CRPIX3F,CRVAL3F

  // The original spectral axis is lost in the current workflow
  // the --slide keyword in the HISTORY does not record the original channel numbers


  sprintf(comment,"deg   (%g arcsec)  ", C->resolution_size);
  bmaj = C->resolution_size / 3600.;
  if((retval=fits_update_key(fptr, TFLOAT,  "BMAJ    ", &bmaj, comment, &status)) != 0)
    {
      printf("BMAJ %f\n",bmaj);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "BMIN    ", &bmaj, comment, &status)) != 0)
    {
      printf("BMIN %f\n",bmaj);
      print_fits_error(status);
    }
  bpa = 0.0;
  if((retval=fits_update_key(fptr, TFLOAT,  "BPA     ", &bpa, cunit, &status)) != 0)
    {
      printf("BPA %f\n",bmaj);
      print_fits_error(status);
    }
      
  strcpy(radesys,"FK5     ");
  equinox = 2000.;
  if((retval=fits_update_key(fptr, TFLOAT,  "EQUINOX ", &equinox, " ", &status)) != 0)
    {
      printf("EQUINOX %f\n", equinox);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TSTRING, "RADESYS ", radesys, " ", &status)) != 0)
    {
      printf("RADESYS %s\n", radesys);
      print_fits_error(status);
    }


  // Might be nice if these are needed for any back-tracking in the future:
  // Header.Sky.ObsVel
  // Header.Sky.BaryVel
  
  strcpy(comment,"VLSR, VSOURCE and ZSOURCE need to be properly looked at ");
  fits_write_comment(fptr, comment, &status);

  strcpy(comment, "Header.LineData.VSource");
  if((retval=fits_update_key(fptr, TFLOAT,  "VLSR", &C->vlsr, comment, &status)) != 0)
    {
      printf("VLSR %f\n", C->vlsr);
      print_fits_error(status);
    }
  if((retval=fits_update_key(fptr, TFLOAT,  "VSOURCE", &C->vlsr, comment, &status)) != 0)
    {
      printf("VSOURCE %f\n", C->vlsr);
      print_fits_error(status);
    }
  // alma: ZSOURCE  (even astropy claims this is the preferred one)
  // alma: LINTRN

  strcpy(comment, "Header.LineData.LineRestFrequency");
  if((retval=fits_update_key(fptr, TDOUBLE,  "RESTFRQ", &C->restfreq, comment, &status)) != 0)
    {
      printf("RESTFRQ %f\n", C->restfreq);
      print_fits_error(status);
    }

  int velref = 257;
  strcpy(comment,"greisen?");
  if((retval=fits_update_key(fptr, TINT,  "VELREF", &velref, comment, &status)) != 0)
    {
      printf("VELREF %d\n", velref);
      print_fits_error(status);
    }
  

  strcpy(comment,"number of channels of original cube");
  if((retval=fits_update_key(fptr, TINT,  "NCHAN0", &C->nchan0, comment, &status)) != 0)
    {
      printf("NCHAN0 %d\n", C->nchan0);
      print_fits_error(status);
    }
  strcpy(comment,"first channel from original cube");
  if((retval=fits_update_key(fptr, TINT,  "CHAN0", &C->chan0, comment, &status)) != 0)
    {
      printf("NCHAN0 %d\n", C->chan0);
      print_fits_error(status);
    }
  

  // @todo this is assumed here, but should be fixed if upstream not VELO-LSR was choosen
  // LSRK vs. LSRD
  strcpy(cunit,"LSRK");
  strcpy(comment,"depends on --x_axis ?");
  if((retval=fits_update_key(fptr, TSTRING,  "SPECSYS", cunit, comment, &status)) != 0)
    {
      printf("SPECSYS %s\n", cunit);
      print_fits_error(status);
    }

  // can't hurt, right?
  strcpy(cunit,"TOPOCENT");
  strcpy(comment,"[Greisen 2009]");
  if((retval=fits_update_key(fptr, TSTRING,  "SSYSOBS", cunit, comment, &status)) != 0)
    {
      printf("SSYSOBS %s\n", cunit);
      print_fits_error(status);
    }

  if (1)
    {
      char toutstr[200];
      time_t t;
      struct tm *tmp;
      const char *fmt = "%Y-%m-%dT%H:%M:%S";

      t = time(NULL);
      tmp = localtime(&t);
      // if (tmp==NULL) error("localtime failed");
      strftime(toutstr, sizeof(toutstr), fmt, tmp);
      // error("strftime returned 0");
      strcpy(comment,"Date of file written");
      if((retval=fits_update_key(fptr, TSTRING,  "DATE", toutstr, comment, &status)) != 0)
	{
	  printf("DATE\n");
	  print_fits_error(status);
	}
    }

  strcpy(comment,"LMT observing numbers in this data: (Header.Obs.ObsNum)");
  fits_write_comment(fptr, comment, &status);
  for (i=0; i<C->nobsnum; i++)
    {
      sprintf(comment,"OBSNUM %d", C->obsnum[i]);   // or should this be %06d
      fits_write_comment(fptr, comment, &status);
    }

  strcpy(comment,"Convert RAW to SpecFile");
  fits_write_comment(fptr, comment, &status);
  fits_write_history(fptr, C->history1, &status);

  strcpy(comment,"Convert SpecFile to FITS");
  fits_write_comment(fptr, comment, &status);
  fits_write_history(fptr, C->history2, &status);

  // write the data cube
  if((retval=fits_write_img(fptr, TFLOAT, 1, C->ncube, buffer, &status)) != 0)
    print_fits_error(status);

  // close the file
  if((retval=fits_close_file(fptr, &status)) != 0)
    print_fits_error(status);

} 

void print_fits_error(int status)
{
  char e[81];
  if(status)
    {
      fits_get_errstatus(status, e);
      fprintf(stderr, "%s\n",e);
      exit(status);
    }
}
