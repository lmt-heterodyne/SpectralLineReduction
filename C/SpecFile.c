#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <netcdf.h>
#include <math.h>

#include "OTFParameters.h"
#include "SpecFile.h"

static int nfiles = 0;

int read_spec_file(SpecFile *S, char *filename)
{
  int i,j,k;
  /* This will be the netCDF ID for the file and data variable. */
  int ncid;
  /* for error handling. */
  int retval;
  /* dimensions and varid's */
  size_t nspec, nchan, nchan0, chan0, ncal, npix;
  int nspec_id, nchan_id, ncal_id, npix_id, npixl_id, data_id, tsys_id, x_id, y_id, pix_id, seq_id, rms_id;
  int nchan0_id, chan0_id;
  char version[20];
  char history[512];

  if (access( filename, F_OK ) == 0 ) {
    printf("Opening SpecFile file \"%s\"\n",filename);
  } else {
    printf("SpecFile file \"%s\" does not exist\n",filename);    
    return -1;
  }
  
  nfiles++;
  printf("specfile %d: %s\n",nfiles,filename);
  S->nfiles = nfiles;
 
  /* Open the file. NC_NOWRITE tells netCDF we want read-only access to the file.*/
  if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
    ERR(retval);

  

  /* Get the dimensions */
  if ((retval = nc_inq_dimid(ncid, "nspec", &nspec_id)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "nchan", &nchan_id)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "nchan0", &nchan0_id)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "chan0", &chan0_id)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "ncal",  &ncal_id)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncid, "npix",  &npix_id)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, nspec_id, &nspec)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, nchan0_id, &nchan0)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, chan0_id, &chan0)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, nchan_id, &nchan)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, ncal_id,  &ncal)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, npix_id,  &npix)))
    ERR(retval);

  int obsnum_id, mapcoord_id, source_id,source_x,source_y,crval_id,crpix_id,cdelt_id,ctype_id,caxis_id;
  int deltaf_id, deltat_id;
  int rf_id, vlsr_id, do_id, do_rc, do_version, do_history;
  /* Get the varids of the observation header */
  if ((retval = nc_inq_varid(ncid, "Header.Obs.ObsNum", &obsnum_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.MapCoord", &mapcoord_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.SourceName", &source_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.XPosition", &source_x)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.YPosition", &source_y)))
    ERR(retval);
  // PJT hack
  if ((retval = nc_inq_varid(ncid, "Header.LineData.LineRestFrequency", &rf_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.LineData.VSource", &vlsr_id)))
    ERR(retval);
  // Header.LineData.ZSource = 0
  
  if ((retval = nc_inq_varid(ncid, "Header.LineData.DeltaFrequency", &deltaf_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.DumpTime", &deltat_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.DateObs", &do_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Obs.Receiver", &do_rc)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.Version", &do_version)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "Header.History", &do_history)))
    ERR(retval);
  

  // printf("Header.Obs complete\n");

  /* Get the varids of the spectrum axis header variables */
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CRVAL", &crval_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CRVAL\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CRPIX", &crpix_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CRPIX\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CDELT", &cdelt_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CDELT\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CTYPE", &ctype_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.CTYPE\n");
  if ((retval = nc_inq_varid(ncid, "Header.SpectrumAxis.CAXIS", &caxis_id)))
    ERR(retval);
  //printf("Header.SpectrumAxis.XAxis\n");

  /* Get the varid of the data variables, based on its name. */
  //printf("Data.Spectra\n");
  if ((retval = nc_inq_varid(ncid, "Data.Spectra", &data_id)))
    ERR(retval);
  //printf("Data.Spectra completed\n");
  if ((retval = nc_inq_varid(ncid, "Data.XPos", &x_id)))
    ERR(retval);
  //printf("Data.XPos\n");
  if ((retval = nc_inq_varid(ncid, "Data.YPos", &y_id)))
    ERR(retval);
  //printf("Data.YPos\n");
  if ((retval = nc_inq_varid(ncid, "Data.Pixel", &pix_id)))
    ERR(retval);
  //printf("Data.Pixel\n");
  if ((retval = nc_inq_varid(ncid, "Data.Sequence", &seq_id)))
    ERR(retval);
  //printf("Data.Sequence\n");
  if ((retval = nc_inq_varid(ncid, "Data.RMS", &rms_id)))
    ERR(retval);
  //printf("Data.RMS\n");
  if ((retval = nc_inq_varid(ncid, "Data.Tsys", &tsys_id)))
    ERR(retval);
  // printf("Data.Tsys\n");

  //  if ((retval = nc_inq_varid(ncid, "Data.PixelList", &pixl_id)))
  //    ERR(retval);

  /* Read the data and load the SpecFile struct */
  printf("file: %s nspec= %zu nchan= %zu npix=%zu ncal=%zu\n",filename,nspec,nchan,npix,ncal);
  S->nspec = nspec;
  S->nchan = nchan;
  S->nchan0 = nchan0;
  S->chan0 = chan0;
  S->npix  = npix;
  S->ncal  = ncal;

  if((retval = nc_get_var_int(ncid,obsnum_id, &S->obsnum)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_int(ncid,mapcoord_id, &S->map_coord)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid,source_id, S->source)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, source_x, &S->x_positionN)) != NC_NOERR)
    ERR(retval);
  if (nfiles == 1) S->x_position = S->x_positionN;
  if((retval = nc_get_var_double(ncid, source_y, &S->y_positionN)) != NC_NOERR)
    ERR(retval);
  if (nfiles == 1) S->y_position = S->y_positionN;
  if((retval = nc_get_var_double(ncid, rf_id, &S->restfreq)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_float(ncid, vlsr_id, &S->vlsr)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_float(ncid, deltaf_id, &S->deltaf)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_float(ncid, deltat_id, &S->deltat)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid,do_id, S->date_obs)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid,do_rc, S->receiver)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid,do_version, version)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid,do_history, S->history)) != NC_NOERR)
    ERR(retval);
  printf("SpecFile version %s\n",version);
  

  if((retval = nc_get_var_double(ncid, crval_id, &S->CRVAL)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, crpix_id, &S->CRPIX)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var_double(ncid, cdelt_id, &S->CDELT)) != NC_NOERR)
    ERR(retval);
  if((retval = nc_get_var(ncid, ctype_id, S->CTYPE)) != NC_NOERR)
    ERR(retval);
  S->CAXIS = (double*)malloc(nchan*sizeof(double));
  if((retval = nc_get_var(ncid, caxis_id, S->CAXIS)) != NC_NOERR)
    ERR(retval);

  S->theData = (float *)malloc(nspec*nchan*sizeof(float));
  if(S->theData == NULL) {
    fprintf(stderr,"SpecFile: Error allocating data array\n");
    exit(1);
  }
  S->XPos = (float *)malloc(nspec*sizeof(float));
  S->YPos = (float *)malloc(nspec*sizeof(float));
  S->RMS  = (float *)malloc(nspec*sizeof(float));
  S->Tsys = (float *)malloc(ncal*npix*nchan*sizeof(float));
  if (S->Tsys == NULL) {
    fprintf(stderr,"SpecFile: Error allocating tsys array\n");
    exit(1);
  }
  S->Pixel = (int *)malloc(nspec*sizeof(int));
  S->Sequence = (int *)malloc(nspec*sizeof(int));
  S->RMS_cut = (float *)malloc(MAXPIXEL*sizeof(float));   // MAXPIXEL
  S->use = (int *)malloc(nspec*sizeof(int));

  if ((retval = nc_get_var_float(ncid, data_id, S->theData)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, x_id, S->XPos)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, y_id, S->YPos)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, rms_id, S->RMS)))
    ERR(retval);
  if ((retval = nc_get_var_float(ncid, tsys_id, S->Tsys)))
    ERR(retval);
  if ((retval = nc_get_var_int(ncid, pix_id, S->Pixel)))
    ERR(retval);
  if ((retval = nc_get_var_int(ncid, seq_id, S->Sequence)))
    ERR(retval);
  if (S->x_positionN != S->x_position ||
      S->y_positionN != S->y_position) {
    // @todo  correct the XPos and YPos for a new center, using SFL projection rule
    double dx = -(S->x_positionN - S->x_position) * cos(S->y_position/57.2958) * 3600;
    double dy = +(S->y_positionN - S->y_position) * 3600;
    printf("Fixing offsets by adding dx=%f dy=%f arcsec for specfile %d\n", dx, dy, S->nfiles);
    for (i=0; i<nspec; i++) {
      S->XPos[i] += dx;
      S->YPos[i] += dy;
    }
  }

  for (i=0; i<nspec; i++)
    S->use[i] = 1;

#if 1
  float tsys_min;
  float tsys_max;
  float tsys, tsys_1, tsys_2;
  tsys_min = tsys_max = S->Tsys[0];
  for (i=1; i<ncal*npix*nchan; i++) {
    if (S->Tsys[i] < tsys_min) tsys_min = S->Tsys[i]; 
    if (S->Tsys[i] > tsys_max) tsys_max = S->Tsys[i]; 
  }
  printf("Global Tsys min/max = %g %g\n", tsys_min, tsys_max);
  // take an average of all nchan values and stuff if in the first channel
  // for get_tsys()
  int idx;
  for (i=0; i<ncal; i++) {
    for (j=0; j<npix; j++) {
      tsys_1 = tsys_2 = 0.0;
      idx = i*npix*nchan + j*nchan;
      for (k=0; k<nchan; k++) {
	tsys = S->Tsys[idx+k];
	tsys_1 += tsys;
	tsys_2 += tsys*tsys;
      }
      tsys = tsys_1/nchan;
      printf("Tsys[%d] = %g +/- %g\n", j, tsys, sqrt(tsys_2/nchan - tsys*tsys));
      S->Tsys[idx] = tsys;
    }
  }
  
#endif

  // compute and report on the extreme positions the array has seen

  double xmin = S->XPos[0],
         xmax = S->XPos[0],
         ymin = S->YPos[0],
         ymax = S->YPos[0];
  double bs = 27.8;      // spacings between beams in the 4x4 grid
         
  for (i=1; i<nspec; i++) {
    if (S->XPos[i] < xmin) xmin = S->XPos[i];
    if (S->XPos[i] > xmax) xmax = S->XPos[i];
    if (S->YPos[i] < ymin) ymin = S->YPos[i];
    if (S->YPos[i] > ymax) ymax = S->YPos[i];
  }
  printf("First/Last beam: [mapcoord=%d] %g %g  %g %g\n",S->map_coord, S->XPos[0], S->YPos[0], S->XPos[nspec-1], S->YPos[nspec-1]);
  printf("X-range: %g %g   Y-range: %g %g arcsec\n",xmin,xmax,ymin,ymax);
  printf("X-ramp:  %g %g   Y-ramp:  %g %g arcsec\n",xmin+3*bs,xmax-3*bs,ymin+3*bs,ymax-3*bs);
  printf("MapSize: %g x %g arcsec\n", xmax-xmin, ymax-ymin);

  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    ERR(retval);
  return 0;
}

void free_spec_file(SpecFile *S)
{
  free(S->Sequence);
  free(S->Pixel);
  free(S->RMS);
  free(S->YPos);
  free(S->XPos);
  free(S->theData);
  free(S->CAXIS);
}

//  S    Spectrum
//  i    spectrum index (0..nspec)
//  returns the address where a spectrum[nchan] is located
float *get_spectrum(SpecFile *S, int i)
{
  int index;

  index = i*S->nchan;
  return(&(S->theData[index]));
}

//   S     spectrum
//   ical  which cal
//   jpix  which pix
//   returns the Tsys in the first element of the Tsys[chan] array
float get_tsys(SpecFile *S, int ical, int jpix)
{
  int index = ical*S->npix*S->nchan + jpix*S->nchan;
  return S->Tsys[index];
}
