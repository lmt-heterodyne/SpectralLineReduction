/*
 *  process an SDFITS file
 *
 * there will be a subtle difference if the SDFITS file is from LMT or not
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <math.h>

#include "OTFParameters.h"
#include "SpecFile.h"
#include "fitsio.h"

typedef struct  {
  fitsfile *fptr;
  int nrows;
  int ncols;
  char **colnames;
  int *tdim;
} bintable, *bintableptr;

void iferror(fitsfile *fptr, int fstatus)
{
  if (fstatus) {
    printf("status=%d\n",fstatus);
    fits_close_file(fptr, &fstatus);
    fits_report_error(stderr, fstatus);
    exit(fstatus);
  }
}

char *fkeyword(char *card, char *keyword)
{
  
}

int keyindex(int ncols, char **colnames, char *keyword)
{
  for (int i=0; i<ncols; i++)
    if (strcmp(colnames[i],keyword)==0) return i;
  return -1;
}


void minmaxi(int n, int *data, int *data_min, int *data_max)
{
  *data_min = *data_max = data[0];
  for (int i=1; i<n; i++) {
    if (data[i] < *data_min) *data_min = data[i];
    if (data[i] > *data_max) *data_max = data[i];
  }
}


/*
 *   double *restfreq_data = get_column_dbl(ftp, "RESTFREQ", nrows, ncols, colnames);
 */


double *get_column_dbl(fitsfile *fptr, char *colname, int nrows, int ncols, char **colnames)
{
  int col = keyindex(ncols, colnames, colname) + 1;
  if (col < 1) return NULL;
  double *data = (double *) malloc(nrows * sizeof(double));
  double nulval = 0.0;
  int anynul;
  int fstatus = 0;
  fits_read_col(fptr, TDOUBLE, col, 1, 1, nrows, &nulval, data, &anynul, &fstatus);
  if (anynul) printf("get_column %s -> %d null's\n", colname, anynul);  
  return data;
}

int *get_column_int(fitsfile *fptr, char *colname, int nrows, int ncols, char **colnames)
{
  int col = keyindex(ncols, colnames, colname) + 1;
  if (col < 1) return NULL;
  int *data = (int *) malloc(nrows * sizeof(int));
  int nulval = 0;
  int anynul;
  int fstatus = 0;
  fits_read_col(fptr, TINT, col, 1, 1, nrows, &nulval, data, &anynul, &fstatus);
  if (anynul) printf("get_column %s -> %d null's\n", colname, anynul);
  return data;
}

char **get_column_str(fitsfile *fptr, char *colname, int nrows, int ncols, char **colnames)
{
  char keyword[FLEN_KEYWORD], data_fmt[FLEN_VALUE];
  int fstatus = 0;
  int col = keyindex(ncols, colnames, colname) + 1;
  if (col < 1) return NULL;
  char **data = (char **) malloc(nrows * sizeof(char *));
  // find the string length of this keyword
  fits_make_keyn("TFORM", col, keyword, &fstatus);
  fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &fstatus);	
  int slen = atoi(data_fmt);  // 1 extra for terminating 0
  printf("DATA in column %d  slen=%d\n",col,slen);
  // allocate the full block of chars for slen*nrows 
  char *vals = (char *) malloc(nrows*(slen+1)*sizeof(char));
  for (int i=0; i<nrows; i++)
    data[i] = &vals[(slen+1)*i];
  char *nulval = "\0";
  int anynul;
  fits_read_col_str(fptr,       col, 1, 1, nrows, nulval, data, &anynul, &fstatus);  
  if (anynul) printf("get_column %s -> %d null's\n", colname, anynul);  
  return data;
}


int read_sdfits_file(SpecFile *S, char *filename)
{
  int i,j,ii;
  int retval;
  int nchan;
  size_t nspec;
  int nspec_id, nchan_id, data_id, x_id, y_id, pix_id, seq_id, rms_id;
  char version[20];
  char history[512];
  fitsfile *fptr;
  int fstatus = 0;
  int fmode = READONLY;
  size_t nrows;
  int ncols;
  int nkeys;
  char keyword[FLEN_KEYWORD], colname[FLEN_VALUE], data_fmt[FLEN_VALUE], card[FLEN_VALUE];
  char date_obs[FLEN_KEYWORD];
  char date_obs2[FLEN_KEYWORD];  
  char telescop[FLEN_KEYWORD];
  char instrume[FLEN_KEYWORD];

  printf("WARNING:   the SDFITS reader is still experimental\n");

  fits_open_file(&fptr, filename, fmode, &fstatus);
  fits_get_hdrspace(fptr, &nkeys, NULL, &fstatus);
  printf("nkeys: %d\n", nkeys);
  for (ii = 1; ii <= nkeys; ii++)  { 
    fits_read_record(fptr, ii, card, &fstatus);
    printf("%s\n", card);
  }
  // DATE-OBS should really be a column
  fits_read_keyword(fptr,"DATE-OBS", date_obs, NULL, &fstatus);
  if (fstatus) {
    strcpy(date_obs,"");
    fstatus = 0;
  }
  fits_read_keyword(fptr,"TELESCOP", telescop, NULL, &fstatus);
  fits_read_keyword(fptr,"INSTRUME", instrume, NULL, &fstatus);
  printf("DATE-OBS: %s\n", date_obs);
  iferror(fptr,fstatus);  
  fits_close_file(fptr, &fstatus);
  iferror(fptr,fstatus);

  // SDFITVER= 'sdfits ver1.2'      / this file was created by sdfits                
  // FITSVER = '1.2     '           / FITS definition version                       

  // EXTNAME = 'SINGLE DISH'
  
  printf("Opening SDfits file %s\n",filename);
  fits_open_table(&fptr, filename, fmode, &fstatus);
  iferror(fptr,fstatus);
  fits_get_hdrspace(fptr, &nkeys, NULL, &fstatus);
  printf("nkeys: %d\n", nkeys);
  fits_get_num_rows(fptr, &nrows, &fstatus);
  iferror(fptr,fstatus);
  fits_get_num_cols(fptr, &ncols, &fstatus);
  iferror(fptr,fstatus);
  printf("%s : Nrows: %ld   Ncols: %d\n",filename,nrows,ncols);


  // gather all column names, but collect which column "DATA" or "SPECTRUM" is
  // also collect 'nchan' from the dimension of the DATA
  // In a BINTABLE each column name is under the FITS keyword TTYPEn
  char **colnames = (char **) malloc(ncols * sizeof(char *));
  int col_data = -1;

  for (i=0; i<ncols; i++) {    // loop over all columns to find the DATA column, awkward
    ii = i + 1;
    fits_make_keyn("TTYPE", ii, keyword, &fstatus);
    fits_read_key(fptr, TSTRING, keyword, colname, NULL, &fstatus);
    colnames[i] = strdup(colname);
    if (strcmp(colname,"DATA")==0 || strcmp(colname,"SPECTRUM")==0) {  // classic SDFITS or CLASS FITS
      col_data = ii;
      fits_make_keyn("TFORM", ii, keyword, &fstatus);
      fits_read_key(fptr, TSTRING, keyword, data_fmt, NULL, &fstatus);	
      nchan = atoi(data_fmt);
      printf("DATA in column %d  nchan=%d\n",col_data,nchan);
    }
  } //for(i)

  int col_crval2 =  keyindex(ncols, colnames, "CRVAL2") + 1;
  int col_crval3 =  keyindex(ncols, colnames, "CRVAL3") + 1;
  int col_tcal   =  keyindex(ncols, colnames, "TCAL") + 1;

  int col_crval1 =  keyindex(ncols, colnames, "CRVAL1") + 1;
  int col_crpix1 =  keyindex(ncols, colnames, "CRPIX1") + 1;
  int col_cdelt1 =  keyindex(ncols, colnames, "CDELT1") + 1;
  int col_ctype1 =  keyindex(ncols, colnames, "CTYPE1") + 1;
  printf("crval1: %d %d %d %d\n", col_crval1, col_crpix1, col_cdelt1, col_ctype1);
  int col_date_obs = keyindex(ncols, colnames, "DATE-OBS") + 1;
  printf("date-obs: %d\n", col_date_obs);

  if (col_tcal < 0)
    printf("Warning: no TCAL, cannot calibrate, assumed all are 1.0\n");
  printf("CRVAL2,3: cols=%d %d\n",col_crval2, col_crval3);


  nspec = nrows;

  S->theData = (float *)malloc(nspec*nchan*sizeof(float));      // the big chunk N-dim waterfall if you wish
  S->XPos = (float *)malloc(nspec*sizeof(float));
  S->YPos = (float *)malloc(nspec*sizeof(float));
  S->RMS = (float *)malloc(nspec*sizeof(float));
  S->Pixel = (int *)malloc(nspec*sizeof(int));
  S->Sequence = (int *)malloc(nspec*sizeof(int));
  S->RMS_cut = (float *)malloc(MAXPIXEL*sizeof(float));   // MAXPIXEL
  S->use = (int *)malloc(nspec*sizeof(int));

  // spectral axis
  double *crval1_data   = get_column_dbl(fptr, "CRVAL1", nrows, ncols, colnames);    
  double *crpix1_data   = get_column_dbl(fptr, "CRPIX1", nrows, ncols, colnames);    
  double *cdelt1_data   = get_column_dbl(fptr, "CDELT1", nrows, ncols, colnames);
  double *restfreq_data = get_column_dbl(fptr, "RESTFREQ", nrows, ncols, colnames);  
  printf("Spectral Axis: %g @ %g %g\n", crval1_data[0], crpix1_data[0], cdelt1_data[0]);
  printf("RESTFREQ: %.8f GHz\n", restfreq_data[0]/1e9);
  
  // RA,DEC (or whatever spatial system)
  double *crval2_data = get_column_dbl(fptr, "CRVAL2", nrows, ncols, colnames);  
  double *crval3_data = get_column_dbl(fptr, "CRVAL3", nrows, ncols, colnames);  
  double crval2_center=0.0, crval3_center=0.0;
  


#if 1
  int anynul = 0;
  char nulvalc = '\0';
  int  nulvali = 0;
  float nulvalf = 0.0;
  double nulvald = 0.0;
  char *nulvals = "\0";

  int col_cal  =    keyindex(ncols, colnames, "CAL") + 1;
  int col_sig  =    keyindex(ncols, colnames, "SIG") + 1;
  printf("tcal: %d   cal: %d sig: %d \n", col_tcal, col_cal, col_sig);

  char *sig_data = (char *) malloc(nspec*sizeof(char));
  fits_read_col(fptr, TBYTE, col_sig, 1, 1, nrows, &nulvalc, sig_data, &anynul, &fstatus);
  iferror(fptr,fstatus);

  char *cal_data = (char *) malloc(nspec*sizeof(char));
  fits_read_col(fptr, TBYTE, col_cal, 1, 1, nrows, &nulvalc, cal_data, &anynul, &fstatus);
  iferror(fptr,fstatus);
#endif

  int *fdnum_data = get_column_int(fptr, "FDNUM", nrows, ncols, colnames);
  int *ifnum_data = get_column_int(fptr, "IFNUM", nrows, ncols, colnames);
  int *plnum_data = get_column_int(fptr, "PLNUM", nrows, ncols, colnames);

  int fd_min, fd_max, if_min, if_max, pl_min, pl_max;
  minmaxi(nrows, fdnum_data, &fd_min, &fd_max);
  minmaxi(nrows, ifnum_data, &if_min, &if_max);
  minmaxi(nrows, plnum_data, &pl_min, &pl_max);
  printf("FDNUM: %d %d\n", fd_min, fd_max);
  printf("IFNUM: %d %d\n", if_min, if_max);
  printf("PLNUM: %d %d\n", pl_min, pl_max);

  for (i=0; i<nrows; i++) {
      printf("%d  %c %c  %d %d %d  %.8g \n", i,
             sig_data[i],  cal_data[i],  
             fdnum_data[i], ifnum_data[i], plnum_data[i], crval1_data[i]/1e9);
  }
  

  // warning: only one IFNUM needs to be selected
  //          two PLNUM's can be combined, if they are the summation type
  //          different FDNUM can be selected


#if 1
  char **date_obs_data = get_column_str(fptr, "DATE-OBS", nrows, ncols, colnames);
#else  
  char **date_obs_data = (char **) malloc(nspec*sizeof(char **));
  date_obs_data[0] = date_obs;
  date_obs_data[1] = &date_obs[23];   //  with 22 it fails on 2nd one, 23 is ok
  //fits_read_col(fptr, TSTRING, col_date_obs, 1, 1, 1, &nulvals, &date_obs, &anynul, &fstatus);
  fits_read_col_str(fptr, col_date_obs, 1, 1, 2, nulvals, date_obs_data, &anynul, &fstatus);
#endif  
  printf("DATE-OBS2: '%s' '%s'\n",date_obs_data[0],date_obs_data[1]);

  // these data should be read in sections if plnum/fdnum/ifnum dictate
  fits_read_col(fptr, TFLOAT,  col_data,   1, 1, nchan*nrows, &nulvalf, S->theData, &anynul, &fstatus);




  // Fill Pixel and Sequence with fake values for now
  // also calculate average of grid position @todo find alternative ways
  for (i=0; i<nrows; i++) {
    S->Pixel[i] = 0;
    S->use[i] = 1;
    S->Sequence[i] = i;
    S->RMS[i] = 1.0;
    crval2_center += crval2_data[i];
    crval3_center += crval3_data[i];
  }
  crval2_center /= nrows;
  crval3_center /= nrows;
  printf("Guessed center: ra=%g dec=%g\n",crval2_center, crval3_center);
  for (i=0; i<nrows; i++) {
    //S->XPos[i] = (crval2_center - crval2_data[i]) * cos(crval2_data[i] / 57.29577951308232087679815) * 3600.0;
    S->XPos[i] = (crval2_data[i] - crval2_center) * cos(crval3_data[i] / 57.29577951308232087679815) * 3600.0;
    S->YPos[i] = (crval3_data[i] - crval3_center) * 3600.0;
  }
  free(crval2_data);
  free(crval3_data);


  printf("file: %s nspec= %ld nchan= %d\n",filename,nspec,nchan);
  S->nspec = nspec;
  S->nchan = nchan;

  S->obsnum = -1;   // OBSID   = 'Skynet_62079'     
  S->x_position = crval2_center;
  S->y_position = crval3_center;

  //S->restfreq = 1420405751.786;
  S->vlsr = 0;
  strcpy(S->source, "LMT_TEST");
  strcpy(S->date_obs, date_obs);
  // S->version =
  // S->history
  
#if 0
  S->CRVAL = crval1_data[0];
  S->CRPIX = crpix1_data[0];
  S->CDELT = cdelt1_data[0];
  strcpy(S->CTYPE, "FREQ-OBS");   // for GBO 20m and 100m 
#else
  S->CRVAL = -700.0;  // km/s
  S->CRPIX = 1.0;
  S->CDELT = 3.0;
  strcpy(S->CTYPE, "VRAD");
#endif  
  
  // For SDFITS there is no map center, thus all RA,DEC positions are averaged and taken as center
  // After this all RA,DEC need to be referenced to this and stored in the S->XPos and S->YPos
  


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
  printf("X-range: %g %g   Y-range: %g %g arcsec\n",xmin,xmax,ymin,ymax);
  printf("X-ramp:  %g %g   Y-ramp:  %g %g arcsec\n",xmin+3*bs,xmax-3*bs,ymin+3*bs,ymax-3*bs);
  printf("MapSize: %g x %g arcsec\n", xmax-xmin, ymax-ymin);


  printf("SDfits file has been read\n");
  return 0;
}


