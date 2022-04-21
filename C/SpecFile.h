#ifndef _SPECFILE_H_
#define _SPECFILE_H_

#include <netcdf.h>

/** error handling for netcdf reading and writing
 */

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

typedef struct
{
  int nspec;                           // number of spectra (see sequence)
  int nchan;                           // number of channels
  int npix;                            // number of pixels (fixed at 16)
  int ncal;                            // number of Tsys cals
  // data from obs header
  int obsnum;
  int map_coord;                       // 0=azel 1=radec 2=latlon
  char source[40];                     // what is enough?
  double x_position, y_position;       // degrees ?
  double restfreq;                     //
  float vlsr;                          // km/s
  float zsource;
  char date_obs[40];                   // e.g. 2021-09-14T07:25:23.370
  char receiver[40];                   // SEQ,OMA,1MM
  // axis parameters from line header
  double CRPIX, CRVAL, CDELT;          // km/s
  double *CAXIS;                       // km/s
  char CTYPE[20];
  // the data from the file
  float *theData;                      // data[nspec][nchan]
  int *Pixel;                          // 0,1,2,...
  float *RMS_cut;  
  int *Sequence;
  int *use;
  float *XPos;                         // arcsec
  float *YPos;                         // arcsec
  float *RMS;                          // K
  float deltaf;                        // deltaFreq * deltaTime for radiometer equation
  float deltat;                        // deltaFreq * deltaTime for radiometer equation
  float *Tsys;                         // tsys[ncal][npix][nchan]   unused pixels have Tsys=0
  char *telescop;
  char *instrume;
  double *Date;        // not used yet
  char history[512];   // cmdline how the specfile was made ( > 10-jan-2021)
} SpecFile;
  

int read_spec_file(SpecFile *S, char *filename);
void make_spec_beam(SpecFile *S);
void free_spec_file(SpecFile *S);
float *get_spectrum(SpecFile *S, int i);
float get_tsys(SpecFile *S, int i, int j);

#endif
