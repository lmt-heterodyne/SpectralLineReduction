#ifndef _OTF_PARAMETERS_H_
#define _OTF_PARAMETERS_H_

#define OTF_FILE_NAME_LENGTH  200
#define N_INPUT_FILES        2000
#define MAXPIXEL               16
#define MAXSEGMENT            512
#define MAXSTRLEN            1024
typedef struct
{
  // also set via command line options
  char i_filename[N_INPUT_FILES][OTF_FILE_NAME_LENGTH];
  int nfiles;
  char o_filename[OTF_FILE_NAME_LENGTH];
  char w_filename[OTF_FILE_NAME_LENGTH];
  char a_filename[OTF_FILE_NAME_LENGTH];
  float resolution_size;
  float cell_size;
  float rms_cutoff;
  float noise_sigma;
  float x_extent,y_extent;
  int otf_select;
  float rmax;
  int nsamples;
  float otf_jinc_a, otf_jinc_b, otf_jinc_c;
  int model;  // set if a_filename was given
  int fuzzy_edge;

  // fixed in code
  int n_cell;
  int n_subcell;
  float model_spectrum_hpw, model_source_amp, model_source_hpw, model_source_x, model_source_y;
  int model_nchan;
  float sample_step;
  float scan_step;
  
  // derived in code
  int nx_samples, ny_samples;
  int use_pixels[MAXPIXEL];
  int nsegment;
  int samples[3*MAXSEGMENT];  // triples of (pixel,sample0,sample1)
} OTFParameters;

void initialize_otf_parameters(OTFParameters *OTF, int argc, char *argv[]);
char **str_split(char*, const char);
void decode_pix_list(OTFParameters*, char*);
void decode_file_list(OTFParameters*, char*);

#endif
