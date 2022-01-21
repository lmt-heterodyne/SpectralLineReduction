#ifndef _OTF_PARAMETERS_H_
#define _OTF_PARAMETERS_H_

#define OTF_FILE_NAME_LENGTH 200
#define N_INPUT_FILES 2000
typedef struct
{
  char i_filename[N_INPUT_FILES][OTF_FILE_NAME_LENGTH];
  int nfiles;
  char o_filename[OTF_FILE_NAME_LENGTH];
  float rms_cutoff;
  float noise_sigma;
  float resolution_size;
  float cell_size;
  int nsamples;
  int n_cell;
  int n_subcell;
  float model_spectrum_hpw, model_source_amp, model_source_hpw, model_source_x, model_source_y;
  int model_nchan;
  float otf_jinc_a, otf_jinc_b, otf_jinc_c;
  int otf_select;
  float rmax;
  float x_extent,y_extent;
  float sample_step;
  float scan_step;
  int nx_samples, ny_samples;
  int use_pixels[16];
} OTFParameters;

void initialize_otf_parameters(OTFParameters *OTF, int argc, char *argv[]);
char **str_split(char*, const char);
void decode_pix_list(OTFParameters*, char*);
void decode_file_list(OTFParameters*, char*);

#endif
