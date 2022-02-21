#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "OTFParameters.h"

char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

void decode_pix_list(OTFParameters *OTF, char *the_list)
{
    int pixel,i,len;
    char p[MAXSTRLEN];
    char** tokens;

    len = strlen(the_list);

    //printf("PJT pix_list=%s\n\n", the_list);

    if (the_list[0] == '[') {
      strncpy(p,&the_list[1],len-2);
      p[len-1] = '\0';
    } else
      strncpy(p,&the_list[0],len);
    tokens = str_split(p, ',');

    if (tokens)
    {
        for (i = 0; *(tokens + i); i++)
        {
            pixel = atoi(*(tokens + i));
            OTF->use_pixels[pixel] = 1;
            //printf("pixel=%d\n", pixel);
            free(*(tokens + i));
        }
	// printf("\n");
        free(tokens);
    }
#if 0    
    for(i=0;i<16;i++)
	printf("pixel %d use_it %d\n",i,OTF->use_pixels[i]);
#endif
}

void decode_sample_list(OTFParameters *OTF, char *the_list)
{
    int id,i,len,nsam;
    char p[MAXSTRLEN];
    char** tokens;
    
    len = strlen(the_list);
    
    printf("sample_list=%s\n\n", the_list);
    
    if (the_list[0] == '[') {
      strncpy(p,&the_list[1],len-2);
      p[len-1] = '\0';
    } else
      strncpy(p,&the_list[0],len);
    tokens = str_split(p, ',');
    nsam = 0;
    

    if (tokens)
    {
        for (i = 0; *(tokens + i); i++)
        {
            id = atoi(*(tokens + i));
            OTF->samples[nsam] = id;
	    nsam++;
            free(*(tokens + i));
        }
        free(tokens);
    }
    if (nsam%3 == 0)
      OTF->nsegment = nsam/3;
    else {
      OTF->nsegment = 0;
      printf("Warning: sample list is not a multiple of 3\n");
    }
#if 1    
    printf("nsamples=%d\n",nsam);
    for (i=0; i<nsam/3; i++)
      printf("Pixel %2d :  Masking %d - %d\n",OTF->samples[3*i], OTF->samples[3*i+1], OTF->samples[3*i+2]);
#endif    
}


void decode_file_list(OTFParameters *OTF, char *the_list)
{
    int i;
    char** tokens;

    printf("file list=%s\n", the_list);

    tokens = str_split(the_list, ',');

    OTF->nfiles = 0;
    if (tokens)
    {
        for (i = 0; *(tokens + i); i++)
        {
	  OTF->nfiles++;
          strcpy(OTF->i_filename[i],*(tokens + i));
	  free(*(tokens + i));
        }
        free(tokens);
    }
}


void initialize_otf_parameters(OTFParameters *OTF, int argc, char *argv[])
{
  int i;
  int coption;

  // debug: show the commandline
#if 1
  printf("CMD: ");
  for (i=0; i<argc; i++)
    printf("%s ", argv[i]);
  printf("\n");
#endif	   


  //  Default Parameters  (grid_data only exposes 15, with the '=' symbol)  
  //                                 // -i=
  strcpy(OTF->o_filename,"test.nc"); // -o=
  strcpy(OTF->w_filename,"");        // -w=   default is no weight file
  strcpy(OTF->a_filename,"");        // -w=   default is no model file  
  OTF->resolution_size    = 14.;     // -l=   would be better to use the default for 115 GHz???
  OTF->cell_size          = 7.;      // -c=
  for(i=0;i<MAXPIXEL;i++) 
    OTF->use_pixels[i]    = 0;       // -u=   what's the deal with or without []
  OTF->rms_cutoff         = 10.;     // -z=
  OTF->noise_sigma        = 0.0;     // -s=
  OTF->x_extent           = 160.0;   // -x=
  OTF->y_extent           = 160.0;   // -y=
  OTF->otf_select         = 1;       // -f=        1=jinc other=gauss
  OTF->rmax               = 3.0;     // -r=
  OTF->nsamples           = 256;               // error??? should this not be -n ???
  OTF->otf_jinc_a         = 1.1;     // -0=
  OTF->otf_jinc_b         = 4.75;    // -1=
  OTF->otf_jinc_c         = 2.0;     // -2=
  
  OTF->n_cell             = 5;       // -n=   // error???
  OTF->n_subcell          = 2;       // -m
#if 0  
  OTF->model_spectrum_hpw = 5.;
  OTF->model_source_amp   = 1.;
  OTF->model_source_hpw   = 16.1;
  OTF->model_source_x     = 0.;
  OTF->model_source_y     = 0.;
  OTF->model_nchan        = 64;
#endif

  OTF->sample_step        = 1.0;     // -p   ???     Xstep 
  OTF->scan_step          = 6.65;    // -q   ???     Ystep

  OTF->model              = 0;       //
  OTF->fuzzy_edge         = 0;       // -e
 
  // parse command line arguments
  opterr = 0;

  

  while(1)
    {
      
      static struct option long_options[] =
	{
	  {"help",             no_argument,       0, 'h'},
	  
	  {"input",            required_argument, 0, 'i'},  // 20 real options here
	  {"output",           required_argument, 0, 'o'},  //
	  {"weight",           optional_argument, 0, 'w'},  //
	  {"model",            optional_argument, 0, 'a'},  //
	  {"resolution_size",  required_argument, 0, 'l'},  // --resolution
	  {"cell_size",        required_argument, 0, 'c'},  // --cell
	  {"pix_list",         required_argument, 0, 'u'},  // --pix_list 
	  {"rms_cutoff",       required_argument, 0, 'z'},  // --rms_cut
	  {"noise_sigma",      required_argument, 0, 's'},  // --noise_sigma
	  {"x_extent",         required_argument, 0, 'x'},  // --x_extent
	  {"y_extent",         required_argument, 0, 'y'},  // --y_extent
	  {"filter",           required_argument, 0, 'f'},  // --otf_select
	  {"rmax",             required_argument, 0, 'r'},  // --rmax
	  {"n_cell",           required_argument, 0, 'n'},  // --n_samples      // wrong one?
	  {"jinc_a",           required_argument, 0, '0'},  // --otf_a
	  {"jinc_b",           required_argument, 0, '1'},  // --otf_b
	  {"jinc_c",           required_argument, 0, '2'},  // --otf_c
	  {"sample",           optional_argument, 0, 'b'},  // --sample
	  {"edge",             optional_argument, 0, 'e'},  // --edge

	  {"n_subcell",        required_argument, 0, 'm'},   // not passed
	  {"sample_step",      required_argument, 0, 'p'},   // not passed
	  {"scan_step",        required_argument, 0, 'q'},   // not passed
	  // undocumented -b flag
	  // undocumented -j flag
	  {0,0,0,0}
	};
      
      int option_index=0;
      //const char *optstring = "hi:o:bjl:z:c:n:f:m:s:r:0:1:2:x:y:p:q:u:";   // original w/ -b,-j
      static const char *optstring = "i:o:w:l:c:u:z:s:x:y:f:r:n:e:0:1:2:b:w:a:p:q:h";
      coption = getopt_long(argc, argv, optstring, long_options,&option_index);
      
      if(coption == -1)
	break;

      //printf("%c\n",coption);
      switch(coption)
	{
	case 'h':
	  for(i=0;i<19;i++)
	    printf("%c %s\n",long_options[i].val,long_options[i].name);
	  exit(0);
	case 'i':
	  decode_file_list(OTF, optarg);
	  //strcpy(OTF->i_filename,optarg);
	  break;
	case 'o':
	  strcpy(OTF->o_filename,optarg);
	  break;
	case 'w':
	  strcpy(OTF->w_filename,optarg);
	  break;
	case 'a':
	  strcpy(OTF->a_filename,optarg);
	  OTF->model = 1;	  
	  break;
	case 'l':
	  OTF->resolution_size = atof(optarg);
	  break;
	case 'c':
	  OTF->cell_size = atof(optarg);
	  printf("DEBUG:%g\n",	  OTF->cell_size);
	  break;
	case 'u':
	  decode_pix_list(OTF,optarg);
	  break;
	case 'b':
	  decode_sample_list(OTF,optarg);
	  break;
	case 'z':
	  OTF->rms_cutoff = atof(optarg);
	  break;
	case 's':
	  OTF->noise_sigma = atof(optarg);
	  break;
	case 'x':
	  OTF->y_extent = atof(optarg);       // X and Y axes are reversed
	  break;
	case 'y':
	  OTF->x_extent = atof(optarg);
	  break;
	case 'f':
	  OTF->otf_select = atoi(optarg);
	  break;
	case 'r':
	  OTF->rmax = atof(optarg);
	  break;
	case 'n':
#if 1
	  OTF->n_cell = atoi(optarg);        // original code
#else	  
  	  OTF->nsamples = atoi(optarg);      // PJT
#endif	  
	  //printf("PJT  n_cell=%d\n",OTF->n_cell);
	  //printf("PJT  nsample=%d\n",OTF->nsamples);	  
	  break;

	  
	case 'e':
	  OTF->fuzzy_edge= atoi(optarg);
	  break;
	case 'm':
	  OTF->n_subcell = atoi(optarg);
	  break;
	case '0':
	  OTF->otf_jinc_a = atof(optarg);
	  break;
	case '1':
	  OTF->otf_jinc_b = atof(optarg);
	  break;
	case '2':
	  OTF->otf_jinc_c = atof(optarg);
	  break;
	case 'p':
	  OTF->sample_step = atof(optarg);
	  break;
	case 'q':
	  OTF->scan_step = atof(optarg);
	  break;
	default:
	  printf("Command Line Argument Problem for %c\n",coption);
	  abort();
	}
    }
  // Derived Values
  OTF->nx_samples = (int)(OTF->x_extent/OTF->sample_step);
  OTF->ny_samples = (int)(OTF->y_extent/OTF->scan_step);
}
