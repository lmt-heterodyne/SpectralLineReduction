# Gridding Programs

## lmtgridder

grid_data.py is the python front-end that calls the gridder.  We have two gridders:

1. spec_driver_fits  - reads SpecFile (netCDF) format
2. lmtgridder - reads SDFITS, and optionally the SpecFile


   grid_data.py               G   details
   ------------              ---  -------
p PP --program_path PP        -   Executable [Default: spec_driver_fits]
-i INPUT --input INPUT        i   Input SpecFile (no default)
-o OUTPUT --output OUTPUT     o   Output map (no default)
-w WEIGHT --weight WEIGHT     w   Output weight map (no default)
--resolution RESOLUTION       f   Resolution in arcsec [Default: 14]
--cell CELL                   c   Cell size in arcsec [Default: 7]
--pix_list PIX_LIST           u   Comma separated list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
--rms_cut RMS_CUT             z   RMS threshold for data, negative allowed for robust MAD method  [Default: 10.0]
--noise_sigma NOISE_SIGMA     s   noise weighting - apply if > 0 [default: 1]
--x_extent X_EXTENT           x   x extent of cube (arcsec) note: cube will go to +/- x_extent [Default: 400]
--y_extent Y_EXTENT           y   y extent of cube (arcsec) note: cube will go to +/- y_extent [Default: 400]
--otf_select OTF_SELECT       f   otf filter code one of (0=box,1=jinc,2=gaussian,3=triangle) [default: 1]
--rmax RMAX                   r   maximum radius of convolution (units lambda/D) [default: 3.0]
--n_samples N_SAMPLES         n   number of samples in convolution filter [default: 256]
--otf_a OTF_A                 1   OTF A parameter [default: 1.1]
--otf_b OTF_B                 2   OTF B parameter [default: 4.75]
--otf_c OTF_C                 3   OTF C parameter [default: 2.0]
--sample P,S0,S1,P,...        b   Blank sample S0 to S1 for pixel P, etc. [Default: -1,0,0]

-h --help                     h   show some help



## gbtgridder

For good measure, here are the options for gbtgridder (development branch)

-c  --channels         Optional channel range to use.  '<start>:<end>' counting from 0.

-a  --average

                       Optionally average channels, keeping only number of channels/naverage channels

-s  --scans            Only use data from these scans.  comma separated list or <start>:<end> range syntax or combination of both


-m  --maxtsys          max Tsys value to use"


-z  --mintsys          min Tsys value to use

    --clobber          Overwrites existing output files if set

-k  --kernel           gridding kernel, default is gauss. options are: gauss, gaussbessel, nearest

    --diameter         Diameter of the telescope the observations were taken on [100]

-o  --output           root output name, instead of source and rest frequency

    --mapcenter        Map center in longitude and latitude of coordinate type used in data (RA/DEC, Galactic, etc) (degrees)",
    
    --size             Image X,Y size (pixels)

    --pixelwidth       Image pixel width on sky (arcsec)
	
-p  --proj             Projection to use for the spatial axes, default is SFL (choises are SFL and TAN)

    --clonecube        A FITS cube to use to set the image size and WCS parameters
    		       in the spatial dimensions.  The cube must have the same axes
        	       produced here, the spatial axes must be of the same type as
        	       found in the data to be gridded, and the projection used in the
        	       cube must be either TAN, SFL, or GLS [which is equivalent to SFL].
        	       Default is to construct the output cube using values appropriate for
        	       gridding all of the input data.  Use of --clonecube overrides any use
        	       of --size, --pixelwidth, --mapcenter and --proj arguments.

    --autoConfirm      Set this to True if you'd like to auto-confirm the program stop and move straight into gridding.
                       Default is False.

    --nocov            Set this to turn off production of the output coverage map
                       Default is False.
		       
    --noweight         Set this to turn off production of the output weight cube
                       Default is False.

    --noline           Set this to turn off prodution of the output line cube
                       Default is False.

    --nocont           Set this to turn off prodution of the output 'cont' image
                       Default is False.

-v  --verbose          set the verbosity level-- 0-1:none,2:errors only, 3:+warnings, 4(default):+user info, 5:+debug"
	
-V  --version          show gbtgridder version

## gbtgridder - old python2 version


      gbtgridder.py [optional arguments] SDFITSfiles [SDFITSfiles ...]


optional arguments:

  -h, --help            show this help message and exit
  -c CHANNELS, --channels CHANNELS
                        Optional channel range to use. '<start>:<end>'
                        counting from 0.
  -a AVERAGE, --average AVERAGE
                        Optionally average channels, keeping only
                        nchan/naverage channels
  -s SCANS, --scans SCANS
                        Only use data from these scans. comma separated list
                        or <start>:<end> range syntax or combination of both
  -m MAXTSYS, --maxtsys MAXTSYS
                        max Tsys value to use

  -z MINTSYS, --mintsys MINTSYS
                        min Tsys value to use
  --clobber             Overwrites existing output files if set.
  -k {gauss,gaussbessel,nearest}, --kernel {gauss,gaussbessel,nearest}
                        gridding kernel, default is gauss
  -o OUTPUT, --output OUTPUT
                        root output name, instead of source and rest frequency
  --mapcenter LONG LAT  Map center in longitude and latitude of coordinate
                        type used in data (RA/DEC, Galactic, etc) (degrees)
  --size X Y            Image X,Y size (pixels)
  --dish DISH           Effective dish size (m)
  --pixelwidth PIXELWIDTH
                        Image pixel width on sky (arcsec)
  --restfreq RESTFREQ   Rest frequency (MHz)
  -p {SFL,TAN}, --proj {SFL,TAN}
                        Projection to use for the spatial axes, default is SFL
  --clonecube CLONECUBE

                        A FITS cube to use to set the image size and WCS
                        parameters in the spatial dimensions. The cube must
                        have the same axes produced here, the spatial axes
                        must be of the same type as found in the data to be
                        gridded, and the projection used in the cube must be
                        either TAN, SFL, or GLS [which is equivalent to SFL].
                        Default is to construct the output cube using values
                        appropriate for gridding all of the input data. Use of
                        --clonecube overrides any use of --size, --pixelwidth,
                        --mapcenter and --proj arguments.
  --eqweight            Set this to use equal data weights for all spectra
  --noweight            Set this to turn off production of the output weight
                        cube
  --noline              Set this to turn off prodution of the output line cube
  --nocont              Set this to turn off prodution of the output 'cont'
                        image
  -v VERBOSE, --verbose VERBOSE
                        set the verbosity level-- 0-1:none, 2:errors only,
                        3:+warnings, 4(default):+user info, 5:+debug
  -V, --version         show program's version number and exit



## todo

there will be more selections needed

band:      select one band, or allow band-merge
pol:       select one pol, or allow averaging, or diff to see problems?
beam:      presumably for OTF, but the -u flag already selects them.
