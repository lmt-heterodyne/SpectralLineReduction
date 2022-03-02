#!/usr/bin/env python

"""Usage: grid_data.py  -i INPUT -o OUTPUT -w WEIGHT [options]

-p PP --program_path PP       Executable [Default: spec_driver_fits]
-i INPUT --input INPUT        Input SpecFile (no default)
-o OUTPUT --output OUTPUT     Output map (no default)
-w WEIGHT --weight WEIGHT     Output weight map (no default)
--resolution RESOLUTION       Resolution in arcsec [Default: 14]
--cell CELL                   Cell size in arcsec [Default: 7]
--pix_list PIX_LIST           Comma separated list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
--rms_cut RMS_CUT             RMS threshold for data, negative allowed for robust MAD method  [Default: 10.0]
--noise_sigma NOISE_SIGMA     noise weighting - apply if > 0 [default: 1]
--x_extent X_EXTENT           x extent of cube (arcsec) note: cube will go to +/- x_extent [Default: 400]
--y_extent Y_EXTENT           y extent of cube (arcsec) note: cube will go to +/- y_extent [Default: 400]
--otf_select OTF_SELECT       otf filter code one of (0=box,1=jinc,2=gaussian,3=triangle) [default: 1]
--rmax RMAX                   maximum radius of convolution (units lambda/D) [default: 3.0]
--n_samples N_SAMPLES         number of samples in convolution filter [default: 256]
--otf_a OTF_A                 OTF A parameter [default: 1.1]
--otf_b OTF_B                 OTF B parameter [default: 4.75]
--otf_c OTF_C                 OTF C parameter [default: 2.0]
--sample P,S0,S1,P,...        Blank sample S0 to S1 for pixel P, etc. [Default: -1,0,0]
--edge EDGE                   Fuzzy edge?  [default: 1]

-h --help                     show this help

 
GRID_DATA Reads one or more SpecFiles and creates a data cube and
weight map in FITS format

Currently the weight map designates how often a pixel has been
weighted, but is not weighted by the RMS of the pixel yet. We will get
another option for this.

Another option in the future will be the treatment of the edge, or
more generally, cells with no pixels positions. This
parameter would designate how many neighboring cells are needed with
pixels to allow interpolation or extrapolation, the latter being
the controversial one.  A value of 4 is suggested

"""

# Python Imports	
import numpy as np		
import matplotlib.pyplot as pl
import subprocess		 
import netCDF4

# command line parsing
from docopt import docopt
import lmtslr.utils.convert as acv


# Line Data Reduction Imports
from lmtslr.spec.spec import *
from lmtslr.ifproc.ifproc import *

from lmtslr.reduction.line_reduction import *

#from lmtslr.utils.parser import HandleGridOptions
from lmtslr.utils.argparser import HandleGridOptions


def main(argv):
    av = docopt(__doc__,options_first=True, version='0.2')
    print(av)   # debug

    # check to see whether output files exists and remove it if it does
    output_file_name = av['--output']
    weight_file_name = av['--weight']
    for f in [output_file_name, weight_file_name]:
        if os.path.isfile(f) == True:
            print("Removing ",f)
            os.remove(f)
    # pix_list is a different beast.  For now, we cheat and add the odd looking [] list
    # but this should be a @todo, why do we need these?? 
    pix_list = '[' + av['--pix_list'] + ']'

    with open('out.txt','w+') as outputfile:
        with open('err.txt','w+') as errorfile:

            #exit_code = subprocess.call(['./test_otf','-i',Opts.input_file_name,'-u',Opts.pix_list],stdout=outputfile,stderr=errorfile)
            exit_code=subprocess.call([av['--program_path'],
                                       '-i',av['--input'],
                                       '-o',output_file_name,
                                       '-w',weight_file_name,
                                       '-l',av['--resolution'],
                                       '-c',av['--cell'],
                                       '-u',pix_list,
                                       '-z',av['--rms_cut'],
                                       '-s',av['--noise_sigma'],
                                       '-x',av['--x_extent'],
                                       '-y',av['--y_extent'],
                                       '-f',av['--otf_select'],
                                       '-r',av['--rmax'],
                                       '-n',av['--n_samples'],
                                       '-0',av['--otf_a'],
                                       '-1',av['--otf_b'],
                                       '-2',av['--otf_c'],
                                       '-b',av['--sample'],
                                       '-e',av['--edge'],
                                       ], 
                                      stdout=outputfile,
                                      stderr=errorfile)
            
            # reset stdout file to read from it
            outputfile.seek(0)
            # save output (if any) in variable
            standard_output=outputfile.read()
            print('STDOUT ****************************')
            print(standard_output)

            # reset stderr file to read from it
            errorfile.seek(0) 
            # save errors (if any) in variable
            standard_error = errorfile.read()
            print('STDERR ****************************')
            print(standard_error)

    print('Exit Code: %d'%(exit_code))


if __name__ == '__main__':
    main(sys.argv[1:])

