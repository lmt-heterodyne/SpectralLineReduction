#!/usr/bin/env python

"""Usage: process_otf_map2.py -p PATH -O OBSNUM -o OUTPUT [options]

-p PATH --path PATH                Path where ifproc and spectrometer/roach* files are
-o OUTPUT --output OUTPUT          Output SpecFile  [test.nc]
-O OBSNUM --obsnum OBSNUM          The obsnum, something like 79448. 
-b BANK --bank BANK                Spectral Bank for processing [default: 0]
--pix_list PIX_LIST                Comma separated list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
--eliminate_list ELIMINATE_LIST    Comma separated list of channels to be blanked
--use_cal                          Use Calibration scan
--tsys TSYS                        If use_cal is False, value of Tsys to use [default: 250.0] ** not used **
--use_otf_cal                      Use calibration within OTF scan (default: False)
--save_tsys                        Should tsys (from CAL) be saved in specfile?
--stype STYPE                      type of spectral line reduction;
                                   0 - median; 1 - single ref spectra; 2 - bracketed ref [Default: 2]
--map_coord MAP_COORD              Override map_coord for output coordinate system. [Default: -1]
                                   -1 - default; 0 - Az/El;  1 - Ra/Dec; 2 - L/B
--x_axis X_AXIS                    select spectral x axis.
                                   options one of VLSR, VSKY, VBARY, VSRC, FLSR, FSKY, FBARY, FSRC [default: VLSR]
--b_order B_ORDER                  set polynomial baseline order [default: 0]
--b_regions B_REGIONS              enter list of lists for baseline regions (default: [[],[]])
--l_regions L_REGIONS              enter list of lists for line fit regions (default: [[],[]])
--slice SLICE                      enter list to specify slice from spectrum for processing
--sample PIXEL,S0,S1               Series of sample sections per pixel to be removed from SpecFile (not implemented yet)
--restfreq RESTFREQ                Override the rest frequency (in GHz) for this bank. [Default: -1]

-h --help                          show this help

Creates a SpecFile from a single OTF mapping observation

Not a few options are still listed here, but don't seem to be doing
anything, as underlying code has changed.

We also list the --sample keyword, which is a proposed way to cull
time-based sections of a selected pixel. Use with care, as the sample
counter is based on the original RAW ON data. Note that the max number
of samples can differ per pixel due to their connection to the roach
board.

"""

# Python Imports	
import sys
import os
import numpy as np		
import matplotlib.pyplot as pl
import netCDF4
#from pypapi import events, papi_high as high

# command line parsing
from docopt import docopt
import lmtslr.utils.convert as acv


# Line Data Reduction Imports
from lmtslr.spec.specfile import SpecFile

#from lmtslr.spec.spec import *
#from lmtslr.ifproc.ifproc import *

#from lmtslr.reduction.line_reduction import *
#from lmtslr.grid.grid import *

from lmtslr.utils.reader import read_obsnum_otf #, count_otf_spectra
#from lmtslr.utils.parser import HandleProcessOptions
from lmtslr.utils.argparser import HandleOTFProcessOptions

import time

def main(argv):
    av = docopt(__doc__,options_first=True, version='0.2')
    print(av)   # debug
    # vslice = acv.listf(av['--slice'], 2)

    # set the command line for HISTORY in SpecFile and beyond
    history = sys.argv[0].split('/')[-1]
    for arg in sys.argv[1:]:
        history = history + " " + arg
        
    Opts = HandleOTFProcessOptions()
    Opts.parse_options(argv, 'process_otf_map', 1, True)
    save_tsys = True     # until it's a real option
    
    # check to see whether output file exists and remove it if it does
    if os.path.isfile(Opts.output_file_name) == True:
        os.remove(Opts.output_file_name) 

    I, S = read_obsnum_otf(Opts.obsnum,
                           Opts.pix_list,
                           Opts.bank,
                           Opts.use_cal,
                           tsys=Opts.tsys,
                           stype=Opts.stype,
                           use_otf_cal=Opts.use_otf_cal,
                           save_tsys=save_tsys,
                           restfreq=Opts.restfreq,
                           map_coord=Opts.map_coord,
                           path=Opts.data_path)

    specfile = SpecFile(I, S, Opts.pix_list)
    specfile.set_history(history)
    specfile.set_line_parameters(vslice=[Opts.slice[0], Opts.slice[1]],
                                 b_order=Opts.b_order,
                                 b_regions=Opts.b_regions,
                                 l_regions=Opts.l_regions,
                                 eliminate_list=Opts.eliminate_list)
    specfile.open_output_netcdf(Opts.output_file_name)
    specfile.write_ncdata()
    
    print('Written %s'% Opts.output_file_name)

if __name__ == '__main__':
    main(sys.argv[1:])

