#!/usr/bin/env python
'''
Creates a SpecFile from a single OTF mapping observation
'''

# Python Imports	
import sys
import os
import numpy as np		
import matplotlib.pyplot as pl
import netCDF4			 

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
    print(time.time(), time.clock())
    Opts = HandleOTFProcessOptions()
    Opts.parse_options(argv, 'process_otf_map', 1, True)

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
                           path=Opts.data_path)

    specfile = SpecFile(I, S, Opts.pix_list)
    specfile.set_line_parameters(vslice=[Opts.slice[0], Opts.slice[1]],
                                 b_order=Opts.b_order,
                                 b_regions=Opts.b_regions,
                                 l_regions=Opts.l_regions,
                                 eliminate_list=Opts.eliminate_list)
    specfile.open_output_netcdf(Opts.output_file_name)
    specfile.write_ncdata()
    
    print('netCDF %s Done'%(Opts.output_file_name))
    print(time.time(), time.clock())

if __name__ == '__main__':
    main(sys.argv[1:])

