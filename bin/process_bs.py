#!/usr/bin/env python
'''
Processes a set of BS spectra
'''

# Python Imports
import numpy as np
import matplotlib.pyplot as pl
import sys

# Line Data Reduction Imports
from lmtslr.spec.spec import *
from lmtslr.utils.roach_file_utils import *
from lmtslr.ifproc.ifproc import *
from lmtslr.utils.ifproc_file_utils import *
from lmtslr.viewer.spec_viewer import *

from lmtslr.reduction.line_reduction import *

from lmtslr.utils.reader import read_obsnum_bs
#from lmtslr.utils.parser import HandleProcessOptions
from lmtslr.utils.argparser import HandlePSProcessOptions

def main(argv):
    #  Header.Bs.Beam = 10, 8 ;
    Opts = HandlePSProcessOptions()
    result = Opts.parse_options(argv, 'process_bs', 1, True)
    print("D:",Opts.data_path)
    #if result == 0:
    # this will be a list of processed spectral lines
    LineList = []
    # here are the number of spectra to be reduced
    nscans = len(Opts.obs_list)
    npixels = len(Opts.pix_list)

    for obs in Opts.obs_list:

        I,S = read_obsnum_bs(obs,
                             Opts.pix_list,
                             Opts.bank,
                             Opts.use_cal,
                             tsys=Opts.tsys,
                             path=Opts.data_path)

        # create a LineData object for each "pixel" in pix_list
        # we could do processing here...
        for i in range(npixels):
            LD = LineData(I,Opts.bank,S.nchan,S.bandwidth,S.roach[i].ps_spectrum)
            LD.set_x_axis(Opts.x_axis)
            LineList.append(LD)

    # show all the plots just to illustrate reduction
    # this will be replaced by write out of spectra to FITS file.
    if True:
        if len(LineList) == 2:
            pl.plot(LineList[0].xarray, LineList[1].yarray - LineList[0].yarray, label='Diff')
    for i in range(len(LineList)):
        pl.plot(LineList[i].xarray,LineList[i].yarray, label='%s' % Opts.pix_list[i])
        #pl.axis([-20,20,-1,1])
        pl.xlabel('VSRC')
    pl.legend()
    pl.show()


if __name__ == '__main__':
    main(sys.argv[1:])
