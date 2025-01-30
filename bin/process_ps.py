#!/usr/bin/env python
'''
Processes a set of PS spectra
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

from lmtslr.utils.reader import read_obsnum_ps
#from lmtslr.utils.parser import HandleProcessOptions
from lmtslr.utils.argparser import HandlePSProcessOptions

def main(argv):
    Opts = HandlePSProcessOptions()
    result = Opts.parse_options(argv, 'process_ps', 1, True)
    #if result == 0:
    # this will be a list of processed spectral lines
    LineList = []
    # here are the number of spectra to be reduced
    nscans = len(Opts.obs_list)
    npixels = len(Opts.pix_list)

    import matplotlib
    if Opts.show:
        matplotlib.use('qt5agg')
    else:
        matplotlib.use('agg')
    import matplotlib.pyplot as plt
    print('mpl backend spectra',matplotlib.get_backend())

    # for 1MM:   0 for pix_list=0,2 and 1 for pix_list=1,3
    bank_hack = -1

    for obs in Opts.obs_list:

        I,S = read_obsnum_ps(obs,
                             Opts.pix_list,
                             Opts.bank,
                             Opts.use_cal,
                             tsys=Opts.tsys,
                             path=Opts.data_path)
        tsys = S.roach[0].tsys_spectrum
        print("PJT TSYS",tsys)
        #   set bank_hack on the first time
        if bank_hack < 0:
            if I.receiver == 'Msip1mm':
                if Opts.pix_list[0] == 0:
                    bank_hack = 0
                elif Opts.pix_list[0] == 1:
                    bank_hack = 1
                else:
                    print("Something odd with your 1MM pix_list",Opts.pix_list)
                    bank_hack = Opts.bank
            else:
                bank_hack = Opts.bank
            print("PJT:  bank_hack=:",bank_hack)
            
        # create a LineData object for each "pixel" in pix_list
        # we could do processing here...
        for i in range(npixels):
            LD = LineData(I,bank_hack,S.nchan,S.bandwidth,S.roach[i].ps_spectrum)
            LD.set_x_axis(Opts.x_axis)
            print("PJT",LD.xarray[0],LD.xarray[-1],S.nchan,S.bandwidth)
            LineList.append(LD)

    # prepare for subsetting the line and do a baseline subtraction
    #print("Line/Baseline:",Opts.b_order, Opts.b_regions, Opts.l_regions)
    Opts.b_order = 3
    print("Line/Baseline:",Opts.b_order,Opts.slice)

    # show all the plots just to illustrate reduction
    # this will be replaced by write out of spectra to FITS file.
    for i in range(len(LineList)):
        pl.plot(LineList[i].xarray,LineList[i].yarray)
        #pl.axis([-20,20,-1,1])
        #pl.xlabel('VSRC')
        pl.xlabel(Opts.x_axis)
    if Opts.show:
        pl.show()
    else:
        pout = 'ps.png'
        pl.savefig(pout)
        print("%s written" % pout)

    # show all the plots just to illustrate reduction
    # this will be replaced by write out of spectra to ASCII and FITS file.
    edge = 64
    if len(LineList) >  0:
        fp = open(Opts.output,"w")
        fp.write("# vlsr  TA(K) Tsys(K)     (average of %d beams  %d channels)\n" % (len(LineList), len(LineList[0])))
        for i in range(edge,len(LineList[0].xarray)-edge):
            ysum = 0.0
            for j in range(len(LineList)):
                ysum = ysum + LineList[j].yarray[i]
            ysum = ysum / len(LineList);
            fp.write("%g %g %g\n" %
                     (LineList[0].xarray[i],ysum,tsys[i]))
        fp.close()

if __name__ == '__main__':
    main(sys.argv[1:])
