#!/usr/bin/env python
'''
Processes a set of BS spectra
'''
import sys
import numpy as np
import matplotlib.pyplot as pl

from lmtslr.reduction.line_reduction import LineData
from lmtslr.utils.reader import read_obsnum_bs
from lmtslr.utils.argparser import HandlePSProcessOptions

def main(argv):
    #  Header.Bs.Beam = 10, 8 ;
    Opts = HandlePSProcessOptions()
    result = Opts.parse_options(argv, 'process_bs', 1, True)
    print("D:",Opts.data_path)
    #if result == 0:
    # this will be a list of processed spectral lines
    LineList = []
    label    = []
    # here are the number of spectra to be reduced
    nscans = len(Opts.obs_list)
    try:
        pix_list = Opts.pix_list
    except:
        pix_list = None
    block = Opts.block
    edge = 64
    
    for obs in Opts.obs_list:

            if block == -2:
                # special loop, creates waterfall
                I,S = read_obsnum_bs(obs,
                                     pix_list,
                                     Opts.bank,
                                     Opts.use_cal,
                                     tsys=Opts.tsys,
                                     block=0,
                                     stype=Opts.stype,                                     
                                     path=Opts.data_path)
                nblocks = S.roach[0].nons
                npixels = len(S.roach)
                for i in range(npixels):
                    data = S.roach[i].ps_spectrum
                    LD = LineData(I,Opts.bank,S.nchan,S.bandwidth,S.roach[i].ps_spectrum)
                    LD.set_x_axis(Opts.x_axis)
                    LineList.append(LD)
                    label.append("%d/%d/%d" % (obs,Opts.pix_list[i],0))
                for iblock in range(1,nblocks):
                    I,S = read_obsnum_bs(obs,
                                         Opts.pix_list,
                                         Opts.bank,
                                         Opts.use_cal,
                                         tsys=Opts.tsys,
                                         block=iblock,
                                         stype=Opts.stype,
                                         path=Opts.data_path)
                    for i in range(npixels):
                        data = S.roach[i].ps_spectrum
                        LD = LineData(I,Opts.bank,S.nchan,S.bandwidth,S.roach[i].ps_spectrum)
                        LD.set_x_axis(Opts.x_axis)
                        LineList.append(LD)
                        label.append("%d/%d/%d" % (obs,Opts.pix_list[i],iblock))
                    
                
            else:
                # creates spectrum
                nblocks = 1
                I,S = read_obsnum_bs(obs,
                                     pix_list,
                                     Opts.bank,
                                     Opts.use_cal,
                                     tsys=Opts.tsys,
                                     block=block,
                                     stype=Opts.stype,                                     
                                     path=Opts.data_path)

                if pix_list is None:
                    pix_list = I.bs_beams
                print("NOTE: for obsnum=%d bs_beams=%s, you specified %s" % (obs,I.bs_beams,pix_list))
                print("Number of blocks: %d" % S.roach[0].nons)

                # create a LineData object for each "pixel" in pix_list
                # we could do processing here...
                npixels = len(S.roach)
                for i in range(npixels):
                    data = S.roach[i].ps_spectrum
                    LD = LineData(I,Opts.bank,S.nchan,S.bandwidth,S.roach[i].ps_spectrum)
                    LD.set_x_axis(Opts.x_axis)
                    LineList.append(LD)
                    label.append("%d/%d" % (obs,pix_list[i]))

    # show all the plots just to illustrate reduction
    # this will be replaced by write out of spectra to FITS file.
    if len(LineList) == 2:
        fp = open(Opts.output,"w")
        fp.write("# vlsr  TA(K)   Tsys1  Tsys2\n")
        for i in range(edge,len(LineList[0].xarray)-edge):
            fp.write("%g %g  %g %g\n" %
                     (LineList[0].xarray[i],
                      0.5*(LineList[0].yarray[i] - LineList[1].yarray[i]),    # @todo   should do Tsys weighted
                      S.roach[0].tsys_spectrum[i] if Opts.use_cal else S.roach[0].tsys_spectrum,
                      S.roach[1].tsys_spectrum[i] if Opts.use_cal else S.roach[1].tsys_spectrum))
        fp.close()
        # pl.plot(LineList[0].xarray[edge:-edge], LineList[1].yarray[edge:-edge] - LineList[0].yarray[edge:-edge], label='Diff')
    pl.figure()
    offset = 0.0
    lscale = -1.0
    for i in range(len(LineList)):
        if nblocks > 1:
            ax = pl.subplot(nblocks,2,i+1)
        lscale = -lscale
        pl.plot(LineList[i].xarray[edge:-edge],
                LineList[i].yarray[edge:-edge]*lscale+offset,
                label='%s' % label[i])
        # offset = offset + 0.5
    #pl.axis([-20,20,-1,1])
    pl.xlabel('VSRC')
    pl.legend()
    pl.suptitle("%s : block=%d bank=%d" % (I.source,block, Opts.bank))
    if Opts.show:
        pl.show()
    else:
        pout = 'bs%d.png' % block
        pl.savefig(pout)
        print("%s written" % pout)


if __name__ == '__main__':
    main(sys.argv[1:])
