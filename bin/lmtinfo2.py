#! /usr/bin/env python
#
#  lmtinfo2.py:    provide some info on the SpecFile in rc format
#

_version="11-jul-2023"

_help = """
Usage: lmtinfo.py specfile

-h --help      This help
-v --version   Script version

version = %s

""" % _version


import os
import sys
import math
import glob
import numpy as np		
import datetime
import netCDF4
from docopt import docopt

# arguments = docopt(_help,options_first=True, version=_version)


def specfile_summary(specfile):
    """   summary of a specfile in rc format
    """
    
    nc = netCDF4.Dataset(specfile)
    
    nchan0 = nc.dimensions['nchan0'].size
    chan0 =  nc.dimensions['chan0'].size

    
    nchan = nc.variables['Header.Line.NChannels'][0]
    cdelt = nc.variables['Header.SpectrumAxis.CDELT'][0]
    crpix = nc.variables['Header.SpectrumAxis.CRPIX'][0]
    crval = nc.variables['Header.SpectrumAxis.CRVAL'][0]
    ctype = b''.join(nc.variables['Header.SpectrumAxis.CTYPE'][:]).decode().strip()
    bank  = nc.variables['Header.LineData.Bank'][0]

    nc.close()
    
    vmin = crval - chan0*cdelt
    vmax = vmin + (nchan0-1)*cdelt
        
    print('# <lmtinfo2>')
    print('# version=%s' % _version)
    print('bank=%d' % bank)
    print('nchan0=%d' % nchan0)
    print('chan0=%d' % chan0)
    print('vmin=%g' % vmin)
    print('vmax=%g' % vmax)
    print("# </lmtinfo2>")

    # -end specfile_summary() 

# ==================================================================================================================            

specfile = sys.argv[1]
specfile_summary(specfile)
