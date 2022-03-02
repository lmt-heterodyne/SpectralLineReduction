#! /usr/bin/env python
#
#    lmtar:   find files for an obsnum
#

"""
Usage: lmtar.py OBSNUM

-h --help  This help


This routine finds all LMT raw files for given OBSNUM

Example 1: 

     cd $DATA_LMT
     rsync -avR `lmtar.py 79447 79448` lma:/lma/data_lmt

Example 2: Giving my collaborator the raw data for an observation

     cd $DATA_LMT
     tar cf ~/public_html/secret.tar `lmtar.py 79447 79448`

"""

import os
import sys
import math
import numpy as np		
import matplotlib.pyplot as pl
import datetime

import sys, os
import glob
import netCDF4
from  docopt import docopt


#arguments = docopt(__doc__,options_first=True, version='0.1')

data_lmt = os.environ['DATA_LMT']

for obsnum in sys.argv[1:]:
    obsnum  = int(obsnum)
    obsnum5 = '%d' % obsnum
    obsnum6 = '%06d' % obsnum


    os.chdir(data_lmt)
    fn = glob.glob('ifproc/ifproc_*_%s*.nc' % obsnum6)
    for f in fn:
        print(f)

    fn = glob.glob('spectrometer/roach?/roach?_%s_*nc' % obsnum5)
    for f in fn:
        if f.find('allantest') > 0: continue
        print(f)

    fn = glob.glob('RedshiftChassis?/RedshiftChassis*_%s_*.nc'  % obsnum6)
    for f in fn:
        print(f)
