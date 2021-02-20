#! /usr/bin/env python
#
#    lmtinfo:    provide some info, in terms of a list of "keyword=value", that could be used
#                by an external pipeline that needs input.
#                This list is both bash and python friendly: we only allow integer/float/string
#
#  To run for all the RSR and SLR in Feb 2021 took 9 mins on "cln"

"""
Usage: lmtinfo.py OBSNUM
       lmtinfo.py IFPROCFILE
       lmtinfo.py PATH OBSNUM
       lmtinfo.py PATH

-h --help  This help


This routine grabs some useful summary information from the ifproc file, ignoring
the roach files. If one unique OBSNUM is given, it will show this information
in a "rc" style for the pipeline. If more OBSNUM are possible, for example by only
giving a PATH, all possible OBSNUMs are listed with a short summary, one OBSNUM
per line. Example of output:

      #     DATE  OBSNUM   OBSPGM SOURCE      RESTFRQ VLSR INTTIME
      2018-11-16  079447   Cal    IRC+10216   115.271  -20       8
      2018-11-16  079448   Map    IRC+10216   115.271  -20     686
      2020-02-18  090910   Cal    NGC5194     115.271  463       7
      2020-02-18  090911   Map    NGC5194     115.271  463    3986
      2020-02-20  091111   Cal    NGC5194     115.271  463       7
      2020-02-20  091112   Map    NGC5194     115.271  463    6940

OBSNUM for early SLR (testing?) are 99nnnnn,
but after 2018-04-14 back to the normal nnnnnn, where 074686
seems to be the first.

"""

import sys
import math
import numpy as np		
import matplotlib.pyplot as pl
import datetime

import sys, os
import glob
import netCDF4

from docopt import docopt

#  ifproc/ifproc_2018-06-29_078085_00_0001.nc
#  spectrometer/roach0/roach0_78085_0_1_CHI-Cyg_2018-06-29_041713.nc
#  RedshiftChassis0/RedshiftChassis0_2015-01-22_033551_00_0001.nc

def slr_summary(ifproc, rc=False):
    """   summary a procnum in a single line
    """
    #   e.g.  M51_data/ifproc/ifproc_2020-02-20_091111_00_0001.nc
    fn  = ifproc.split("/")[-1].replace('.nc','').split('_')
    #   e.g. ['ifproc', '2020-02-20', '091111', '00', '0001']
    
    nc = netCDF4.Dataset(ifproc)
    vlsr = nc.variables['Header.Source.Velocity'][0]
    src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
    skyfreq  = nc.variables['Header.Sequoia.SkyFreq'][0]
    restfreq = nc.variables['Header.Sequoia.LineFreq'][0]
    bbtime = nc.variables['Data.IfProc.BasebandTime']
    # Header.Dcs.ObsNum 
    obspgm = b''.join(nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
    # the following Map only if obspgm=='Map'
    if obspgm=='Map':
        xlen = nc.variables['Header.Map.XLength'][0] * 206264.806
        ylen = nc.variables['Header.Map.YLength'][0] * 206264.806
        hpbw = nc.variables['Header.Map.HPBW'][0]

    date_obs = nc.variables['Data.TelescopeBackend.TelTime'][0].tolist()
    date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')
    
    
    t0 = float(bbtime[0].data)
    t1 = float(bbtime[-1].data)
    dt = t1-t0
    nc.close()

    if rc:
        print('# <lmtinfo>')
        print('# ifproc="%s"' % ifproc)
        print('# date-obs="%s"' % date_obs)
        print('# inttime=%g sec' % dt)
        print('# obspgm="%s"' % obspgm)
        print('vlsr=%g' % vlsr)
        print('skyfreq=%g' % skyfreq)
        print('restfreq=%g' % restfreq)
        print('src="%s"' % src)
        resolution = math.ceil(1.0 * 299792458 / skyfreq / 1e9 / 50.0 * 206264.806)
        print('resolution=%g' % resolution)
        print('cell=%g' % (resolution/2.0))
        print('x_extent=%g' % xlen)
        print('y_extent=%g' % ylen)
        
        print("# </lmtinfo>")
    else:    
        print("%s %s  %-5s %-20s %g %g %g" % (date_obs, fn[2], obspgm, src, restfreq, vlsr, dt))


def rsr_summary(rsr_file, rc=False):
    def new_date_obs(date):
        """
        date_obs from RSR have a few common non-ISO formats:
        date = '30/03/2016 03:50:08'   case-1
        date = '2013-12-16 21:10:07'   case-2
        date = '05-03-2020 02:19:06'   case-3
        date = '01:57:07 20/05/18'     case-4
        """
        d  = date.split()
        nd = len(d)
        if nd == 1:
            return date
        if nd == 2:
            if date[2]==':':                      # case-4
                dmy = d[1].split('/')
                return '20%s-%s-%sT%s' % (dmy[2],dmy[1],dmy[0],d[0])
            if date[2]=='-':                      # case-3
                dmy = d[0].split('-')
                return '%s-%s-%sT%s' % (dmy[2],dmy[1],dmy[0],d[1])                
            if date[4]=='-':                      # case-2
                return '%sT%s' % (d[0],d[1])
            if date[2]=='/':                      # case-1
                dmy = d[0].split('/')
                return '%s-%s-%sT%s' % (dmy[2],dmy[1],dmy[0],d[1])
        # uncaught cases
        return date
        
                
    # RedshiftChassis2/RedshiftChassis2_2015-01-22_033551_00_0001.nc
    nc = netCDF4.Dataset(rsr_file)

    # Header.Source.SourceName
    src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
    
    # Header.Dcs.ObsNum = 33551 ;
    obsnum = nc.variables['Header.Dcs.ObsNum'][0] 
    
    # Header.Radiometer.UpdateDate = "21/01/2015 23:12:07
    date_obs = b''.join(nc.variables['Header.Radiometer.UpdateDate'][:]).decode().strip()
    date_obs = new_date_obs(date_obs)
    
    # Header.Weather.UpdateDate = "22/01/15 0:39:48
    # Header.Source.Ra
    # Header.Source.Dec
    ra  = nc.variables['Header.Source.Ra'][0]  * 57.2957795131
    dec = nc.variables['Header.Source.Dec'][0] * 57.2957795131
    
        
    nc.close()

    # no rc mode, only one line summary
    print("%s  %d  RSR   %-30s  %.6f %.6f" %   (date_obs, obsnum, src, ra, dec))


#  although we grab the command line arguments here, they are actually not
#  used in the way most scripts use them. Below there is a more hardcoded
#  parsing of arguments based on how many there are, which are files, and
#  which are directories.
arguments = docopt(__doc__,options_first=True, version='0.1')
#print(arguments)

if len(sys.argv) == 2:
                                                     # mode 1: obsnum or nc_file or path
    obsnum = sys.argv[1]
    fn = glob.glob('*/ifproc/ifproc_*%s*.nc' % obsnum)
    if len(fn) > 0:
        ifproc = fn[0]
    else:
        ifproc = sys.argv[1]

    if os.path.isdir(ifproc):
        path = ifproc

        # pick one, but they all seem to have different # data, 1 has the most
        #RedshiftChassis0_2011-05-08_001809_00_0001.nc - RedshiftChassis0_2020-03-05_092087_00_0001.nc
        #RedshiftChassis1_2011-05-09_001819_00_0001.nc - RedshiftChassis1_2020-03-11_092345_00_0001.nc
        #RedshiftChassis2_2011-05-09_001819_00_0001.nc - RedshiftChassis2_2020-03-11_092345_00_0001.nc
        #RedshiftChassis3_2013-05-04_007484_00_0001.nc - RedshiftChassis3_2020-03-11_092345_00_0001.nc

        chassis = 1
        if chassis < 0:
            globs = '%s/RedshiftChassis?/RedshiftChassis?_*.nc'  % (path)
        else:
            globs = '%s/RedshiftChassis%d/RedshiftChassis%d_*.nc'  % (path,chassis,chassis)            
        fn = glob.glob(globs)
        for f in fn:
            # print('RSR',f)
            try:
                rsr_summary(f)
            except:
                print("Failing on ",f)
            

        globs = '%s/ifproc/ifproc*.nc' % path
        fn = glob.glob(globs)
        for f in fn:
            try:
                slr_summary(f)
            except:
                try:
                    yyyymmdd = f.split('/')[-1].split('_')[1]
                    obsnum   = f.split('/')[-1].split('_')[2]
                    print("%s          %s   failed on %s" % (yyyymmdd,obsnum,f))
                except:
                    print("1900-00-00       failed on %s" % f)

        sys.exit(0)
    elif os.path.exists(ifproc):
        try:
            slr_summary(ifproc,rc=True)
        except:
            print("%s: failed" % ifproc)
               
elif len(sys.argv) == 3:
                                                     # mode 2: path and obsnum : differentiate between SLR and RSR
    path = sys.argv[1]
    obsnum = sys.argv[2]
    globs = '%s/ifproc/ifproc*%s*.nc' % (path,obsnum)
    fn = glob.glob(globs)
    if len(fn) > 0:
        ifproc = fn[0]
        slr_summary(ifproc,True)
    else:
        globs = '%s/RedshiftChassis?/RedshiftChassis?_*%s*.nc'  % (path,obsnum)
        print("Trying RSR %s" % globs)
        fn = glob.glob(globs)
        if len(fn) > 0:
            for f in fn:
                print(f)
        else:
            print("Warning - no RSR files found")
else:
                                                     # no other modes
    print("Usage : %s [path] obsnum" % sys.argv[0])
    sys.exit(0)


