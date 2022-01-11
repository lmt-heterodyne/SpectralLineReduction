#! /usr/bin/env python
#
#    lmtinfo:    provide some info, in terms of a list of "keyword=value", that could be used
#                by an external pipeline that needs input.
#                This list is both bash and python friendly: we only allow integer/float/string
#
#  To run for all the RSR and SLR in Feb 2021 took 9 mins on "cln"
#
#  Examples of using the online version:
#      http://187.248.54.232/lmtmc/notes/LmtShiftReportByDate.html
#      http://187.248.54.232/cgi-bin/lmtmc/mc_sql.cgi?-s=2018-12-15&-e=2018-12-20&-project=&-obsGoal=all&-format=html&-instrument=all
#
#  @todo       show with "0123" that a roach/chassis is present. If not, put a * , e.g. "01*3" means roach2 is missing.
#              ifproc            uses %06d_%02d_%04d for ObsNum,SubObsNum,ScanNum
#              spectrometer      uses %d_%d_%d
#              RedshiftChassisN  uses %06d_%02d_%04d 
#
#
"""
Usage: lmtinfo.py OBSNUM
       lmtinfo.py IFPROCFILE
       lmtinfo.py PATH OBSNUM
       lmtinfo.py PATH
       lmtinfo.py grep TERM1 [TERM2 ...]
       lmtinfo.py build

-h --help  This help


This routine grabs some useful summary information from the ifproc file, ignoring
the roach files. If one unique OBSNUM is given, it will show this information
in a "rc" style for the pipeline. If more OBSNUM are possible, for example by only
giving a PATH, all possible OBSNUMs are listed with a short summary, one OBSNUM
per line. Example of early output:


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

grep:     search in database, terms are logically AND-ed

build:    rebuild the database

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

version="11-jan-2022"

def grep(terms):
    """
    search a predefined $DATA_LMT/data_lmt.log file for terms
    @todo check if the log file exists
    """

    if len(terms) == 1:
        cmd = "grep -i %s $DATA_LMT/data_lmt.log" % (terms[0])
    elif len(terms) == 2:
        cmd = "grep -i %s $DATA_LMT/data_lmt.log | grep -i %s" % (terms[0],terms[1])
    elif len(terms) == 3:
        cmd = "grep -i %s $DATA_LMT/data_lmt.log | grep -i %s | grep -i %s" % (terms[0],terms[1],terms[2])
    elif len(terms) == 4:
        cmd = "grep -i %s $DATA_LMT/data_lmt.log | grep -i %s | grep -i %s | grep -i %s" % (terms[0],terms[1],terms[2],terms[3])
    os.system(cmd)


def build():
    """
    search a predefined $DATA_LMT/data_lmt.log file for terms
    @todo check if the log file exists
    """
    cmd = "cd $DATA_LMT; make new"
    os.system(cmd)


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
    if 'Header.Sequoia.SkyFreq' in nc.variables:
        skyfreq  = nc.variables['Header.Sequoia.SkyFreq'][0]
        restfreq = nc.variables['Header.Sequoia.LineFreq'][0]
        instrument = 'SEQ'
        bbtime = nc.variables['Data.IfProc.BasebandTime'][:]
    else:
        # @todo is 1MM the only option
        skyfreq  = nc.variables['Header.Msip1mm.SkyFreq'][0]
        restfreq = nc.variables['Header.Msip1mm.LineFreq'][0]
        instrument = '1MM'        
        bbtime = nc.variables['Data.IfProc.BasebandTime'][:]
        if len(bbtime.shape) > 1:
            print('Warning: PJT1',bbtime,'pjt2',bbtime[0],'pjt3',bbtime.shape)
            bbtime = bbtime[:,0]
        
    bufpos = nc.variables['Data.TelescopeBackend.BufPos'][:]
    ubufpos = np.unique(bufpos)
    # Header.Dcs.ObsNum

    # Header.Sequoia.NumBands or Header.IfProc.NumBands
    numbands = nc.variables['Header.Sequoia.NumBands'][0]

    obsnum = nc.variables['Header.Dcs.ObsNum'][0]
    try:
        calobsnum = nc.variables['Header.IfProc.CalObsNum'][0]
    except:
        calobsnum = -1
        
    obspgm = b''.join(nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
    try:
        pid = b''.join(nc.variables['Header.Dcs.ProjectId'][:]).decode().strip()
    except:
        pid = "Unknown"
    
    # the following Map only if obspgm=='Map'
    if obspgm=='Map':
        xlen = nc.variables['Header.Map.XLength'][0] * 206264.806
        ylen = nc.variables['Header.Map.YLength'][0] * 206264.806
        xoff = nc.variables['Header.Map.XOffset'][0] * 206264.806
        yoff = nc.variables['Header.Map.YOffset'][0] * 206264.806
        xram = nc.variables['Header.Map.XRamp'][0]   * 206264.806
        yram = nc.variables['Header.Map.YRamp'][0]   * 206264.806
        hpbw = nc.variables['Header.Map.HPBW'][0]    * 206264.806
        xstep= nc.variables['Header.Map.XStep'][0] 
        ystep= nc.variables['Header.Map.YStep'][0] 
    else:
        xlen = 0
        ylen = 0
        xoff = 0
        yoff = 0
        xram = 0
        yram = 0
        hpbw = 0
        xstep= 0
        ystep= 0

    date_obs = nc.variables['Data.TelescopeBackend.TelTime'][0].tolist()
    date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')

    ra  = nc.variables['Header.Source.Ra'][0]  * 57.2957795131
    dec = nc.variables['Header.Source.Dec'][0] * 57.2957795131
    az  = nc.variables['Header.Sky.AzReq'][0]  * 57.2957795131
    el  = nc.variables['Header.Sky.ElReq'][0] * 57.2957795131

    az1 = nc.variables['Header.Sky.AzOff'][1] * 206264.81
    el1 = nc.variables['Header.Sky.ElOff'][1] * 206264.81

    tint = nc.variables['Header.Dcs.IntegrationTime'][0]
    
    # Header.Dcs.ObsGoal
    # Header.ScanFile.Valid = 1 ;

    t0 = float(bbtime[0])
    t1 = float(bbtime[-1])
    t2 = float(bbtime[-2])
    tsky = t1-t0 + (t1-t2)
    
    nc.close()
        
    if rc:
        print('# <lmtinfo>')
        print('# version=%s' % version)
        print('# ifproc="%s"' % ifproc)
        print('date_obs="%s"' % date_obs)
        print('# skytime=%g sec' % tsky)
        print('# inttime=%g sec' % tint)
        print('obsnum="%s"' % obsnum)
        print('calobsnum="%s"' % calobsnum)
        print('obspgm="%s"' % obspgm)
        print('ProjectId="%s"' % pid)
        print('# SkyOff=%g %g' % (az1,el1))
        #print('# bufpos=%s' % str(ubufpos))
        print('# HPBW=%g arcsec' % hpbw)
        print('# XYLength=%g %g arcsec' % (xlen,ylen))
        print('# XYRamp=%g %g arcsec' % (xram,yram))
        print('# XYoff=%g %g arcsec' % (xoff,yoff))
        print('# XYstep=%g %g' % (xstep,ystep))
        print('numbands=%d' % numbands)
        print('vlsr=%g        # km/s' % vlsr)
        print('skyfreq=%g     # GHz' % skyfreq)
        print('restfreq=%g    # Ghz' % restfreq)
        print('src="%s"' % src)
        resolution = math.ceil(1.0 * 299792458 / skyfreq / 1e9 / 50.0 * 206264.806)
        print('resolution=%g  # arcsec' % resolution)
        print('cell=%g   # arcsec' % (resolution/2.0))
        # @todo https://github.com/astroumd/lmtoy/issues/9     xlen needs to be equal to ylen
        print('x_extent=%g   # arcsec' % xlen)
        print('y_extent=%g   # arcsec' % ylen)
        print('instrument="%s"' % instrument)
        print("# </lmtinfo>")
    else:
        print("%-20s %7s  %-5s %-30s %8.4f %5.f    %6.1f  %10.6f %10.6f  %5.1f %5.1f  %g %g" % (date_obs, obsnum, obspgm, src, restfreq, vlsr, tint, ra, dec, az, el, az1,el1))

#       print("%-20s %7d  %-5s %-30s RSR  0      %5.1f  %10.6f %10.6f  %5.1f %5.1f" %   (date_obs, obsnum, obspgm, src, tint, ra, dec, az, el))

    # -end slr_summary() 

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

    try:
        pid = b''.join(nc.variables['Header.Dcs.ProjectId'][:]).decode().strip()
    except:
        pid = "Unknown"
        
    # Header.Dcs.ObsNum = 33551 ;
    obsnum = nc.variables['Header.Dcs.ObsNum'][0]
    # yuck, with the RSR filenameconvention this is the trick to find the chassic
    chassis   = rsr_file[ rsr_file.rfind('/') + 16 ]
    con_name  = 'Header.RedshiftChassis_%s_.CalObsNum' % chassis
    calobsnum = nc.variables[con_name][0]
    
    # Bs, Cal
    obspgm = b''.join(nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()

    # Header.Radiometer.UpdateDate = "21/01/2015 23:12:07
    date_obs = b''.join(nc.variables['Header.Radiometer.UpdateDate'][:]).decode().strip()
    date_obs = new_date_obs(date_obs)
    
    # Header.Weather.UpdateDate = "22/01/15 0:39:48
    # Header.Source.Ra
    # Header.Source.Dec
    ra  = nc.variables['Header.Source.Ra'][0]  * 57.2957795131
    dec = nc.variables['Header.Source.Dec'][0] * 57.2957795131

    az  = nc.variables['Header.Sky.AzReq'][0]  * 57.2957795131
    el  = nc.variables['Header.Sky.ElReq'][0] * 57.2957795131

    t = nc.variables['Data.Sky.Time'][:].tolist()
    tint = t[-1]-t[0] + (t[-1]-t[-2])

    nc.close()

    if rc:
        print('# <lmtinfo>')
        print('# %s' % rsr_file)
        print('date_obs="%s"' % date_obs)
        # print('# skytime=%g sec' % tsky)
        print('inttime=%g # sec' % tint)
        print('obsnum="%s"' % obsnum)
        print('calobsnum="%s"' % calobsnum)
        print('obspgm="%s"' % obspgm)
        print('ProjectId="%s"' % pid)
        #print('# SkyOff=%g %g' % (az1,el1))
        #print('# bufpos=%s' % str(ubufpos))
        #print('# HPBW=%g arcsec' % hpbw)
        #print('# XYLength=%g %g arcsec' % (xlen,ylen))
        #print('# XYRamp=%g %g arcsec' % (xram,yram))
        print('src="%s"' % src)
        #resolution = math.ceil(1.0 * 299792458 / skyfreq / 1e9 / 50.0 * 206264.806)
        #print('resolution=%g  # arcsec' % resolution)
        #print('cell=%g   # arcsec' % (resolution/2.0))
        # @todo https://github.com/astroumd/lmtoy/issues/9     xlen needs to be equal to ylen
        #print('x_extent=%g   # arcsec' % xlen)
        #print('y_extent=%g   # arcsec' % ylen)
        print('instrument="RSR"')
        print("# </lmtinfo>")

    else:     # one line summary
        print("%-20s %7d  %-5s %-30s     RSR      0    %6.1f  %10.6f %10.6f  %5.1f %5.1f" %   (date_obs, obsnum, obspgm, src, tint, ra, dec, az, el))

    # -end  rsr_summary()
    

#   SLR
#   print("%-20s %7s  %-5s %-30s %g %g %g" % (date_obs, fn[2], obspgm, src, restfreq, vlsr, dt))

#  although we grab the command line arguments here, they are actually not
#  used in the way most scripts use them. Below there is a more hardcoded
#  parsing of arguments based on how many there are, which are files, and
#  which are directories.
arguments = docopt(__doc__,options_first=True, version='0.1')
#print(arguments)

if len(sys.argv) == 2:

    if sys.argv[1] == "build":
        print("Special rebuild of a new")
        build()
        sys.exit(0)        

    print("# Y-M-D   T H:M:S     ObsNum ObsPgm SourceName                     RestFreq  VLSR   TINT     RA        DEC          AZ    EL")
    
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
                # Failing on  /home/teuben/LMT/data_lmt/RedshiftChassis1/RedshiftChassis1_2013-04-18_006779_00_0004.nc
                try:
                    yyyymmdd = f.split('/')[-1].split('_')[1]
                    obsnum   = f.split('/')[-1].split('_')[2]
                except:
                    yyyymmdd = "1900-00-00"
                    obsnum   = " "
                print("%-20s %7s  failed for %s" % (yyyymmdd,obsnum,f))                    

        globs = '%s/ifproc/ifproc*.nc' % path
        fn = glob.glob(globs)
        for f in fn:
            try:
                slr_summary(f)
            except:
                try:
                    yyyymmdd = f.split('/')[-1].split('_')[1]
                    obsnum   = f.split('/')[-1].split('_')[2]
                except:
                    yyyymmdd = "1900-00-00"
                    obsnum   = " "
                print("%-20s %7s  failed for %s" % (yyyymmdd,obsnum,f))
        sys.exit(0)
    elif os.path.exists(ifproc):
        try:
            slr_summary(ifproc,rc=True)
        except:
            print("%s: failed" % ifproc)
               
elif len(sys.argv) == 3:

    # special case:
    if sys.argv[1] == "grep":
        grep(sys.argv[2:])
        sys.exit(0)

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
        # print("Trying RSR %s" % globs)
        fn = glob.glob(globs)
        if len(fn) > 0:
            rsr_summary(fn[0], rc=True)
        else:
            print("# Warning - no RSR files found, possibly unknown obsnum")
            print("obsnum=0")
else:
    # special case:
    if sys.argv[1] == "grep":
        grep(sys.argv[2:])
        sys.exit(0)
    
                                                     # no other modes
    print("Usage : %s [path] obsnum" % sys.argv[0])
    sys.exit(0)
