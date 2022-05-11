#! /usr/bin/env python
#
#  lmtinfo.py:    provide some info on the RAW data
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
       lmtinfo.py data
       lmtinfo.py build
       lmtinfo.py grep TERM1 [TERM2 ...]
       lmtinfo.py find TERM1 [TERM2 ...]

-h --help  This help


This routine grabs some useful summary information for LMT raw data.
For SLR they are taken from the ifproc file, ignoring the roach files.

If an OBSNUM (5 or 6 digits) is given, it will show this information
in a "rc" style for the pipeline. All OBSNUM summaries are given in
tabular format.
Example of early output:


      #     DATE  OBSNUM   OBSPGM SOURCE      RESTFRQ VLSR INTTIME
      2018-11-16  079447   Cal    IRC+10216   115.271  -20       8
      2018-11-16  079448   Map    IRC+10216   115.271  -20     686
      2020-02-18  090910   Cal    NGC5194     115.271  463       7
      2020-02-18  090911   Map    NGC5194     115.271  463    3986
      2020-02-20  091111   Cal    NGC5194     115.271  463       7
      2020-02-20  091112   Map    NGC5194     115.271  463    6940

OBSNUM for early SLR (testing?) are 99nnnnn,
but after 2018-04-14 back to the normal nnnnnn, where 074686
seems to be the first data here.

data:     show the database, no sorting and culling
build:    rebuild the sorted database (needs write permission in $DATA_LMT)
grep:     search in database, terms are logically AND-ed
find:     search in database, terms are logically AND-ed

"""

import os
import sys
import math
import glob
import numpy as np		
import datetime
import netCDF4
from docopt import docopt

version="10-may-2022"

if "DATA_LMT" in os.environ:
    data_lmt = os.environ["DATA_LMT"]
else:
    data_lmt = "/data_lmt/"

arguments = docopt(__doc__,options_first=True, version='0.3')

header = "# Y-M-D   T H:M:S     ObsNum ObsGoal       ObgPgm    SourceName                ProjectId                 RestFreq  VLSR   TINT     RA        DEC          AZ    EL"

def grep(terms):
    """
    search a predefined $DATA_LMT/data_lmt.log file for matching terms
    """
    logfile = "%s/%s" % (data_lmt, "data_lmt.log")
    if not os.path.exists(logfile):
        print("Logfile %s does not exist, use the build option to create it" % logfile)
        sys.exit(1)

    if len(terms) == 1:
        cmd = "grep -i %s %s" % (terms[0],logfile)
    elif len(terms) == 2:
        cmd = "grep -i %s %s | grep -i %s" % (terms[0],logfile,terms[1])
    elif len(terms) == 3:
        cmd = "grep -i %s %s | grep -i %s | grep -i %s" % (terms[0],logfile,terms[1],terms[2])
    elif len(terms) == 4:
        cmd = "grep -i %s %s | grep -i %s | grep -i %s | grep -i %s" % (terms[0],logfile,terms[1],terms[2],terms[3])
    else:
        print("too many arguments")
        sys.exit(1)        
    os.system(cmd)


def build():
    """
    build the greppable database $DATA_LMT/data_lmt.log 
    """
    cmd = "cd $DATA_LMT; make -f $LMTOY/data_lmt/Makefile new"
    os.system(cmd)

def alist(x):
    """
    print a python array or list as a comma separated string of values
    [1,2,3] -> "1,2,3"
    """
    n = len(x)
    s = repr(x[0])
    if n==1:
        return s
    
    for i in range(1,n):
        s = s + ",%s" % repr(x[i])
    return s


#  Examples:
#  ifproc/ifproc_2018-06-29_078085_00_0001.nc
#         ifproc_2018-02-26_9901395_00_0001.nc        (older)
#  spectrometer/roach0/roach0_78085_0_1_CHI-Cyg_2018-06-29_041713.nc
#  RedshiftChassis0/RedshiftChassis0_2015-01-22_033551_00_0001.nc

def slr_summary(ifproc, rc=False):
    """   summary of a procnum in a single line, or rc format
    """
    #   e.g.  M51_data/ifproc/ifproc_2020-02-20_091111_00_0001.nc
    fn  = ifproc.split("/")[-1].replace('.nc','').split('_')
    #   e.g. ['ifproc', '2020-02-20', '091111', '00', '0001']
    
    nc = netCDF4.Dataset(ifproc)
    obsnum = nc.variables['Header.Dcs.ObsNum'][0]
    receiver = b''.join(nc.variables['Header.Dcs.Receiver'][:]).decode().strip()
    tau = nc.variables['Header.Radiometer.Tau'][0]
    
    vlsr = nc.variables['Header.Source.Velocity'][0]
    src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
    if receiver == 'Sequoia':
        instrument = 'SEQ'
        numbands = nc.variables['Header.Sequoia.NumBands'][0]        
        skyfreq  = nc.variables['Header.Sequoia.SkyFreq'][:numbands]
        restfreq = nc.variables['Header.Sequoia.LineFreq'][:numbands]
        bbtime = nc.variables['Data.IfProc.BasebandTime'][:]
    elif receiver == 'Msip1mm':
        instrument = '1MM'        
        numbands = nc.variables['Header.Msip1mm.NumBands'][0]
        skyfreq  = nc.variables['Header.Msip1mm.SkyFreq'][:numbands]
        restfreq = nc.variables['Header.Msip1mm.LineFreq'][:numbands]
        bbtime = nc.variables['Data.IfProc.BasebandTime'][:]
        if len(bbtime.shape) > 1:
            # For the 1mm receiver, some signals are sampled 5 times faster to better sample the chopped beam.
            # print('# Warning: PJT bbtime',bbtime.shape,'for obsnum',obsnum)
            bbtime = bbtime[:,0]
    else:
        print('receiver %s not supported yet' % receiver)
        sys.exit(1)

    bufpos = nc.variables['Data.TelescopeBackend.BufPos'][:]
    ubufpos = np.unique(bufpos)

    try:
        calobsnum = nc.variables['Header.IfProc.CalObsNum'][0]
    except:
        calobsnum = -1
        
    obspgm  = b''.join(nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
    #  ObsGoal is in Dcs as well as IfProc
    obsgoal = b''.join(nc.variables['Header.Dcs.ObsGoal'][:]).decode().strip()    
    try:
        pid = b''.join(nc.variables['Header.Dcs.ProjectId'][:]).decode().strip()
    except:
        pid = "Unknown"
    
    # the following Map only if obspgm=='Map'
    if obspgm=='Map':
        map_coord = b''.join(nc.variables['Header.Map.MapCoord'][:]).decode().strip()
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
        print('skytime=%g' % tsky)
        print('inttime=%g' % tint)
        print('obsnum="%s"' % obsnum)
        print('calobsnum="%s"' % calobsnum)
        print('obspgm="%s"' % obspgm)
        print('obsgoal="%s"' % obsgoal)
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
        print('skyfreq=%s     # GHz' % alist(skyfreq))
        print('restfreq=%s    # Ghz' % alist(restfreq))
        print('src="%s"' % src)
        resolution = math.ceil(1.0 * 299792458 / skyfreq[0] / 1e9 / 50.0 * 206264.806)
        print('resolution=%g  # arcsec' % resolution)
        print('cell=%g   # arcsec' % (resolution/2.0))
        # @todo https://github.com/astroumd/lmtoy/issues/9     xlen needs to be equal to ylen
        print('x_extent=%g   # arcsec' % xlen)
        print('y_extent=%g   # arcsec' % ylen)
        print('instrument="%s"' % instrument)
        print('tau=%g' % tau)
        print("# </lmtinfo>")
    else:
        print("%-20s %7s  %-12s %-9s %-25s %-30s %8.4f %5.f    %6.1f  %10.6f %10.6f  %5.1f %5.1f  %g %g" %
              (date_obs, obsnum, obsgoal, obspgm +(('('+map_coord+')') if obspgm=='Map' else ''), src, pid, restfreq[0], vlsr, tint, ra, dec, az, el, az1, el1))

    # -end slr_summary() 

def rsr_summary(rsr_file, rc=False):
    """
    summary of a RSR in a single line, or rc formatted
    """
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
    # ObsGoal
    obsgoal = b''.join(nc.variables['Header.Dcs.ObsGoal'][:]).decode().strip()

    # Header.Radiometer.UpdateDate = "21/01/2015 23:12:07
    #date_obs = b''.join(nc.variables['Header.Radiometer.UpdateDate'][:]).decode().strip()
    #date_obs = new_date_obs(date_obs)
    date_obs = nc.variables['Data.Sky.Time'][0].tolist()
    date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')

    tau = nc.variables['Header.Radiometer.Tau'][0]    
    
    # Header.Weather.UpdateDate = "22/01/15 0:39:48
    # Header.Source.Ra
    # Header.Source.Dec
    ra  = nc.variables['Header.Source.Ra'][0]  * 57.2957795131
    dec = nc.variables['Header.Source.Dec'][0] * 57.2957795131

    az  = nc.variables['Header.Sky.AzReq'][0]  * 57.2957795131
    el  = nc.variables['Header.Sky.ElReq'][0] * 57.2957795131

    t = nc.variables['Data.Sky.Time'][:].tolist()
    if len(t) > 2:
        tint = t[-1]-t[0] + (t[-1]-t[-2])
    elif len(t) > 1:
        tint = t[-1]-t[0]
    elif len(t) > 0:
        tint = 1
    else:
        tint = 0

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
        print('obsgoal="%s"' % obsgoal)
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
        print('tau=%g' % tau)
        print("# </lmtinfo>")

    else:     # one line summary
        #import pdb; pdb.set_trace()
        #print(date_obs, obsnum, obsgoal, obspgm, src,  pid,         tint,   ra,    dec,   az,   el)
        if True:
            print("%-20s %7d  %-12s %-9s %-25s %-30s RSR      0    %6.1f  %10.6f %10.6f  %5.1f %5.1f" %
              (date_obs, obsnum, obsgoal, obspgm, src,  pid,         tint,   ra,    dec,   az,   el))

    # -end  rsr_summary()
    

#   SLR
#   print("%-20s %7s  %-5s %-30s %g %g %g" % (date_obs, fn[2], obspgm, src, restfreq, vlsr, dt))

#  although we grab the command line arguments here, they are actually not
#  used in the way most scripts use them. Below there is a hardcoded
#  parsing of arguments based on how many there are, which are files, and
#  which are directories.

if len(sys.argv) == 2:

    # build
    if sys.argv[1] == "build":
        print("Rebuilding $DATA_LMT/data_lmt.log for grep")
        build()
        sys.exit(0)

    # single obsnum (for SLR or RSR)
    obsnum = sys.argv[1]
    if obsnum.isnumeric():
        # obsnum=int(obsnum)
        # globs = '%s/ifproc/ifproc_*_%06d_*.nc' % (data_lmt,obsnum)
        globs = '%s/ifproc/ifproc_*_*%s_*.nc' % (data_lmt,obsnum)
        #print("# GLOBS slr:",globs)
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        if len(fn) == 1:
            slr_summary(fn[0],rc=True)
            sys.exit(0)
        elif len(fn) > 0:
            print("Multiple ifproc? ",fn)
            sys.exit(0)

        # since no SLR found, try RSR
        globs = '%s/RedshiftChassis?/RedshiftChassis?_*%s*.nc'  % (data_lmt,obsnum)
        # print("# GLOBS rsr:" % globs)
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        if len(fn) > 0 and len(fn) < 5:
            rsr_summary(fn[0], rc=True)
            sys.exit(0)
        elif len(fn) > 4:
            print("Multiple RSR ",fn)
            sys.exit(0)

        # since no RSR found, give up
        print("# No matching OBSNUM %s" % obsnum)
        sys.exit(0)

    # replacement for $DATA_LMT
    if sys.argv[1] == 'data':
        print(header)
    elif os.path.isdir(sys.argv[1]):
        data_lmt = sys.argv[1]
        print(header)
    else:
        print("no more valid options")
        sys.exit(0)

    if True:

        # pick one, but they all seem to have different # data, 1 has the most
        #RedshiftChassis0_2011-05-08_001809_00_0001.nc - RedshiftChassis0_2020-03-05_092087_00_0001.nc
        #RedshiftChassis1_2011-05-09_001819_00_0001.nc - RedshiftChassis1_2020-03-11_092345_00_0001.nc
        #RedshiftChassis2_2011-05-09_001819_00_0001.nc - RedshiftChassis2_2020-03-11_092345_00_0001.nc
        #RedshiftChassis3_2013-05-04_007484_00_0001.nc - RedshiftChassis3_2020-03-11_092345_00_0001.nc

        chassis = 1
        if chassis < 0:
            globs = '%s/RedshiftChassis?/RedshiftChassis?_*.nc'  % (data_lmt)
        else:
            globs = '%s/RedshiftChassis%d/RedshiftChassis%d_*.nc'  % (data_lmt,chassis,chassis)            
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        print("# Found %d RSR with %s" % (len(fn),globs))        
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
                print("%-20s %7s  failed for rsr %s" % (yyyymmdd,obsnum,f))                    

        globs = '%s/ifproc/ifproc_*.nc' % data_lmt
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        print("# Found %d SLR with %s" % (len(fn),globs))
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
                print("%-20s %7s  failed for slr %s" % (yyyymmdd,obsnum,f))
        sys.exit(0)
               
elif len(sys.argv) == 3:

    # special case:
    if sys.argv[1] == "grep" or sys.argv[1] == "find":
        print(header)
        grep(sys.argv[2:])
        sys.exit(0)

    # there should be no more options now
    print("Illegal option ",sys.argv[1])


else:
    # grep allows more terms
    if sys.argv[1] == "grep":
        print(header)
        grep(sys.argv[2:])
        sys.exit(0)

    # otherwise illegal options, so give help
    
                                                     # no other modes
    print("Usage : %s [path] obsnum" % sys.argv[0])
    sys.exit(0)
