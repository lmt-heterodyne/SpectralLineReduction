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

_version="9-mar-2024"

_help = """
Usage: lmtinfo.py OBSNUM
       lmtinfo.py data
       lmtinfo.py build
       lmtinfo.py last
       lmtinfo.py new OBSNUM
       lmtinfo.py lmtot OBSNUM
       lmtinfo.py grep  [TERM1 [TERM2 ...]]
       lmtinfo.py grepw [TERM1 [TERM2 ...]]
       lmtinfo.py find  [TERM1 [TERM2 ...]]

-h --help      This help
-v --version   Script version

This routine grabs some useful summary information for LMT raw data.
For SLR they are taken from the ifproc file, ignoring the roach files.

If an OBSNUM (5 or 6 digits) is given, it will show this information
in a "rc" style for the pipeline. All OBSNUM summaries are given in
tabular format.

Example of output (some columns not shown here)


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
last:     report the last known obsnum
new:      build the database with only new obsnums since the last build
lmtot:    return the LMTOT (observing) file to stdout
grep:     search in database, terms are logically AND-ed - partial matches allowed
grepw:    search in database, terms are logically AND-ed and words need to match exactly
find:     search in database, terms are logically AND-ed

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

if "DATA_LMT" in os.environ:
    data_lmt = os.environ["DATA_LMT"]
else:
    data_lmt = "/data_lmt/"

arguments = docopt(_help,options_first=True, version=_version)

header = "# Y-M-D   T H:M:S     ObsNum  Receiver   ObsGoal      ObgPgm      SourceName                ProjectId                   RestFreq      VLSR    TINT    RA          DEC        AZ     EL"

def grep(terms, flags=""):
    """
    search a predefined $DATA_LMT/data_lmt.log file for matching terms
    """
    logfile = "%s/%s" % (data_lmt, "data_lmt.log")
    if not os.path.exists(logfile):
        print("Logfile %s does not exist, use the build option to create it" % logfile)
        sys.exit(1)

    if len(terms) == 0:
        cmd = "grep %s -v ^# %s" %\
        (flags,logfile)
    elif len(terms) == 1:
        cmd = "grep %s -i %s %s" %\
        (flags,terms[0],logfile)
    elif len(terms) == 2:
        cmd = "grep %s -i %s %s | grep %s -i %s" %\
        (flags,terms[0],logfile,flags,terms[1])
    elif len(terms) == 3:
        cmd = "grep %s -i %s %s | grep %s -i %s | grep %s -i %s" %\
        (flags,terms[0],logfile,flags,terms[1],flags,terms[2])
    elif len(terms) == 4:
        cmd = "grep %s -i %s %s | grep %s -i %s | grep %s -i %s | grep %s -i %s" %\
        (flags,terms[0],logfile,flags,terms[1],flags,terms[2],flags,terms[3])
    else:
        print("too many arguments for now")
        sys.exit(1)        
    os.system(cmd)

def last():
    """
    report the last known obsnum
    """
    fn = data_lmt + '/last.obsnum'
    if os.path.exists(fn):
        lines = open(fn).readlines()
        return int(lines[0])
    print("Warning: no %s found" % fn)
    return -1

def build():
    """
    build the greppable database $DATA_LMT/data_lmt.log
    this does nothing more than "lmtinfo.py $DATA_LMT > data_lmt.log"
    """
    cmd = "cd $DATA_LMT; make -f $LMTOY/data_lmt/Makefile new"
    os.system(cmd)

def new(obsnum):
    """
    update the database from a given obsnum up
    """
    cmd = "cd $DATA_LMT; make -f $LMTOY/data_lmt/Makefile new2 OBSNUM0="
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

def iau(src):
    """
    sanitize a source name to prevent spaces and meta characters?
    e.g. "NGC6946_(CO)"  ->  "NGC6946"
    """
    # for now case by case basis
    if src=="NGC6946_(CO)":   return "NGC6946"
    return src

def pid_sanitize(pid):
    """
    sanitize badly formatted ProjectID's
    For regular projects,   YYYY-S1-XX-YY
    For commisioning:       YYYY{S1}{Instrument}Commissioning
    """
    if pid == "2018S1-MU-8":
        return "2018-S1-MU-8"      # 90648..90666 were mis-labeled
    if pid == "2022MSIP1mmCommissioning":
        return "2022S1MSIP1mmCommissioning"
    if pid == "2024MSIP1mmCommissioning":
        return "2024S1MSIP1mmCommissioning"

    return pid

def dataverse_old(pid):
    """
    input:    Projectid (string)
    output:   dictionary of dataverse key/val pairs
    
    for DataVerse ingestion we need a few new parameters in the rc file
      PIName
      projectTitle
    The LMTOY environment needs to be present for this
    
    """
    # PID,PI,Title
    dbname = os.environ['LMTOY'] + '/etc/ProjectId.tab'
    print("# project info:",dbname)
    fp = open(dbname)
    lines = fp.readlines()
    fp.close()
    #
    db = {}
    for line in lines:
        w = line.strip().split()
        if len(w) < 3:
            continue
        if w[0] == pid:
            db['PIName'] = w[1]
            db['projectTitle'] = w[2]
            return db
    db['PIName']       = 'Unknown'
    db['projectTitle'] = 'Unknown'
    return db

def dataverse(pid):
    """
    input:    Projectid (string)
    output:   dictionary of dataverse key/val pairs
    
    for DataVerse ingestion we need a few new parameters in the rc file
      PIName
      projectTitle
    The LMTOY environment needs to be present for this
    
    """
    db = {}
    if pid.find('Commissioning') >= 0:
        db['PIName'] = 'LMT'
        db['projectTitle'] = 'Commissioning'
        return db
    # PID,Title,PI
    dbname = os.environ['LMTOY'] + '/etc/ProjectId.csv'
    print("# project info:",dbname)
    fp = open(dbname,encoding='utf-8', errors='ignore')
    lines = fp.readlines()
    fp.close()
    #
    for line in lines:
        w = line.strip().split(',')
        if len(w) < 3:
            continue
        if w[0].strip() == pid:
            db['PIName'] = w[2].strip()
            db['projectTitle'] = w[1].strip()
            return db
    # last resort (should not happen)
    db['PIName']       = 'Unknown'
    db['projectTitle'] = 'Unknown'
    return db


def date_obs_utdate(date):
    """ convert from a Header.TimePlace.UTDate fractional year 
            this is a hack until we know how they got this UTDate
    """
    y1 = int(date)
    e0=datetime.datetime(1970,1,1,0,0,0)        
    e1=datetime.datetime(y1,1,1,0,0,0)
    e2=datetime.datetime(y1+1,1,1,0,0,0)
    ys=(e2-e1).total_seconds()
    yd=(date-y1)*ys
    y70 = (e1-e0).total_seconds() + yd
    #print('PJT0 secs',ys,ys/24/3600)
    #print('PJT0 secs1970',(e1-e0).total_seconds())
    #print('PJT0 secs-now',yd)
    #print('PJT0 secs-all',y70)
    date_obs = datetime.datetime.fromtimestamp(y70).strftime('%Y-%m-%dT%H:%M:%S')
    #print('PJT0 date_obs',date_obs)
    return date_obs


def seq_bandwidth(nchan):
    """  older data don't store the 'Header.SpecBackend.Bandwidth' variable
    so we need to retrieve is via nchan.... but neither do we have nchan,
    so this routine isn't used (yet) at the moment
    """
    if nchan==2048:  return 0.8
    if nchan==4096:  return 0.4
    if nchan==8192:  return 0.2
    bw = 0.7999
    print("# WARNING: unknown nchan=%d for SEQ - probably old data and assuming bandwidth=%g" % (nchan,bw))
    return bw
    
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
    subobsnum = nc.variables['Header.Dcs.SubObsNum'][0]
    scannum = nc.variables['Header.Dcs.ScanNum'][0]    
    receiver = b''.join(nc.variables['Header.Dcs.Receiver'][:]).decode().strip()
    tau = nc.variables['Header.Radiometer.Tau'][0]
    
    vlsr = nc.variables['Header.Source.Velocity'][0]
    src = b''.join(nc.variables['Header.Source.SourceName'][:]).decode().strip()
    src = iau(src)
    if receiver == 'Sequoia':
        instrument = 'SEQ'
        numbands = nc.variables['Header.Sequoia.NumBands'][0]        
        skyfreq  = nc.variables['Header.Sequoia.SkyFreq'][:numbands]
        restfreq = nc.variables['Header.Sequoia.LineFreq'][:numbands]
        try:
            bandwidth = nc.variables['Header.SpecBackend.Bandwidth'][:numbands]
        except:
            bandwidth = [seq_bandwidth(1)]
            if numbands > 1:
                print("Warning: numbands=%d and Header.SpecBackend.Bandwidth missing" % numbands)
                bandwidth.append(bandwidth[0])
            
        bbtime = nc.variables['Data.IfProc.BasebandTime'][:]
    elif receiver == 'Msip1mm':
        instrument = '1MM'        
        numbands = nc.variables['Header.Msip1mm.NumBands'][0]
        skyfreq  = nc.variables['Header.Msip1mm.SkyFreq'][:numbands]
        restfreq = nc.variables['Header.Msip1mm.LineFreq'][:numbands]
        bandwidth = nc.variables['Header.SpecBackend.Bandwidth'][:numbands]        
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
    #  ObsGoal is in Dcs as well as IfProc, but only after 2014-10-14
    try:
        obsgoal = b''.join(nc.variables['Header.Dcs.ObsGoal'][:]).decode().strip()
    except:
        obsgoal = "Unknown"
    # ProjectId was not used before some date ?
    try:
        pid = b''.join(nc.variables['Header.Dcs.ProjectId'][:]).decode().strip()
    except:
        pid = "Unknown"
    pid = pid_sanitize(pid)
    
    # the following Map only if obspgm=='Map'
    if obspgm=='Map':
        map_coord = b''.join(nc.variables['Header.Map.MapCoord'][:]).decode().strip()
        map_motion = b''.join(nc.variables['Header.Map.MapMotion'][:]).decode().strip() 
        xlen = nc.variables['Header.Map.XLength'][0] * 206264.806
        ylen = nc.variables['Header.Map.YLength'][0] * 206264.806
        xoff = nc.variables['Header.Map.XOffset'][0] * xlen
        yoff = nc.variables['Header.Map.YOffset'][0] * ylen
        xram = nc.variables['Header.Map.XRamp'][0]   * 206264.806
        yram = nc.variables['Header.Map.YRamp'][0]   * 206264.806
        hpbw = nc.variables['Header.Map.HPBW'][0]    * 206264.806
        srate= nc.variables['Header.Map.ScanRate'][0]* 206264.806
        xstep= nc.variables['Header.Map.XStep'][0] 
        ystep= nc.variables['Header.Map.YStep'][0]
        tsamp= nc.variables['Header.Map.TSamp'][0]        
    else:
        xlen = 0
        ylen = 0
        xoff = 0
        yoff = 0
        xram = 0
        yram = 0
        hpbw = 0
        srate= 0
        xstep= 0
        ystep= 0
        tsamp= 0

    try:
        # @todo   this variable can be missing, e.g. 108778
        date_obs = nc.variables['Data.TelescopeBackend.TelTime'][0].tolist()
        date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')
    except:
        pass

    # Header.TimePlace.UTDate
    date_obs = nc.variables['Header.TimePlace.UTDate'][0].tolist()
    #print('PJT Header.TimePlace.UTDate',date_obs)
    date_obs = date_obs_utdate(date_obs)
    #print('PJT0 date_obs',date_obs)        
    

    ra  = nc.variables['Header.Source.Ra'][0]  * 57.2957795131
    dec = nc.variables['Header.Source.Dec'][0] * 57.2957795131
    az  = nc.variables['Header.Sky.AzReq'][0]  * 57.2957795131
    el  = nc.variables['Header.Sky.ElReq'][0] * 57.2957795131

    az1 = nc.variables['Header.Sky.AzOff'][1] * 206264.81
    el1 = nc.variables['Header.Sky.ElOff'][1] * 206264.81

    tint = nc.variables['Header.Dcs.IntegrationTime'][0]

    # this estimate only works for rectilinear maps
    if xlen > 0 and ylen>0:
        nrows = math.ceil((ylen+2*yram)/(ystep*hpbw))
        rowtint = (xlen+2*xram)/(srate)
        # rowtint = (xlen+2*xram)/(xstep*hpbw)*tsamp    # two ways to get the same
        tint_on = nrows * rowtint
    else:
        nrows = 0
        rowtint = 0
        tint_on = -1
    
    # Header.Dcs.ObsGoal
    # Header.ScanFile.Valid = 1 ;

    t0 = float(bbtime[0])
    t1 = float(bbtime[-1])
    t2 = float(bbtime[-2])
    tsky = t1-t0 + (t1-t2)
    
    nc.close()
        
    if rc:
        # SEQ
        nppb = 2.0
        print('# <lmtinfo>')
        print('# version=%s' % _version)
        print('rawnc="%s"' % ifproc)
        print('date_obs="%s"' % date_obs)
        print('skytime=%g' % tsky)
        print('inttime=%g' % tint)
        print('obsnum=%s' % obsnum)
        print('subobsnum=%s' % subobsnum)
        print('scannum=%s' % scannum)
        print('calobsnum=%s' % calobsnum)
        print('obspgm="%s"' % obspgm)
        if obspgm=='Map':
            print('map_coord="%s"' % map_coord)
            print('map_motion="%s"' % map_motion)
        print('ProjectId="%s"' % pid)
        print("ra=%f  # deg" % ra)
        print("dec=%f # deg" % dec)
        print('# SkyOff=%g %g' % (az1,el1))
        #print('# bufpos=%s' % str(ubufpos))
        print('# HPBW=%g arcsec' % hpbw)
        print('# XYLength=%g %g arcsec' % (xlen,ylen))
        print('# XYRamp=%g %g arcsec' % (xram,yram))
        print('# XYoff=%g %g arcsec' % (xoff,yoff))
        print('# XYstep=%g %g' % (xstep,ystep))
        print('nrows=%d' % nrows)
        print('rowtint=%g' % rowtint)
        print('tint_on=%g  # excluding 1.6sec (?) switchback' % tint_on)
        print('numbands=%d' % numbands)
        print('vlsr=%g        # km/s' % vlsr)
        print('skyfreq=%s     # GHz' % alist(skyfreq))
        print('restfreq=%s    # Ghz' % alist(restfreq))
        print('bandwidth=%s   # Ghz' % alist(bandwidth))        
        # Header.SpecBackend.Bandwidth
        if numbands == 2 and restfreq[1] == 0.0:
            print('numbands=1   # overriding')
            print('restfreq=%s' % repr(restfreq[0]))
        print('src="%s"' % src)
        resolution = 1.0 * 299792458 / skyfreq[0] / 1e9 / 50.0 * 206264.806
        # why is this an integer again?
        resolution = math.ceil(resolution)
        print('resolution=%g  # arcsec' % resolution)
        print('nppb=%g   # number of points per beam' % nppb)
        print('cell=%g   # arcsec' % (resolution/nppb))
        # @todo https://github.com/astroumd/lmtoy/issues/9     xlen needs to be equal to ylen
        print('x_extent=%g   # arcsec' % xlen)
        print('y_extent=%g   # arcsec' % ylen)
        print('instrument="%s"' % instrument)
        print('tau=%g' % tau)
        dv = dataverse(pid)
        if dv != None:
            for k in dv.keys():
                print('%s="%s"' % (k,dv[k]))
        else:
            print("# no dataverse info")
        print("config=%s__SEQ__%s__%s__%s" % (pid,obspgm,src,alist(restfreq)))
        print("# </lmtinfo>")
    else:
        print("%-20s %7s  %-10s %-12s %-11s %-25s %-30s %8.4f %5.f    %6.1f  %10.6f %10.6f  %5.1f %5.1f  %g %g" %
              (date_obs, obsnum, instrument, obsgoal, obspgm +(('('+map_coord+'/'+map_motion[0]+')') if obspgm=='Map' else ''), src, pid, restfreq[0], vlsr, tint, ra, dec, az, el, az1, el1))

    # -end slr_summary() 

def rsr_summary(rsr_file, rc=False):
    """
    summary of a RSR in a single line, or rc formatted
    """
    def new_date_obs(date):
        """
        - not used anymore - 
        date_obs from RSR have a few common non-ISO formats:   (see Data.Sky.Time)
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
    src = iau(src)

    try:
        pid = b''.join(nc.variables['Header.Dcs.ProjectId'][:]).decode().strip()
    except:
        pid = "Unknown"
        
    # Header.Dcs.ObsNum = 33551 ;
    obsnum = nc.variables['Header.Dcs.ObsNum'][0]
    subobsnum = nc.variables['Header.Dcs.SubObsNum'][0]
    scannum = nc.variables['Header.Dcs.ScanNum'][0]    
    receiver = b''.join(nc.variables['Header.Dcs.Receiver'][:]).decode().strip()
    instrument = "RSR"
    # yuck, with the RSR filenameconvention this is the trick to find the chassic
    chassis   = rsr_file[ rsr_file.rfind('/') + 16 ]
    # this is how you'd do it.
    for k in nc.variables.keys():
        if k.startswith('Header.RedshiftChassis') and k.endswith('.ChassisNumber'):
            chassis = nc.variables[k][0]
    con_name  = 'Header.RedshiftChassis_%s_.CalObsNum' % chassis
    calobsnum = nc.variables[con_name][0]
    
    # Bs, Cal
    obspgm = b''.join(nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
    if obspgm=='Map':
        map_coord = b''.join(nc.variables['Header.Map.MapCoord'][:]).decode().strip()
    # ObsGoal
    try:
        obsgoal = b''.join(nc.variables['Header.Dcs.ObsGoal'][:]).decode().strip()
    except:
        obsgoal = "Unknown"
    # Header.Bs.Beam
    try:
        #bsbeam = b''.join(nc.variables['Header.Bs.Beam'][:]).decode().strip()
        bsbeam1 =  nc.variables['Header.Bs.Beam'][0]
        bsbeam2 =  nc.variables['Header.Bs.Beam'][1]
        bsbeam = [bsbeam1,bsbeam2]
        print("PJT",bsbeam)
    except:
        bsbeam = [-1,-1]

    # this variable isn't present in old data (e.g. 11654)
    # Header.Radiometer.UpdateDate = "21/01/2015 23:12:07
    date_obs = b''.join(nc.variables['Header.Radiometer.UpdateDate'][:]).decode().strip()
    date_obs = new_date_obs(date_obs)
    #print('# Header.Radiometer.UpdateDate',date_obs)

    try:
        # @todo   e.g. 108991 misses this variable
        date_obs = nc.variables['Data.Sky.Time'][0].tolist()
        #print('# Data.Sky.Time',date_obs)
        if date_obs < 1000000:
            #print('PJT  bad time, this was < 2015 when it was "since boot"')
            pass
        date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')
        #print('# date_obs',date_obs)
    except:
        pass

    # Header.TimePlace.UTDate
    date_obs = nc.variables['Header.TimePlace.UTDate'][0].tolist()
    #print('PJT Header.TimePlace.UTDate',date_obs)
    date_obs = date_obs_utdate(date_obs)
    #print('PJT0 date_obs',date_obs)        

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

    try:
        tint = nc.variables['Header.Dcs.IntegrationTime'][0]
    except:
        # older data (e.g. 28190) missing this??? - mark it with 300.1 so we know it's "fake"
        tint = 300.1

    nc.close()

    if rc:
        # RSR
        print('# <lmtinfo>')
        print('# version=%s' % _version)
        print('rawnc="%s"' % rsr_file)
        print('date_obs="%s"' % date_obs)
        # print('# skytime=%g sec' % tsky)
        print('inttime=%g # sec' % tint)
        print('obsnum=%s' % obsnum)
        print('subobsnum=%s' % subobsnum)
        print('scannum=%s' % scannum)
        print('calobsnum=%s' % calobsnum)
        print('obspgm="%s"' % obspgm)
        print('obsgoal="%s"' % obsgoal)
        print('ProjectId="%s"' % pid)
        print("ra=%f  # deg" % ra)
        print("dec=%f # deg" % dec)
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
        print('bsbeam=%d,%d' % (bsbeam[0],bsbeam[1]))
        print('tau=%g' % tau)
        dv = dataverse(pid)
        if dv != None:
            for k in dv.keys():
                print('%s="%s"' % (k,dv[k]))
        else:
            print("# no dataverse info")
        print("config=%s__RSR__%s__%s" % (pid,obspgm,src))
        print("# </lmtinfo>")

    else:     # one line summary
        #import pdb; pdb.set_trace()
        #print(date_obs, obsnum, obsgoal, obspgm, src,  pid,         tint,   ra,    dec,   az,   el)
        if True:
            print("%-20s %7d  %-10s %-12s %-11s %-25s %-30s              0    %6.1f  %10.6f %10.6f  %5.1f %5.1f" %
              (date_obs, obsnum, instrument, obsgoal, obspgm+(('('+map_coord+')') if obspgm=='Map' else ''), src,  pid,         tint,   ra,    dec,   az,   el))

    # -end  rsr_summary()


def nc_find(obsnum, rawnc=False):
    """ find the raw NC file that belongs to an obsnum
    """
    if obsnum.isnumeric():
        # obsnum=int(obsnum)
        # globs = '%s/ifproc/ifproc_*_%06d_*.nc' % (data_lmt,obsnum)
        globs = '%s/ifproc/ifproc_*_*%s_*.nc' % (data_lmt,obsnum)
        #print("# GLOBS slr:",globs)
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        if len(fn) == 1:
            if rawnc:
                return fn[0]
            slr_summary(fn[0],rc=True)
            sys.exit(0)
        elif len(fn) > 0:
            print("# lmtinfo: Multiple ifproc? ",fn)
            sys.exit(1)

        # since no SLR found, try RSR
        globs = '%s/RedshiftChassis?/RedshiftChassis?_*%s*.nc'  % (data_lmt,obsnum)
        # print("# GLOBS rsr:" % globs)
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        if len(fn) > 0 and len(fn) < 5:
            if rawnc:
                return fn[0]
            rsr_summary(fn[0], rc=True)
            sys.exit(0)
        elif len(fn) > 4:
            print("# lmtinfo: Multiple RSR ",fn)
            sys.exit(1)

        # since no RSR found, give up
        print("# lmtinfo: No matching OBSNUM %s in %s" % (obsnum,data_lmt))
        sys.exit(1)

def find_newer(root,newer):
    """ find files newer than a given file using the unix find command
    """
    # add a '/' since it may be a symlink
    cmd = 'find %s/ -name \*.nc -newer %s' % (root,newer)
    pipe = os.popen(cmd,'r')
    lines = pipe.readlines()
    pipe.close()
    fn = []
    for line in lines:
        fn.append(line.strip())
    return fn
        
def rsr_find(newer=None):
    """
    find all RSR
    """
    if newer == None:
        chassis = 1
        if chassis < 0:
            globs = '%s/RedshiftChassis?/RedshiftChassis?_*.nc'  % (data_lmt)
        else:
            globs = '%s/RedshiftChassis%d/RedshiftChassis%d_*.nc'  % (data_lmt,chassis,chassis)            
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        print("# Found %d RSR with %s" % (len(fn),globs))
    else:
        fn = find_newer('RedshiftChassis1',newer)
        
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
            print("# %-20s %7s  failed for rsr %s" % (yyyymmdd,obsnum,f))                    

    
def seq_find(newer = None):
    """
    find all RSR
    """
    if newer == None:
        globs = '%s/ifproc/ifproc_*.nc' % data_lmt
        fn = glob.glob(globs)
        fn.sort(key=os.path.getmtime)
        print("# Found %d SLR with %s" % (len(fn),globs))
    else:
        fn = find_newer('ifproc',          newer)
        
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
            print("# %-20s %7s  failed for slr %s" % (yyyymmdd,obsnum,f))


# ==================================================================================================================            


if len(sys.argv) == 2:

    # build
    if sys.argv[1] == "build":
        print("Rebuilding $DATA_LMT/data_lmt.log")
        build()
        sys.exit(0)

    # last
    if sys.argv[1] == "last":
        print(last())
        sys.exit(0)

    # new
    if sys.argv[1] == "old":
        print("Updating $DATA_LMT/data_lmt.log")
        new()
        sys.exit(0)

    # single obsnum (for SLR or RSR)
    obsnum = sys.argv[1]
    nc_find(obsnum)
    


    # replacement for $DATA_LMT
    if sys.argv[1] == 'data':
        print(header)
    elif sys.argv[1] == 'grep' or sys.argv[1] == 'find':
        print(header)
        grep([])
    elif sys.argv[1] == 'grepw' or sys.argv[1] == 'findw':
        print(header)
        grep([],"-w")
    elif os.path.isdir(sys.argv[1]):
        data_lmt = sys.argv[1]
        print(header)
    else:
        print("no more valid options")
        sys.exit(1)

    rsr_find()
    seq_find()
    sys.exit(0)
               
elif len(sys.argv) == 3:

    # special cases:
    if sys.argv[1] == "grep" or sys.argv[1] == "find":
        print(header)
        grep(sys.argv[2:])
        sys.exit(0)
        
    if sys.argv[1] == "grepw" or sys.argv[1] == "findw":
        print(header)
        grep(sys.argv[2:],"-w")
        sys.exit(0)

    # newer than an obsnum for incremental build
    if sys.argv[1] == "new":
        obsnum = sys.argv[2]
        rawnc = nc_find(obsnum, rawnc=True)
        rsr_find(newer=rawnc)
        seq_find(newer=rawnc)        
        sys.exit(0)

    if sys.argv[1] == "lmtot":
        obsnum = sys.argv[2]
        cmd = 'wget -q http://taps.lmtgtm.org/cgi-bin/script/x.cgi?-obsnum=%s -O -' % obsnum
        os.system(cmd)
        sys.exit(0)
    # there should be no more options now
    print("Illegal option 3",sys.argv[1])


else:
    # grep allows more terms
    if sys.argv[1] == "grep":
        print(header)
        grep(sys.argv[2:])
        sys.exit(0)
    if sys.argv[1] == "grepw":
        print(header)
        grep(sys.argv[2:],"-w")
        sys.exit(0)

    # otherwise illegal options, so give help
    
    print("Usage : %s [path] obsnum" % sys.argv[0])
