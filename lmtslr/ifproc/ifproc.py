"""
Module for reading and operating on IFPROC files 

classes: IFProc, IFProcData, IFProcCal
function: 
uses: numpy, netCDF4, os, fnmatch, RSRUtilities.TempSens
author: FPS
date: May 2018
changes: 
KS changes for online system
FPS added automatic calibration step
python 3
PJT changes to handle any of AE/RD/LL coordinates
"""

import numpy as np
import datetime
import netCDF4
import os
import fnmatch
import ast
from scipy.signal import detrend
from scipy import interpolate
import traceback

from lmtslr.ifproc.RSRUtilities import TempSens # move into utils folder?
from lmtslr.utils.ifproc_file_utils import lookup_ifproc_file
"""
def lookup_ifproc_file(obsnum,path='/data_lmt/ifproc/'):
    filename = ''
    for file in os.listdir(path):
        if fnmatch.fnmatch(file,'*_%06d_*.nc'%(obsnum)):
            print('found %s'%(file))
            filename = path+file
    if filename == '':
        print('lookup_ifproc_file: no file for obsnum ', obsnum)
        if 'lmttpm' not in path:
            print('look in lmttpm')
            return lookup_ifproc_file(obsnum,path='/data_lmt/lmttpm/')
    return(filename)
"""

def MapCoord(map_coord, obsgoal, source_coord_sys, obsnum):
    """   translate ascii MapCoord to an index  (0,1,2)
    
          Exceptions are for Science (only Ra/Dec or L/B are returned)
          and for Pointing (Az/El is returned) obsgoal's.
          The rest returns whatever map_coord implies.
    """
    import sys
    print('MapCoord:',map_coord, obsgoal, source_coord_sys, obsnum)
    if obsgoal == "Science":
        if source_coord_sys == 2:
            return 2
        if source_coord_sys != 1:
            print("Warning: science data, assuming Ra/Dec")
        return 1
    if obsgoal == "Pointing":
        print("Warning: pointing data, assumed Az/El")
        return 0
    if "Az"  in map_coord: return 0
    if "El"  in map_coord: return 0
    if "Ra"  in map_coord: return 1
    if "Dec" in map_coord: return 1
    if "L"   in map_coord: return 2
    if "B"   in map_coord: return 2
    # illegal MapCoord
    return -1
    

class IFProcQuick():
    """
    Base class for reading quick information from IFPROC
    """
    def __init__(self, filename, instrument='Sequoia'):
        """
        Constructor for IFProcQuick class.
        Args:
            filename (str): name of target NetCDF data file
            instrument (str): target instrument (default is Sequoia)
        Returns:
            none
        """
        self.filename = filename
        if os.path.isfile(self.filename):
            self.nc = netCDF4.Dataset(self.filename)
            self.obspgm = b''.join(self.nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
            self.obsgoal = b''.join(self.nc.variables['Header.Dcs.ObsGoal'][:]).decode().strip()
            self.obsnum = self.nc.variables['Header.Dcs.ObsNum'][0]
            self.receiver = b''.join(self.nc.variables['Header.Dcs.Receiver'][:]).decode().strip()
            self.nc.close()
        else:
            print('IFProcQuick: file \'%s\' is not found'%(self.filename))

class IFProc():
    """
    Base class for reading generic header information from IFPROC.
    """
    def __init__(self, filename, instrument='Sequoia'):
        """
        Constructor for IfProc class.
        Args:
            filename (str): name of target NetCDF data file
            instrument (str): target instrument (default is Sequoia)
        Returns:
            none
        """
        self.filename = filename
        if os.path.isfile(self.filename):
            self.nc = netCDF4.Dataset(self.filename)

            # header information
            self.source = b''.join(self.nc.variables['Header.Source.SourceName'][:]).decode().strip()
            self.source_coord_sys = self.nc.variables['Header.Source.CoordSys'][0]
            self.vlsr = self.nc.variables['Header.Source.Velocity'][0]

            date_obs = self.nc.variables['Data.TelescopeBackend.TelTime'][0].tolist()
            self.date_obs = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%dT%H:%M:%S')
            self.date_ymd = datetime.datetime.fromtimestamp(date_obs).strftime('%Y-%m-%d')
            print("%s obs-start %s" % (self.date_obs, self.filename))

            date_obs2 = self.nc.variables['Data.TelescopeBackend.TelTime'][-1:].tolist()[0]
            date_obs2 = datetime.datetime.fromtimestamp(date_obs2).strftime('%Y-%m-%dT%H:%M:%S')
            print("%s obs-stop  %s" % (date_obs2, self.filename))
            delta1 = self.nc.variables['Data.TelescopeBackend.TelTime'][-1:].tolist()[0] - \
                     self.nc.variables['Data.TelescopeBackend.TelTime'][0].tolist()
            print("delta1 %6.1f sec" % delta1)
            if True:
                # report on the stats of BufPos
                delta2 = 0
                on = self.nc.variables['Header.Dcs.ObsNum'][0]
                tt = self.nc.variables['Data.TelescopeBackend.TelTime'][:]
                bp = self.nc.variables['Data.TelescopeBackend.BufPos'][:]            
                nt  = len(tt)
                tt0 = tt[0]
                bp0 = bp[0]
                for i in range(1,nt):
                    if bp[i] != bp0:
                        print("BufPos %3d  %6.1f sec %s" % (bp0, tt[i] - tt0, on))
                        delta2 = delta2 + tt[i] - tt0
                        bp0 = bp[i]
                        tt0 = tt[i]
                print("BufPos %3d  %6.1f sec %s" % (bp[-1], tt[-1] - tt0, on))
                delta2 = delta2 + tt[-1] - tt0
                print("delta2 %6.1f sec" % delta2)
            try:
                # only in newer data (@todo after when?)
                self.dumptime = self.nc.variables['Header.SpecBackend.DumpTime'][0]
            except:
                print("Old data, assuming Header.SpecBackend.DumpTime = 0.1")
                self.dumptime = 0.1
                
            self.source_RA = self.nc.variables['Header.Source.Ra'][0]
            self.source_Dec = self.nc.variables['Header.Source.Dec'][0]
            # PJT: mapcoords add L,B
            self.source_L = self.nc.variables['Header.Source.L'][0]
            self.source_B = self.nc.variables['Header.Source.B'][0]            
            self.obspgm = b''.join(self.nc.variables['Header.Dcs.ObsPgm'][:]).decode().strip()
            self.obsgoal = b''.join(self.nc.variables['Header.Dcs.ObsGoal'][:]).decode().strip()
            if 'ifproc' in filename:
                self.calobsnum = self.nc.variables['Header.IfProc.CalObsNum'][0]
            elif 'lmttpm' in filename:
                self.calobsnum = self.nc.variables['Header.LmtTpm.CalObsNum'][0]
            else:
                self.calobsnum = 0
                    
            self.obsnum = self.nc.variables['Header.Dcs.ObsNum'][0]
            self.utdate = self.nc.variables['Header.TimePlace.UTDate'][0]
            self.ut1_h = self.nc.variables['Header.TimePlace.UT1'][0] / 2 / np.pi * 24
            self.azim = self.nc.variables['Header.Telescope.AzDesPos'][0] * 180 / np.pi
            self.elev = self.nc.variables['Header.Telescope.ElDesPos'][0] * 180 / np.pi
            self.m1ZernikeC0 = self.nc.variables['Header.M1.ZernikeC'][0]

            key = 'Header.M1.ReqPos'
            if key in self.nc.variables:
                self.m1ReqPos = self.nc.variables[key][:]
            else:
                self.m1ReqPos = np.zeros(720)

            self.m2x = self.nc.variables['Header.M2.XReq'][0]
            self.m2y = self.nc.variables['Header.M2.YReq'][0]
            self.m2z = self.nc.variables['Header.M2.ZReq'][0]
            self.m2xPcor = self.nc.variables['Header.M2.XPcor'][0]
            self.m2yPcor = self.nc.variables['Header.M2.YPcor'][0]
            self.m2zPcor = self.nc.variables['Header.M2.ZPcor'][0]

            # rotation about X
            self.m2tip = self.nc.variables['Header.M2.TipCmd'][0]
            # rotation about Y
            self.m2tilt = self.nc.variables['Header.M2.TiltCmd'][0]
            self.zc0 = self.nc.variables['Header.M1.ZernikeC'][0]
            self.zc_enabled = self.nc.variables['Header.M1.ZernikeEnabled'][0]

            # sometimes the Receiver designation is wrong; check and warn but don't stop
            self.receiver = b''.join(self.nc.variables['Header.Dcs.Receiver'][:]).decode().strip()
            try:
                print('before read npix')
                self.npix = int(self.nc.variables['Header.' + self.receiver + '.NumPixels'][0])
                print('from pixels npix =', self.npix)
                if 'ifproc' in filename:
                    if 'Data.IfProc.BasebandLevel_ylen' in self.nc.dimensions:
                        self.npix = len(self.nc.dimensions['Data.IfProc.BasebandLevel_ylen'])
                    else:
                        self.npix = len(self.nc.dimensions['Data.IfProc.BasebandLevel_xlen'])
                    if 'Data.IfProc.DetectorLevel_ylen' in self.nc.dimensions:
                        self.npix += len(self.nc.dimensions['Data.IfProc.DetectorLevel_ylen'])
                elif 'lmttpm' in filename:
                    self.npix = len(self.nc.dimensions['Data.LmtTpm.Signal_xlen'])
                    if True or self.receiver == 'B4r':
                        self.npix = 1
                else:
                        self.npix = 1
                print('from xlen npix =', self.npix)
                self.tracking_beam = self.nc.variables['Header.' + self.receiver + '.BeamSelected'][0]
                if self.tracking_beam != -1:
                    print('TRACKING ' + self.receiver + ' PIXEL ', self.tracking_beam)
            except Exception as e:
                print(e)
                print('WARNING: NOT AN HETERODYNE FILE')
                self.tracking_beam = -1

            # sideband information
            self.sideband = np.zeros(2)
            try:
                self.sideband[0] = self.nc.variables[
                    'Header.' + self.receiver + '.SideBand1Lo'][0]
                self.sideband[1] = self.nc.variables[
                    'Header.' + self.receiver + '.SideBand1Lo'][1]
            except Exception as e:
                self.sideband[0] = 0
                self.sideband[1] = 0
                print(e)
                print('WARNING: NOT AN HETERODYNE FILE')

            # Pointing Variables
            self.modrev = self.nc.variables['Header.PointModel.ModRev'][0]
            self.az_user = self.nc.variables['Header.PointModel.AzUserOff'][0] * 206264.8
            self.el_user = self.nc.variables['Header.PointModel.ElUserOff'][0] * 206264.8
            self.az_paddle = self.nc.variables['Header.PointModel.AzPaddleOff'][0] * 206264.8
            self.el_paddle = self.nc.variables['Header.PointModel.ElPaddleOff'][0] * 206264.8
            self.az_total = self.nc.variables['Header.PointModel.AzTotalCor'][0] * 206264.8
            self.el_total = self.nc.variables['Header.PointModel.ElTotalCor'][0] * 206264.8
            self.az_receiver = self.nc.variables['Header.PointModel.AzReceiverOff'][0] * 206264.8
            self.el_receiver = self.nc.variables['Header.PointModel.ElReceiverOff'][0] * 206264.8
            self.az_m2 = self.nc.variables['Header.PointModel.AzM2Cor'][0] * 206264.8
            self.el_m2 = self.nc.variables['Header.PointModel.ElM2Cor'][0] * 206264.8
            self.az_point_model_cor = self.nc.variables['Header.PointModel.AzPointModelCor'][0] * 206264.8
            self.el_point_model_cor = self.nc.variables['Header.PointModel.ElPointModelCor'][0] * 206264.8

            # TILTMETER Information                                                                                  
            self.tilt0_x = self.nc.variables['Header.Tiltmeter_0_.TiltX'][0] * 206264.8
            self.tilt0_y = self.nc.variables['Header.Tiltmeter_0_.TiltY'][0] * 206264.8
            self.tilt1_x = self.nc.variables['Header.Tiltmeter_1_.TiltX'][0] * 206264.8
            self.tilt1_y = self.nc.variables['Header.Tiltmeter_1_.TiltY'][0] * 206264.8

            # TEMPERATURE SENSOR Information
            self.T = TempSens(self.nc.variables['Header.TempSens.TempSens'][:] / 100)

            # WEATHER
            self.weather_temperature = self.nc.variables['Header.Weather.Temperature'][0]
            

            # map parameters 
            try:
                self.hpbw = self.nc.variables['Header.Map.HPBW'][0] * 206264.8
                self.xlength = self.nc.variables['Header.Map.XLength'][0] * 206264.8
                self.ylength = self.nc.variables['Header.Map.YLength'][0] * 206264.8
                self.xstep = self.nc.variables['Header.Map.XStep'][0]
                self.ystep = self.nc.variables['Header.Map.YStep'][0]
                self.xoffset = self.nc.variables['Header.Map.XOffset'][0]
                self.yoffset = self.nc.variables['Header.Map.YOffset'][0]
                self.rows = self.nc.variables['Header.Map.RowsPerScan'][0]
                # check the coordinate system AzEl = 0; RaDec = 1; LatLon = 2; default =0
                test_map_coord = b''.join(self.nc.variables['Header.Map.MapCoord'][:]).decode().strip()
                self.map_coord = MapCoord(test_map_coord, self.obsgoal, self.source_coord_sys, self.obsnum)
                if self.map_coord < 0:
                    print("Warning: unknown map_coord ",test_map_coord)
                    self.map_coord = 0

                self.map_motion = b''.join(self.nc.variables['Header.Map.MapMotion'][:]).decode().strip()
                self.scanang = self.nc.variables['Header.Map.ScanAngle'][0] * 206264.8/3600.
                print('Map Parameters: %s %s'%(test_map_coord, self.map_motion))
                print('HPBW=%5.1f XLength=%8.1f YLength=%8.1f XStep=%6.2f YStep=%6.2f ScanAngle=%6.2f'
                      %(self.hpbw, self.xlength, self.ylength, self.xstep, self.ystep, self.scanang))
            except Exception as e:
                print(e)
                self.map_motion = None
                print('%s does not have map parameters'%(self.filename))

            # bs parameters 
            try:
                self.bs_beams = self.nc.variables['Header.Bs.Beam'][:]
            except:
                self.bs_beams = []
                print('%s does not have bs parameters'%(self.filename))
                
            # Spectral Information
            self.velocity = self.nc.variables['Header.Source.Velocity'][0]
            self.velocity_system = self.nc.variables['Header.Source.VelSys'][0]

            try:
                self.line_list = ast.literal_eval(str(netCDF4.chartostring(
                    self.nc.variables['Header.Source.LineList'][:])
                    ).decode().strip())
                self.baseline_list = ast.literal_eval(str(
                    netCDF4.chartostring(self.nc.variables[
                        'Header.Source.BaselineList'][:])).decode().strip())
            except Exception as e:
                self.line_list = []
                self.baseline_list = []
            try:
                self.line_rest_frequency = self.nc.variables['Header.' + self.receiver + '.LineFreq'][0:2]
                print("PJT: line_rest_freq ", self.line_rest_frequency)
                # @todo   reset self.line_rest_frequency[board]
                self.doppler_track = self.nc.variables['Header.' + self.receiver + '.DopplerTrack'][0]
                self.observatory_velocity = self.nc.variables['Header.Sky.ObsVel'][0]
                self.barycenter_velocity = self.nc.variables['Header.Sky.BaryVel'][0]
                self.sky_frequency = self.nc.variables['Header.' + self.receiver + '.SkyFreq'][0:2]
                self.lo_1_frequency = self.nc.variables['Header.' + self.receiver + '.Lo1Freq'][0]
                self.lo_2_frequency = self.nc.variables['Header.' + self.receiver + '.Lo2Freq'][0:2]
                self.if_1_frequency = self.nc.variables['Header.' + self.receiver + '.If1Freq'][0:2]
                self.if_2_frequency = self.nc.variables['Header.' + self.receiver + '.If2Freq'][0:2]
                self.synthesizer_harmonic = self.nc.variables['Header.' + self.receiver + '.SynthHarm'][0:2]
                self.synthesizer_frequency = self.nc.variables['Header.' + self.receiver + '.SynthFreq'][0:2]
                self.sideband_1_lo_type = self.nc.variables['Header.' + self.receiver + '.SideBand1LoType'][0:2]
                self.sideband_2_lo_type = self.nc.variables['Header.' + self.receiver + '.SideBand2LoType'][0:2]
                self.sideband_1_lo = self.nc.variables['Header.' + self.receiver + '.SideBand1Lo'][0:2]
                self.sideband_2_lo = self.nc.variables['Header.' + self.receiver + '.SideBand2Lo'][0:2]
                self.velocity_definition = self.nc.variables['Header.' + self.receiver + '.VelocityDefinition'][0]
                self.frequency_offset = self.nc.variables['Header.' + self.receiver + '.LineOffset'][0:2]
                self.line_redshift = self.nc.variables['Header.' + self.receiver + '.LineRedshift'][0:2]
                
            except Exception as e:
                self.line_rest_frequency = 0
                self.doppler_track = 0
                print(e)
                print('WARNING: NOT AN HETERODYNE FILE')
        else:
            print('ifproc: file "%s" is not found'%(self.filename))

    def close_nc(self):
        """
        Closes open NetCDF file.
        Args:
            none
        Returns:
            none
        """
        self.nc.close()

    def process_chopped_encoder(self, chop, chan,
                                thresholds=[[[15,45,181],[105,135]],
                                            [[15,45,181],[105,135]],
                                            [[15,45,181],[105,135]],
                                            [[15,45,181],[105,135]],
                                            [[0,45,155],[65,135]],
                                            [[0,45,155],[65,135]]]):
        # create array of indices for main and ref based on chop array
        ang = (chop/8000*360)%180

        midx = np.where(np.logical_or(np.logical_and(ang > thresholds[chan][0][0], ang < thresholds[chan][0][1]),np.logical_and(ang>thresholds[chan][0][2],ang<=180)))[0]
        ridx = np.where(np.logical_and(ang > thresholds[chan][1][0], ang < thresholds[chan][1][1]))[0]
        return midx, ridx

    def process_chopped_signal(self, bb_level, chop, chop_option, ww=25,
                               thresholds=[[[15,45,181],[105,135]],
                                           [[15,45,181],[105,135]],
                                           [[15,45,181],[105,135]],
                                           [[15,45,181],[105,135]],
                                           [[0,45,155],[65,135]],
                                           [[0,45,155],[65,135]]]):


        '''
        gated chopper signal processor
        inputs:
             bb_level is npts by nchannels 2D array with baseband if data samples
             chop is chopper wheel position 0 to 8000 corresponds to 0 to 360 degrees
             ww defines smoothing window of ww points.  The total smoothing
              window must span at least one chop cycle
             thresholds are positions for including data points in the main and ref
               thresholds[0] elements give main limits in degrees from 0 to 180
                             data are included if between thresholds[0][0] and thresholds[0][1]
                             OR if greater than thresholds[0][2]
               thresholds[1] elements give reference limits
                             data are included if between thresholds[1][0] and thresholds[1][1]
        output:
             result is a 2D array with npts samples to match input arrays and nchannels.
        '''

        # look at the shape of the arrays to determine if super sampled and reshape
        s1 = np.shape(bb_level)
        s2 = np.shape(chop)
        if len(s1) == 3 and len(s2) == 2:
            if s1[0] != s2[0] or s1[1] != s2[1]:
                return None
            bb_level = bb_level.reshape(s1[0]*s1[1], s1[2])
            chop = chop.reshape(s2[0]*s2[1])
            super_sample = s1[1]
            ww = ww*super_sample
        else:
            super_sample = 1

        window = int(int(ww-1)/2)
        print('window =', window)
        npts = len(chop)
        nchannels = np.shape(bb_level)[-1]

        if chop_option == 8 or chop_option == 16:
            print(' chopping')
            
            result = np.zeros((npts,nchannels))
            
            # define the smoothing window
            ww = 2*window+1

            for i in range(nchannels):
                # find indices where encoder value are within a range
                midx, ridx = self.process_chopped_encoder(chop, i, thresholds=thresholds)

                msig = np.zeros(npts)
                msig[midx] = 1
                rsig = np.zeros(npts)
                rsig[ridx] = 1

                channel_level = bb_level[:,i] # gets rid of "masked array

                # create a rolling sum of the main points
                msum = np.cumsum(np.insert(msig*channel_level,0,0))
                mrollsum = msum[ww:]-msum[:-ww]

                # to do this accurately we also need a rolling sum for normalization
                mnorm = np.cumsum(np.insert(msig,0,0))
                mrollnorm = mnorm[ww:]-mnorm[:-ww]

                # same procedure for reference points
                rsum = np.cumsum(np.insert(rsig*channel_level,0,0))
                rrollsum = rsum[ww:]-rsum[:-ww]

                # same normalization procedure for reference points
                rnorm = np.cumsum(np.insert(rsig,0,0))
                rrollnorm = rnorm[ww:]-rnorm[:-ww]

                # now compute difference between main and ref for all points 
                result[window:npts-window,i] = mrollsum/mrollnorm - rrollsum/rrollnorm
                result[:window,i] = result[window,i]*np.ones(window)
                result[npts-window:,i] = result[npts-window-1,i]*np.ones(window)

        else:
            print(' not chopping')
            result = bb_level

        # average the arrays back down if super sampled
        if super_sample > 1:
            result = np.mean(result.reshape(-1, super_sample, nchannels), axis=1)
        return(result)

class IFProcData(IFProc):
    """
    Class for reading IFPROC data file to obtain time sequence of total
    power measurements.
    """
    def __init__(self, filename, npix=16):
        """
        Constructor for IFProcData class.
        Args:
            filename (str): name of target NetCDF data file
            npix (int): number of pixels or beams (default is 16)
        Returns:
            none
        """
        self.npix = npix
        IFProc.__init__(self, filename)

        if not hasattr(self, "obspgm"):
            return

        # identify the obspgm
        self.map_coord = -1
        # PJT  ->  'Az'

        if self.obspgm == 'Bs':
            print('%d is a Bs observation'%(self.obsnum))
            # bs parameters
            try:
                self.nrepeats = self.nc.variables['Header.Bs.NumRepeats'][0]
                self.nscans = self.nc.variables['Header.Bs.NumScans'][0]
                self.tsamp = self.nc.variables['Header.Bs.TSamp'][0]
                self.nsamp = self.nc.variables['Header.Bs.NSamp'][0]
                self.bs_pixel_ids = self.nc.variables['Header.Bs.Beam'][:]
            except:
                print('%s does not have Bs parameters'%(self.filename))

        elif self.obspgm == 'Ps':
            print('%d is a Ps observation'%(self.obsnum))
            # ps parameters
            param = ''
            try:
                param = 'Header.Ps.NumRepeats'
                self.nrepeats = self.nc.variables[param][0]
                param = 'Header.Ps.NumScans'
                self.nscans = self.nc.variables[param][0]
                param = 'Header.Ps.TMain'
                self.tmain = self.nc.variables[param][0]
                param = 'Header.Ps.TRef'
                self.tref = self.nc.variables[param][0]
                param = 'Header.Ps.NSamp'
                self.nsamp = self.nc.variables[param][0]
                param = 'Header.Ps.Mode'
                self.mode = ''.join(self.nc.variables[param][:]).strip()
                param = 'Header.Ps.RefSwitch'
                self.refswitch = ''.join(self.nc.variables[param][:]).strip()
            except:
                print('%s does not have Ps parameters %s'%(self.filename, 
                                                           param))

        elif self.obspgm == 'Map':
            print('%d is a Map observation'%(self.obsnum))
            # map parameters 
            try:
                self.hpbw = self.nc.variables['Header.Map.HPBW'][0] * 206264.8
                self.xlength = self.nc.variables[
                    'Header.Map.XLength'][0]*206264.8
                self.ylength = self.nc.variables[
                    'Header.Map.YLength'][0]*206264.8
                self.xstep = self.nc.variables['Header.Map.XStep'][0]
                self.ystep = self.nc.variables['Header.Map.YStep'][0]
                self.xoffset = self.nc.variables['Header.Map.XOffset'][0]
                self.yoffset = self.nc.variables['Header.Map.YOffset'][0]
                self.rows = self.nc.variables['Header.Map.RowsPerScan'][0]
                # check the coordinate system Az = 0; Ra = 1; default =0
                test_map_coord = b''.join(self.nc.variables['Header.Map.MapCoord'][:]).decode().strip()
                self.map_coord = MapCoord(test_map_coord, self.obsgoal, self.source_coord_sys, self.obsnum)
                if self.map_coord < 0:
                    print("Warning: unknown map_coord ",test_map_coord)
                    self.map_coord = 0

                self.map_motion = b''.join(self.nc.variables['Header.Map.MapMotion'][:]).decode().strip()
            except Exception as e:
                print('e1', e)
                print('%s does not have map parameters'%(self.filename))

        elif self.obspgm == 'Cal':
            print('WARNING: %d is a Cal observation'%(self.obsnum))

        elif self.obspgm == 'On':
            print('WARNING: %d is an On observation'%(self.obsnum))

        else:
            print('WARNING: ObsPgm type %s for Obsum %d is not identified'%(self.obspgm, self.obsnum))

        # data arrays
        self.time = self.nc.variables['Data.TelescopeBackend.TelTime'][:]
        self.bufpos = self.nc.variables['Data.TelescopeBackend.BufPos'][:]
        # AzEl map
        self.azmap = self.nc.variables['Data.TelescopeBackend.TelAzMap'][:]* 206264.8
        self.elmap = self.nc.variables['Data.TelescopeBackend.TelElMap'][:]* 206264.8
        self.parang = self.nc.variables['Data.TelescopeBackend.ActParAng'][:]
        try:
            self.galang = self.nc.variables['Data.TelescopeBackend.ActGalAng'][:]
        except:
            self.galang = np.zeros(len(self.parang))

        # RaDec map
        self.ramap = (self.nc.variables['Data.TelescopeBackend.SourceRaAct'][:] - self.source_RA) * np.cos(self.source_Dec) * 206264.8
        self.decmap = (self.nc.variables['Data.TelescopeBackend.SourceDecAct'][:] - self.source_Dec) * 206264.8

        # set the l/b map
        self.lmap = (self.nc.variables['Data.TelescopeBackend.SourceLAct'][:] - self.source_L) * np.cos(self.source_B) * 206264.8
        self.bmap = (self.nc.variables['Data.TelescopeBackend.SourceBAct'][:] - self.source_B) * 206264.8

        if self.map_coord == 1:
            self.xmap = self.ramap
            self.ymap = self.decmap
        elif self.map_coord == 2:
            self.xmap = self.lmap
            self.ymap = self.bmap
        else:
            self.xmap = self.azmap
            self.ymap = self.elmap
            
        self.chop_option = 0
        if 'ifproc' in filename:
            self.bb_level = self.nc.variables['Data.IfProc.BasebandLevel'][:]
            try:
                print('get chop')
                self.chop = self.nc.variables['Data.Msip1mm.BeamChopperActPos'][:]
                self.chop_option = self.nc.variables['Header.Msip1mm.BeamChopperActState'][0]
                self.level = self.process_chopped_signal(self.bb_level, self.chop, self.chop_option)
                if 'Data.IfProc.DetectorLevel' in self.nc.variables:
                    self.detector_level = self.nc.variables['Data.IfProc.DetectorLevel'][:]
                    self.detector = self.process_chopped_signal(self.detector_level, self.chop, self.chop_option)
                    self.bb_level = np.concatenate((self.bb_level, self.detector_level), axis=2)
                    self.level = np.concatenate((self.level, self.detector), axis=1)
            except Exception as e:
                print(e)
                traceback.print_exc()
                print(' no chop')
                self.level = self.bb_level
        elif 'lmttpm' in filename:
            self.level = detrend(self.nc.variables['Data.LmtTpm.Signal'][:], axis=0)
        else:
            self.level = np.zeros(0)
            
        self.nsamp = len(self.level)

        # initialize calibration flag
        self.cal_flag = False

        self.close_nc()

    def calibrate_data(self, CAL):
        """
        Calibrates data using constants from CAL object, or replicates 
        if CAL object is absent.
        Args:
            CAL (object): Cal object
        Returns:
            none
        """
        self.caldata = np.zeros((self.npix, self.nsamp))
        self.bias = np.zeros(self.npix)
        self.tsys = np.zeros(self.npix)
        for ipix in range(self.npix):
            self.caldata[ipix, :] = (self.level[:, ipix] - CAL.calcons[ipix, 1]) / CAL.calcons[ipix, 0]
            self.bias[ipix] = np.median(self.caldata[ipix, :])
            self.tsys[ipix] = CAL.tsys[ipix]
        self.cal_flag = True

    def dont_calibrate_data(self):# why can't this be an option in calibrate_data with the CAL object absent?
        """
        Sets data. Sets tsys = 0.
        Args:
            none
        Returns:
            none
        """
        self.caldata = np.zeros((self.npix, self.nsamp))
        self.bias = np.zeros(self.npix)
        self.tsys = np.zeros(self.npix)
        for ipix in range(self.npix):
            self.caldata[ipix, :] = self.level[:, ipix]
            self.bias[ipix] = np.median(self.caldata[ipix, :])
            # set to zero for the case of no calibration
            self.tsys[ipix] = 0
        self.cal_flag = False

    def find_map_pixel_index(self, ipixel):
        """
        Returns the target pixel index.
        Args:
            ipixel: target pixel
        Returns:
            ipixel: target pixel
        """
        return(ipixel)

    def create_map_data(self):
        """
        Sets the map data.
        Args:
            none
        Returns:
            none
        """
        self.map_data = []
        self.map_x = []
        self.map_y = []
        self.map_az = []
        self.map_el = []
        self.map_ra = []
        self.map_dec = []
        self.map_l = []
        self.map_b = []
        self.map_n = []
        self.map_p = []
        self.map_g = []
        print('PJT map_coord',self.map_coord)
        for i in range(self.npix):
            self.map_x.append(self.xmap)
            self.map_y.append(self.ymap)
            self.map_az.append(self.azmap)
            self.map_el.append(self.elmap)
            self.map_ra.append(self.ramap)
            self.map_dec.append(self.decmap)
            self.map_l.append(self.lmap)
            self.map_b.append(self.bmap)
            self.map_p.append(self.parang)
            self.map_g.append(self.galang)
            self.map_n.append(self.nsamp)
            self.map_data.append(self.caldata[i,:] - self.bias[i])
        self.map_x = np.array(self.map_x)
        self.map_y = np.array(self.map_y)
        self.map_az = np.array(self.map_az)
        self.map_el = np.array(self.map_el)
        self.map_ra = np.array(self.map_ra)
        self.map_dec = np.array(self.map_dec)
        self.map_l = np.array(self.map_l)
        self.map_b = np.array(self.map_b)
        self.map_p = np.array(self.map_p)
        self.map_g = np.array(self.map_g)
        self.map_n = np.array(self.map_n)
        self.map_data = np.array(self.map_data)

class IFProcCal(IFProc):
    """
    Class for reading IFPROC calibration file, which contains a 
    sequence of observations on Hot and Sky.
    """
    def __init__(self, filename, npix=16):
        """
        Constructor for IFProcCal class.
        Args:
            filename (str): name of target NetCDF data file
            npix (int): number of pixels (default is 16)
        Returns:
            none
        """
        self.npix = npix
        IFProc.__init__(self,filename)

        # check observation program type
        self.map_coord = -1
        if self.obspgm == 'Cal':
            print('%d is a Cal observation'%(self.obsnum))
        else:
            print('WARNING: %d is NOT a Cal observation : %s'%(self.obsnum, 
                                                               self.obspgm))

        # data arrays
        self.time = self.nc.variables['Data.TelescopeBackend.TelTime'][:]
        self.azmap = self.nc.variables['Data.TelescopeBackend.TelAzMap'][:]
        self.elmap = self.nc.variables['Data.TelescopeBackend.TelElMap'][:]
        self.xmap = self.azmap
        self.ymap = self.elmap
        self.ramap =  np.zeros(len(self.azmap))
        self.decmap =  np.zeros(len(self.azmap))
        self.lmap =  np.zeros(len(self.azmap))
        self.bmap =  np.zeros(len(self.azmap))
        self.parang = np.zeros(len(self.azmap))
        self.galang = np.zeros(len(self.azmap))
        self.bufpos = self.nc.variables['Data.TelescopeBackend.BufPos'][:]
        self.chop_option = 0
        if 'ifproc' in filename:
            self.bb_level = self.nc.variables['Data.IfProc.BasebandLevel'][:]
            try:
                print('get chop cal')
                self.chop = self.nc.variables['Data.Msip1mm.BeamChopperActPos'][:]
                self.chop_option = self.nc.variables['Header.Msip1mm.BeamChopperActState'][0]
                self.level = self.process_chopped_signal(self.bb_level, self.chop, self.chop_option)
                if 'Data.IfProc.DetectorLevel' in self.nc.variables:
                    self.detector_level = self.nc.variables['Data.IfProc.DetectorLevel'][:]
                    self.detector = self.process_chopped_signal(self.detector_level, self.chop, self.chop_option)
                    self.bb_level = np.concatenate((self.bb_level, self.detector_level), axis=2)
                    self.level = np.concatenate((self.level, self.detector), axis=1)
            except Exception as e:
                print(e)
                traceback.print_exc()
                print(' no chop cal')
                self.level = self.bb_level
        elif 'lmttpm' in filename:
            self.level = detrend(self.nc.variables['Data.LmtTpm.Signal'][:], axis=0)
        else:
            self.level = np.zeros(0)
        self.nsamp = len(self.level)
        self.tamb = 280.
        self.receiver = b''.join(self.nc.variables['Header.Dcs.Receiver'][:]).decode().strip()
        try:
            self.blank_level = self.nc.variables['Header.' + self.receiver + '.BlankLevel'][0]
        except:
            if self.receiver == 'B4r':
                self.blank_level = -8.9
            else:
                self.blank_level = 0
        #self.tamb = self.nc.variables['Header.'+self.receiver+'.LoadAmbientTemp'][0]

        self.close_nc()

    def compute_calcons(self):
        """
        Computes the calibration constants.
        Args:
            none
        Returns:
            none
        """
        hot_list = np.where(self.bufpos == 3)
        sky_list = np.where(self.bufpos == 2)
        self.calcons = np.zeros((self.npix, 2))
        if self.chop_option == 8 or self.chop_option == 16:
            bb_level = self.bb_level
            if len(np.shape(bb_level)) == 3:
                bb_level = np.mean(bb_level, axis=1)
        for ipix in range(self.npix):
            if self.chop_option == 8 or self.chop_option == 16:
                self.calcons[ipix,0] = (np.median(bb_level[hot_list, ipix]) - 
                                        np.median(bb_level[sky_list, ipix])) / self.tamb
                self.calcons[ipix,1] = np.median(bb_level[sky_list, ipix])
            else:
                self.calcons[ipix,0] = (np.median(self.level[hot_list, ipix]) - 
                                        np.median(self.level[sky_list, ipix])) / self.tamb
                self.calcons[ipix,1] = np.median(self.level[sky_list, ipix])

    def compute_tsys(self):
        """
        Computes the system temperature, based on blank = 0V.
        Args:
            none
        Returns:
            none
        """
        hot_list = np.where(self.bufpos == 3)
        sky_list = np.where(self.bufpos == 2)
        self.tsys = np.zeros((self.npix))
        if self.chop_option == 8 or self.chop_option == 16:
            bb_level = self.bb_level
            chop = self.chop
            if len(np.shape(bb_level)) == 3:
                bb_level = np.mean(bb_level, axis=1)
            if len(np.shape(chop)) == 2:
                chop = np.mean(chop.reshape(-1, np.shape(chop)[1]), axis=1)
            chop_load = chop[hot_list]
        for ipix in range(self.npix):
            if self.chop_option == 8 or self.chop_option == 16:
                midx, ridx = self.process_chopped_encoder(chop_load, ipix)
                level = bb_level[:,ipix]
                level_load = level[hot_list]
                yhot = level_load[midx]
                vhot = np.median(yhot)
                vsky = np.median(level[sky_list])
                vzero = self.blank_level
                self.tsys[ipix] = self.tamb * (vsky - vzero) / (vhot - vsky)
            else:
                vsky = np.median(self.level[sky_list, ipix])
                vhot = np.median(self.level[hot_list, ipix])
                vzero = self.blank_level
                self.tsys[ipix] = self.tamb * (vsky - vzero) / (vhot - vsky)
        
