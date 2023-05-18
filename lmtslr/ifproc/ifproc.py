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
    """
    import sys
    print(map_coord, obsgoal, source_coord_sys, obsnum)
    if obsgoal == "Science":
        if source_coord_sys == 2:
            return 2
        return 1
    if obsgoal == "Pointing":
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
            print("%s begin %s" % (self.date_obs, self.filename))

            date_obs2 = self.nc.variables['Data.TelescopeBackend.TelTime'][-1:].tolist()[0]
            date_obs2 = datetime.datetime.fromtimestamp(date_obs2).strftime('%Y-%m-%dT%H:%M:%S')
            print("%s end   %s" % (date_obs2, self.filename))
            if True:
                # report on the stats of BufPos
                on = self.nc.variables['Header.Dcs.ObsNum'][0]
                tt = self.nc.variables['Data.TelescopeBackend.TelTime'][:]
                bp = self.nc.variables['Data.TelescopeBackend.BufPos'][:]            
                nt  = len(tt)
                tt0 = tt[0]
                bp0 = bp[0]
                for i in range(1,nt):
                    if bp[i] != bp0:
                        print("BufPos %3d  %6.1f sec %s" % (bp0, tt[i] - tt0, on))
                        bp0 = bp[i]
                        tt0 = tt[i]
                print("BufPos %3d  %6.1f sec %s" % (bp[-1], tt[-1] - tt0, on))
            try:
                # only in newer data
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

        self.tel_utc = 180/15/np.pi*self.nc.variables['Data.TelescopeBackend.TelUtc'][:][:]
        utdate = self.utdate
        utdate = datetime.datetime.utcfromtimestamp(self.time[0]).date()
        print(utdate, self.tel_utc[0])
        import dateutil.parser as dparser
        print([dparser.parse(str(utdate)+' '+str(datetime.timedelta(hours=self.tel_utc[0]))+' UTC', fuzzy=True) for i in range(1)])
        self.sky_time = np.array([dparser.parse(str(utdate)+' '+str(datetime.timedelta(hours=self.tel_utc[i]))+' UTC', fuzzy=True).timestamp() for i in range(len(self.tel_utc))])

        # RaDec map
        self.ramap_file = (self.nc.variables['Data.TelescopeBackend.SourceRaAct'][:] - self.source_RA) * np.cos(self.source_Dec) * 206264.8
        self.decmap_file = (self.nc.variables['Data.TelescopeBackend.SourceDecAct'][:] - self.source_Dec) * 206264.8

        # interpolate ra/dec based on tel time
        if True:
            sl = slice(0, len(self.time), 1)
            ra_file = self.nc.variables['Data.TelescopeBackend.SourceRaAct'][:]
            dec_file = self.nc.variables['Data.TelescopeBackend.SourceDecAct'][:]
            self.ra_interpolation_function = interpolate.interp1d(self.sky_time[sl],
                                                                  ra_file[sl],
                                                                  bounds_error=False,
                                                                  kind='previous',
                                                                  fill_value='extrapolate')
            self.dec_interpolation_function = interpolate.interp1d(self.sky_time[sl],
                                                                   dec_file[sl],
                                                                   bounds_error=False,
                                                                   kind='previous',
                                                                   fill_value='extrapolate')
            ra_interp = self.ra_interpolation_function(
                np.ma.getdata(self.time, subok=False))
            dec_interp = self.dec_interpolation_function(
                np.ma.getdata(self.time, subok=False))
            self.ramap_interp = (ra_interp - self.source_RA) * np.cos(self.source_Dec) * 206264.8
            self.decmap_interp = (dec_interp - self.source_Dec) * 206264.8
            if False:
                print('sky_time-time', self.sky_time-self.time)
                print('ramap_file', self.ramap_file)
                print('sky_time', self.sky_time)
                print('time', self.time)
                print('ra_interp', ra_interp)
                print('ramap_interp', self.ramap_interp)

        # compute ra/dec from astropy
        if True:
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            from astropy.time import Time
            from astroplan import Observer
            from astropy import coordinates as coord
            from pytz import timezone

            _lmt_info = {
                'instru': 'lmt',
                'name': 'LMT',
                'name_long': "Large Millimeter Telescope",
                'location': {
                    'lon': '-97d18m52.5s',
                    'lat': '+18d59m9.6s',
                    'height': 4640 << u.m,
                    },
                'timezone_local': 'America/Mexico_City'
            }

            lmt_location = coord.EarthLocation.from_geodetic(**_lmt_info['location'])
            """The local of LMT."""

            lmt_timezone_local = timezone(_lmt_info['timezone_local'])
            """The local time zone of LMT."""

            lmt_observer = Observer(
                name=_lmt_info['name_long'],
                location=lmt_location,
                timezone=lmt_timezone_local,
            )
            
            observer = lmt_observer

            tel_time = Time(self.nc['Data.TelescopeBackend.TelTime'][:], format='unix', scale='utc', location=lmt_observer.location)
            tel_az = self.nc['Data.TelescopeBackend.TelAzAct'][:] << u.rad
            tel_alt = self.nc['Data.TelescopeBackend.TelElAct'][:] << u.rad
            tel_az_cor = self.nc['Data.TelescopeBackend.TelAzCor'][:] << u.rad
            tel_alt_cor = self.nc['Data.TelescopeBackend.TelElCor'][:] << u.rad
            tel_az_tot = tel_az - (tel_az_cor) / np.cos(tel_alt)
            tel_alt_tot = tel_alt - (tel_alt_cor)
            altaz_frame = observer.altaz(time=tel_time)
            tel_icrs_astropy = SkyCoord(tel_az_tot, tel_alt_tot, frame=altaz_frame).transform_to('icrs')
            # update variables and save
            parang = observer.parallactic_angle(time=tel_time, target=tel_icrs_astropy)
            self.parang_astropy = parang.radian
            self.ramap_astropy = (tel_icrs_astropy.ra.radian - self.source_RA) * np.cos(self.source_Dec) * 206264.8
            self.decmap_astropy = (tel_icrs_astropy.dec.radian - self.source_Dec) * 206264.8

            az_ = tel_az_tot.value
            el_ = tel_alt_tot.value
            ra_ = tel_icrs_astropy.ra.radian
            dec_ = tel_icrs_astropy.dec.radian
            ut_ = tel_time.ut1
            lst_ = tel_time.sidereal_time('apparent')


        # set the ra/dec map
        #self.ramap = self.ramap_interp
        #self.decmap = self.decmap_interp
        self.ramap = self.ramap_file
        self.decmap = self.decmap_file

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
            
        if False:
            def stat_change(d, d_orig, unit, name):
                #dd = (d - d_orig).to_value(unit)
                dd = (d - d_orig)
                dd = dd[np.isfinite(dd)]
                print(f"{name} changed with diff ({unit}): min={dd.max()} max={dd.min()} mean={dd.mean()} std={np.std(dd)}")
            stat_change(self.parang_astropy, self.parang, u.deg, 'ActParAng') 
            stat_change(self.ramap_file, self.ramap_interp, u.arcsec, 'file-interp') 
            stat_change(self.decmap_file, self.decmap_interp, u.arcsec, 'file-interp')
            stat_change(self.ramap_file, self.ramap_astropy, u.arcsec, 'file-astropy') 
            stat_change(self.decmap_file, self.decmap_astropy, u.arcsec, 'file-astropy')
            stat_change(self.ramap_interp, self.ramap_astropy, u.arcsec, 'interp-astropy') 
            stat_change(self.decmap_interp, self.decmap_astropy, u.arcsec, 'interp-astropy')
                
        if False:
            import matplotlib.pyplot as pl
            sl = np.where(self.bufpos == 0)
            if False:
                # traces
                ax = pl.subplot()
                ax.plot(self.time[sl],self.ramap_file[sl], 'rx')
                ax.plot(self.time[sl],self.ramap_interp[sl], 'mx')
                ax.plot(self.time[sl],self.ramap_astropy[sl], 'yx')
            pl.show()
            import sys
            #sys.exit(0)


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
        
