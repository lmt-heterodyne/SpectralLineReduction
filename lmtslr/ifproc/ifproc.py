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
"""

import numpy as np
import netCDF4
import os
import fnmatch
import ast
from scipy.signal import detrend
from scipy import interpolate

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
            self.obspgm = b''.join(self.nc.variables['Header.Dcs.ObsPgm'][:]
                                  ).decode().strip()
            self.obsnum = self.nc.variables['Header.Dcs.ObsNum'][0]
            self.receiver = b''.join(self.nc.variables['Header.Dcs.Receiver'
                                                      ][:]).decode().strip()
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
            self.source = b''.join(self.nc.variables['Header.Source.SourceName'
                                                    ][:]).decode().strip()
            self.source_RA = self.nc.variables['Header.Source.Ra'][0]
            self.source_Dec = self.nc.variables['Header.Source.Dec'][0]
            self.obspgm = b''.join(self.nc.variables['Header.Dcs.ObsPgm'][:]
                                  ).decode().strip()
            if 'ifproc' in filename:
                self.calobsnum = self.nc.variables['Header.IfProc.CalObsNum'
                                                  ][0]
            elif 'lmttpm' in filename:
                self.calobsnum = self.nc.variables['Header.LmtTpm.CalObsNum'
                                                  ][0]
            else:
                self.calobsnum = 0
                    
            self.obsnum = self.nc.variables['Header.Dcs.ObsNum'][0]
            self.utdate = self.nc.variables['Header.TimePlace.UTDate'][0]
            self.ut1_h = self.nc.variables['Header.TimePlace.UT1'
                                          ][0] / 2 / np.pi * 24
            self.azim = self.nc.variables['Header.Telescope.AzDesPos'
                                         ][0] * 180 / np.pi
            self.elev = self.nc.variables['Header.Telescope.ElDesPos'
                                         ][0] * 180 / np.pi
            self.m1ZernikeC0 = self.nc.variables['Header.M1.ZernikeC'][0]

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

            # sometimes the Receiver designation is wrong; check and warn but \
            # don't stop
            self.receiver = b''.join(self.nc.variables['Header.Dcs.Receiver'
                                                      ][:]).decode().strip()
            try:
                print('before read npix')
                self.npix = int(self.nc.variables['Header.' + self.receiver + 
                                                  '.NumPixels'][0])
                print('from pixels npix =', self.npix)
                if 'ifproc' in filename:
                    self.npix = len(self.nc.dimensions[
                        'Data.IfProc.BasebandLevel_xlen'])
                elif 'lmttpm' in filename:
                    self.npix = len(self.nc.dimensions[
                        'Data.LmtTpm.Signal_xlen'])
                    if True or self.receiver == 'B4r':
                        self.npix = 1
                else:
                        self.npix = 1
                print('from xlen npix =', self.npix)
                self.tracking_beam = self.nc.variables[
                    'Header.' + self.receiver + '.BeamSelected'][0]
                if self.tracking_beam != -1:
                    print('TRACKING ' + self.receiver + ' PIXEL ', 
                          self.tracking_beam)
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
            self.az_user = self.nc.variables['Header.PointModel.AzUserOff'][0]\
                           * 206264.8
            self.el_user = self.nc.variables['Header.PointModel.ElUserOff'][0]\
                           * 206264.8
            self.az_paddle = self.nc.variables['Header.PointModel.AzPaddleOff'
                                              ][0] * 206264.8
            self.el_paddle = self.nc.variables['Header.PointModel.ElPaddleOff'
                                              ][0] * 206264.8
            self.az_total = self.nc.variables['Header.PointModel.AzTotalCor'
                                             ][0] * 206264.8
            self.el_total = self.nc.variables['Header.PointModel.ElTotalCor'
                                             ][0] * 206264.8
            self.az_receiver = self.nc.variables[
                'Header.PointModel.AzReceiverOff'][0] * 206264.8
            self.el_receiver = self.nc.variables[
                'Header.PointModel.ElReceiverOff'][0] * 206264.8
            self.el_m2 = self.nc.variables['Header.PointModel.ElM2Cor'][0]\
                         * 206264.8
            self.az_m2 = self.nc.variables['Header.PointModel.AzM2Cor'][0]\
                         * 206264.8

            # TILTMETER Information                                                                                  
            self.tilt0_x = self.nc.variables['Header.Tiltmeter_0_.TiltX'][0]\
                           * 206264.8
            self.tilt0_y = self.nc.variables['Header.Tiltmeter_0_.TiltY'][0]\
                           * 206264.8
            self.tilt1_x = self.nc.variables['Header.Tiltmeter_1_.TiltX'][0]\
                           * 206264.8
            self.tilt1_y = self.nc.variables['Header.Tiltmeter_1_.TiltY'][0]\
                           * 206264.8

            # TEMPERATURE SENSOR Information
            self.T = TempSens(self.nc.variables['Header.TempSens.TempSens'][:]
                              / 100)

            # map parameters 
            try:
                self.hpbw = self.nc.variables['Header.Map.HPBW'][0] * 206264.8
                self.xlength = self.nc.variables['Header.Map.XLength'][0]\
                               * 206264.8
                self.ylength = self.nc.variables['Header.Map.YLength'][0]\
                               * 206264.8
                self.xstep = self.nc.variables['Header.Map.XStep'][0]
                self.ystep = self.nc.variables['Header.Map.YStep'][0]
                self.rows = self.nc.variables['Header.Map.RowsPerScan'][0]
                # check the coordinate system Az = 0; Ra = 1; default =0
                test_map_coord = b''.join(self.nc.variables[
                    'Header.Map.MapCoord'][:]).decode().strip()
                if test_map_coord[0] == 'A':
                    self.map_coord = 0
                elif test_map_coord[0] == 'R':
                    self.map_coord = 1
                else:
                    self.map_coord = 0

                self.map_motion = b''.join(self.nc.variables[
                    'Header.Map.MapMotion'][:]).decode().strip()
                print('Map Parameters: %s %s'%(test_map_coord, self.map_motion
                                              ))
                print('HPBW=%5.1f XLength=%8.1f YLength=%8.1f XStep=%6.2f \
                       YStep=%6.2f'%(self.hpbw, self.xlength, self.ylength, 
                                     self.xstep, self.ystep))
            except:
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
            self.velocity_system = self.nc.variables['Header.Source.VelSys'
                                                    ][0]

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
                self.line_rest_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.LineFreq'][0:2]
                self.doppler_track = self.nc.variables[
                    'Header.' + self.receiver + '.DopplerTrack'][0]
                self.observatory_velocity = self.nc.variables[
                    'Header.Sky.ObsVel'][0]
                self.barycenter_velocity = self.nc.variables[
                    'Header.Sky.BaryVel'][0]
                self.sky_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.SkyFreq'][0:2]
                self.lo_1_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.Lo1Freq'][0]
                self.lo_2_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.Lo2Freq'][0:2]
                self.if_1_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.If1Freq'][0:2]
                self.if_2_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.If2Freq'][0:2]
                self.synthesizer_harmonic = self.nc.variables[
                    'Header.' + self.receiver + '.SynthHarm'][0:2]
                self.synthesizer_frequency = self.nc.variables[
                    'Header.' + self.receiver + '.SynthFreq'][0:2]
                self.sideband_1_lo_type = self.nc.variables[
                    'Header.' + self.receiver + '.SideBand1LoType'][0:2]
                self.sideband_2_lo_type = self.nc.variables[
                    'Header.' + self.receiver + '.SideBand2LoType'][0:2]
                self.sideband_1_lo = self.nc.variables[
                    'Header.' + self.receiver + '.SideBand1Lo'][0:2]
                self.sideband_2_lo = self.nc.variables[
                    'Header.' + self.receiver + '.SideBand2Lo'][0:2]
                self.velocity_definition = self.nc.variables[
                    'Header.' + self.receiver + '.VelocityDefinition'][0]
                self.frequency_offset = self.nc.variables[
                    'Header.' + self.receiver + '.LineOffset'][0:2]
                self.line_redshift = self.nc.variables[
                    'Header.' + self.receiver + '.LineRedshift'][0:2]
                
            except Exception as e:
                self.line_rest_frequency = 0
                self.doppler_track = 0
                print(e)
                print('WARNING: NOT AN HETERODYNE FILE')
        else:
            print('ifproc: file \'%s\' is not found'%(self.filename))

    def close_nc(self):
        """
        Closes open NetCDF file.
        Args:
            none
        Returns:
            none
        """
        self.nc.close()

    def process_chopped_signal(self, bb_level, chop, window=6, 
                               thresholds=[[15,45,181],[110,135]]):
        """
        Gated chopper signal processor.
        Args:
            bb_level (array): npts by nchannels 2D array with baseband 
                if data samples
            chop (int): chopper wheel position from 0 to 8000 
                corresponding to 0 to 360 degrees
            window (int): defines smoothing window of 2*window + 1 
                points. The total smoothing window must span at least 
                one chop cycle.
            thresholds (list): positions for including data points in 
                the main and ref.
                thresholds[0] elements give main limits in degrees from
                    0 to 180. data are included if between 
                    thresholds[0][0] and thresholds[0][1] OR if greater
                    than thresholds[0][2].
                thresholds[1] elements give reference limits. data are 
                    included if between thresholds[1][0] and 
                    thresholds[1][1].
        Returns:
            result (array): 2D array with npts samples to match input 
                arrays and nchannels.
        """
        npts = len(chop)
        nchannels = np.shape(bb_level)[1]
        # define the smoothing window
        ww = 2 * window + 1

        # create array of indices for main, ref, blank based on chop array
        ang = (chop / 8000 * 360)%180

        # find indices where cos of encoder value are within a range
        midx = np.where(np.logical_or(np.logical_and(ang > thresholds[0][0], 
            ang < thresholds[0][1]), np.logical_and(ang > thresholds[0][2], 
                                                    ang <= 180)))[0]
        ridx = np.where(np.logical_and(ang > thresholds[1][0], 
                                       ang < thresholds[1][1]))[0]
=======
    def process_chopped_signal(self, bb_level, chop, window=6, thresholds=[[15,45,181],[110,135]]):
        '''
        gated chopper signal processor
        inputs:
             bb_level is npts by nchannels 2D array with baseband if data samples
             chop is chopper wheel position 0 to 8000 corresponds to 0 to 360 degrees
             window defines smoothing window of 2*window + 1 points.  The total smoothing
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

        npts = len(chop)
        nchannels = np.shape(bb_level)[1]
        # define the smoothing window
        ww = 2*window+1

        # create array of indices for main, ref, blank based on chop array
        ang = (chop/8000*360)%180

        # find indices where cos of encoder value are within a range
        midx = np.where(np.logical_or(np.logical_and(ang > thresholds[0][0], ang < thresholds[0][1]),np.logical_and(ang>thresholds[0][2],ang<=180)))[0]
        ridx = np.where(np.logical_and(ang > thresholds[1][0], ang < thresholds[1][1]))[0]
>>>>>>> f13161f003f01c8a8e191eb96afb88b07958c95d

        msig = np.zeros(npts)
        msig[midx] = 1
        rsig = np.zeros(npts)
        rsig[ridx] = 1

        result = np.zeros((npts,nchannels))

        for i in range(nchannels):
            channel_level = bb_level[:,i] # gets rid of "masked array

            # create a rolling sum of the main points
<<<<<<< HEAD
            msum = np.cumsum(np.insert(msig * channel_level, 0, 0))
            mrollsum = msum[ww:] - msum[:-ww]

            # to do this accurately we also need a rolling sum for \
            # normalization
            mnorm = np.cumsum(np.insert(msig, 0, 0))
            mrollnorm = mnorm[ww:] - mnorm[:-ww]

            # same procedure for reference points
            rsum = np.cumsum(np.insert(rsig * channel_level, 0, 0))
            rrollsum = rsum[ww:] - rsum[:-ww]

            # same normalization procedure for reference points
            rnorm = np.cumsum(np.insert(rsig, 0, 0))
            rrollnorm = rnorm[ww:] - rnorm[:-ww]

            # now compute difference between main and ref for all \
            # points
            result[window:npts - window, i] = mrollsum / mrollnorm - \
                                              rrollsum / rrollnorm
            result[:window, i] = result[window, i] * np.ones(window)
            result[npts - window:, i] = result[npts - window - 1, i] * \
                                        np.ones(window)
        return(result)

    def process_chopped_signal_f(self, bb_level, chop, window=62, harm=2):
        """
        Gated chopper signal processor.
        Args:
            bb_level (array): npts by nchannels 2D array with baseband 
                if data samples
            chop (int): chopper wheel position from 0 to 8000 
                corresponding to 0 to 360 degrees
            window (int): defines smoothing window of 2*window + 1 
                points. The total smoothing window must span at least 
                one chop cycle.
            harm (int): number of harmonic
        Returns:
            aresult (array), presult (array): 
        """
        npts = len(chop)
        nchannels = np.shape(bb_level)[1]

        ww = 2 * window + 1

        # create array of second harmonic values
        tcos = np.cos(harm * chop / 8000 * 2 * np.pi)
        tsin = np.sin(harm * chop / 8000 * 2 * np.pi)

        cresult = np.zeros((npts, nchannels))
        sresult = np.zeros((npts, nchannels))
        aresult = np.zeros((npts, nchannels))
        presult = np.zeros((npts, nchannels))

        for i in range(nchannels):
            c_level = bb_level[:,i] * tcos
            s_level = bb_level[:,i] * tsin

            for j in range(window, npts - window):
                cresult[j,i] = np.sum(c_level[j - window : j + window + 1])
                sresult[j,i] = np.sum(s_level[j - window : j + window + 1])

            cresult[:window, i] = cresult[window, i] * np.ones(window)
            cresult[npts - window :, i] = cresult[npts - window - 1, i] * \
                                         np.ones(window)
            sresult[:window, i] = sresult[window, i] * np.ones(window)
            sresult[npts - window :, i] = sresult[npts - window - 1, i] * \
                                          np.ones(window)

            aresult[:,i] = 2 / ww * np.sqrt(cresult[:,i]**2 + sresult[:,i]**2)
            presult[:,i] = np.arctan2(sresult[:,i], cresult[:,i])
=======
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
        return(result)

    def process_chopped_signal_f(self, bb_level, chop, window=62, harm=2):
        npts = len(chop)
        nchannels = np.shape(bb_level)[1]

        ww = 2*window+1

        # create array of second harmonic values
        tcos = np.cos(harm*chop/8000*2.*np.pi)
        tsin = np.sin(harm*chop/8000*2.*np.pi)

        cresult = np.zeros((npts,nchannels))
        sresult = np.zeros((npts,nchannels))
        aresult = np.zeros((npts,nchannels))
        presult = np.zeros((npts,nchannels))

        for i in range(nchannels):
            c_level = bb_level[:,i]*tcos
            s_level = bb_level[:,i]*tsin

            for j in range(window, npts-window):
                cresult[j,i] = np.sum(c_level[j-window:j+window+1])
                sresult[j,i] = np.sum(s_level[j-window:j+window+1])

            cresult[:window,i] = cresult[window,i]*np.ones(window)
            cresult[npts-window:,i] = cresult[npts-window-1,i]*np.ones(window)
            sresult[:window,i] = sresult[window,i]*np.ones(window)
            sresult[npts-window:,i] = sresult[npts-window-1,i]*np.ones(window)

            aresult[:,i] = 2/ww*np.sqrt(cresult[:,i]**2+sresult[:,i]**2)
            presult[:,i] = np.arctan2(sresult[:,i],cresult[:,i])
>>>>>>> f13161f003f01c8a8e191eb96afb88b07958c95d

        return(aresult, presult)


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
            npix (int): number of pixels (default is 16)
        Returns:
            none
        """
        self.npix = npix
        IFProc.__init__(self, filename)

        # identify the obspgm
        self.map_coord = 0 # set this up to be nominal for all cases

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
                test_map_coord = b''.join(self.nc.variables[
                    'Header.Map.MapCoord'][:]).decode().strip()
                if test_map_coord[0] == 'A':
                    self.map_coord = 0
                elif test_map_coord[0] == 'R':
                    self.map_coord = 1
                else:
                    self.map_coord = 0

                self.map_motion = b''.join(self.nc.variables[
                    'Header.Map.MapMotion'][:]).decode().strip()
            except:
                print('%s does not have map parameters'%(self.filename))

        elif self.obspgm == 'Cal':
            print('WARNING: %d is a Cal observation'%(self.obsnum))

        elif self.obspgm == 'On':
            print('WARNING: %d is an On observation'%(self.obsnum))

        else:
            print('WARNING: ObsPgm type %s for Obsum %d is not identified'%(
                self.obspgm, self.obsnum))

        # data arrays
        self.time = self.nc.variables['Data.TelescopeBackend.TelTime'][:]
        if self.map_coord == 0:
            self.azmap = self.nc.variables['Data.TelescopeBackend.TelAzMap'
                                          ][:]* 206264.8
            self.elmap = self.nc.variables['Data.TelescopeBackend.TelElMap'
                                          ][:]* 206264.8
            self.parang = self.nc.variables['Data.TelescopeBackend.ActParAng'][:]
#            self.parang = np.zeros(len(self.azmap))
        elif self.map_coord == 1:
            # OK it is not really azmap and elmap
            self.azmap = (self.nc.variables[
                'Data.TelescopeBackend.SourceRaAct'][:] - self.source_RA)\
                * np.cos(self.source_Dec) * 206264.8
            self.elmap = (self.nc.variables[
                'Data.TelescopeBackend.SourceDecAct'][:] - self.source_Dec)\
                * 206264.8
            self.parang = self.nc.variables['Data.TelescopeBackend.ActParAng'
                                           ][:]
        else:
            self.azmap = self.nc.variables['Data.TelescopeBackend.TelAzMap'
                                          ][:] * 206264.8
            self.elmap = self.nc.variables['Data.TelescopeBackend.TelElMap'
                                          ][:] * 206264.8
            self.parang = np.zeros(len(self.azmap))
            
        self.bufpos = self.nc.variables['Data.TelescopeBackend.BufPos'][:]
        if 'ifproc' in filename:
            self.bb_level = self.nc.variables['Data.IfProc.BasebandLevel'][:]
            try:
                print('get chop')
                chop = self.nc.variables['Data.Msip1mm.BeamChopperActPos'][:]
                chop_option = self.nc.variables['Header.Msip1mm.BeamChopperActState'][0]
                if chop_option == 8 or chop_option == 16:
                    print(' chopping')
                    self.level = self.process_chopped_signal(self.bb_level, chop)
                else:
                    print(' not chopping')
                    self.level = self.bb_level
            except Exception as e:
                print(e)
                print(' no chop')
                self.level = self.bb_level
        elif 'lmttpm' in filename:
            self.level = detrend(self.nc.variables['Data.LmtTpm.Signal'][:],
                                 axis=0)
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
            self.caldata[ipix, :] = (self.level[:, ipix] - 
                    CAL.calcons[ipix, 1]) / CAL.calcons[ipix, 0]
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
        self.map_n = []
        for i in range(self.npix):
            self.map_x.append(self.azmap)
            self.map_y.append(self.elmap)
            self.map_n.append(self.nsamp)
            self.map_data.append(self.caldata[i,:] - self.bias[i])
        self.map_x = np.array(self.map_x)
        self.map_y = np.array(self.map_y)
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
        self.map_coord = 0 # set this up to be nominal for all cases
        if self.obspgm == 'Cal':
            print('%d is a Cal observation'%(self.obsnum))
        else:
            print('WARNING: %d is NOT a Cal observation : %s'%(self.obsnum, 
                                                               self.obspgm))

        # data arrays
        self.time = self.nc.variables['Data.TelescopeBackend.TelTime'][:]
        self.azmap = self.nc.variables['Data.TelescopeBackend.TelAzMap'][:]
        self.elmap = self.nc.variables['Data.TelescopeBackend.TelElMap'][:]
        self.parang = np.zeros(len(self.azmap))
        self.bufpos = self.nc.variables['Data.TelescopeBackend.BufPos'][:]
        if 'ifproc' in filename:
            self.bb_level = self.nc.variables['Data.IfProc.BasebandLevel'][:]
            try:
                print('get chop cal')
                chop = self.nc.variables['Data.Msip1mm.BeamChopperActPos'][:]
                chop_option = self.nc.variables['Header.Msip1mm.BeamChopperActState'][0]
                if chop_option == 8 or chop_option == 16:
                    print(' chopping cal')
                    self.level = self.process_chopped_signal(self.bb_level, chop)
                else:
                    print(' not chopping cal')
                    self.level = self.bb_level
            except Exception as e:
                print(e)
                print(' no chop cal')
                self.level = self.bb_level
        elif 'lmttpm' in filename:
            self.level = detrend(self.nc.variables['Data.LmtTpm.Signal'][:], 
                                 axis=0)
        else:
            self.level = np.zeros(0)
        self.nsamp = len(self.level)
        self.tamb = 280.
        self.receiver = b''.join(self.nc.variables['Header.Dcs.Receiver'][:]
                                ).decode().strip()
        try:
            self.blank_level = self.nc.variables['Header.' + self.receiver + 
                                                 '.BlankLevel'][0]
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
        for ipix in range(self.npix):
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
        for ipix in range(self.npix):
            vsky = np.median(self.level[sky_list, ipix])
            vhot = np.median(self.level[hot_list, ipix])
            vzero = self.blank_level
            self.tsys[ipix] = self.tamb * (vsky - vzero) / (vhot - vsky)
        
