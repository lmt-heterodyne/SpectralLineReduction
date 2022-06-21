"""
Module for dealing with spectral line data
classes: Line, LineDataHeader, LineData, Accum, NetCDFLineHeader
methods: read_obsnum_ps, read_obsnum_otf, count_otf_spectra, 
uses: 
author: FPS
date: September 2019 (documented)
changes:
python 3
moved reader functions to utils.reader
"""

import numpy as np
import math
from scipy.interpolate import interp1d
from lmtslr.spec.spec import SpecBankCal, SpecBankData
from lmtslr.ifproc.ifproc import IFProcData, IFProcCal
from lmtslr.utils.roach_file_utils import create_roach_list, \
    lookup_roach_files
from lmtslr.utils.ifproc_file_utils import lookup_ifproc_file
import netCDF4


class Line(object):
    """
    Class with methods that operate on spectral line data.
        
    The methods here are common to any spectral line reduction. They
    assume a common format with arrays of data for the x axis and y 
    axis of the spectrum.
    """
    def __init__(self, iarray, xarray, yarray, tarray, xlabel):
        """
        Constructor for Line class.
        Args:
            iarray (np array): array of channel numbers
            xarray (np array): array of "x axis" values.  Could be 
                velocity, frequency, etc.
            yarray (np array): spectral line data values corresponding 
                to xarray/iarray.
            tarray (np array or None): array of tsys
            xlabel (string): label for the x axis; used in plots
        """
        self.nchan = len(xarray)
        self.iarray = iarray
        self.xarray = xarray
        self.yarray = yarray
        self.tarray = tarray
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xname = xlabel

    def __len__(self):
        return self.nchan

    def eliminate(self, list, interpolate=True):
        """
        Sets yarray values to nan according to input list so that 
        channels can be ignored in functions.
        Args:
            list (list): list of channels that are to be eliminated
            interpolate (boolean):   if to interpolate instead
        Returns:
            none
        """
        if interpolate:
            for i in list:
                if i<=0: continue
                j = np.where(self.iarray == i)[0][0]                
                self.yarray[j] = (self.yarray[j-1] + self.yarray[j+1])/2.0
        else:
            for i in list:
                self.yarray[np.where(self.iarray == i)] = np.nan

    def baseline(self, list, n, baseline_order=0):
        """
        Removes a baseline from the working spectrum.
        Args:
            list (list): list of lists specifies indices in the 
                spectrum where baseline is to be fit
            n (int): not used
            baseline_order (int): order of polynomial to be fit to 
                baseline (used in np.polyfit)
        Result:
            baseline (np array): polynomial from baseline fit
            rms (float): rms of data about polynomial fit in regions
            yarray (np.array): yarray is modified by subtracting the 
                baseline polynomial
        """
        input_baseline_list = np.array(list)
        baseline_list = input_baseline_list[np.where(np.isfinite(
            self.yarray[input_baseline_list]) == True)[0]]
        n_baseline_list = len(baseline_list)
        params = np.polyfit(baseline_list, self.yarray[baseline_list], 
                            baseline_order)
        self.baseline = np.polyval(params, np.arange(self.nchan))
        residuals = self.yarray[baseline_list] - np.polyval(params, 
                                                            baseline_list)
        self.rms = np.sqrt(np.dot(residuals.transpose(), residuals) 
                           / n_baseline_list)
        self.yarray = self.yarray - self.baseline
    
    # special definition of generic line_statistics method
    def line_stats(self, list, n):
        """
        Computes line statistics for a spectrum.
        Args:
            list (list): indices in spectrum for line integration
            n (int): not used
        Line Statistics are:
            ymax = maximum value in channel_list
            cmax = channel of maximum value
            xmax = velocity of maximum value
            yint = integral of line over channel list in units of K km/s
            yerr = error on YINT
            xmean = first moment of line integrated over channel list
            xwidth = estimate of line width based on Peak value and integrated
                intensity
        """
        input_channel_list = np.array(list)
        channel_list = input_channel_list[np.where(np.isfinite(
            self.yarray[input_channel_list]) == True)[0]]
        n_channel_list = len(channel_list)
        self.ymax = np.max(self.yarray[channel_list])
        self.cmax = channel_list[np.where(
            self.yarray[channel_list] == self.ymax)[0]]
        self.xmax = self.xarray[self.cmax]
        self.yint = np.sum(self.yarray[channel_list]) * self.dxdc
        self.yerr = self.dxdc * self.rms * np.sqrt(n_channel_list)
        self.xmean = np.sum(self.yarray[channel_list] 
            * self.xarray[channel_list]) / np.sum(self.xarray[channel_list])
        self.xwidth = self.yint / self.ymax

    def smo(self,nw,type='box'):
        """
        Smooths a spectrum using a smoothing function
        Args:
            nw (int): number of channels in smoothing function
            type (string): identifies smoothing function
                           'box' (default) - boxcar smoothing function
                           'hanning' - hanning window
                           'hamming' - hamming window
        Result:
            yarray replaced with smoothed version of spectrum sampled 
            at original spacing
        """
        if type == 'box':
            window = np.ones(nw)
        elif type == 'hanning':
            window = np.hanning(nw)
        elif type == 'hamming':
            window = np.hamming(nw)
        else:
            print('smo: window %s not supported'%(type))
        denom = np.sum(window)
        self.yarray = np.convolve(self.yarray,window,mode='same') / denom

    def xlist(self,regions):
        """
        Generates a list of lists with indices corresponding to x axis 
        values.
        Args:
            regions (list): a list of lists of regions in spectrum 
                identified according to x axis values
        Returns:
            channel_list (list): list of lists with channel numbers 
                corresponding to input regions
            nchannels (int): number of channels in the list
        """
        nregions = len(regions)
        if self.xarray[1] > self.xarray[0]:
            channel_list = []
            for i in range(nregions):
                idx = np.where((self.xarray>regions[i][0]) == 
                               (self.xarray < regions[i][1]))[0]
                channel_list = channel_list + idx.tolist()
        else:
            channel_list = []
            for i in range(nregions):
                idx = np.where((self.xarray>regions[i][1]) ==  
                               (self.xarray < regions[i][0]))[0]
                channel_list = channel_list + idx.tolist()

        nchannels = len(channel_list)
        return channel_list, nchannels

    def clist(self, regions):
        """
        Generates a list of lists with indices corresponding to channel
        number.
        Args:
            regions (list): a list of lists of channels in spectrum
        Returns:
            channel_list (list): list of lists with indices 
                corresponding to input regions
            nchannels (int): number of channels in the list
        """
        nregions = len(regions)
        channel_list = []
        for i in range(nregions):
            idx = np.where((self.iarray > regions[i][0]) == 
                           (self.iarray < regions[i][1]))[0]
            channel_list = channel_list + idx.tolist()
        nchannels = len(channel_list)
        return channel_list, nchannels


class LineDataHeader(object):
    """
    Class to create header information for line observation and provide
    methods to transform "x" axis.
    units of velocity: km/s
    units of frequency: Hertz
    """
    def __init__(self, ifproc, bank, nchan, bandwidth):
        """
        Constructor for LineDataHeader.
        Args:
            ifproc (object): ifproc object with ifproc information
            bank (int): selects which bank of the spectrometer is 
                relevant 
            nchan (int): number of channels in the spectrum
            bandwidth (float): bandwidth of the spectrum (MHz)
        Returns:
            none
        """
        self.c = 299792.458 # speed of light in km/s
        self.bank = bank
        self.v_source = ifproc.velocity # km/s
        self.v_system = ifproc.velocity_system
        self.z_source = ifproc.line_redshift[bank] # z
        self.frest = ifproc.line_rest_frequency[bank] * np.float64(1.0e9) # Hz
        self.dop_track  = ifproc.doppler_track
        self.vbary = ifproc.barycenter_velocity # km/s
        self.vobs = ifproc.observatory_velocity # km/s
        self.fsky = ifproc.sky_frequency[bank] * np.float64(1.0e9) # Hz
        self.flo1  = ifproc.lo_1_frequency * np.float64(1.0e9) # Hz
        self.flo2  = ifproc.lo_2_frequency[bank] * np.float64(1.0e9) # Hz
        self.fif1 = ifproc.if_1_frequency[bank] * np.float64(1.0e9) # Hz
        self.fif2 = ifproc.if_2_frequency[bank] * np.float64(1.0e9) # Hz
        self.fsyn = ifproc.synthesizer_frequency[bank] * np.float64(1.0e9) # Hz
        self.foff = ifproc.frequency_offset[bank] * np.float64(1.0e9) # Hz  
        self.voff = -self.foff/self.frest * self.c # km/z
        self.sideband = ifproc.sideband[bank]

        self.nchan = nchan 
        self.ch0 = float(nchan - 1) / 2.
        self.bw    = bandwidth * np.float64(1.0e6) # Hz  
        if self.sideband == 1:
            self.dfdc = -self.bw / self.nchan
        else:
            self.dfdc = self.bw / self.nchan
        self.dvdc = -self.dfdc / self.frest * self.c # radio definition

        self.iarray = np.arange(self.nchan)

        # default is system velocity for the x array
        if self.v_system == 0:
            self.v0_lsr = self.v_source
            self.v0_bary = self.v_source + self.vbary
            self.v0_sky = self.vobs + self.v_source
            self.x_vlsr()
        elif self.v_system == 1:
            self.v0_lsr = self.v_source + self.vbary
            self.v0_bary = self.v_source
            self.v0_sky = self.vobs + self.v_source - self.vbary
            self.x_vbary()


    def set_x_axis(self, option):
        """
        Sets x axis according to option supplied.
        Args:
            option (str): string specifying x axis option.
                          vslr/VLSR: velocity units wrt Local Standard 
                            of Rest
                          vbary/VBARY: velocity units wrt Solar System 
                            Barycenter
                          vsky/VSKY: velocity units wrt topocentric 
                            system
                          vsrc/VSRC: velocity units wrt the source 
                            velocity
                          fsky/FSKY: frequency units accounting for 
                            topocentric doppler shift
                          flsr/FLSR: frequency units accounting for 
                            Local Standard of Rest and doppler shift
                          fbary/FBARY: frequency units accounting for 
                            S.S. Barycenter and doppler shift
                          fsrc/FSRC: frequency units accounting for 
                            source doppler shift
        """
        if option == 'VLSR' or option == 'vlsr':
            self.x_vlsr()
        elif option == 'VBARY' or option == 'vbary':
            self.x_vbary()
        elif option == 'VSKY' or option == 'vsky':
            self.x_vsky()
        elif option == 'VSRC' or option == 'vsrc':
            self.x_vsrc()
        elif option == 'FSKY' or option == 'fsky':
            self.x_fsky()
        elif option == 'FLSR' or option == 'flsr':
            self.x_flsr()
        elif option == 'FBARY' or option == 'fbary':
            self.x_fbary()
        elif option == 'FSRC' or option == 'fsrc':
            self.x_fsrc()
        else:
            print('LineDataHeader:set_x_axis: bad option %s'%(option))
            
            
    def x_vlsr(self):
        """
        Sets x axis to be in velocity units wrt Local Standard of Rest.
        Args:
            none
        Returns:
            none
        """
        self.xv0 = self.v0_lsr + self.voff
        self.xarray = (self.iarray - self.ch0) * self.dvdc + self.xv0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 0
        self.xname = 'VLSR'

    def x_vbary(self):
        """
        Sets x axis to be in velocity units wrt Solar System Barycenter.
        Args:
            none
        Returns:
            none
        """
        self.xv0 = self.v0_bary + self.voff
        self.xarray = (self.iarray - self.ch0) * self.dvdc + self.xv0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 0
        self.xname = 'VBARY'

    def x_vsky(self):
        """
        Sets x axis to be in velocity units wrt topocentric system.
        Args:
            none
        Returns:
            none
        """
        self.xv0 = self.v0_sky + self.voff
        self.xarray = (self.iarray - self.ch0) * self.dvdc + self.xv0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 0  
        self.xname = 'VSKY'

    def x_vsrc(self):
        """
        Sets x axis to be in velocity units wrt the source velocity.
        Args:
            none
        Returns:
            none
        """
        self.xv0 = self.voff
        self.xarray = (self.iarray - self.ch0) * self.dvdc + self.xv0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 0  
        self.xname = 'VSRCE'

    def x_fsky(self):
        """
        Sets x axis to be in frequency units accounting for topocentric
        doppler shift.
        Args:
            none
        Returns:
            none
        """
        self.xarray = (self.iarray - self.ch0) * self.dfdc + self.fsky
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 1  
        self.xname = 'FSKY'

    def x_flsr(self):
        """
        Sets x axis to be in frequency units accounting for Local 
        Standard of Rest and doppler shift.
        Args:
            none
        Returns:
            none
        """
        self.xf0 = self.fsky + self.frest / self.c * (self.v0_sky
                                                      - self.v0_lsr)
        self.xarray = (self.iarray - self.ch0) * self.dfdc + self.xf0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 1  
        self.xname = 'FLSR'

    def x_fbary(self):
        """
        Sets x axis to be in frequency units accounting for S.S. 
        Barycenter and doppler shift.
        Args:
            none
        Returns:
            none
        """
        self.xf0 = self.fsky + self.frest / self.c * (self.v0_sky
                                                      - self.v0_bary)
        self.xarray = (self.iarray - self.ch0) * self.dfdc + self.xf0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 1  
        self.xname = 'FLSR'

    def x_fsrc(self):
        """
        Sets x axis to be in frequency units accounting for source 
        doppler shift.
        Args:
            none
        Returns:
            none
        """
        self.xf0 = self.fsky + self.frest / self.c * self.v0_sky
        self.xarray = (self.iarray - self.ch0) * self.dfdc + self.xf0
        self.dxdc = self.xarray[1] - self.xarray[0]
        self.xtype = 1
        self.xname = 'FSRCE'

    def v2c(self, v):
        """
        Returns velocity converted into channel number.
        Args:
            v (float/np array): velocity for conversion
        Returns:
            c (int/np array): channels corresponding to velocity input
        """
        if self.xtype == 0:
            c = np.array(np.round((v - self.xv0) / self.dvdc + self.ch0), 
                         dtype=int)
        else:
            print('v2c: xarray is in frequency')
        return(c)

    def f2c(self, f):
        """
        Returns frequency converted into channel number.
        Args:
            f (float/np array): frequency for conversion
        Returns:
            c (int/np array): channel numbers corresponding to 
                frequency input
        """
        if self.xtype == 1:
            c = np.array(np.round((f - self.xf0) / self.dfdc + self.ch0), 
                         dtype=int)
        else:
            print('f2c: xarray is in velocity')
        return(c)

    def c2v(self, c):
        """
        Returns channel number converted into velocity.
        Args:
            c (int): input channel number
        Returns:
            v (float): velocity of channel
        """
        if self.xtype == 0:
            v = (c - self.ch0) * self.dvdc + self.xv0
        else:
            print('c2v: xarray is in frequency')
        return(v)

    def c2f(self, c):
        """
        Returns channel number converted into frequency.
        Args:
            c (int): input channel number
        Returns:
            f (float): frequency of channel
        """
        if self.xtype == 1:
            f = (c - self.ch0) * self.dfdc + self.xf0
        else:
            print('c2f: xarray is in velocity')
        return(v)

    def make_velocity_list(self, regions):
        """
        Returns the list and number of velocities corresponding to 
        target regions.
        Args:
            regions (list): list of target regions
        Returns:
            channel_list (list): list of corresponding velocity 
                channels
            nchannels (int): number of corresponding velocity channels
        """
        if self.xtype == 0:
            nregions = len(regions)
            channel_list = []
            for i in range(nregions):
                c0 = self.v2c(regions[i][0])
                c1 = self.v2c(regions[i][1]) + 1
                if c1>c0:
                    channel_list = channel_list + range(c0, c1)
                else:
                    channel_list = channel_list + range(c1, c0)
            nchannels = len(channel_list)
        else:
            print('make_velocity_list: xarray not in velocity units')
        return channel_list, nchannels

    def make_frequency_list(self, regions):
        """
        Returns the list and number of frequencies corresponding to 
        target regions.
        Args:
            regions (list): list of target regions
        Returns:
            channel_list (list): list of corresponding frequency 
                channels
            nchannels (int): number of corresponding frequency channels
        """
        if self.xtype == 1:
            nregions = len(regions)
            channel_list = []
            for i in range(nregions):
                c0 = self.f2c(regions[i][0])
                c1 = self.f2c(regions[i][1]) + 1
                if c1>c0:
                    channel_list = channel_list + range(c0, c1)
                else:
                    channel_list = channel_list + range(c1, c0)
            nchannels = len(channel_list)
        else:
            print('make_frequency_list: xarray not in frequency units')
        return channel_list, nchannels


    def make_channel_list(self, regions):
        """
        Returns list of channels given a list of lists of intervals to 
        be processed.
        for example: regions = [[1,3],[5,8]] creates [1,2,3,5,6,7,8]
        Args:
            regions (list): list of target regions
        Returns:
            channel_list (list): list of corresponding channels
            nchannels (int): number of corresponding channels
        """
        nregions = len(regions)
        channel_list = []
        for i in range(nregions):
            channel_list = channel_list + range(regions[i][0], 
                                                regions[i][1] + 1)
        nchannels = len(channel_list)
        return channel_list, nchannels



class LineData(LineDataHeader):
    """
    Class to provide methods to manipulate line data.
    units of velocity: km/s
    units of frequency: Hertz
    """
    def __init__(self, ifproc, bank, nchan, bandwidth, data, tsys=None):
        """
        Constructor for LineData class.
        Args:
            ifproc (object): ifproc object with ifproc information
            bank (int): selects which bank of the spectrometer is 
                relevant 
            nchan (int): number of channels in the spectrum
            bandwidth (float): bandwidth of the spectrum (MHz)
            data (array): the spectral line data
            tsys (array): [optional] tsys array
        Returns:
            none
        """
        LineDataHeader.__init__(self, ifproc, bank, nchan, bandwidth)
        self.data  = data
        # yarray is the working spectrum - initialized as copy of data
        self.yarray = data.copy()
        if type(tsys) == np.ndarray:
            self.tarray = tsys.copy()
        else:
            #self.tarray = data * 0 + 1
            self.tarray = None

    def __len__(self):
        return len(self.data)

    def line(self):
        """
        Returns a Line object using data.
        Args:
            none
        Returns:
            Line(...) (object): Line object generated from data in 
                iarray, xarray, yarray, and xname.
        """
        return Line(self.iarray, self.xarray, self.yarray, self.tarray, self.xname)

    def revert(self):
        """
        Resets the spectrum to the original values.
        Args:
            none
        Returns:
            none
        """
        self.yarray = self.data.copy()

    def cgen(self, c0, c1, dc=1):
        """
        Generate a Line object with data covering range from channel c0
        to channel c1.
        Args:
            c0 (int): beginning channel number for new spectrum
            c1 (int): ending channel number for new spectrum
            dc (int): spacing of channels for new specrum
        Returns:
            result (Line object): Line object with resampled spectrum
        """
        new_iarray = np.arange(c0, c1, dc)
        fx = interp1d(self.iarray, self.xarray)
        fy = interp1d(self.iarray, self.yarray)
        new_xarray = fx(new_iarray)
        new_yarray = fy(new_iarray)
        if type(self.tarray) == np.ndarray:
            ft = interp1d(self.iarray, self.tarray)
            new_tarray = ft(new_iarray)
        else:
            new_tarray = None
        result = Line(new_iarray, new_xarray, new_yarray, new_tarray, self.xname)
        return result

    def xgen(self, x0, x1, dx):
        """
        Generate a Line object with data covering range of x axis 
        values from x0 to x1.
        Args:
            x0 (int): beginning x value for new spectrum
            x1 (int): ending x value for new spectrum
            dx (int): spacing of samples for new specrum
        Returns:
            result (Line object): Line object with resampled spectrum
        """
        new_xarray = np.arange(x0, x1, dx)
        fi = interp1d(self.xarray, self.iarray)
        fy = interp1d(self.xarray, self.yarray)
        new_iarray = fi(new_xarray)
        new_yarray = fy(new_xarray)
        if type(self.tarray) == np.ndarray:
            ft = interp1d(self.xarray, self.tarray)
            new_tarray = ft(new_xarray)
        else:
            new_tarray = None

        result = Line(new_iarray, new_xarray, new_yarray, new_tarray, self.xname)
        return result

    def cslice(self, c0, c1):
        """
        Generate a Line object by taking a slice from original array.
        Args:
            c0 (int): starting channel number of slice
            c1 (int): ending channel number of slice   @todo  +1 ???
        Returns:
            result (Line object): Line object with the slice
        """
        if c0>c1:
            c0,c1 = c1,c0
        new_iarray = np.arange(c0, c1)
        new_xarray = self.xarray[c0:c1]
        new_yarray = self.yarray[c0:c1]
        if type(self.tarray) == np.ndarray:
            new_tarray = self.tarray[c0:c1]
        else:
            new_tarray = None
        result = Line(new_iarray, new_xarray, new_yarray, new_tarray, self.xname)
        return result

    def vslice(self, v0, v1):
        """
        Generate a Line object by taking a slice from original array; 
        slice obtained by giving velocity limits.
        Args:
            v0 (int): starting velocity of slice
            v1 (int): ending velocity of slice
        Returns:
            result (Line object): Line object with the slice
        """
        c0 = self.v2c(v0)
        c1 = self.v2c(v1)
        if c0>c1:
            c0,c1 = c1,c0 
        new_iarray = self.iarray[c0:c1]
        new_xarray = self.xarray[c0:c1]
        new_yarray = self.yarray[c0:c1]
        if type(self.tarray) == np.ndarray:
            new_tarray = self.tarray[c0:c1]
        else:
            new_tarray = None

        result = Line(new_iarray, new_xarray, new_yarray, new_tarray, self.xname)
        return result

    def fslice(self, f0, f1):
        """
        Generate a Line object by taking a slice from original array; 
        slice obtained by giving frequency limits.
        Args:
            f0 (int): starting freqeuncy of slice
            f1 (int): ending frequency of slice
        Returns:
            result (Line object): Line object with the slice
        """
        c0 = self.f2c(f0)
        c1 = self.f2c(f1)
        if c0>c1:
            c0,c1 = c1,c0
        new_iarray = self.iarray[c0:c1]
        new_xarray = self.xarray[c0:c1]
        new_yarray = self.yarray[c0:c1]
        if type(self.tarray) == np.ndarray:
            new_tarray = self.tarray[c0:c1]
        else:
            new_tarray = None

        result = Line(new_iarray, new_xarray, new_yarray, new_tarray, self.xname)
        return result

    
class Accum():
    """
    Class to manage spectral averaging.
    All spectra to be averaged must be aligned in channel/index space.
    Weighted averages are possible using weight provided when spectrum 
    is added to the stack.
    """
    def __init__(self):
        """
        Constructor for Accum class.
        Args:
            none
        Returns:
            none
        """
        self.count = 0
        self.alist = []
        self.wlist = []

    def load(self, data, wt=1):
        """
        Loads a set of data to be averaged.
        Args:
            data (array): data array
            wt (float): weighting of data (default is 1)
        Returns:
            none
        """
        ndata = len(data)
        if self.count == 0:
            self.nchannels = ndata
            self.count += 1 
            self.alist.append(data)
            self.wlist.append(wt)
        elif self.count>0 and ndata!=self.nchannels:
            print('accum:load - WARNING - number of points in spectrum doesnt \
                   match')
        else:
            self.count += 1 
            self.alist.append(data)
            self.wlist.append(wt)

    def ave(self, type=0):
        """
        Averages data sets.
        Args:
            type (int): type of averaging (not used)
        Returns:
            none
        """
        self.average = np.zeros(self.nchannels)
        a = np.array(self.alist)
        w = np.array(self.wlist)
        for i in range(self.nchannels):
            slist = np.where(np.isfinite(a[:, i]))[0]
            self.average[i] = np.sum(a[slist, i] * w[slist]) / np.sum(w[slist])


class NetCDFLineHeader():
    """
    Class to manage the header information to write to netcdf file.
    Assumes an open NetCDF file : nc
    """
    def __init__(self, nc):
        """
        Constructor for NetCDFLineHeader class.
        Args:
            nc (object): NetCDF file opened with netCDF4 module
        Returns:
            none
        """
        self.nc = nc

    def write_line_data_header_variables(self, L):
        """
        Writes creates line header variables for a NetCDF file.
        Args:
            L (object): LineData object
        Returns:
            none
        """
        t = self.nc.createVariable('Header.LineData.Bank', 'i4')
        self.nc.variables['Header.LineData.Bank'][0] = L.bank

        t = self.nc.createVariable('Header.LineData.VSource', 'f8')
        self.nc.variables['Header.LineData.VSource'][0] = L.v_source
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.VSystem', 'f8')
        self.nc.variables['Header.LineData.VSystem'][0] = L.v_system
        t.units = 'km/s'
        
        t = self.nc.createVariable('Header.LineData.ZSource', 'f8')
        self.nc.variables['Header.LineData.ZSource'][0] = L.z_source
        t.units = 'z'

        t = self.nc.createVariable('Header.LineData.LineRestFrequency', 'f8')
        self.nc.variables['Header.LineData.LineRestFrequency'][0] = L.frest
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.DopplerTrack', 'i4')
        self.nc.variables['Header.LineData.DopplerTrack'][0] = L.dop_track
        t.units = 'T/F'

        t = self.nc.createVariable('Header.LineData.VBarycenter', 'f8')
        self.nc.variables['Header.LineData.VBarycenter'][0] = L.vbary
        t.units = 'km/s'
 
        t = self.nc.createVariable('Header.LineData.VObservatory', 'f8')
        self.nc.variables['Header.LineData.VObservatory'][0] = L.vobs
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.SkyFrequency', 'f8')
        self.nc.variables['Header.LineData.SkyFrequency'][0] = L.fsky
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.LO1Frequency', 'f8')
        self.nc.variables['Header.LineData.LO1Frequency'][0] = L.flo1
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.LO2Frequency', 'f8')
        self.nc.variables['Header.LineData.LO2Frequency'][0] = L.flo2
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.IF1Frequency', 'f8')
        self.nc.variables['Header.LineData.IF1Frequency'][0] = L.fif1
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.IF2Frequency', 'f8')
        self.nc.variables['Header.LineData.IF2Frequency'][0] = L.fif2
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.SynthesizerFrequency', 'f8')
        self.nc.variables['Header.LineData.SynthesizerFrequency'][0] = L.fsyn
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.FrequencyOffset', 'f8')
        self.nc.variables['Header.LineData.FrequencyOffset'][0] = L.foff
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.VelocityOffset', 'f8')
        self.nc.variables['Header.LineData.VelocityOffset'][0] = L.voff
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.Sideband', 'i4')
        self.nc.variables['Header.LineData.Sideband'][0] = L.sideband
        t.units = 'USB/LSB'

        t = self.nc.createVariable('Header.LineData.NChannels', 'i4')
        self.nc.variables['Header.LineData.NChannels'][0] = L.nchan
        
        t = self.nc.createVariable('Header.LineData.Bandwidth', 'f8')
        self.nc.variables['Header.LineData.Bandwidth'][0] = L.bw
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.DeltaFrequency', 'f8')
        self.nc.variables['Header.LineData.DeltaFrequency'][0] = L.dfdc
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.DeltaVelocity', 'f8')
        self.nc.variables['Header.LineData.DeltaVelocity'][0] = L.dvdc
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.LSRFrameV0', 'f8')
        self.nc.variables['Header.LineData.LSRFrameV0'][0] = L.v0_lsr
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.BarycenterFrameV0', 'f8')
        self.nc.variables['Header.LineData.BarycenterFrameV0'][0] = L.v0_bary
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.SkyFrameV0', 'f8')
        self.nc.variables['Header.LineData.SkyFrameV0'][0] = L.v0_sky
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.XType', 'i4')
        self.nc.variables['Header.LineData.XType'][0] = L.xtype
        t.long_name = 'Axis Type'

        if False:
            # already written in:    Header.Line.ChannelNumbe
            iarray = self.nc.createVariable('Header.LineData.ChannelNumber', 'i4',
                                        ('nchan',))
            iarray[0:L.nchan] = L.iarray[0:L.nchan]
            iarray.long_name = 'Channel Number'

        if True:
            # what's below was done elsewhere
            return

        # the following are the generic spectrum axis parameters

        t = self.nc.createVariable('Header.SpectrumAxis.CRPIX', 'f8')
        self.nc.variables['Header.SpectrumAxis.CRPIX'][0] = L.ch0
        t.long_name = 'CRPIX'

        t = self.nc.createVariable('Header.SpectrumAxis.CRVAL', 'f8')
        self.nc.variables['Header.SpectrumAxis.CRVAL'][0] = L.xv0
        t.long_name = 'CRVAL'
        if L.xtype == 0:
            t.units = 'km/s'
        else:
            t.units = 'Hz'

        t = self.nc.createVariable('Header.SpectrumAxis.CDELT', 'f8')
        self.nc.variables['Header.SpectrumAxis.CDELT'][0] = L.dxdc
        t.long_name = 'CDELT'
        if L.xtype == 0:
            t.units = 'km/s'
        else:
            t.units = 'Hz'

        ctype = self.nc.createVariable('Header.SpectrumAxis.CTYPE', 'c', 
                                       ('nlabel',))
        ctype[0:len(L.xname)] = L.xname[0:len(L.xname)]
        ctype.long_name = 'CTYPE'

        xarray = self.nc.createVariable('Header.SpectrumAxis.CAXIS', 'f8', 
                                        ('nchan',))
        xarray[0:L.nchan] = L.xarray[0:L.nchan]
        if L.xtype == 0:
            xarray.units = 'km/s'
        else:
            xarray.units = 'Hz'
        xarray.long_name = L.xname

    def read_line_data_header_variables(self, L):
        """
        Reads the header data from a NetCDF file.
        Args:
            L (object): LineData object
        Returns:
            none
        """
        L.bank = self.nc.variables['Header.LineData.Bank'][0]
        L.v_source = self.nc.variables['Header.LineData.VSource'][0]
        L.v_system = self.nc.variables['Header.LineData.VSystem'][0]
        L.z_source = self.nc.variables['Header.LineData.ZSource'][0]
        L.frest = self.nc.variables['Header.LineData.LineRestFrequency'][0]
        L.dop_track = self.nc.variables['Header.LineData.DopplerTrack'][0]
        L.vbary = self.nc.variables['Header.LineData.VBarycenter'][0]
        L.vobs = self.nc.variables['Header.LineData.VObservatory'][0]
        L.fsky = self.nc.variables['Header.LineData.SkyFrequency'][0]
        L.flo1 = self.nc.variables['Header.LineData.LO1Frequency'][0]
        L.flo2 = self.nc.variables['Header.LineData.LO2Frequency'][0]
        L.fif1 = self.nc.variables['Header.LineData.IF1Frequency'][0]
        L.fif2 = self.nc.variables['Header.LineData.IF2Frequency'][0]
        L.fsyn = self.nc.variables['Header.LineData.SynthesizerFrequency'][0]
        L.foff = self.nc.variables['Header.LineData.FrequencyOffset'][0]
        L.voff = self.nc.variables['Header.LineData.VelocityOffset'][0]
        L.sideband = self.nc.variables['Header.LineData.Sideband'][0]
        L.nchan = self.nc.variables['Header.LineData.NChannels'][0]
        L.bw = self.nc.variables['Header.LineData.Bandwidth'][0]
        L.dfdc = self.nc.variables['Header.LineData.DeltaFrequency'][0]
        L.dvdc = self.nc.variables['Header.LineData.DeltaVelocity'][0]
        L.v0_lsr = self.nc.variables['Header.LineData.LSRFrameV0'][0]
        L.v0_bary = self.nc.variables['Header.LineData.BarycenterFrameV0'][0]
        L.v0_sky = self.nc.variables['Header.LineData.SkyFrameV0'][0]
        L.xtype = self.nc.variables['Header.LineData.XType'][0]
        L.iarray = self.nc.variables['Header.LineData.ChannelNumber'][:]
        # the following are the generic spectrum axis parameters
        L.ch0 = self.nc.variables['Header.SpectrumAxis.CRPIX'][0]
        L.xv0 = self.nc.variables['Header.SpectrumAxis.CRVAL'][0]
        L.dxdc = self.nc.variables['Header.SpectrumAxis.CDELT'][0]
        L.xname = netCDF4.chartostring(self.nc.variables['Header.SpectrumAxis.CTYPE'][:])
        L.xarray = self.nc.variables['Header.SpectrumAxis.CAXIS'][:]


    def write_line_header_variables(self, L):
        """
        Writes line header variables to L (Line) object.
        Args:
            L (object): Line object
        Returns:
            none
        """
        t = self.nc.createVariable('Header.Line.NChannels', 'i4')
        self.nc.variables['Header.Line.NChannels'][0] = L.nchan

        iarray = self.nc.createVariable('Header.Line.ChannelNumber', 'i4', 
                                        ('nchan',))
        iarray[0:L.nchan] = L.iarray[0:L.nchan]
        iarray.long_name = 'Channel Number'
        
        # the following are the generic spectrum axis parameters
        t = self.nc.createVariable('Header.SpectrumAxis.CDELT', 'f8')
        self.nc.variables['Header.SpectrumAxis.CDELT'][0] = L.dxdc
        
        t = self.nc.createVariable('Header.SpectrumAxis.CRPIX', 'f8')
        self.nc.variables['Header.SpectrumAxis.CRPIX'][0] = 0.0

        t = self.nc.createVariable('Header.SpectrumAxis.CRVAL', 'f8')
        self.nc.variables['Header.SpectrumAxis.CRVAL'][0] = L.xarray[0]

        ctype = self.nc.createVariable('Header.SpectrumAxis.CTYPE', 'c', 
                                       ('nlabel',))
        ctype[0:len(L.xname)] = L.xname[0:len(L.xname)]

        xarray = self.nc.createVariable('Header.SpectrumAxis.CAXIS', 'f8', 
                                        ('nchan',))
        xarray[0:L.nchan] = L.xarray[0:L.nchan]

    def read_line_header_variables(self, L):
        """
        Reads line header variables of L (Line) object.
        Args:
            L (object): Line object
        Returns:
            none
        """
        L.nchan = self.nc.variables['Header.Line.NChannels'][0]
        L.iarray = self.nc.variables['Header.Line.ChannelNumber'][:]
        # the following are the generic spectrum axis parameters
        L.dxdc = self.nc.variables['Header.SpectrumAxis.CDELT'][0]
        L.xname = netCDF4.chartostring(self.nc.variables['Header.SpectrumAxis\
                                                          .CTYPE'][:])
        L.xarray = self.nc.variables['Header.SpectrumAxis.CAXIS'][:]


