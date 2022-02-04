"""
Module to read a set of spectra from the spectrometer

classes: RoachSpec, SpecBank, SpecBankData, SpecBankCal
methods: lookup_roach_files, find_roach_from_pixel, create_roach_list
uses:    IFProc, Grid
author:  FPS
date:    May 2018
changes:  
python 3
"""
import numpy as np
import math
import os
import sys
import fnmatch
import netCDF4
from scipy import interpolate
from operator import itemgetter
from itertools import groupby

from lmtslr.grid.grid import Grid
from lmtslr.ifproc.ifproc import IFProc

# define all the pixels in the roach boards they appear in
roach_pixels_all = [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]]



class RoachSpec():
    """
    Base class to deal with a single time series of spectra.
    """
    def __init__(self, obsnum, roach_id, roach_input, nchan, bandwidth, nspec,
                 raw_spec, spec_time, xmap, ymap, pmap, bufpos):
        """
        Constructor for RoachSpec class.
        Args:
           obsnum (int): target observation number
           roach_id (int): Roach ID number
           roach_input (): not used 
           nchan (int): number of channels
           bandwidth (): channel bandwidth
           nspec (int): number of spectra
           raw_spec (array): raw spectrum data
           spec_time (array): masked array of data timestamps
           xmap (array): not used
           ymap (array): not used
           pmap (array): not used
           bufpos (int): type of observation
        Returns:
            none
        """
        self.obsnum = obsnum
        self.roach_id = roach_id
        self.roach_input = roach_input
        self.nchan = nchan
        self.bandwidth = bandwidth
        self.nspec = nspec
        self.raw_spec = raw_spec
        self.spec_time = spec_time
        self.xmap = xmap
        self.ymap = ymap
        self.pmap = pmap
        self.bufpos = bufpos
        self.tsys_aver = False
        self.tsyscal = None
        self.ncal = 0

        self.ons,  self.on_ranges,  self.nons  = self.get_ranges(0)
        self.refs, self.ref_ranges, self.nrefs = self.get_ranges(1)
        self.skys, self.sky_ranges, self.nskys = self.get_ranges(2)
        self.hots, self.hot_ranges, self.nhots = self.get_ranges(3)

        if False:
            # some debugging on times when bufpos changed and how many
            nt = len(spec_time)
            print("RoachSpec: %d %d - %d %d %d %d" % (roach_id,nt,self.nons,self.nrefs,self.nskys,self.nhots))
            time0 = spec_time[0]
            print("TB",0,time0,bufpos[0],"first")
            for i in range(1,nt):
                if bufpos[i] != bufpos[i-1]:
                    print("TB",i,spec_time[i]-time0,bufpos[i])
            i = nt-1
            print("TB",i,spec_time[i]-time0,bufpos[i],"last")
            

    def get_ranges(self, value=1):
        """
        Searches bufpos and finds indices with values corresponding to 
        value. 
        Args:
            value (int): type of observation for search.
                         0 is ons
                         1 is refs
                         2 is sky
                         3 is hot
        Returns:
            (idx (list), ranges (list), len(ranges) (int)): idx is the 
                indices, ranges is the list of (begin,end), len(ranges) is the 
                length of ranges.
        """
        if value == 0:
            idx = np.where((self.bufpos == value) | (self.bufpos >= 100))[0]
        else:
            idx = np.where(self.bufpos == value)[0]
        ranges = []
        for k, g in groupby(enumerate(idx), lambda x:x[0] - x[1]):  
            group = list(map(itemgetter(1), g))
            ranges.append((group[0], group[-1]))
        return(idx, ranges, len(ranges))


    def compute_reference_spectra(self):
        """
        Computes the average reference spectrum for each block of refs,
        which are returned as a list.
        Args:
            none
        Returns:
            none
        """
        self.reference_spectra = np.zeros((self.nrefs, self.nchan))
        for iref in range(self.nrefs):
            the_refs = range(self.ref_ranges[iref][0], 
                             self.ref_ranges[iref][1] + 1)
            self.reference_spectra[iref, :] = np.mean(
                   self.raw_spec[the_refs, :], axis=0)

    def compute_main_spectra(self):
        """
        Computes the average main (on) spectrum for each block of on-
        source observations, which are returned as a list.
        Args:
            none
        Returns:
            none
        """
        self.main_spectra = np.zeros((self.nons, self.nchan))
        for ion in range(self.nons):
            the_ons = range(self.on_ranges[ion][0], self.on_ranges[ion][1] + 1)
            self.main_spectra[ion, :] = np.mean(self.raw_spec[the_ons, :], axis=0)

    def compute_median_spectrum(self):
        """
        Computes the median spectrum from ALL spectra in the file.
        Args:
            none
        Returns:
            none
        """
        self.median_spectrum = np.median(self.raw_spec[:,:],axis=0)

    def compute_reference_spectrum(self):
        """
        Computes the average reference spectrum using all ref 
        observations.
        Args:
            none
        Returns:
            none
        """
        self.reference_spectrum = np.mean(self.raw_spec[self.refs, :],axis=0)
        
    def compute_main_spectrum(self):
        """
        Computes the average main spectrum using all main (on) 
        observations.
        Args:
            none
        Returns:
            none
        """
        self.main_spectrum = np.mean(self.raw_spec[self.ons, :], axis=0)

    def compute_tsys_spectra(self, bdrop=100, edrop=100):
        """
        Computes the TSYS, typically for an otf_cal=1
        """
        self.tsys_spectra = np.zeros((self.nhots, self.nchan))
        tsys = np.zeros(self.nhots)
        for ihot in range(self.nhots):
            hot_spectrum = np.median(self.raw_spec[self.hot_ranges[ihot][0]:self.hot_ranges[ihot][1], :], axis=0)
            sky_spectrum = np.median(self.raw_spec[self.sky_ranges[ihot][0]:self.sky_ranges[ihot][1], :], axis=0)
            # @todo why 280, and not a measured ambient?            
            tsys_spec = 280.0 * sky_spectrum/(hot_spectrum - sky_spectrum)
            # find the index where tsys_spec is finite
            indx_fin = np.where(np.isfinite(tsys_spec))
            # compute tsys as mean of tsys_spec
            tsys[ihot] = np.mean(tsys_spec[indx_fin][bdrop:self.nchan-edrop])
            # find index where tsys_spec is infinte
            indx_inf = np.where(np.isinf(tsys_spec))
            # replace with mean
            tsys_spec[indx_inf] = tsys[ihot]
            if self.tsys_aver:
                self.tsys_spectra[ihot, :] = tsys[ihot]
            else:
                self.tsys_spectra[ihot, :] = tsys_spec
            #print("TSYS: ",self.tsys_spectra[ihot, :],tsys[ihot])
        # can we use pixel= this way?
        pixel = 4*self.roach_id + self.roach_input
        print("TSYS[%d] otf_cal %s" % (pixel,repr(tsys)))
            
        
    def get_nearest_reference(self, index, left=True):
        """
        Given a dump index and if it is to the left of ONs (default)
        or right of ONs returns the correct reference index
        """
        if left:
            arr = [abs(index-r[1]) for r in self.ref_ranges]
            return arr.index(min(arr))
        else:
            arr = [abs(index-r[0]) for r in self.ref_ranges]
            return arr.index(min(arr))            

    def get_previous_hot(self, index, nearest=False):
        """
        Given a dump index finds the previous HOT
        and returns the correct hot index
        
        Using nearest=True an older version of the code
        where "accidentally?" the nearest was used, can be returned.
        """
        if nearest:
            # nearest_hot (original code, arguably wrong)
            arr = [abs(index-r[1]) for r in self.hot_ranges]
            return arr.index(min(arr))
        else:            
            # true previous_hot (after Heyer)
            arr = [(index-r[1]) for r in self.hot_ranges]
            arr=np.array(arr)
            idx=np.where(arr > 0)
            return np.argmin(arr[idx])

    
    def reduce_on_spectrum(self, calibrate=False, tsys_spectrum=0,
                           tsys_no_cal=1):
        """
        Creates a ON spectrum returned as self.on_spectrum. Reduction 
        procedure depends on the stype parameter.
        Args:
            calibrate (bool): True when we want to calibrate the 
                ps_spectrum (default is False)
            tsys_spectrum (float): tsys values for the calibration used
                when calibrate=True (default is 0)
            tsys_no_cal (float): tsys value for the calibration used 
                when calibrate=False (default is 1)
        """
        self.compute_main_spectrum()
        self.on_spectrum = self.main_spectrum
        # calibrate it if requested
        if calibrate == True:
            self.on_spectrum = self.on_spectrum * tsys_spectrum
        else:
            self.on_spectrum = self.on_spectrum * tsys_no_cal

    def reduce_ps_spectrum(self, stype=2, normal_ps=True, calibrate=False, 
                           tsys_spectrum=0, tsys_no_cal=1):
        """
        Creates a PS spectrum returned as self.ps_spectrum. Reduction 
        procedure depends on the stype parameter.
        Args:
            stype (int): type of reduction to run (default is 2)
                        0: not defined
                        1: compute average of all refs and mains and 
                           take the difference
                        2: compute individual on's and off's in blocks;
                           difference each pair and average differences
            normal_ps (bool): True for (Main-Ref)/Ref 
                              False for (Ref-Main)/Main (beam switch)
                              (default is True)
            calibrate (bool): option to calibrate ps_spectrum (default 
                is False)
            tsys_spectrum (array): spectrum of tsys values for 
                calibration when calibrate=True
            tsys_no_cal (array): spectrum of tsys values for 
                calibration when calibrate=False
        Returns:
            none
        """
        if stype == 1:
            self.compute_main_spectrum()
            self.compute_reference_spectrum()
            if normal_ps == True:
                self.ps_spectrum = (self.main_spectrum - 
                                    self.reference_spectrum) \
                                   / self.reference_spectrum
            else:
                self.ps_spectrum = (self.reference_spectrum - 
                                    self.main_spectrum) / self.main_spectrum
        else: # type == 2:
            self.compute_main_spectra()
            self.compute_reference_spectra()
            if self.nons == self.nrefs:
                ps_list = np.zeros((self.nons,self.nchan))
                if normal_ps == True:
                    for i in range(self.nons):
                        ps_list[i,:] = (self.main_spectra[i,:] - 
                                        self.reference_spectra[i,:]) \
                                       / self.reference_spectra[i,:]
                else:
                    for i in range(self.nons):
                        ps_list[i,:] = (self.reference_spectra[i,:] - 
                                        self.main_spectra[i,:]) \
                                       / self.main_spectra[i,:]
                self.ps_spectrum = np.mean(ps_list[:,:],axis=0)
            else:
                print('check number of ons and refs %d %d'%(self.nons, 
                                                            self.nrefs))
            
        # calibrate it if requested
        if calibrate == True:
            self.ps_spectrum = self.ps_spectrum * tsys_spectrum
        else:
            self.ps_spectrum = self.ps_spectrum * tsys_no_cal

    def reduce_spectra(self, stype=0, calibrate=False,
                       tsys_spectrum=0, tsys_no_cal=1,
                       use_otf_cal=False):
        """
        Creates a list of all the "on" spectra in the file. The on 
        spectra are reduced according to the value of "type".
        Args:
            stype (int): type of reduction to run (default is 0)
                        0: use the median spectrum for the whole thing
                        1: use a single reference spectra which is 
                           average of all refs
                        2: use the average of refs which bracket the 
                           ons
                        3: use weighted average of refs which bracket the
                           ons
            calibrate (bool): option to calibrate ps_spectrum (default 
                is False)
            tsys_spectrum (float): tsys value for calibration when 
                calibrate=True
            tsys_no_cal (float): tsys value for calibration when 
                calibrate=False
            use_otf_cal (bool): option to use in-scan hot and sky to derive cals
                (default is False)
        Returns:
            none
        """
        spectra = []

        if self.nrefs == 0:
            type = 0

        if use_otf_cal:
            self.compute_tsys_spectra()
            
        if stype == 0:
            self.compute_median_spectrum()
            for i in self.ons:
                spectra.append((self.raw_spec[i,:] - self.median_spectrum[:]) 
                                / self.median_spectrum[:])
            #if use_otf_cal:
            #    self.compute_tsys_spectrum()
            #    self.reduced_spectra = np.array(spectra) * self.tsys_spectrum
        elif stype == 1:
            if self.nrefs != 0:
                self.compute_reference_spectrum()
                for i in self.ons:
                    spectra.append((self.raw_spec[i,:] - 
                                    self.reference_spectrum[:]) / 
                                   self.reference_spectrum[:])
            else:
                for i in self.ons:
                    spectra.append((self.raw_spec[i,:]))
            #if use_otf_cal:
            #    self.compute_tsys_spectrun()
            #    self.reduced_spectra = np.array(spectra) * self.tsys_spectrum
        elif stype == 2: # type == 2:
            if self.nrefs != 0:
                self.compute_reference_spectra()
                #nbins = self.nrefs - 1
                nbins = self.nons
                for ibin in range(nbins):
                    istart = self.on_ranges[ibin][0]
                    start_ref_bin = self.get_nearest_reference(istart-1, left=True)
                    istop = self.on_ranges[ibin][1]
                    stop_ref_bin = self.get_nearest_reference(istop+1, left=False)
                    if use_otf_cal:
                        tsys_bin = self.get_previous_hot(istart-1)
                        # print("OTF_CAL %d %d %d %d %d -> %d" % (ibin,istart,istop,start_ref_bin,stop_ref_bin,tsys_bin))
                    for i in range(istart, istop + 1):
                        ref = (self.reference_spectra[start_ref_bin] + 
                               self.reference_spectra[stop_ref_bin]) / 2
                        if use_otf_cal:
                            spectra.append( ((self.raw_spec[i, :] - ref) / ref) * self.tsys_spectra[tsys_bin] )
                        else:
                            spectra.append((self.raw_spec[i,:] - ref) / ref)
            else:
                for i in self.ons:
                    spectra.append((self.raw_spec[i,:]))
        elif stype == 3: # type == 3:
            if self.nrefs != 0:
                self.compute_reference_spectra()
                #nbins = self.nrefs - 1
                nbins = self.nons
                for ibin in range(nbins):
                    istart = self.on_ranges[ibin][0]
                    start_ref_bin = self.get_nearest_reference(istart-1, left=True)                    
                    istop = self.on_ranges[ibin][1]
                    stop_ref_bin = self.get_nearest_reference(istop+1, left=False)                    
                    mapwidth = istop - istart
                    if use_otf_cal:
                        tsys_bin = self.get_previous_hot(istart-1)                    
                    for i in range(istart, istop + 1):
                        wt1 = float(mapwidth - (i - istart))/float(mapwidth)
                        wt2 = 1 - wt1
                        ref = (wt1 * self.reference_spectra[start_ref_bin] + 
                               wt2 * self.reference_spectra[stop_ref_bin])
                        if use_otf_cal:
                            spectra.append( ((self.raw_spec[i,:] - ref) / ref) * self.tsys_spectra[tsys_bin] )
                        else:
                            spectra.append((self.raw_spec[i,:] - ref) / ref)
            else:
                for i in self.ons:
                    spectra.append((self.raw_spec[i,:]))
                        
        # save reduced spectra as a 2D numpy array
        # calibrate it if requested
        # do special things if use_otf_cal is True
        if use_otf_cal:
            if stype in (0, 1):
                self.reduced_spectra = np.array(spectra) * self.tsys_spectrum
            else:
                self.reduced_spectra = np.array(spectra)
        else:
            if calibrate == False:
                self.reduced_spectra = np.array(spectra) * tsys_no_cal
            else:
                self.reduced_spectra = np.array(spectra) * tsys_spectrum

    def baseline(self, spectrum, baseline_list, n_baseline_list,
                 baseline_order=0):
        """
        Computes a baseline given a spectrum, using a list of channels.
        Args:
            spectrum (array): array of spectrum data
            baseline_list (list): list of channels to use
            n_baseline_list (int): number of channels
            baseline_order (int): order of fitting function (default is
                0)
        Returns:
            (baseline (array), rms (array)): baseline is the array of 
                computed baseline values, rms is the array of root-
                mean-square error
        """
        baseline_list = [i for i in baseline_list if i < len(spectrum)]
        params = np.polyfit(baseline_list, spectrum[baseline_list], 
                            baseline_order)
        x = np.arange(0, self.nchan)
        baseline = np.polyval(params, x)
        residuals = spectrum[baseline_list] - np.polyval(params, 
                                                         baseline_list)
        rms = np.sqrt(np.dot(residuals.transpose(), residuals) /  
                      n_baseline_list)
        return(baseline, rms)

    def integrate_spectra(self, channel_list, n_channel_list, baseline_list, 
                          n_baseline_list, baseline_order=0, type=0):
        """
        Computes the integral of the spectrum over some range of 
        channels (channel_list) with baseline removed using 
        baseline_list.
        Args:
            channel_list (list): list of channels to use
            n_channel_list (int): number of channels to use
            baseline_list (list): list of baselines to use
            n_baseline_list (int): number of baselines to use
            baseline_order (int): order of fitting function (default is
                0)
            type (int): type of integration to use. 0 is YINT and 1 is
                YMAX.
        Returns:
            np.array(spectra) (array): array of integrated spectra 
                data
        """
        spectra = []
        for i,isample in enumerate(self.ons):
            baseline,rms = self.baseline(self.reduced_spectra[i], 
                                         baseline_list, n_baseline_list, 
                                         baseline_order)
            baselined_spectrum = self.reduced_spectra[i] - baseline
            if type == 0:
                spectra.append(np.sum(baselined_spectrum[channel_list]))
            else:
                spectra.append(np.max(baselined_spectrum[channel_list]))

        return(np.array(spectra))

    def integrate_spectrum(self, on_list, channel_list, n_channel_list,
                           baseline_list, n_baseline_list, baseline_order=0, 
                           type=0):
        """
        Averages reduced spectra over list of indices (on_list) and 
        then computes the integral of the spectrum over some range of 
        channels (channel_list) with baseline removed using 
        baseline_list.
        Args:
            on_list (list): list of indices to use
            channel_list (list): list of channels to use
            n_channel_list (int): number of channels to use
            baseline_list (list): list of baselines to use
            n_baseline_list (int): number of baselines to use
            baseline_order (int): order of fitting function (default is
                0)
            type (int): type of integration to use. 0 is YINT and 1 is
                YMAX.
        Returns:
            (baselined_spectrum (array), result (float)): 
                baselined_spectrum is the spectrum with baseline 
                subtracted, result is the integrated spectrum.
        """
        spectrum = np.mean(self.reduced_spectra[on_list,:], axis=1)[0]
        #print(on_list)
        #print(np.mean(self.xmap[on_list]),np.mean(self.ymap[on_list]))
        #print(spectrum)
        baseline,rms = self.baseline(spectrum, baseline_list, n_baseline_list,
                                     baseline_order)
        baselined_spectrum = spectrum - baseline 
        if type == 0:
            result = np.sum(baselined_spectrum[channel_list]) # YINT
        else:
            result = np.max(baselined_spectrum[channel_list]) # YMAX
            
        return(baselined_spectrum, result)

    def compute_tsys_spectrum(self, bdrop=100, edrop=100):
        """
        Computes a system temperature spectrum and average value of 
        tsys.
        Args:
            bdrop (int): number of channels excluded from beginning of 
                spectrum
            edrop (int): number of channels excluded from end of 
                spectrum
        Returns:
            none
        """
        if ((self.nhots>0) and (self.nskys>0)):
            hot_spectrum = np.median(self.raw_spec[self.hots,:], axis=0)
            sky_spectrum = np.median(self.raw_spec[self.skys,:], axis=0)
            # @todo why 280, and not a measured ambient?
            self.tsys_spectrum = 280 * sky_spectrum / (hot_spectrum - 
                                                       sky_spectrum)
            # find the index where tsys_spectrum in finite
            indx_fin = np.where(np.isfinite(self.tsys_spectrum))
            # compute tsys as the mean of finite tsys_spectrum
            self.tsys = np.mean(self.tsys_spectrum[indx_fin][bdrop:self.nchan - edrop])
            tsysstd = np.std(self.tsys_spectrum[indx_fin][bdrop:self.nchan - edrop])
            # find the index where tsys_spectrum in infinite
            indx_inf = np.where(np.isinf(self.tsys_spectrum))
            # replace infinite tsys_spectrum with the mean
            self.tsys_spectrum[indx_inf] = self.tsys
            if self.tsys_aver:
                self.tsys_spectrum[:] = self.tsys
                
            #print("TSYS: ",self.tsys_spectrum[:])
            # @todo can we use pixel= this way?
            pixel = 4*self.roach_id + self.roach_input
            print("TSYS[%d] = %g +/- %g (%d channels)" % (pixel,self.tsys, tsysstd, len(self.tsys_spectrum)))
        else:
            print('ObsNum %d Roach %d does not have calibration data'%(
                self.obsnum, self.roach_id))
            self.tsys_spectrum = np.zeros(self.nchan)
            self.tsys = 0.

    class LineStatistics():
        """
        Class to handle line statistics.
        """
        def __init__(self, parent, v, spectrum, channel_list, n_channel_list,
                     baseline_list, n_baseline_list, baseline_order=0, 
                     pixel_id_0=0, pixel_id_1=1, obspgm=None):
            """
            Constructor for LineStatistics class inside RoachSpec class.
            Removes a baseline and computes line statistics for a 
            spectrum.
            Line Statistics are:
            YMAX = maximum value in channel_list
            CMAX = channel of maximum value
            XMAX = velocity of maximum value
            YINT = integral of line over channel list in units of K km/s
            YERR = error on YINT
            XMEAN = first moment of line integrated over channel list
            XWIDTH = estimate of line width based on Peak value and 
                integrated intensity
            RMS = rms of baseline fit

            Args:
                v (array): array of velocities
                spectrum (array): array of spectrum data
                channel_list (list): list of channels to use
                n_channel_list (int): number of channels to use
                baseline_list (list): list of baselines to use
                n_baseline_list (int): number of baselines to use
                baseline_order (int): order of fitting function 
                    (default is 0)
                pixel_id_0 (int): ID of 0th target pixel (default is 0)
                pixel_id_1 (int): ID of 1st target pixel (default is 1)
                obspgm (str): obspgm from IFProc data
            Returns:
                none
            """
            self.pixel_id_0 = pixel_id_0
            self.pixel_id_1 = pixel_id_1
            self.obspgm = obspgm
            delta_v = np.abs(v[1] - v[0])
            baseline, rms = parent.baseline(spectrum, baseline_list, n_baseline_list, 
                                    baseline_order)
            baselined_spectrum = spectrum - baseline
            self.v = v
            self.spectrum = baselined_spectrum
            self.rms = rms
            self.ymax = np.max(baselined_spectrum[channel_list])
            x = np.arange(0, parent.nchan)
            self.cmax = x[channel_list[np.where(
                baselined_spectrum[channel_list] == self.ymax)[0][0]]]
            self.xmax = v[self.cmax]
            self.yint = np.sum(baselined_spectrum[channel_list]) * delta_v
            self.yerr = delta_v * self.rms * np.sqrt(n_channel_list)
            self.xmean = np.sum(baselined_spectrum[channel_list] 
                                * v[channel_list]) / n_channel_list
            self.xwidth = self.yint / self.ymax

        def to_string(self):
            """
            Returns the line statistics represented in a string.
            Args:
                none
            Returns:
                str (str): string representing the line statistics 
                    values
            """
            str = '%s Pix=%02d'%(self.obspgm, self.pixel_id_0)
            if self.pixel_id_0 != self.pixel_id_1:
                str += '/%02d'%self.pixel_id_1
            str += ' ymax=%.3f cmax=%d xmax=%.3f yint=%.3f yerr=%.3f\
                     xmean=%.3f xwidth=%.3f rms=%.3f'%(self.ymax, self.cmax, 
                        self.xmax, self.yint, self.yerr, self.xmean, 
                        self.xwidth, self.rms)
            return str


class SpecBank():
    """
    Base class for dealing with a complete "bank" of spectra.
    """
    def __init__(self, roach_files, ifproc_data, 
                 pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], 
                 time_offset=[-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,
                              -0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03], 
                 bank=0, save_tsys=False):
        """
        Constructor for SpecBank class.
        Args: 
            roach files (list): list of files for ROACH boards 
                (nominally 4)
            ifproc_data (object): data from corresponding ifproc file
            pixel_list (list): list of pixels to process (default is 
                all)
            bank (int): bank number for receiver
            time_offset (float): time lag for each ROACH board in units
                of seconds (s)
        Returns:
            none
        """
        #self.ifproc = ifproc_data(ifproc_file)
        self.ifproc = ifproc_data
        # most header variables in ifproc object; copy these two since often used
        self.obsnum = self.ifproc.obsnum
        self.calobsnum = self.ifproc.calobsnum
        self.elev = self.ifproc.elev
        self.obspgm = self.ifproc.obspgm
        self.source = self.ifproc.source
        self.vlsr = self.ifproc.vlsr              # PJT
        self.date_obs = self.ifproc.date_obs      # PJT
        self.map_coord = self.ifproc.map_coord

        # timing offsets for each roach board
        self.time_offset = time_offset

        self.bank = bank
        
        self.x_interpolation_function = interpolate.interp1d(self.ifproc.time,
            self.ifproc.azmap, bounds_error=False)
        self.y_interpolation_function = interpolate.interp1d(self.ifproc.time,
            self.ifproc.elmap, bounds_error=False)
        self.p_interpolation_function = interpolate.interp1d(self.ifproc.time,
            self.ifproc.parang, bounds_error=False)
        self.b_interpolation_function = interpolate.interp1d(self.ifproc.time,
            self.ifproc.bufpos, kind='nearest', bounds_error=False)

        self.gaps = np.where(self.ifproc.bufpos[:-1] != \
                    self.ifproc.bufpos[1:])[0]
        self.time_gaps_0 = np.append([self.ifproc.time[0]], \
                           self.ifproc.time[self.gaps + 1])
        self.time_gaps_1 = np.append(self.ifproc.time[self.gaps], 
                                     [self.ifproc.time[-1]])

        self.nroach = len(roach_files)
        # No need to report here.
        # if self.nroach != len(roach_pixels_all):
        #    print("Warning: only %d/%d roach files found" % (self.nroach,len(roach_pixels_all)))

        self.roach_pixel_ids = []
        self.roach = []
        for i in range(self.nroach):
            self.roach_pixel_ids = self.roach_pixel_ids + \
                                   self.read_roach(roach_files[i], i, 
                                                   pixel_list)
        self.roach_pixel_ids = np.array(self.roach_pixel_ids)
        self.npix = len(self.roach)
        self.nchan = self.roach[0].nchan
        self.bandwidth = self.roach[0].bandwidth
        self.channel_0 = (self.nchan - 1) / 2
        self.velocity_0 = self.ifproc.velocity
        self.line_rest_frequency = self.ifproc.line_rest_frequency[bank]
        self.receiver = self.ifproc.receiver
        self.sideband = self.ifproc.sideband[bank]
        self.tracking_beam = self.ifproc.tracking_beam

        # dfdc sign depends on the sideband (USB:- and LSB+)
        if self.sideband == 1:
            self.dfdc = -self.bandwidth*np.float64(1e6) / self.nchan
        else:
            self.dfdc = +self.bandwidth*np.float64(1e6) / self.nchan
        self.dvdc = -self.dfdc / (self.line_rest_frequency * np.float64(1e9)) \
                    * np.float64(299792.458)

        self.cal_flag = False
        self.save_tsys = save_tsys
        
        # define the windows for line integrations and baseline fits
        self.nc = 0
        self.clist = []
        self.baseline_order = 0
        self.nb = 0
        self.blist = []

        
    def make_pixel_id_list(self, roach_files, pixel_list):
        """
        Creates a list of pixel id numbers corresponding to list of roach 
        spectra that is created.
        Args:
            roach files (list): list of files for ROACH boards 
                (nominally 4)
            pixel_list (list): list of pixels to process
        Returns:
            np.array(l (list)): l is the list of pixel ids
        """
        l = []
        for f in roach_files:
            for i in range(4):
                roach_name = 'roach%d'%i
                if roach_name in f:
                    for p in roach_pixels_all[i]:
                        if p in pixel_list:
                            l = l + [p]
                    break
        return(np.array(l))

    def make_channel_list(self, regions):
        """
        Creates a list of channels given a list of lists of intervals 
        to be processed.
        for example: regions = [[1,3],[5,8]] creates [1,2,3,5,6,7,8]
        Args:
            regions (list): list of lists of intervals to be processed
        Returns:
            (channel_list (list), nchannels (int)): channel_list is the
                list of channels, nchannels is the number of channels
        """
        nregions = len(regions)
        channel_list = []
        for i in range(nregions):
            # python 3 requires range to be converted to list
            channel_list = channel_list + list(range(regions[i][0], 
                                               regions[i][1] + 1))
        nchannels = len(channel_list)
        return(channel_list, nchannels)

    def make_velocity_list(self, velocity_regions, id='line'):
        """
        Creates a channel list from a set of velocity intervals.
        Args:
            velocity_regions (list): set of velocity intervals.
            id (str): type of line
                      'line' designates line (default)
                      'baseline' designates baseline
        Returns:
            (channel_list (list), nchannels (int)): channel_list is 
                the list of channels, nchannels is the number of 
                channels
        """
        nregions = len(velocity_regions)
        channel_list = []
        for i in range(nregions):
            c0 = self.v2c(velocity_regions[i][0])
            c1 = self.v2c(velocity_regions[i][1])+1
            # python 3 requires range to be converted to list
            if c1>c0:
                channel_list = channel_list + list(range(c0, c1))
            else:
                channel_list = channel_list + list(range(c1, c0))
        nchannels = len(channel_list)
        if id == 'line':
            self.clist = channel_list
            self.nc = nchannels
        elif id == 'baseline':
            self.blist = channel_list
            self.nb = nchannels
        return(channel_list, nchannels)

    def v2c(self, v):
        """
        Converts velocity into channel number.
        Args:
            v (float): velocity in units of 
        Returns:
            c (array): array of channel numbers corresponding to 
                velocity
        """
        c = np.array(np.round((v - self.velocity_0) / self.dvdc 
                              + self.channel_0), dtype=int)
        return(c)
    
    def c2v(self, c):
        """
        Converts a channel number into a velocity.
        Args:
            c (int or array): channel number(s)
        Returns:
            v (float or array): velocity(ies) corresponding to channel 
                in units of 
        """
        v = (c - self.channel_0) * self.dvdc + self.velocity_0
        return(v)

    def create_velocity_scale(self):
        """
        Computes the velocity scale corresponding to the channels in 
        the spectrum.
        Args:
            none
        Returns:
            self.c2v(range(self.nchan)) (array): velocity scale 
                corresponding to channels
        """
        return(self.c2v(range(self.nchan)))

    def find_pixel_index(self, ipixel):
        """
        Looks up the index in our array of spectra corresponding to a 
        specific pixel.
        Args:
            ipixel (int): index of pixel
        Returns:
            (int(index[0]) (int)) if len(index[0]) else 0: index of 
                pixel in spectra. If not found, returns 0.
        """
        index = np.where(self.roach_pixel_ids == ipixel)
        return (int(index[0])) if len(index[0]) else 0

    def find_map_pixel_index(self, ipixel):
        """
        Looks up the index in our map of spectra corresponding to a 
        specific pixel.
        Args:
            ipixel (int): index of pixel
        Returns:
            (int(index[0]) (int)) if len(index[0]) else 0: index of 
                pixel in map. If not found, returns 0.
        """
        index = np.where(self.map_pixel_list == ipixel)
        return (int(index[0])) if len(index[0]) else 0

    def read_roach(self, filename, roach_id, pixel_list, as_float=True):
        """
        Reads a roach file and creates spectrum (roach class) for each 
        input.
        Args:
            filename (str): name of the roach file
            roach_id (int): number of the roach. e.g. roach0 would be 0
            pixel_list (list): list of pixels we want to process
            as_float (bool): if True, uses float instead of double.
        """
        pixel_ids = []
        if os.path.isfile(filename):
            nc = netCDF4.Dataset(filename)

            # header information
            obsnum = nc.variables['Header.Telescope.ObsNum'][0]
            nchan = nc.variables['Header.Mode.numchannels'][0]
            bandwidth = nc.variables['Header.Mode.Bandwidth'][0]
            ninputs = 4

            datatime = nc.variables['Data.Integrate.time'][:]
            if as_float:
                rawdata = np.asarray(nc.variables['Data.Integrate.Data'][:],np.float32)
            else:
                rawdata = nc.variables['Data.Integrate.Data'][:]
            input_list = nc.variables['Data.Integrate.Inputs'][:]

            acc_n = nc.variables['Data.Integrate.acc_n'][:]
            sync_time = nc.variables['Data.Integrate.sync_time'][:]
            read_time = nc.variables['Data.Integrate.read_time'][:]
            datatime = datatime - read_time

            print('read_roach %s     nspec,nchan=%d,%d' %
                  (filename, rawdata.shape[0], rawdata.shape[1]))
            
            # get roach index
            roach_index = roach_id
            for i in range(4):
                roach_name = 'roach%d'%i
                if roach_name in filename:
                    roach_index = i
                    break

            for input_chan in range(ninputs):
                # check to see whether this input matches one of the \
                # pixels we would like
                # if no match, then don't process it into the list
                if roach_pixels_all[roach_index][input_chan] not in pixel_list:
                    continue
                # if there is a match, then add the pixel id to our list and process
                pixel_ids = pixel_ids + [roach_pixels_all[roach_index][input_chan]]
                ilist  = np.where(input_list == input_chan)
                ilist0  = np.where(input_list == 0)
                raw_spec = rawdata[ilist,:][0]

                print('r:%d inp:%d pix:%d time_offset:%f'%(roach_index, input_chan, 
                    roach_pixels_all[roach_index][input_chan], 
                    self.time_offset[roach_pixels_all[roach_index][input_chan]]))

                spec_time = datatime[ilist0] + \
                    self.time_offset[roach_pixels_all[roach_index][input_chan]]
                nspec = len(spec_time)
                # use the interpolation functions to find map positions,\
                # paralactic angle, and bufpos
                # in python 3 we must fix spec_time which was a masked \
                # array
                xmap = self.x_interpolation_function(
                    np.ma.getdata(spec_time, subok=False))
                ymap = self.y_interpolation_function(
                    np.ma.getdata(spec_time, subok=False))
                pmap = self.p_interpolation_function(
                    np.ma.getdata(spec_time, subok=False))
                bufpos = self.b_interpolation_function(
                    np.ma.getdata(spec_time, subok=False))

                # correct the interpolated arrays to remove points not\
                # in range for interpolation
                l = len(self.time_gaps_1)
                cond = False
                for i in range(l):
                    cond = cond | ((spec_time > self.time_gaps_0[i]) & 
                                   (spec_time < self.time_gaps_1[i]))
                spec_time = spec_time[cond]
                xmap = xmap[cond]
                ymap = ymap[cond]
                pmap = pmap[cond]
                bufpos = bufpos[cond].astype(int)
                raw_spec = raw_spec[cond]
                
                nspec = len(spec_time)
                indx_small, = np.where(raw_spec[:,1024] < 100)
                if indx_small.any() > 0:
                    print('low_power', obsnum, roach_index, input_chan, 
                          indx_small)
                # now we append the individual roach spectrum object to\
                # our list
                self.roach.append(RoachSpec(obsnum, roach_id, input_chan, 
                                            nchan, bandwidth, nspec, raw_spec,
                                            spec_time, xmap, ymap, pmap, 
                                            bufpos))
            nc.close()
        else:
            print('%s does not exist for roach_id=%d'%(self.filename, 
                                                       roach_id))
        return pixel_ids                            

class SpecBankData(SpecBank):
    """
    Base class to deal with a complete "bank" of spectra.
    """
    def __init__(self,roach_files, ifproc_data,
                 pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
                 time_offset=[-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,
                              -0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03],
                 bank=0, save_tsys=False):
        """
        Constructor for SpecBankData class.
        Args:
            roach files (list): list of files for ROACH boards 
                (nominally 4)
            ifproc_data (object): data from corresponding ifproc file
            pixel_list (list): list of pixels to process (default is 
                all)
            time_offset (list): list of time offsets for pixels
        Returns:
            none
        """
        SpecBank.__init__(self, roach_files, ifproc_data, 
                          pixel_list=pixel_list, time_offset=time_offset, 
                          bank=bank, save_tsys=save_tsys)
    
    def create_map_data(self, channel_list, n_channel_list, baseline_list, 
                        n_baseline_list, baseline_order=0, 
                        pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], 
                        type=0):
        """
        Processes a list of pixel ids to make a set of integrated 
        spectra for mapping (and fitting).
        Args:
            channel_list (list): list of channels to use
            n_channel_list (int): number of channels to use
            baseline_list (list): list of baselines to use
            n_baseline_list (int): number of baselines to use
            baseline_order (int): order of fitting function (default is
                0)
            pixel_list (list): list of pixels to process (default is 
                all)
            type (int): type of integration to use. 0 is YINT and 1 is
                YMAX.
        Returns:
            none
        """
        t_list = []
        data_list = []
        x_list = []
        y_list = []
        p_list = []
        n_list = []
        mp_list = []
        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            mp_list.append(ipix)
            t_list.append(self.roach[i].spec_time[self.roach[i].ons])
            x_list.append(self.roach[i].xmap[self.roach[i].ons])
            y_list.append(self.roach[i].ymap[self.roach[i].ons])
            p_list.append(self.roach[i].pmap[self.roach[i].ons])
            n_list.append(len(self.roach[i].xmap[self.roach[i].ons]))
            data_list.append(self.roach[i].integrate_spectra(channel_list, 
                n_channel_list, baseline_list, n_baseline_list, 
                baseline_order, type=type))
        self.map_pixel_list = np.array(mp_list)
        self.map_t = np.array(t_list)
        self.map_x = np.array(x_list)
        self.map_y = np.array(y_list)
        self.map_p = np.array(p_list)
        self.map_n = np.array(n_list)
        self.map_bufpos = np.array(n_list)
        self.map_data = np.array(data_list)

    def create_map_grid_bufpos_from_grid(self, xgrid, ygrid, tole, 
                                         channel_list, n_channel_list, 
                                         baseline_list, n_baseline_list, 
                                         baseline_order=0, 
                                         pixel_list=[0,1,2,3,4,5,6,7,8,9,10,
                                                     11,12,13,14,15], 
                                         type=0):
        """
        Processes a list of pixel ids to make a set of bufpos grids for
        mapping (and fitting).
        Args:
            xgrid (array): 
            ygrid (array): 
            tole (): 
            channel_list (list): list of channels to use
            n_channel_list (int): number of channels to use
            baseline_list (list): list of baselines to use
            n_baseline_list (int): number of baselines to use
            baseline_order (int): order of fitting function (default is
                0)
            pixel_list (list): list of pixels to process (default is 
                all)
            type (int): type of integration to use. 0 is YINT and 1 is
                YMAX.
        Returns:
            none
        """
        ngridx = len(xgrid)
        ngridy = len(ygrid)

        bufpos = 100
        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            for iypt in range(ngridy):
                ypt = ygrid[iypt]
                for ixpt in range(ngridx):
                    xpt = xgrid[ixpt]
                    #grid_list = np.where( np.sqrt( (self.roach[i].xmap[self.roach[i].ons]-xpt)**2+(self.roach[i].ymap[self.roach[i].ons]-ypt)**2) < tole)
                    grid_list = np.where(np.sqrt((self.roach[i].xmap - xpt)**2
                                                  + (self.roach[i].ymap - ypt)**2)
                                         < tole)
                    self.roach[i].bufpos[grid_list] = bufpos
                    #print(xpt,ypt,bufpos)
                    bufpos = bufpos + 1
        #np.set_printoptions(threshold=np.nan)
        #print(self.roach[i].bufpos)
                    
    def create_map_grid_data_from_grid(self, xgrid, ygrid, tole, channel_list,
                                       n_channel_list, baseline_list, 
                                       n_baseline_list, baseline_order=0, 
                                       pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,
                                                   12,13,14,15],
                                       type=0):
        """
        Processes a list of pixel ids to make a set of grid data for 
        mapping (and fitting).
        Args:
            xgrid (array): 
            ygrid (array): 
            tole (): 
            channel_list (list): list of channels to use
            n_channel_list (int): number of channels to use
            baseline_list (list): list of baselines to use
            n_baseline_list (int): number of baselines to use
            baseline_order (int): order of fitting function (default is
                0)
            pixel_list (list): list of pixels to process (default is 
                all)
            type (int): type of integration to use. 0 is YINT and 1 is
                YMAX.
        Returns:
            none
        """
        ngridx = len(xgrid)
        ngridy = len(ygrid)

        t_list = []
        data_list = []
        test_x_list = []
        test_y_list = []
        bufpos_list = []
        x_list = []
        y_list = []
        n_list = []
        mp_list = []
        spectrum_list = []

        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            mp_list.append(ipix)
            pt_list = []
            pn_list = []
            pbufpos_list = []
            px_list = []
            py_list = []
            ptest_x_list = []
            ptest_y_list = []
            pdata_list = []
            pspectrum_list = []
            for iypt in range(ngridy):
                ypt = ygrid[iypt]
                for ixpt in range(ngridx):
                    xpt = xgrid[ixpt]
                    #print(xpt,ypt)
                    grid_list = np.where(np.sqrt(
                        (self.roach[i].xmap[self.roach[i].ons] - xpt)**2 
                        + (self.roach[i].ymap[self.roach[i].ons] - ypt)**2) 
                        < tole)
                    pt_list.append(np.mean(
                        self.roach[i].spec_time[self.roach[i].ons][grid_list])
                                  )
                    pn_list.append(
                        len(self.roach[i].xmap[self.roach[i].ons][grid_list]))
                    pbufpos_list.append(np.mean(
                        self.roach[i].bufpos[self.roach[i].ons][grid_list]))
                    px_list.append(np.mean(
                        self.roach[i].xmap[self.roach[i].ons][grid_list]))
                    py_list.append(np.mean(
                        self.roach[i].ymap[self.roach[i].ons][grid_list]))

                    ptest_x_list.append(
                        self.roach[i].xmap[self.roach[i].ons][grid_list])
                    ptest_y_list.append(
                        self.roach[i].ymap[self.roach[i].ons][grid_list])


                    # z is the integrated intensity value for the averaged spectrum
                    zspec,z = self.roach[i].integrate_spectrum(grid_list, 
                        channel_list, n_channel_list, baseline_list, 
                        n_baseline_list, baseline_order, type)
                    pdata_list.append(z)
                    pspectrum_list.append(zspec)

            t_list.append(pt_list)
            n_list.append(len(xgrid)*len(ygrid))
            bufpos_list.append(pbufpos_list)
            x_list.append(px_list)
            y_list.append(py_list)
            test_x_list.append(ptest_x_list)
            test_y_list.append(ptest_y_list)

            
            data_list.append(pdata_list)
            spectrum_list.append(pspectrum_list)

        self.map_pixel_list = np.array(mp_list)
        self.map_t = np.array(t_list)
        self.map_bufpos = np.array(bufpos_list)
        self.map_x = np.array(x_list)
        self.map_y = np.array(y_list)
        self.map_test_x = np.array(test_x_list)
        self.map_test_y = np.array(test_y_list)
        self.map_n = np.array(n_list)
        self.map_data = np.array(data_list)
        self.map_spectra = np.array(spectrum_list)

    def create_map_grid_data(self, channel_list, n_channel_list, 
                             baseline_list, n_baseline_list, baseline_order, 
                             pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
                             type=0):
        """
        Processes a list of pixel ids to make a set of grid data for 
        mapping (and fitting).
        Args:
            channel_list (list): list of channels to use
            n_channel_list (int): number of channels to use
            baseline_list (list): list of baselines to use
            n_baseline_list (int): number of baselines to use
            baseline_order (int): order of fitting function (default is
                0)
            pixel_list (list): list of pixels to process (default is 
                all)
            type (int): type of integration to use. 0 is YINT and 1 is
                YMAX.
        Returns:
            none
        """
        t_list = []
        data_list = []
        bufpos_list = []
        x_list = []
        y_list = []
        n_list = []
        mp_list = []
        spectrum_list = []

        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            mp_list.append(ipix)
            pt_list = []
            pn_list = []
            pbufpos_list = []
            px_list = []
            py_list = []
            pdata_list = []
            pspectrum_list = []
            # find the grid point data
            b = np.bincount(self.roach[i].bufpos)
            e = np.where(b>0)[0]
            h = np.where(e>=100)[0]
            hlen = len(h)
            for ih in range(hlen):
                grid_list = np.where(
                    self.roach[i].bufpos[self.roach[i].ons] == e[h[ih]])
                pt_list.append(np.mean(
                    self.roach[i].spec_time[self.roach[i].ons][grid_list]))
                pn_list.append(len(
                    self.roach[i].xmap[self.roach[i].ons][grid_list]))
                pbufpos_list.append(np.mean(
                    self.roach[i].bufpos[self.roach[i].ons][grid_list]))
                px_list.append(np.mean(
                    self.roach[i].xmap[self.roach[i].ons][grid_list]))
                py_list.append(np.mean(
                    self.roach[i].ymap[self.roach[i].ons][grid_list]))
                # z is the integrated intensity value for the averaged spectrum
                zspec,z = self.roach[i].integrate_spectrum(grid_list, 
                    channel_list, n_channel_list, baseline_list, 
                    n_baseline_list, baseline_order, type)
                pdata_list.append(z)
                pspectrum_list.append(zspec)

            t_list.append(pt_list)
            n_list.append(hlen)
            bufpos_list.append(pbufpos_list)
            x_list.append(px_list)
            y_list.append(py_list)
            data_list.append(pdata_list)
            spectrum_list.append(pspectrum_list)

        self.map_pixel_list = np.array(mp_list)
        self.map_t = np.array(t_list)
        self.map_bufpos = np.array(bufpos_list)
        self.map_x = np.array(x_list)
        self.map_y = np.array(y_list)
        self.map_n = np.array(n_list)
        self.map_data = np.array(data_list)
        self.map_spectra = np.array(spectrum_list)
    

    def create_grid(self):
        """
        Creates the map grid from map parameters.
        Args:
            none
        Returns:
            none
        """
        # to allow arange to find the last point
        xepsilon = self.ifproc.xlength * 0.001
        xstep = self.xstep * self.hpbw
        self.xgrid = np.arange(-self.ifproc.xlength / 2, 
            self.ifproc.xlength / 2 + xepsilon, xstep)
        self.nx = len(xgrid)
        # to allow arange to find the last point
        yepsilon = self.ifproc.ylength * 0.001
        ystep = self.ystep * self.hpbw
        self.ygrid = np.arange(-self.ifproc.ylength / 2, 
            self.ifproc.ylength / 2 + yepsilon, ystep)
        self.ny = len(ygrid)
        self.ngrid = self.nx * self.ny

    def find_grid_location(self, ix, iy):
        """
        Returns grid location of grid point ix, iy.
        Args:
            ix (int): x index
            iy (int): y index
        Returns:
            self.xgrid[ix] (): x location
            self.ygrid[iy] (): y location
        """
        return(self.xgrid[ix], self.ygrid[iy])

    def find_grid_point_id(self, ix, iy):
        """
        Returns bufpos number of grid point ix, iy.
        Args:
            ix (int): x index
            iy (int): y index
        Returns:
            index (float) + 100: index is the bufpos number before 
                offset is added.
        """
        index = iy * self.nx + ix
        return(index + 100)



class SpecBankCal(SpecBank):
    """
    Base class to deal with a complete "bank" of spectra
    """
    def __init__(self, roach_files, ifproc_data, 
                 pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], 
                 bdrop=50, edrop=100):
        """
        Args: 
            roach files (list): list of files for ROACH boards 
                (nominally 4)
            ifproc_data (object): data from corresponding ifproc file
            pixel_list (list): list of pixels to process (default is 
                all)
            bdrop (int): number of channels excluded from beginning of 
                spectrum
            edrop (int): number of channels excluded from end of 
                spectrum
        Returns:
            none
        """
        # @todo no bank?
        SpecBank.__init__(self, roach_files, ifproc_data, 
                          pixel_list=pixel_list)
    
        # go ahead and process the calibration files
        for ipix in range(self.npix):
            self.roach[ipix].compute_tsys_spectrum(bdrop, edrop)

    def test_cal(self, DATA_FILE):
        """
        Checks correspondance of observation parameters against map.
        Args:
            DATA_FILE (object): data file object
        Returns:
            result (float): value indicating number of parameters that 
                do not correspond
        """
        result = 0
        if DATA_FILE.line_rest_frequency != self.line_rest_frequency:
            print('Warning: Line Rest Frequency of Cal observation does not \
                   match map')
            result = result + 1
        if DATA_FILE.nchan != self.nchan:
            print('Warning: number of channels for Cal observation does not \
                   match map')
            result = result + 1
        if DATA_FILE.bandwidth != self.bandwidth:
            print('Warning: spectrometer bandwidth of Cal observation does \
                   not match map')
            result = result + 1
        return(result)
