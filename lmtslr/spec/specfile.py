"""
A class for intermediate Spectrometer Files

classes: SpecFile
author:  GN
date:    Feb 2020

- mods for writing AzEl/RaDec/LatLong style files
"""

import os
import sys
import time
import numpy as np
import netCDF4
from astropy.stats import mad_std

from lmtslr.spec.spec import SpecBankData
from lmtslr.reduction.line_reduction import LineData, NetCDFLineHeader
from lmtslr.utils.reader import count_otf_spectra
from lmtslr.grid.grid import Grid 

class SpecFile():
    def __init__(self, ifproc, specbank, pix_list):
        self.version = "18-apr-2022"         # modify this if anything in the output SpecFile has been changed
        self.ifproc = ifproc
        self.specbank = specbank
        self.pix_list = pix_list
        self.outnc_filename = None
        self.vslice = []
        self.b_order = 0
        self.b_regions, self.l_regions = [], []
        self.eliminate_list = []

    def set_history(self, history):
        self.history = history
        
        
    def set_line_parameters(self, vslice=[], b_order=0,
                            b_regions=[], l_regions=[],
                            eliminate_list=[]):
        self.vslice = vslice
        self.b_order = b_order
        self.b_regions = b_regions
        self.l_regions = l_regions
        self.eliminate_list = eliminate_list
        self.velocity_slice()
        
    def velocity_slice(self):
        # make a dummy spectrum to count the channels after processing steps
        # also reports the vellocity range in the data
        LD = LineData(self.ifproc, self.specbank.bank, self.specbank.nchan,
                      self.specbank.bandwidth, np.zeros(self.specbank.nchan))
        if len(self.vslice) == 2:
            self.L = LD.vslice(self.vslice[0], self.vslice[1])
            self.nchan_to_save = self.L.nchan
        else:
            self.nchan_to_save = self.specbank.nchan
            self.L = LD.vslice(-10000, 10000) # extreme limits to include whole spectrum
            # @todo   a wrong VLSR/RESTFREQ could cause this not to work
        vmin = self.specbank.c2v(self.specbank.nchan-1)
        vmax = self.specbank.c2v(0)
        print("Spectral Band velocity range: %g  %g km/s" % (vmin,vmax))
            
    def _create_nc_dimensions(self):
        # count the total number of spectra that will be processed and written to file
        total_spectra = count_otf_spectra(self.specbank, self.pix_list)
        
        # dimension of number of spectra is from total number count
        nc_dimension_nspec = self.ncout.createDimension('nspec', total_spectra)
    
        # dimension of number of channels in spectrum is from trial reduction step
        nc_dimension_nchan = self.ncout.createDimension('nchan', self.nchan_to_save)
    
        # just doing 20 characters in string
        nc_dimension_nlabel = self.ncout.createDimension('nlabel', 20)

        # history is long... why is netcdf so irrationally complicated to handle strings
        nc_dimension_nhist = self.ncout.createDimension('nhist', 512)

        # number of pixels in this NC
        #nc_dimension_npix = self.ncout.createDimension('npix', len(self.pix_list))
        nc_dimension_npix = self.ncout.createDimension('npix', 16)

        # number of tsyscals
        nc_dimension_ncal = self.ncout.createDimension('ncal', self.specbank.ncal)

    def _create_nc_header(self):
        # a version header
        nc_version = self.ncout.createVariable('Header.Version', 'c', ('nlabel',))
        nc_version[0:len(self.version)] = self.version

        # history how the SpecFile was created
        nc_history = self.ncout.createVariable('Header.History', 'c', ('nhist',))
        nc_history[0:len(self.history)] = self.history
        
        # the Observation Header
        nc_obsnum = self.ncout.createVariable('Header.Obs.ObsNum', 'i4')
        self.ncout.variables['Header.Obs.ObsNum'][0] = self.specbank.obsnum

        # the MapCoord Header
        nc_mapcoord = self.ncout.createVariable('Header.Obs.MapCoord', 'i4')
        self.ncout.variables['Header.Obs.MapCoord'][0] = self.specbank.map_coord

        # copy the source name into netCDF header
        nc_source = self.ncout.createVariable('Header.Obs.SourceName', 'c', ('nlabel',))
        if len(self.specbank.source) > 19:
            nc_source[0:19] = self.specbank.source[0:19]
        else:
            nc_source[0:len(self.specbank.source)] = self.specbank.source[0:len(self.specbank.source)]

        nc_x_position = self.ncout.createVariable('Header.Obs.XPosition', 'f4')
        nc_y_position = self.ncout.createVariable('Header.Obs.YPosition', 'f4')
        # PJT:  1=RA/DEC   2=Lat/Lon 
        if self.specbank.map_coord == 1:
            self.ncout.variables['Header.Obs.XPosition'][0] = self.specbank.ifproc.source_RA/np.pi*180.0
            self.ncout.variables['Header.Obs.YPosition'][0] = self.specbank.ifproc.source_Dec/np.pi*180.0
        else:
            self.ncout.variables['Header.Obs.XPosition'][0] = self.specbank.ifproc.source_L/np.pi*180.0
            self.ncout.variables['Header.Obs.YPosition'][0] = self.specbank.ifproc.source_B/np.pi*180.0

        # PJT new DATE-OBS
        nc_do = self.ncout.createVariable('Header.Obs.DateObs', 'c', ('nlabel',))
        nc_do[0:len(self.specbank.date_obs)] = self.specbank.date_obs[0:len(self.specbank.date_obs)]

        # PJT new RECEIVER
        nc_rc = self.ncout.createVariable('Header.Obs.Receiver', 'c', ('nlabel',))
        nc_rc[0:len(self.specbank.receiver)] = self.specbank.receiver[0:len(self.specbank.receiver)]
         
        # using line header information derived from spec bank

        ncl = NetCDFLineHeader(self.ncout)
        ncl.write_line_header_variables(self.L) # write using the result of trial run

        # PJT write the additional Header.LineData (hacked)
        LD = LineData(self.ifproc, self.specbank.bank, self.specbank.nchan,
                      self.specbank.bandwidth, np.zeros(self.specbank.nchan))
        ncl.write_line_data_header_variables(LD)

    def _create_nc_data(self):
        # set up the grid geometry
        theGrid = Grid(self.specbank.receiver)

        nc_pix = self.ncout.createVariable('Data.Pixel', 'i4', ('nspec',))
        nc_seq = self.ncout.createVariable('Data.Sequence', 'i4', ('nspec',))
        #nc_pix_list = self.ncout.createVariable('Data.PixelList', 'f4', ('npix'))
        nc_x = self.ncout.createVariable('Data.XPos', 'f4', ('nspec',))
        nc_x.units = 'arcsec'
        nc_y = self.ncout.createVariable('Data.YPos', 'f4', ('nspec',))
        nc_y.units = 'arcsec'
        nc_rms = self.ncout.createVariable('Data.RMS', 'f4', ('nspec',))
        nc_rms.units = 'K'
        nc_data = self.ncout.createVariable('Data.Spectra', 'f4', ('nspec','nchan'))
        nc_data.units = 'K'
        nc_tsys = self.ncout.createVariable('Data.Tsys', 'f4', ('ncal','npix'))
        nc_tsys.units = 'K'

        ncal = self.specbank.ncal
        npix = 16
        
        time0 = time.time()

        fast_nc = True
        if fast_nc:
            # find total nspec to allocate arrays to make
            count = 0
            for ipix in self.pix_list:
                i = self.specbank.find_pixel_index(ipix)
                if count == 0:
                    L = LineData(self.ifproc, self.specbank.bank,
                                 self.specbank.nchan, self.specbank.bandwidth,
                                 self.specbank.roach[i].reduced_spectra[0], None)
                    LL = L.vslice(self.vslice[0], self.vslice[1])
                count = count + len(self.specbank.roach[i].xmap[self.specbank.roach[i].ons])
            nspec = count
            nchan = len(LL)
            print("NSPEC: %d" % nspec)
            print("NCHAN: %d" % nchan)
            tmp_pix = np.zeros(nspec, dtype=int)
            tmp_seq = np.zeros(nspec, dtype=int)
            tmp_pixl= np.zeros(npix,  dtype=int)
            tmp_rms = np.zeros(nspec)
            tmp_x   = np.zeros(nspec)
            tmp_y   = np.zeros(nspec)
            tmp_data= np.zeros(nspec*nchan).reshape(nspec,nchan)
            tmp_tsys= np.zeros(ncal*npix).reshape(ncal,npix)
        print("FAST_NC:", fast_nc)
        #nc_pix_list[:] = self.pix_list[:]
        
        count = 0
        
        print("Looping over %d pixel list %s: " % (len(self.pix_list),str(self.pix_list)))
        print("Processing %d CAL's for Tsys" % ncal)
        print("Pix Nspec  Mean Std    MAD_std Min  Max      <RMS> RMS_max    Warnings")

        # MH: calculate for all 3 coords: gx_ra,gy_ra, gx_az,gx_el, gx_glon,gy_glon or make gx,gy arrays with dim=3

        # @todo ensure the pix_list is sorted
        for ipix in self.pix_list:
            count0 = count
            i = self.specbank.find_pixel_index(ipix)
            n_spectra = len(self.specbank.roach[i].xmap[self.specbank.roach[i].ons])
            x_spectra =     self.specbank.roach[i].xmap[self.specbank.roach[i].ons] # x coordinate
            y_spectra =     self.specbank.roach[i].ymap[self.specbank.roach[i].ons] # y coordinate
            if self.ifproc.map_coord == 0:
                gx,gy = theGrid.azel(self.specbank.elev/180. * np.pi,
                                     self.ifproc.tracking_beam)
            elif self.ifproc.map_coord == 1:
                parang = np.mean(self.specbank.roach[i].pmap[self.specbank.roach[i].ons]) # average parang
                gx,gy = theGrid.radec(self.specbank.elev/180. * np.pi, parang,
                                      self.ifproc.tracking_beam)
            elif self.ifproc.map_coord == 2:
                parang = np.mean(self.specbank.roach[i].pmap[self.specbank.roach[i].ons]) # average parang
                galang = np.mean(self.specbank.roach[i].gmap[self.specbank.roach[i].ons]) # average galang
                gx,gy = theGrid.latlon(self.specbank.elev/180. * np.pi, parang, galang,
                                       self.ifproc.tracking_beam)
            else:
                print("Illegal value self.ifproc.map_coord =",self.ifproc.map_coord)
                
            for j in range(n_spectra):
                # process each spectrum
                if j < ncal:
                    # tsyscal is allowed to be None, in which case it's never written
                    # but since they don't change, only on the first spectrum for this pixel it's needed
                    ### S.roach[0].tsys_spectra
                    if ncal > 1:
                        tsys = self.specbank.roach[i].tsys_spectra[j]
                    elif ncal == 1:
                        tsys = self.specbank.roach[i].tsyscal
                    else:
                        tsys = None
                else:
                    tsys = None
                
                L = LineData(self.ifproc, self.specbank.bank,
                             self.specbank.nchan, self.specbank.bandwidth,
                             self.specbank.roach[i].reduced_spectra[j],
                             tsys)
                LL = L.vslice(self.vslice[0], self.vslice[1])
                LL.eliminate(self.eliminate_list)
                bbase, nbase = LL.xlist(self.b_regions)
                LL.baseline(bbase, nbase, baseline_order=self.b_order)

                # write the reduced line into the NetCDF file
                if fast_nc:
                    tmp_data[count,:] = LL.yarray
                else:
                    nc_data[count,:] = LL.yarray

                # tricked: only the first ncal tsys are for real 
                if type(LL.tarray) == np.ndarray:
                    idx = self.pix_list.index(ipix)
                    if fast_nc:
                        tmp_tsys[j,ipix] = LL.tarray.mean()
                    else:
                        nc_tsys[j,ipix] = LL.tarray.mean()
                    #nc_tsys[j,ipix,:] = LL.tarray
                    t = LL.tarray
                    print("TSYS[%d] slice: %g (%g)  minmax: %g %g" % (ipix,t.mean(),t.std(),t.min(),t.max()))


                # MH: calculate for all 3 coords: nc_x_az,nc_y_el,nc_x_ra,nc_y_dec,nc_x_glon, nc_y_glon or make nc_x,nc_y as arrays with dim=3
                # PJT    - need 3 arrays for AE/RD/LB
                if not fast_nc:
                    nc_rms[count]  = LL.rms
                    nc_pix[count]  = ipix
                    nc_seq[count]  = j
                    nc_x[count]    = x_spectra[j]-gx[ipix]
                    nc_y[count]    = y_spectra[j]-gy[ipix]
                else:
                    tmp_rms[count] = LL.rms
                    tmp_pix[count] = ipix
                    tmp_seq[count] = j
                    tmp_x[count]   = x_spectra[j]-gx[ipix]
                    tmp_y[count]   = y_spectra[j]-gy[ipix]
                    
                count = count + 1                
            if fast_nc:
                prms  = tmp_rms[count0:count]
                pdata = tmp_data[count0:count,:]                
            else:
                prms  = nc_rms[count0:count]
                pdata = nc_data[count0:count,:]
            # Attempt to flag data that mis-behave
            s1 = pdata.mean()
            s2 = pdata.std()
            s3 = mad_std(pdata)
            s4 = pdata.min()            
            s5 = pdata.max()
            s6 = prms.mean()
            s7 = prms.max()
            msg = ""
            if s2/s3 > 1.2:  msg = msg + " *P %.1f" % (s2/s3)
            s8 = prms.std()
            s9 = mad_std(prms)
            if s8/s9 > 1.2:  msg = msg + " *M %.1f" % (s8/s9)

            print("%d %d   %.3f %.3f %.3f %.3f %.3f   %.3f %.3f     %s" %
                  (ipix,count-count0,s1,s2,s3,s4,s5,s6,s7,msg))

        if fast_nc:
            # another braindead netcdf feature? can't use without []
            print("CPU TIME: %g sec" % (time.time()-time0))
            nc_rms[:]      = tmp_rms[:]
            nc_pix[:]      = tmp_pix[:]
            nc_seq[:]      = tmp_seq[:]
            nc_x[:]        = tmp_x[:]
            nc_y[:]        = tmp_y[:]
            nc_data[:,:]   = tmp_data[:,:]
            nc_tsys[:,:]   = tmp_tsys[:,:]
        print("CPU TIME: %g sec" % (time.time()-time0))
        
        print("Warnings:")
        print("  *P ratio:      ratio of std/mad too high for data")
        print("  *M ratio:      ratio of std/mad too high for RMS")
        if fast_nc:
            # @todo can't do this using the nc_; why netcdf, why?
            xmin = tmp_x.min()
            xmax = tmp_x.max()
            ymin = tmp_y.min()
            ymax = tmp_y.max()
        else:
            xd = np.zeros(count)
            yd = np.zeros(count)
            for i in range(count):
                xd[i] = nc_x[i]
                yd[i] = nc_y[i]
            xmin = xd.min()
            xmax = xd.max()
            ymin = yd.min()
            ymax = yd.max()
        print("X-range: %g %g   Y-range: %g %g arcsec\n" %(xmin, xmax, ymin, ymax))
                  
            
            
    def open_output_netcdf(self, output_file_name):
        if os.path.exists(output_file_name):
            print("Error: Filename %s already exists. Please rename or delete it." % output_file_name)
            raise
        self.outnc_filename = os.path.abspath(output_file_name)
        self.ncout = netCDF4.Dataset(output_file_name, 'w', format='NETCDF4')
        
    def write_ncdata(self):
        if not hasattr(self, 'ncout'):
            print("First open the output netcdf file")
            raise
        self._create_nc_dimensions()
        self._create_nc_header()
        self._create_nc_data()
        self.ncout.close()
        
class OTFSpecFile(SpecFile):
    def __init__(self, ifproc, specbank, pix_list):
        SpecFile.__init__(self, ifproc, specbank, pix_list)
        
