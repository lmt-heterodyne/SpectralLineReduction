"""
A class for intermediate Spectrometer Files

classes: SpecFile
author: GN
date: Feb 2020
"""

import netCDF4
import os
from lmtslr.spec.spec import SpecBankData
from lmtslr.reduction.line_reduction import LineData, NetCDFLineHeader
from lmtslr.utils.reader import count_otf_spectra
from lmtslr.grid.grid import Grid 
import numpy as np

class SpecFile():
    def __init__(self, ifproc, specbank, pix_list):
        self.ifproc = ifproc
        self.specbank = specbank
        self.pix_list = pix_list
        self.outnc_filename = None
        self.vslice = []
        self.b_order = 0
        self.b_regions, self.l_regions = [], []
        self.eliminate_list = []
        
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
        LD = LineData(self.ifproc, self.specbank.bank, self.specbank.nchan,
                      self.specbank.bandwidth, np.zeros(self.specbank.nchan))
        if len(self.vslice) == 2:
            self.L = LD.vslice(self.vslice[0], self.vslice[1])
            self.nchan_to_save = self.L.nchan
        else:
            self.nchan_to_save = self.specbank.nchan
            self.L = LD.vslice(-10000, 10000) # extreme limits to include whole spectrum
            
    def _create_nc_dimensions(self):
        # count the total number of spectra that will be processed and written to file
        total_spectra = count_otf_spectra(self.specbank, self.pix_list)
        
        # dimension of number of spectra is from total number count
        nc_dimension_nspec = self.ncout.createDimension('nspec', total_spectra)
    
        # dimension of number of channels in spectrum is from trial reduction step
        nc_dimension_nchan = self.ncout.createDimension('nchan', self.nchan_to_save)
    
        # just doing 20 characters in string
        nc_dimension_nlabel = self.ncout.createDimension('nlabel', 20)

    def _create_nc_header(self):
        # the Observation Header
        nc_obsnum = self.ncout.createVariable('Header.Obs.ObsNum', 'i4')
        self.ncout.variables['Header.Obs.ObsNum'][0] = self.specbank.obsnum

        # copy the source name into netCDF header
        nc_source = self.ncout.createVariable('Header.Obs.SourceName', 'c', ('nlabel',))
        if len(self.specbank.source) > 19:
            nc_source[0:19] = self.specbank.source[0:19]
        else:
            nc_source[0:len(self.specbank.source)] = self.specbank.source[0:len(self.specbank.source)]

        nc_x_position = self.ncout.createVariable('Header.Obs.XPosition', 'f4')
        nc_y_position = self.ncout.createVariable('Header.Obs.YPosition', 'f4')
        if self.specbank.map_coord == 1:
            self.ncout.variables['Header.Obs.XPosition'][0] = \
                                                              self.specbank.ifproc.source_RA/np.pi*180.0
            self.ncout.variables['Header.Obs.YPosition'][0] = \
                                                              self.specbank.ifproc.source_Dec/np.pi*180.0
        else:
            self.ncout.variables['Header.Obs.XPosition'][0] = 0.0
            self.ncout.variables['Header.Obs.YPosition'][0] = 0.0
        # using line header information derived from spec bank
        ncl = NetCDFLineHeader(self.ncout)
        ncl.write_line_header_variables(self.L) # write using the result of trial run             

    def _create_nc_data(self):
        # set up the grid geometry
        theGrid = Grid()

        nc_pix = self.ncout.createVariable('Data.Pixel', 'i4', ('nspec',))
        nc_seq = self.ncout.createVariable('Data.Sequence', 'i4', ('nspec',))
        nc_x = self.ncout.createVariable('Data.XPos', 'f4', ('nspec',))
        nc_x.units = 'arcsec'
        nc_y = self.ncout.createVariable('Data.YPos', 'f4', ('nspec',))
        nc_y.units = 'arcsec'
        nc_rms = self.ncout.createVariable('Data.RMS', 'f4', ('nspec',))
        nc_rms.units = 'K'
        nc_data = self.ncout.createVariable('Data.Spectra', 'f4', ('nspec','nchan'))
        nc_data.units = 'K'

        count = 0
        for ipix in self.pix_list:
            i = self.specbank.find_pixel_index(ipix)
            n_spectra = len(self.specbank.roach[i].xmap[self.specbank.roach[i].ons])
            x_spectra = self.specbank.roach[i].xmap[self.specbank.roach[i].ons] # x coordinate
            y_spectra = self.specbank.roach[i].ymap[self.specbank.roach[i].ons] # y coordinate
            if self.ifproc.map_coord == 0:
                gx,gy = theGrid.azel(self.specbank.elev/180. * np.pi,
                                     self.ifproc.tracking_beam)
            else:
                parang = np.mean(self.specbank.roach[i].pmap[self.specbank.roach[i].ons]) # average parang
                gx,gy = theGrid.radec(self.specbank.elev/180. * np.pi, parang,
                                      self.ifproc.tracking_beam)

            for j in range(n_spectra):
                # process each spectrum
                L = LineData(self.ifproc, self.specbank.bank,
                             self.specbank.nchan, self.specbank.bandwidth,
                             self.specbank.roach[i].reduced_spectra[j])
                LL = L.vslice(self.vslice[0], self.vslice[1])
                LL.eliminate(self.eliminate_list)
                bbase, nbase = LL.xlist(self.b_regions)
                LL.baseline(bbase, nbase, baseline_order=self.b_order)

                # write the reduced line into the NetCDF file
                nc_data[count,:] = LL.yarray
                nc_rms[count] = LL.rms
                nc_pix[count] = ipix
                nc_seq[count] = j
                nc_x[count] = x_spectra[j]-gx[ipix]
                nc_y[count] = y_spectra[j]-gy[ipix]
                count = count + 1

            
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
        
