#!/usr/bin/env python
'''
Creates a SpecFile from a single OTF mapping observation
'''

# Python Imports
import numpy as np
import matplotlib.pyplot as pl
import subprocess
import netCDF4

# Line Data Reduction Imports
from lmtslr.spec.spec import *
from lmtslr.ifproc.ifproc import *

from lmtslr.reduction.line_reduction import *
from lmtslr.grid.grid import *

# set up the grid geometry
theGrid = Grid()

try:
    from  pipeline_controller import *
except:
    print("Warning: no pipeline_controller")
    


def main(argv):
    Opts = HandlePipelineOptions()
    Opts.parse_options(argv,'process_map',1)
    
    # check to see whether files exist, and remove if they do
    if os.path.isfile(Opts.spec_file_name) == True: 
        os.remove(Opts.spec_file_name) 

    # read the data from the relevant files for this obsnum

    # note use of "data_path" to set the path to the data.  Default path is '/data_lmt/'
    I,S = read_obsnum_otf(Opts.obsnum,
                          Opts.list_of_pixels,
                          Opts.bank,
                          Opts.use_cal,
                          tsys=Opts.tsys,
                          path=Opts.data_path)

    # count the total number of spectra that will be processed and written to file
    total_spectra = count_otf_spectra(S,Opts.list_of_pixels)

    # make a dummy spectrum to count the channels after processing steps
    LD = LineData(I,Opts.bank,S.nchan,S.bandwidth,np.zeros(S.nchan))
    L = LD.vslice(Opts.slice_list[0],Opts.slice_list[1])
    nchan_to_save = L.nchan

    # write the netCDF file

    # open Dataset.  If the file exists, we stop here!
    nc = netCDF4.Dataset(Opts.spec_file_name, 'w', format='NETCDF4')

    # dimension of number of spectra is from total number count
    nc_dimension_nspec = nc.createDimension('nspec',total_spectra)

    # dimension of number of channels in spectrum is from trial reduction step
    nc_dimension_nchan = nc.createDimension('nchan',nchan_to_save)

    # just doing 20 characters in string
    nc_dimension_nlabel = nc.createDimension('nlabel',20)

    # the Observation Header
    nc_obsnum = nc.createVariable('Header.Obs.ObsNum','i4')
    nc.variables['Header.Obs.ObsNum'][0] = S.obsnum

    # copy the source name into netCDF header
    nc_source = nc.createVariable('Header.Obs.SourceName','c',('nlabel',))
    if len(S.source) > 19:
        nc_source[0:19] = S.source[0:19]
    else:
        nc_source[0:len(S.source)] = S.source[0:len(S.source)]

    nc_x_position = nc.createVariable('Header.Obs.XPosition','f4')
    nc_y_position = nc.createVariable('Header.Obs.YPosition','f4')
    if S.map_coord == 1:
        nc.variables['Header.Obs.XPosition'][0] = S.ifproc.source_RA/np.pi*180.0
        nc.variables['Header.Obs.YPosition'][0] = S.ifproc.source_Dec/np.pi*180.0
    else:
        nc.variables['Header.Obs.XPosition'][0] = 0.0
        nc.variables['Header.Obs.YPosition'][0] = 0.0

    # using line header information derived from spec bank
    ncl = NetCDFLineHeader(nc)
    ncl.write_line_header_variables(L) # write using the result of trial run 
          
    nc_pix = nc.createVariable('Data.Pixel','i4',('nspec',))
    nc_seq = nc.createVariable('Data.Sequence','i4',('nspec',))
    nc_x = nc.createVariable('Data.XPos','f4',('nspec',))
    nc_x.units = 'arcsec'
    nc_y = nc.createVariable('Data.YPos','f4',('nspec',))
    nc_y.units = 'arcsec'
    nc_rms = nc.createVariable('Data.RMS','f4',('nspec',))
    nc_rms.units = 'K'
    nc_data = nc.createVariable('Data.Spectra','f4',('nspec','nchan'))
    nc_data.units = 'K'

    count = 0
    for ipix in Opts.list_of_pixels:
        i = S.find_pixel_index(ipix)
        n_spectra = len(S.roach[i].xmap[S.roach[i].ons])
        x_spectra = S.roach[i].xmap[S.roach[i].ons] # x coordinate
        y_spectra = S.roach[i].ymap[S.roach[i].ons] # y coordinate
        if I.map_coord == 0:
            gx,gy = theGrid.azel(S.elev/180.*np.pi,I.tracking_beam)
        else:
            parang = np.mean(S.roach[i].pmap[S.roach[i].ons]) # average parang
            gx,gy = theGrid.radec(S.elev/180.*np.pi,parang,I.tracking_beam)

        for j in range(n_spectra):
            # process each spectrum
            L = LineData(I,Opts.bank,S.nchan,S.bandwidth,S.roach[i].reduced_spectra[j])
            LL = L.vslice(Opts.slice_list[0],Opts.slice_list[1])
            bbase,nbase = LL.xlist(Opts.baseline_regions)
            LL.baseline(bbase,nbase,baseline_order=Opts.baseline_order)
        
            # write the reduced line into the NetCDF file
            nc_data[count,:] = LL.yarray
            nc_rms[count] = LL.rms
            nc_pix[count] = ipix
            nc_seq[count] = j
            nc_x[count] = x_spectra[j]-gx[ipix]
            nc_y[count] = y_spectra[j]-gy[ipix]
            count = count + 1


    nc.close()        
    print('netCDF Done')

if __name__ == '__main__':
    main(sys.argv[1:])

