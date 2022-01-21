#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits
import astropy.wcs as wcs

#from lmtslr.utils.parser import HandleViewCubeOptions
from lmtslr.utils.argparser import HandleViewCubeOptions

def channel_range(h,vlo,vhi):
    ''' 
    vlo and vhi are velocities of the end of a range 
    '''
    plo = np.round((vlo-h['CRVAL3'])/h['CDELT3'] + h['CRPIX3'])
    phi = np.round((vhi-h['CRVAL3'])/h['CDELT3'] + h['CRPIX3'])
    if h['CDELT3']<0.0:
        temp = plo
        plo = phi
        phi = temp
    return int(plo),int(phi)

def pixel_location(h,dx,dy):
    ''' 
    dx,dy are offsets from reference pixel
    opx, opy are pixels for dx,dy offsets
    px, py are locations for dx,dy offsets
    '''
    opx = np.round(dx/h['CDELT1'] + h['CRPIX1'])
    opy = np.round(dy/h['CDELT2'] + h['CRPIX2'])
    
    px = (opx-h['CRPIX1'])*h['CDELT1'] + h['CRVAL1']
    py = (opy-h['CRPIX2'])*h['CDELT2'] + h['CRVAL2']
    
    return int(opx),int(opy),px,py

def pixel_range(h,d):
    ''' 
    d is list of pixel coordinates [xlo,xhi,ylo,yhi]
    xlo,xhi,ylo,yhi are offsets from reference pixel
    '''
    pxlo = np.round(d[0]/h['CDELT1'] + h['CRPIX1'])
    pxhi = np.round(d[1]/h['CDELT1'] + h['CRPIX1'])
    if h['CDELT1']<0.:
        temp = pxlo
        pxlo = pxhi
        pxhi = temp
    pylo = np.round(d[2]/h['CDELT2'] + h['CRPIX2'])
    pyhi = np.round(d[3]/h['CDELT2'] + h['CRPIX2'])
    if h['CDELT2']<0.:
        temp = pylo
        pylo = pyhi
        pyhi = temp
    return [int(pxlo),int(pxhi),int(pylo),int(pyhi)]

def create_caxis_arrays(h):
    ''' creates axes for FITS data cube given FITS Header information
        INPUT:
            h - FITS header
        RETURNS:
            CAXIS1 - array with values along FITS axis 1 - fast axis
            CAXIS2 - array with values along FITS axis 2 - intermediate axis
            CAXIS3 - array with values along FITS axis 3 - slow axis
    '''
    CAXIS1 = (np.arange(h['NAXIS1'])-h['CRPIX1'])*h['CDELT1'] + h['CRVAL1']
    CAXIS2 = (np.arange(h['NAXIS2'])-h['CRPIX2'])*h['CDELT2'] + h['CRVAL2']
    CAXIS3 = (np.arange(h['NAXIS3'])-h['CRPIX3'])*h['CDELT3'] + h['CRVAL3']
    return CAXIS1,CAXIS2,CAXIS3



def main(argv):

    Opts = HandleViewCubeOptions()
    Opts.parse_options(argv,'view_cube', 1, True)
    

    pl.ion()
    
    # open the file
    hdulist = fits.open(Opts.input_file_name)
    hdu = hdulist[0]
    
    # T is the data cube
    T = hdu.data

    # set channels
    lo_channel, hi_channel = channel_range(hdu.header,
                                           Opts.v_range[0]*Opts.v_scale,
                                           Opts.v_range[1]*Opts.v_scale)
    print(Opts.v_range[0],lo_channel,Opts.v_range[1],hi_channel)

    # set pixel for display
    x_pixel, y_pixel, x_loc, y_loc = pixel_location(hdu.header,
                                                    Opts.location[0]*Opts.scale,
                                                    Opts.location[1]*Opts.scale)
    print(Opts.location[0], x_pixel, Opts.location[1], y_pixel)

    # map limits
    map_limits = [Opts.limits[0]*Opts.scale,
                  Opts.limits[1]*Opts.scale,
                  Opts.limits[2]*Opts.scale,
                  Opts.limits[3]*Opts.scale] 
    
    pixel_limits = pixel_range(hdu.header, map_limits)
    print(Opts.limits)
    print(pixel_limits)
    
    # make the integrated intensity for plotting
    TINT = np.sum(T,axis=0)*np.abs(hdu.header['CDELT3']/Opts.v_scale)
    
    # find the peak intensity for plotting
    TMAX = np.max(T,axis=0)

    # create axes for the heck of it
    CAXIS1,CAXIS2,CAXIS3 = create_caxis_arrays(hdu.header)

    # print out the FITS header
    print(repr(hdu.header))

    TINT = np.sum(T[lo_channel:hi_channel,:,:],axis=0)
    #int_plot_levels = [1.,2.,3.,4.,5.]
    
    # display slices of the data cube
    fig,ax = pl.subplots(figsize=[12,8])
    ax.imshow(T[:,:,x_pixel],origin='lower',extent=[CAXIS1[0],CAXIS1[-1],CAXIS3[0]/Opts.v_scale,CAXIS3[-1]/Opts.v_scale],clim=Opts.tmax_range,aspect='auto',cmap=pl.cm.jet)
    pl.xlabel(hdu.header['CTYPE1'])
    pl.ylabel(hdu.header['CTYPE3'])
    pl.title('%s DEC=%f'%(hdu.header['OBJECT'],y_loc))
    
    fig,ax = pl.subplots(figsize=[12,8])
    ax.imshow(T[:,y_pixel,:],
              origin='lower',
              extent=[CAXIS2[0],CAXIS2[-1],CAXIS3[0]/Opts.v_scale,CAXIS3[-1]/Opts.v_scale],
              clim=Opts.tmax_range,
              aspect='auto',
              cmap=pl.cm.jet)
    pl.xlabel(hdu.header['CTYPE2'])
    pl.ylabel(hdu.header['CTYPE3'])
    pl.title('%s RA=%f'%(hdu.header['OBJECT'],x_loc))
    
    fig,ax = pl.subplots(figsize=[12,8])
    if Opts.plot_type == 'TINT':
        ax.imshow(TINT,
                  origin='lower',
                  extent=[CAXIS1[0],CAXIS1[-1],CAXIS2[0],CAXIS2[-1]],
                  clim=Opts.tint_range,
                  aspect='auto',
                  cmap=pl.cm.jet)
        print(x_loc,y_loc,TINT[y_pixel,x_pixel])
        ax.plot(x_loc,y_loc,'w+')
        pl.title('%s TINT [%f, %f]'%(hdu.header['OBJECT'],Opts.v_range[0],Opts.v_range[1]))

    else:
        ax.imshow(TMAX,
                  origin='lower',
                  extent=[CAXIS1[0],CAXIS1[-1],CAXIS2[0],CAXIS2[-1]],
                  clim=Opts.tmax_range,
                  aspect='auto',
                  cmap=pl.cm.jet)
        print(x_loc,y_loc,TMAX[y_pixel,x_pixel])
        ax.plot(x_loc,y_loc,'w+')
        pl.title('%s TMAX'%(hdu.header['OBJECT']))
        
    pl.xlabel(hdu.header['CTYPE1'])
    pl.ylabel(hdu.header['CTYPE2'])
    pl.axis('equal')

    pl.figure(figsize=(12,8))
    pl.plot(CAXIS3/Opts.v_scale,T[:,y_pixel,x_pixel],'k')
    pl.xlabel(hdu.header['CTYPE3'])
    pl.ylabel('TA*')
    pl.title('%s (%5.1f, %5.1f)'%(hdu.header['OBJECT'],Opts.location[0],Opts.location[1]))

    w = wcs.WCS(hdu.header)

    # print the WCS geometry information
    #print('WCS INFORMATION - All Axes')
    #print(w)
    
    # we must drop the "v" axis from consideration to make maps
    # in this example (and according to most conventions) this is FITS axis 3 (index 2)
    ww = w.dropaxis(2)
    
    # print the keywords for map display
    #print('\n\n')
    #print('WCS INFORMATION - Map Axes')
    #print(ww)

    pl.figure(figsize=[12,8])

    # map the individual channel identified as "my_channel"
    ax = pl.subplot(projection=ww)

    # Set up "Longitude" Axis
    ax.set_xlim(pixel_limits[0],pixel_limits[1]) # note that axis limits are given in pixels
    lon = ax.coords[0]
    lon.set_major_formatter('hh:mm:ss')
    lon.set_axislabel('RA (J2000)')

    # Set up "Latitude" Axis
    ax.set_ylim(pixel_limits[2],pixel_limits[3]) # note that axis limits are given in pixels
    lat = ax.coords[1]
    lat.set_axislabel('Dec (J2000)')

    # plot the image and the contours
    if Opts.plot_type == 'TINT':
        ax.imshow(TINT,
                  interpolation=Opts.interp,
                  origin='lower',
                  clim=Opts.tint_range,
                  cmap=pl.cm.jet)
        pl.title('%s TINT [%f, %f]'%(hdu.header['OBJECT'],Opts.v_range[0],Opts.v_range[1]))

    else:
        ax.imshow(TMAX,
                  interpolation=Opts.interp,
                  origin='lower',
                  clim=Opts.tmax_range,
                  cmap=pl.cm.jet)
        pl.title('%s TMAX'%(hdu.header['OBJECT']))
        
        #ax.contour(TINT,levels=int_plot_levels,colors='white',alpha=0.5)

        # here is the command to invert the x axis if it was not stored correctly in the first place.
        #ax.invert_xaxis()
    pl.ioff()
    pl.show()

if __name__ == '__main__':    
    main(sys.argv[1:])
