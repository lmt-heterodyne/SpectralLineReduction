import numpy as np
import matplotlib.pyplot as pl
import netCDF4
from astropy.stats import mad_std
from astropy.io import fits
from lmtslr.viewer.plots import Plots

class SpecFileViewer():
    """
    Class for viewing spectrum files and conversion to a waterfall 
    """
    def __init__(self, netcdf_filename):
        """
        Constructor for SpecFileViewer class.
        Args:
            netcdf_filename
        Returns:
            none
        """
        nc = netCDF4.Dataset(netcdf_filename, 'r', format='NETCDF4')
        self.obsnum = nc.variables['Header.Obs.ObsNum'][0]
        self.source_name = netCDF4.chartostring(nc.variables['Header.Obs.SourceName'][:])

        self.x_position =  nc.variables['Header.Obs.XPosition'][0]
        self.y_position =  nc.variables['Header.Obs.YPosition'][0]

        self.nchan = nc.variables['Header.Line.NChannels'][0]
        self.chan = nc.variables['Header.Line.ChannelNumber'][:]
        self.cdelt = nc.variables['Header.SpectrumAxis.CDELT'][0]
        self.crpix = nc.variables['Header.SpectrumAxis.CRPIX'][0]
        self.crval = nc.variables['Header.SpectrumAxis.CRVAL'][0]
        self.ctype = netCDF4.chartostring(nc.variables['Header.SpectrumAxis.CTYPE'][:])
        self.caxis = nc.variables['Header.SpectrumAxis.CAXIS'][:]

        self.pixel = nc.variables['Data.Pixel'][:]
        self.sequence = nc.variables['Data.Sequence'][:]
        self.xpos = nc.variables['Data.XPos'][:]
        self.ypos = nc.variables['Data.YPos'][:]
        self.rms = nc.variables['Data.RMS'][:]
        self.data = nc.variables['Data.Spectra'][:]
        try:
            # code version 6-mar-2021 and beyond
            self.have_tsys = True
            self.tsys = nc.variables['Data.Tsys'][:]
        except:
            self.have_tsys = False

        nc.close()

    def sequoia_waterfall_plot(self, pixel_list, rms_cut, plot_range=[-1,1],figsize=8):
        """
        Makes waterfall plots of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [-1,1])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        fig1, ax1 = pl.subplots(4, 4, sharex='col', sharey='row', 
                                gridspec_kw={'hspace': 0, 'wspace': 0}, 
                                figsize=(figsize, figsize))
        fig1.text(0.02, 0.5, self.ctype, va='center', rotation='vertical')
        fig1.text(0.5, 0.1, 'Sample', ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            if rms_cut < 0:
                rindex = np.where(self.rms[pindex] >= 0)[0]
            else:
                rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax1[np.mod(the_pixel, 4), the_pixel // 4].imshow(
                self.data[pindex[rindex]].transpose(), origin='lower', 
                extent=[0, float(len(rindex)), self.caxis[0], self.caxis[-1]],
                clim=plot_range, aspect='auto')
            ax1[np.mod(the_pixel, 4), the_pixel // 4].text(0.05 * len(rindex),
                self.caxis[0] + 0.85 * (self.caxis[-1] - self.caxis[0]), 
                '%d'%(the_pixel))
        Plots.savefig()

    def pixel_waterfall_plot(self, the_pixel, rms_cut, plot_range=[-1,1]):
        """
        Makes waterfall plot of spectrum from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [-1,1])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        if rms_cut < 0:
            rindex = np.where(self.rms[pindex] >= 0)[0]
        else:
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.imshow(self.data[pindex[rindex]].transpose(), origin='lower', 
                  extent=[0, float(len(rindex)), self.caxis[0], 
                  self.caxis[-1]], clim=plot_range, aspect='auto')
        pl.title('PIXEL: %d'%(the_pixel))
        pl.ylabel(self.ctype)
        pl.xlabel('Sample')
        pl.colorbar()
        Plots.savefig()

    def sequoia_rms_plot(self, pixel_list, rms_cut, plot_range=[0,10], figsize=8):
        """
        Makes rms plots of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [0,10])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        fig2, ax2 = pl.subplots(4, 4, sharex='col', sharey='row', 
                                gridspec_kw={'hspace': 0, 'wspace': 0}, 
                                figsize=(figsize,figsize))
        fig2.text(0.02, 0.5, 'RMS', va='center', rotation='vertical')
        fig2.text(0.5, -0.1, 'Sample', ha='center')

        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            if rms_cut < 0:
                rindex = np.where(self.rms[pindex] >= 0)[0]
            else:
                rindex = np.where(self.rms[pindex] < rms_cut)[0]
            
            ax2[np.mod(the_pixel,4), the_pixel // 4].plot(self.rms[pindex[rindex]], 'k.')
            #ax2[np.mod(the_pixel,4), the_pixel//4].text(0.05*len(rindex),plot_range[0] + 0.9*(plot_range[-1]-plot_range[0]), '%d'%(the_pixel))
            #ax2[np.mod(the_pixel,4), the_pixel//4].ylim(plot_range)
        Plots.savefig()

    def pixel_rms_plot(self, the_pixel, rms_cut, plot_range=[0,10]):
        """
        Makes rms plot of spectrum from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            plot_range (list): range of plot (default is [0,10])
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        if rms_cut < 0:
            rindex = np.where(self.rms[pindex] >= 0)[0]
        else:
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.plot(self.rms[pindex[rindex]], 'k.')
        pl.ylim(plot_range)
        pl.ylabel('RMS')
        pl.xlabel('Sample')
        pl.title('PIXEL: %d'%(the_pixel))
        Plots.savefig()
        
    def xy_position_plot(self, all=True, title='Sky Coverage'):
        """
        Makes x-y position plot.   xpos > 0 means larger RA, to the left
                                   ypos > 0 means larger DEC, to the top
        Args:
            none
        Returns:
            none
        """
        Plots.figure()
        s0 = 0
        s1 = max(self.sequence)+1
        print("Max sequence=%d" % s1)
        if all:
            pl.plot(-self.xpos,self.ypos, 'k.', markersize=0.2)
        else:
            pl.plot(-self.xpos[s0:s1],self.ypos[s0:s1], 'k.', markersize=0.2)
        xlim = pl.xlim()
        ylim = pl.ylim()
        pmin = min(xlim[0],ylim[0])
        pmax = max(xlim[1],ylim[1])
        pmax = max(abs(pmax),abs(pmin))
        pl.xlim([-pmax,pmax])
        pl.ylim([-pmax,pmax])
        axes=pl.gca()
        axes.set_aspect("equal")
        
        pl.title(title)
        pl.xlabel('X offset [arcsec]')
        pl.ylabel('Y offset [arcsec]')
        Plots.savefig()        

    def sx_position_plot(self, all=True):
        """
        Makes sequence-x plot.
        Args:
            none
        Returns:
            none
        """
        Plots.figure()
        s0 = 0
        s1 = max(self.sequence)+1
        if all:
            pl.plot(self.sequence,self.xpos, 'k.')
        else:
            pl.plot(self.sequence[s0:s1],self.xpos[s0:s1], 'k.')
        pl.xlabel('Sequence')
        pl.ylabel('X')
        Plots.savefig()        

    def sy_position_plot(self, all=True):
        """
        Makes sequence-y position plot.
        Args:
            none
        Returns:
            none
        """
        Plots.figure()
        s0 = 0
        s1 = max(self.sequence)+1
        if all:
            pl.plot(self.sequence,self.ypos, 'k.')
        else:
            pl.plot(self.sequence[s0:s1],self.ypos[s0:s1], 'k.')
        pl.xlabel('Sequence')
        pl.ylabel('Y')
        Plots.savefig()        
        
    def sequoia_rms_histogram(self, pixel_list, rms_cut, figsize=8):
        """
        Makes rms histogram of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        fig3, ax3 = pl.subplots(4, 4, sharex='col', sharey='row', 
            gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(figsize,figsize))
        fig3.text(0.5, -0.1, 'RMS', ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            if rms_cut < 0:
                rindex = np.where(self.rms[pindex] >= 0)[0]
            else:
                rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax3[np.mod(the_pixel, 4), the_pixel // 4].hist(
                self.rms[pindex[rindex]], bins=np.arange(0,3.02,.02))
        Plots.savefig()

    def pixel_rms_histogram(self, the_pixel, rms_cut):
        """
        Makes rms histogram of spectra from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        if rms_cut < 0:
            rindex = np.where(self.rms[pindex] >= 0)[0]
        else:
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.hist(self.rms[pindex[rindex]],bins = np.arange(0,3.02,.02))
        pl.xlabel('RMS')
        pl.ylabel('N')
        pl.title('PIXEL: %d'%(the_pixel))
        Plots.savefig()        
        
    def sequoia_tsys_spectra_plot(self, pixel_list, plot_range=[50,250], figsize=8):
        """
        Makes mean spectra plot of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        if not self.have_tsys: return
        print("NEW TSYS PIC")
        Plots.figure()
        fig4, ax4 = pl.subplots(4, 4, sharex='col', sharey='row', 
                                gridspec_kw={'hspace': 0, 'wspace': 0},
                                figsize=(figsize,figsize))
        fig4.text(0.02, 0.5, 'TSYS [K]', va='center', rotation='vertical')
        fig4.text(0.5, -0.1, self.ctype, ha='center')
        (ncal,npix,nchan) = self.tsys.shape
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            if len(pindex) == 0:
                continue
            for ical in range(ncal):
                ax4[np.mod(the_pixel, 4), the_pixel // 4].plot(self.caxis,self.tsys[ical][the_pixel])
                #ax4[np.mod(the_pixel, 4), the_pixel // 4]   #  .set_xlim([-300,300])
                #ax4[np.mod(the_pixel, 4), the_pixel // 4]   #  .set_ylim(plot_range)
        Plots.savefig()

    def pixel_tsys_spectra_plot(self, the_pixel, plot_range=[50,250]):
        """
        Makes mean spectra plot of spectra from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        if not self.have_tsys: return
        print("NEW TSYS PIC")        
        Plots.figure()
        #pindex = np.where(self.pixel == the_pixel)[0]
        #rindex = np.where(self.rms[pindex] < rms_cut)[0]
        (ncal,npix,nchan) = self.tsys.shape
        for ical in range(ncal):
            pl.plot(self.caxis, self.tsys[ical][the_pixel])
        pl.xlabel(self.ctype)
        pl.ylabel('Tsys')
        pl.title('PIXEL: %d'%(the_pixel))
        Plots.savefig()                

    def sequoia_mean_spectra_plot(self, pixel_list, rms_cut, figsize=8):
        """
        Makes mean spectra plot of spectra from pixels in pixel_list.
        Args:
            pixel_list (list): list of pixel IDs to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        fig4, ax4 = pl.subplots(4, 4, sharex='col', sharey='row', 
            gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(figsize,figsize))
        fig4.text(0.02, 0.5, 'TA*', va='center', rotation='vertical')
        fig4.text(0.5, -0.1, self.ctype, ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel == the_pixel)[0]
            if rms_cut < 0:
                rindex = np.where(self.rms[pindex] >= 0)[0]
            else:
                rindex = np.where(self.rms[pindex] < rms_cut)[0]
            ax4[np.mod(the_pixel, 4), the_pixel // 4].plot(self.caxis, 
                np.mean(self.data[pindex[rindex]], axis=0))
            ax4[np.mod(the_pixel, 4), the_pixel // 4].text(0.05 * len(rindex),
                self.caxis[0] + 0.85 * (self.caxis[-1] - self.caxis[0]), 
                '%d'%(the_pixel))
        Plots.savefig()

    def pixel_mean_spectrum_plot(self, the_pixel, rms_cut):
        """
        Makes mean spectra plot of spectra from pixel the_pixel.
        Args:
            the_pixel (int): pixel ID to plot
            rms_cut (float): rms cutoff value for plot
            figsize (float): size of figure in inches (default is 8)
        Returns:
            none
        """
        Plots.figure()
        pindex = np.where(self.pixel == the_pixel)[0]
        if rms_cut < 0:
            rindex = np.where(self.rms[pindex] >= 0)[0]
        else:
            rindex = np.where(self.rms[pindex] < rms_cut)[0]
        pl.plot(self.caxis, np.mean(self.data[pindex[rindex]], axis=0))
        pl.xlabel(self.ctype)
        pl.ylabel('TA*')
        pl.title('PIXEL: %d'%(the_pixel))
        Plots.savefig()                
                
    def pixel_mean_spectrum_plot2(self, pixel_list, rms_cut, location, radius, use_mean=False, use_diff=False):
        """
        Overplots mean spectra plot of spectra from a pixel_list
        Args:
            pixel_list (list of int): pixel IDs to plot
            rms_cut (float): rms cutoff value for plot (negative allowed)
            location (list of 2 floats):   location in grid
            radius (float); radius of circle within which points selected
        Returns:
            none
        """

        Plots.figure()
        
        if radius > 0:
            dx = self.xpos - location[0]
            dy = self.ypos - location[1]
            r2 = dx*dx+dy*dy
            rad2 = radius*radius

        npix  =  len(pixel_list)
        nchan =  len(self.caxis)
        sp    = np.zeros(npix*nchan).reshape(npix,nchan)

        for (i,the_pixel) in zip(range(npix),pixel_list):
            pindex = np.where(self.pixel == the_pixel)[0]

            # print("Warning: MAD test",pindex)
            med1 = np.median(self.rms[pindex])
            std1 = mad_std(self.rms[pindex])
            if rms_cut  < 0:
                cut1 = med1 - rms_cut*std1
            else:
                cut1 = rms_cut
            print("MAD:",the_pixel,med1,std1,cut1)
        
            rindex = np.where(self.rms[pindex] < cut1)[0]
            if radius > 0:
                cindex = np.where(r2[pindex[rindex]] < rad2)[0]
                # print('Pixel %d %d' % (the_pixel, len(cindex)))
                sp[i,:] = np.mean(self.data[pindex[rindex[cindex]]], axis=0)
                #pl.plot(self.caxis, sp_i, label="%d" % the_pixel)
            else:
                # print('Pixel %d %d' % (the_pixel, len(rindex)))
                sp[i,:] = np.mean(self.data[pindex[rindex]], axis=0)
                #pl.plot(self.caxis, sp_i, label="%d" % the_pixel)
        sp_mean = np.mean(sp, axis=0)
        for (i,the_pixel) in zip(range(npix),pixel_list):
            if use_diff:
                pl.plot(self.caxis, sp[i,:]-sp_mean, label="%d" % the_pixel)
            else:
                pl.plot(self.caxis, sp[i,:], label="%d" % the_pixel)            
        if use_mean:
            pl.plot(self.caxis,sp_mean,label="M", color='black', linewidth=3)
        
        pl.xlabel(self.ctype)
        pl.ylabel('TA*')
        pl.title('PIXEL: %s   rms_cut: %g'%(str(pixel_list),cut1))
        pl.legend()
        Plots.savefig()
                
    def write_fits(self, pixel_list, fits_file, binning=1):
        """
        Write selected pixels in a waterfall fits cube, each
        pixel is a different plane.
        Along X is the sample (time)
        Along Y in the channel
        Along Z is the pixel
        Args:
            pixel_list (list of int): pixel IDs to plot
            fits_file (string):
            binning (int):  binning along sample, not implemented
        Returns:
            none
        """
        def rebin(arr, new_shape):
            """Rebin 3D array arr to shape new_shape by averaging."""
            shape = (new_shape[0], arr.shape[0] // new_shape[0],
                     new_shape[1], arr.shape[1] // new_shape[1],
                     new_shape[2], arr.shape[2] // new_shape[2])
            return arr.reshape(shape).mean(-1).mean(-2).mean(-3)

        npix  =  len(pixel_list)
        nchan =  len(self.caxis)

        nsamp = 0
        for (i,the_pixel) in zip(range(npix),pixel_list):
            pindex = np.where(self.pixel == the_pixel)[0]
            if len(pindex)==0:
                print("Warning: no data for pixel %d" % the_pixel)
            if nsamp == 0:
                nsamp = len(pindex)
            elif len(pindex) < nsamp:
                nsamp = len(pindex)

        print('pix',pixel_list)
        print('Min nsamp',nsamp)
        print('Binning',binning)

        nsamp1 = (nsamp//binning[0])*binning[0]
        nchan1 = (nchan//binning[1])*binning[1]

        #print("OLD: ",nsamp, nchan)
        #print("NEW: ",nsamp1,nchan1)

        sp    = np.zeros(npix*nchan1*nsamp1, dtype=np.float32).reshape(npix,nsamp1*nchan1)
        print("FITS file will be %d x %d x %d" % (nsamp//binning[0], nchan1//binning[1],npix))

        for (i,the_pixel) in zip(range(npix),pixel_list):
            pindex = np.where(self.pixel == the_pixel)[0]
            if len(pindex) > nsamp1:
                pindex = pindex[:nsamp1]
            if nchan == nchan1:
                sp[i,:] = self.data[pindex].ravel()
            else:
                x = self.data[pindex]
                x = x[:nsamp1,:nchan1]
                sp[i,:] = x.ravel()
                
        # an expensive operation: swapping last two axes
        sp = sp.reshape(npix,nsamp1,nchan1)
        sp = np.moveaxis(sp,1,2)

        if binning[0]*binning[1] > 1:
            sp = rebin(sp,(npix,nchan1//binning[1],nsamp1//binning[0]))
            sp = sp.squeeze()
            print("New shape:",sp.shape)

        hdu = fits.PrimaryHDU(sp)
        if True:
            # ds9 isn't able to display the WCS with this.carta
            hdr = fits.Header()
            # Sample (time)
            hdr['CRPIX1'] = 0.5 + 0.5/binning[0]
            hdr['CRVAL1'] = 1.0
            hdr['CDELT1'] = 1.0 * binning[0]
            hdr['CTYPE1'] = 'Sample'
            # VLSR
            hdr['CRPIX2'] = 0.5 + 0.5/binning[1]
            hdr['CRVAL2'] = float(self.crval.data) 
            hdr['CDELT2'] = float(self.cdelt.data) * binning[1]  
            hdr['CTYPE2'] = 'VELO-LSR'
            hdr['CUNIT2'] = 'km/s'
            # Beam
            hdr['CRPIX3'] = 1.0
            hdr['CRVAL3'] = 0.0
            hdr['CDELT3'] = 1.0
            hdr['CTYPE3'] = 'Beam'
            # List beams
            hdr['COMMENT'] = 'Waterfall plot (Sample-Velocity-Beam)'
            hdr['COMMENT'] = 'pix_list: %s' % str(pixel_list)
            hdr['COMMENT'] = 'binning: %s' % str(binning)
            hdr['COMMENT'] = 'lmtslr 15-mar-2022'
            
        fits.writeto(fits_file, hdu.data, hdr)
        print("Written waterfall cube to %s" % fits_file)


            
