""" Module for viewing SpecBank data

classes: SpecViewer, SpecBankViewer, SpecCalViewer
uses: SpecBank, roach, IFProc, Grid
author: FPS
date: May 2018
changes:
python 3
"""
import matplotlib.pyplot as pl
import matplotlib.mlab as mlab

from lmtslr.ifproc.ifproc import IFProc
from lmtslr.spec.spec import SpecBank
from lmtslr.grid.grid import Grid

class SpecViewer():
    ''' Base class for viewer.  This handles some plotting basics.
    '''
    def __init__(self,figure=1):
        """ constructor
        """
        self.figure = figure

    def set_figure(self,figure):
        """ sets the figure for the plots
        """
        self.figure = figure

    def open_figure(self):
        """ opens the figure window
        """
        pl.figure(self.figure)
        pl.clf()

    def close_figure(self):
        """ closes the figure window
        """
        pl.close(self.figure)


class SpecBankViewer(SpecViewer):
    def __init__(self,figure=1):
        """ SpecBankViewer provides methods to view data in a SpecBank
            S is the SpecBank to be viewed
        """
        SpecViewer.__init__(self,figure)
        
    def raw_waterfall(self,S,pixel,window,plot_range):
        """ raw_waterfall makes a waterfall plot of the reduced data; no baseline removed
            window indicates start and stop indices in the spectra time series: [start,stop]
            plot_range is the range of the intensity scale: [ta_min, ta_max]
        """
        index = S.find_pixel_index(pixel)
        r = range(window[0],window[1])
        pl.imshow(S.roach[index].reduced_spectra[r],clim=plot_range)
        pl.xlabel('Channel')
        pl.ylabel('Spectrum')
        pl.suptitle('Obsnum:%d Pixel:%d'%(S.ifproc.obsnum,pixel))

    def plot_spectrum(self,S,pixel,ispec,plot_axis,baseline_list,n_baseline_list):
        """ plots a specific spectrum in the time series for a specific pixel
            S is the SpecBank object
            pixel is the pixel id for the spectra in the plot
            ispec identifies the specific spectrum in the time series to be plotted
            plot_axis gives the desired axis for the plot: [xlow, xhigh, ylow, yhigh]
            baseline_list is a list of channels which will be averaged to provide a constant baseline
            n_baseline_list is number of channels in the baseline_list
        """
        index = S.find_pixel_index(pixel)
        baseline = S.roach[index].baseline(ispec,baseline_list,n_baseline_list)
        pl.plot((S.roach[pixel].reduced_spectra[ispec]-baseline))
        pl.axis(plot_axis)

    def find_peak_spectrum(self,S,pixel):
        """ locates spectrum in time series which gives the maximum value in map
        """
        map_index = S.find_map_pixel_index(pixel)
        mx = np.max(S.map_data[map_index])
        ispec = np.where(S.map_data[map_index] == mx)
        return(int(ispec[0]))
        
    def plot_peak_spectrum(self,S,pixel,plot_axis,baseline_list,n_baseline_list):
        """ plots the spectrum which gives maximum value in map
            S is the SpecBank object
            pixel is the pixel id for the spectra in the plot
            plot_axis gives the desired axis for the plot: [xlow, xhigh, ylow, yhigh]
            baseline_list is a list of channels which will be averaged to provide a constant baseline
            n_baseline_list is number of channels in the baseline_list
        """
        index = S.find_pixel_index(pixel)
        ispec = self.find_peak_spectrum(S,pixel)
        x = S.roach[index].reduced_spectra[ispec]
        baseline = np.sum(x[baseline_list])/n_baseline_list
        v = S.create_velocity_scale()
        prange = np.where(np.logical_and(v >= plot_axis[0], v <= plot_axis[1]))
        plot_axis[2] = (x-baseline)[prange].min()*1.1
        plot_axis[3] = (x-baseline)[prange].max()*1.1
        legend_label = 'Pixel %d\nPeak Spectrum %d'%(pixel,ispec)
        pl.plot(v,(x-baseline),label=legend_label)
        legend = pl.legend(fontsize='x-small')
        pl.xlabel('Velocity (km/s)')
        pl.suptitle('ObsNum %d: %s %s %sGHz\n Pixel %d Peak Spectrum %d'%(S.obsnum,S.receiver,S.source,S.line_rest_frequency,pixel,ispec)) 
        pl.axis(plot_axis)

    def plot_all_spectra(self,S,pixel,plot_axis,baseline_list,n_baseline_list):
        pixel_index = S.find_pixel_index(pixel)
        v = S.create_velocity_scale()
        peak_index = np.argmax(S.map_data[pixel_index])
        prange = np.where(np.logical_and(v >= plot_axis[0], v <= plot_axis[1]))
        plot_axis[2] = S.map_spectra[pixel_index][:,prange].min()*1.1
        plot_axis[3] = S.map_spectra[pixel_index][:,prange].max()*1.1
        plot_axis2 = np.zeros(4)
        plot_axis2[0] = plot_axis[0]
        plot_axis2[1] = plot_axis[1]
        plot_axis2[2] = 0#S.map_data.min()*1.1
        plot_axis2[3] = S.map_data.max()*1.1
        plen = len(S.map_spectra[pixel_index])
        point_list = []
        xlen = ylen = int(math.sqrt(plen))

        # create an 2d list of array indices
        a = [[i+j*ylen for i in range(xlen)] for j in range(ylen)]
        # change into a numpy array to manipulate
        a = np.array(a)
        # flip the order of every other row
        a[1::2, :] = a[1::2, ::-1]
        # flip the whole array
        a = np.flipud(a)
        # flatten the array to get a 1d list
        a = a.flatten()

        for index in range(plen):
            plot_index = a[index]
            ax = pl.subplot(xlen,ylen,plot_index+1)
            ax.tick_params(axis='both',which='major',labelsize=6)
            ax.tick_params(axis='both',which='minor',labelsize=6)
            ax.plot(v,(S.map_spectra[pixel_index][index]-np.sum(S.map_spectra[pixel_index][index][S.blist])/S.nb))
            ax.axis(plot_axis)
            ax.text(plot_axis[0]+0.1*(plot_axis[1]-plot_axis[0]),plot_axis[3]-0.2*(plot_axis[3]-plot_axis[2]),'%5.2f %5.2f %5.2f'%(S.map_x[pixel_index][index],S.map_y[pixel_index][index], S.map_data[pixel_index][index]),size='6')
            ax2 = ax.twinx()
            ax2.tick_params(axis='both',which='major',labelsize=6)
            ax2.tick_params(axis='both',which='minor',labelsize=6)
            ax2.plot(0.5*(plot_axis[0]+plot_axis[1]),S.map_data[pixel_index][index], 'or')
            ax2.axis(plot_axis2)
        pl.tight_layout(rect=[0, 0.03, 1, 0.9])
        pl.suptitle('ObsNum %d: %s %s %sGHz\n Pixel %d'%(S.obsnum,S.receiver,S.source,S.line_rest_frequency,pixel)) 


    def waterfall(self,S,pixel,window,plot_range,baseline_list,n_baseline_list):
        """ waterfall makes a waterfall plot of the reduced data with baseline removed
            S is the SpecBank object
            pixel is the pixel id for the spectra in the plot
            window indicates start and stop indices in the spectra time series: [start,stop]
            plot_range is the range of the intensity scale: [ta_min, ta_max]
            baseline_list is a list of channels which will be averaged to provide a constant baseline
            n_baseline_list is number of channels in the baseline_list
        """
        index = S.find_pixel_index(pixel)
        ispec = self.find_peak_spectrum(S,pixel)
        start = ispec+window[0]
        stop = ispec+window[1]
        if start < 0:
            start = 0
        if stop > S.roach[index].nspec:
            stop = S.roach[index].nspec
        r = range(start,stop)
        w = np.zeros((len(r),S.nchan))
        for i,ispec in enumerate(r):
            baseline,rms = S.roach[index].baseline(S.roach[index].reduced_spectra[ispec],baseline_list,n_baseline_list)
            for j in range(S.nchan):
                w[i,j] = (S.roach[index].reduced_spectra[ispec][j] - baseline[j])
        pl.imshow(w,clim=plot_range,origin='lower',extent=[S.c2v(0),S.c2v(S.nchan-1),start,stop])
        pl.suptitle('Obsnum: %d Pixel: %d Scans: %d to %d'%(S.obsnum,pixel,start,stop))
        pl.xlabel('Velocity (km/s)')
                

    def sanchez_map(self,S,map_region,grid_spacing,plot_range,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]):
        """ sanchez_map makes a "Sanchez Map" with all pixels displayed on sky in their grid positions
            S is the SpecBank object
            map_region is the extent of the map: [low left, low right, high left, high right] (arcsec)
            grid_spacing is the size of the map cells (arcsec)
            pixel_list identifies the beams to be plotted
            Note: you have to create the map data before running this.
        """
        nx = (map_region[1]-map_region[0])/grid_spacing+1
        ny = (map_region[3]-map_region[2])/grid_spacing+1
        xi = np.linspace(map_region[0],map_region[1],nx)
        yi = np.linspace(map_region[2],map_region[3],ny)
        zi_sum = np.zeros((nx,ny))
        for pixel in pixel_list:
            index = S.find_map_pixel_index(pixel)
            zi = mlab.griddata(S.map_x[index],S.map_y[index],S.map_data[index],xi,yi,interp='linear')
            zi_sum = zi_sum + zi
        pl.imshow(zi_sum,clim=plot_range,interpolation='bicubic',cmap=pl.cm.jet,origin='lower',extent=map_region)
        pl.xlabel('dAz (")')
        pl.ylabel('dEl (")')
        pl.suptitle('Spectral Line Sanchez Map: %d'%(S.obsnum))
        
    def map(self,S,map_region,grid_spacing,plot_range,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]):
        """ map aligns the individual maps for each pixel according to a nominal grid model
            S is the SpecBank object
            map_region is the extent of the map: [low left, low right, high left, high right] (arcsec)
            grid_spacing is the size of the map cells (arcsec)
            pixel_list identifies the beams to be plotted
            Note: you have to create the map data before running this.
        """
        g = grid()
        gx,gy = g.azel(S.elev/180.*np.pi,S.tracking_beam)
        nx = (map_region[1]-map_region[0])/grid_spacing+1
        ny = (map_region[3]-map_region[2])/grid_spacing+1
        xi = np.linspace(map_region[0],map_region[1],nx)
        yi = np.linspace(map_region[2],map_region[3],ny)
        zi_sum = np.zeros((nx,ny))
        wi_sum = np.zeros((nx,ny))
        for pixel in pixel_list:
            index = S.find_map_pixel_index(pixel)
            wdata = np.ones(len(S.map_data[index]))
            zi = mlab.griddata(S.map_x[index]-gx[pixel],S.map_y[index]-gy[pixel],S.map_data[index],xi,yi,interp='linear')
            wi = mlab.griddata(S.map_x[index]-gx[pixel],S.map_y[index]-gy[pixel],wdata,xi,yi,interp='linear')
            zi_sum = zi_sum + zi
            wi_sum = wi_sum + wi
        pl.imshow(zi_sum/wi_sum,clim=plot_range,interpolation='bicubic',cmap=pl.cm.jet,origin='lower',extent=map_region)
        pl.axis(map_region)
        pl.xlabel('dAz (")')
        pl.ylabel('dEl (")')
        pl.suptitle('Spectral Line Aligned Map: %d'%(S.obsnum))

    def pixel_map(self,S,pixel,map_region,grid_spacing,show_points=False):
        """ pixel_map maps results for a single pixel
            S is the SpecBank object
            pixel is the pixel id to be plotted
            map_region is the extent of the map: [low left, low right, high left, high right] (arcsec)
            grid_spacing is the size of the map cells (arcsec)
            show_points=True will plot the locations of the spectra in the map
            Note: you have to create the map data before running this.
        """
        index = S.find_map_pixel_index(pixel)
        nx = (map_region[1]-map_region[0])/grid_spacing+1
        ny = (map_region[3]-map_region[2])/grid_spacing+1
        xi = np.linspace(map_region[0],map_region[1],nx)
        yi = np.linspace(map_region[2],map_region[3],ny)
        zi = mlab.griddata(S.map_x[index],S.map_y[index],S.map_data[index],xi,yi,interp='linear')
        pl.imshow(zi,interpolation='bicubic',cmap=pl.cm.jet,origin='lower',extent=map_region)
        pl.xlabel('X (")')
        pl.ylabel('Y (")')
        pl.suptitle('Spectral Line Map: Obsnum: %d  Pixel: %d'%(S.ifproc.obsnum,pixel))
        pl.axis('equal')
        if show_points:
            pl.plot(S.map_x[index],S.map_y[index],'w.')
            pl.axis(map_region)
        pl.colorbar()

    def plot_on(self,S):
        plot_order = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16]
        for ipix in range(S.npix):
            pixel_id = S.roach_pixel_ids[ipix]
            ax = pl.subplot(4,4,plot_order[pixel_id])
            ax.tick_params(axis='both',which='major',labelsize=6)
            ax.tick_params(axis='both',which='minor',labelsize=6)
            ax.plot(S.roach[ipix].on_spectrum)
            l = len(S.roach[ipix].on_spectrum)
            pl.xticks(np.arange(0, l+1, l/4))
        pl.suptitle('%s: ObsNum %d\n%s %s GHz'%(S.obspgm,S.obsnum,S.receiver,S.line_rest_frequency)) 

    def plot_ps(self,S,baseline_order,plot_axis=[-200,200,-0.5,2.0],line_stats_all=[]):
        plot_order = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16]
        line_stats = line_stats_all[0]
        prange = np.where(np.logical_and(line_stats.v >= plot_axis[0], line_stats.v <= plot_axis[1]))
        plot_axis[2] = line_stats.spectrum[prange].min()*1.1
        plot_axis[3] = line_stats.spectrum[prange].max()*1.4
        for ipix in range(S.npix):
            if ipix == 0: continue
            line_stats = line_stats_all[ipix]
            min_ps = line_stats.spectrum[prange].min()*1.1
            max_ps = line_stats.spectrum[prange].max()*1.4
            if min_ps < plot_axis[2]: plot_axis[2] = min_ps
            if max_ps > plot_axis[3]: plot_axis[3] = max_ps

        for ipix in range(S.npix):
            pixel_id = S.roach_pixel_ids[ipix]
            ax = pl.subplot(4,4,plot_order[pixel_id])
            ax.tick_params(axis='both',which='major',labelsize=6)
            ax.tick_params(axis='both',which='minor',labelsize=6)

            # for each line, fit baseline and compute line statistics
            line_stats = line_stats_all[ipix]
            ax.plot(line_stats.v,line_stats.spectrum)
            ax.axis(plot_axis)
            xtext = plot_axis[0]+0.05*(plot_axis[1]-plot_axis[0])
            ytext = plot_axis[3]-0.05*(plot_axis[3]-plot_axis[2])
            ax.text(xtext,ytext,'%2d I=%8.3f(%8.3f)'%(pixel_id,line_stats.yint,line_stats.yerr),horizontalalignment='left',verticalalignment='top',fontsize=6)
            ytext = plot_axis[3]-0.15*(plot_axis[3]-plot_axis[2])
            pl.text(xtext,ytext,'V=%8.3f RMS=%8.3f'%(line_stats.xmean,line_stats.rms),horizontalalignment='left',verticalalignment='top',fontsize=6)
        pl.suptitle('%s: ObsNum %d\n%s %s GHz'%(S.obspgm,S.obsnum,S.receiver,S.line_rest_frequency)) 

    def plot_bs(self,S,baseline_order,plot_axis=[-200,200,-0.5,2.0],line_stats=None):
        plot_order = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16];
        prange = np.where(np.logical_and(line_stats.v >= plot_axis[0], line_stats.v <= plot_axis[1]))
        plot_axis[2] = line_stats.spectrum[prange].min()*1.1
        plot_axis[3] = line_stats.spectrum[prange].max()*1.4
        ax = pl.subplot(1,1,1)
        ax.tick_params(axis='both',which='major',labelsize=6)
        ax.tick_params(axis='both',which='minor',labelsize=6)

        ax.plot(line_stats.v,line_stats.spectrum)
        ax.axis(plot_axis)
        xtext = plot_axis[0]+0.05*(plot_axis[1]-plot_axis[0])
        ytext = plot_axis[3]-0.05*(plot_axis[3]-plot_axis[2])
        ax.text(xtext,ytext,'%2d/%2d I=%8.3f(%8.3f)'%(S.roach_pixel_ids[0],S.roach_pixel_ids[1],line_stats.yint,line_stats.yerr),horizontalalignment='left',verticalalignment='top',fontsize=10)
        ytext = plot_axis[3]-0.1*(plot_axis[3]-plot_axis[2])
        pl.text(xtext,ytext,'V=%8.3f RMS=%8.3f'%(line_stats.xmean,line_stats.rms),horizontalalignment='left',verticalalignment='top',fontsize=10)
        pl.suptitle('%s: ObsNum %d\n%s %s GHz'%(S.obspgm,S.obsnum,S.receiver,S.line_rest_frequency)) 


class SpecCalViewer(SpecViewer):
    def __init__(self,figure=1):
        """ SpecCalViewer provides methods to view data in a spec_cal
            S is the spec_cal to be viewed
        """
        SpecViewer.__init__(self,figure)
        
    def plot_tsys(self,S):
        plot_scale = 0.
        nscale = 0
        for ipix in range(S.npix):
            indx_fin = np.where(np.isfinite(S.roach[ipix].tsys_spectrum))
            indx_inf = np.where(np.isinf(S.roach[ipix].tsys_spectrum))
            indx_nan = np.where(np.isnan(S.roach[ipix].tsys_spectrum))
            print(ipix,'fin------------',indx_fin[0])
            print(ipix,'inf------------',indx_inf[0])
            print(ipix,'nan------------',indx_nan[0])
            l_fin = len(indx_fin[0])
            if l_fin > 0 and S.roach[ipix].tsys > 0 and S.roach[ipix].tsys < 500:
                plot_scale = plot_scale + S.roach[ipix].tsys
                nscale = nscale + 1
        if nscale > 0:
            plot_scale = plot_scale/nscale
        plot_order = [1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16];
        for ipix in range(S.npix):
            pixel_id = S.roach_pixel_ids[ipix]
            ax = pl.subplot(4,4,plot_order[pixel_id])
            ax.tick_params(axis='both',which='major',labelsize=6)
            ax.tick_params(axis='both',which='minor',labelsize=6)
            indx_fin = np.where(np.isfinite(S.roach[ipix].tsys_spectrum))
            l_fin = len(indx_fin[0])
            if l_fin > 0:
                pl.plot(S.roach[ipix].tsys_spectrum)
                pl.text(S.nchan/2,10,'%d %6.0fK'%(pixel_id,S.roach[ipix].tsys),horizontalalignment='center')
                pl.axis([0,S.nchan,0,plot_scale*1.5])
            else:
                pl.text(0.1,0.5,'%d NaN'%(pixel_id))
        pl.suptitle('TSys: ObsNum %d\n%s %s GHz'%(S.obsnum,S.receiver,S.line_rest_frequency)) 
