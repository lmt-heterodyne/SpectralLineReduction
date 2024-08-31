#!/usr/bin/env python
#

"""Usage: view_spec_file.py -i INPUT [options]

-i INPUT --input INPUT                 Input SpecFile filename (default: None)
--pix_list PIX_LIST                    list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
-p SHOW_PIXEL --show_pixel SHOW_PIXEL  Show one particular pixel (default is all pixels from pix_list)
--rms_cut RMS_CUT                      rms threshold for data [Default: 10]
--plot_range PLOT_RANGE                set plot range for plots [Default: -1,3]
--box BOX                              Fix boxsize -BOX .. BOX in X and Y. Optional.
--plots METHOD                         Plotting method for output plot defaults to on screen.
                                       Optional file type extension after a comma, e.g. png or pdf
--skip_tsys                            Skip the tsys plot [Default: False]

-h --help                              show this help


Reads a SpecFile from a single OTF mapping observation and creates
various visualizations

Tsys is only present for SpecFiles created after 6-mar-2021, data produce before this needs --skip_tsys.

Bug:  tsys pixel number can be wrong is not all pixels are in the SpecFile

The --plots METHOD is a new feature to allow switching between interactive
(the default method) and a batch style. For example
    option:    --plots M51
creates M51.1.png (and .2., .3. etc. as many as there are). But
    option:    --plots M31,pdf,3
 would produce M31.3.pdf (and .4., .5., etc. as many as there are).
Plot files are silently overwritten if they existed before.

"""

# Python Imports
import sys

# command line parsing
from docopt import docopt
import lmtslr.utils.convert as acv


def main(argv):
    av = docopt(__doc__,options_first=True, version='0.1')
    print(av)   # debug

    input_file_name = av['--input']
    pix_list        = acv.listi(av['--pix_list'],  16)
    rms_cut         = acv.listf(av['--rms_cut'],    1)
    plot_range      = acv.listf(av['--plot_range'], 2)
    plots           = av['--plots']
    box             = av['--box']
    skip_tsys       = av['--skip_tsys']

    import matplotlib
    if plots == None:
        matplotlib.use('qt5agg')
    else:
        matplotlib.use('agg')
    import matplotlib.pyplot as pl

    #from lmtslr.utils.parser import HandleViewSpecFileOptions
    from lmtslr.utils.argparser import HandleViewSpecFileOptions
    from lmtslr.viewer.spec_file_viewer import SpecFileViewer
    from lmtslr.viewer.plots import Plots
        
    if av['--show_pixel'] != None:
        show_all_pixels = False        
        show_pixel  = acv.listi(av['--show_pixel'], 1)
    else:
        show_all_pixels = True
    
    Plots.init(plots)
    
    SV = SpecFileViewer(input_file_name)

    if show_all_pixels:
        SV.xy_position_plot(box=box)
        #SV.sx_position_plot()
        #SV.sy_position_plot()        
        SV.sequoia_waterfall_plot(pix_list, rms_cut, plot_range=plot_range)
        SV.sequoia_rms_plot(pix_list, rms_cut, plot_range=[0.,plot_range[1]])
        SV.sequoia_rms_histogram(pix_list, rms_cut)
        SV.sequoia_mean_spectra_plot(pix_list, rms_cut)
        if not skip_tsys: SV.sequoia_tsys_spectra_plot(pix_list)
    else:
        SV.xy_position_plot(False,box=box)
        #SV.sx_position_plot(False)
        #SV.sy_position_plot(False)
        SV.pixel_waterfall_plot(show_pixel, rms_cut, plot_range=plot_range)
        SV.pixel_rms_plot(show_pixel, rms_cut, plot_range=[0.,plot_range[1]])
        SV.pixel_rms_histogram(show_pixel, rms_cut)
        SV.pixel_mean_spectrum_plot(show_pixel, rms_cut)
        if not skip_tsys: SV.pixel_tsys_spectra_plot(show_pixel)
    
    Plots.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])
