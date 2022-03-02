#!/usr/bin/env python


"""Usage: view_spec_point.py  -i INPUT [options]

-i INPUT --input INPUT        Input SpecFile (no default)
--pix_list PIX_LIST           Comma separated list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
--location X,Y                Location [Default: 0,0]
-r RAD --radius RAD           Radius around location in arcsec. Use < 0 for all spectra. [Default: -1] 
-z RMS_CUT --rms_cut RMS_CUT  RMS threshold for data, negative allowed for robust MAD method  [Default: 10.0]
--mean                        Show the mean spectrum as well.
--diff                        Show the difference from the mean instead of the spectrum
--plot_range PLOT_RANGE       Plotting range (deprecated) [Default: -1,4]
-p PIXEL --show_pixel PIXEL   Show a particular pixel flags code (deprecated) [Default: -1]
--plots METHOD                Plotting method for batch, defaults to interactive.
-h --help                     show this help


VIEW_SPEC_POINT will for selected pixels and a selected location and
radius in the field plot their mean spectra, a different color for each
pixel. No weighting scheme is applied.

For convenience it still has an option (using -p) to make it behave
like the old VIEW_SPEC_FILE but only mean_spectra_plot and rms_plot
are plotted. This mode does not support rms_cut<0 or the --plots flag.

The --plots METHOD is a new feature to allow switching between interactive
(the default method) and a batch style. For example
    option:    --plots M51
creates M51.1.png (and .2., .3. etc. as many as there are). But
    option:    --plots M31,pdf,3
 would produce M31.3.pdf (and .4., .5., etc. as many as there are).
Plot files are silently overwritten if they existed before.

"""

import sys
from docopt import docopt
import matplotlib.pyplot as pl
from lmtslr.utils.argparser import HandleViewSpecFileOptions
from lmtslr.viewer.spec_file_viewer import SpecFileViewer
from lmtslr.viewer.plots import Plots
import lmtslr.utils.convert as acv

def main(argv):
    av = docopt(__doc__,options_first=True, version='0.3')

    nc_file    = av['--input']
    pix_list   = acv.listi(av['--pix_list'],  16)
    location   = acv.listf(av['--location'],   2)
    radius     = acv.listf(av['--radius'],     1)
    rms_cut    = acv.listf(av['--rms_cut'],    1)
    show_pixel = acv.listi(av['--show_pixel'], 1)   # deprecating
    plot_range = acv.listf(av['--plot_range'], 2)   # deprecating
    plots      = av['--plots']
    use_mean   = av['--mean']
    use_diff   = av['--diff']

    Plots.init(plots)
    
    SV = SpecFileViewer(nc_file)

    if show_pixel < 0:
        # new routine
        SV.pixel_mean_spectrum_plot2(pix_list, rms_cut, location, radius, use_mean, use_diff)
    else:
        # classic one
        SV.pixel_mean_spectrum_plot(show_pixel, rms_cut)
        SV.pixel_rms_plot(show_pixel, rms_cut, plot_range=[0.,plot_range[1]])

    Plots.show()
    
if __name__ == '__main__':
    main(sys.argv[1:])
