#!/usr/bin/env python


"""Usage: make_spec_fits.py  -i INPUT -o OUTPUT [options]

-i INPUT --input INPUT        Input SpecFile (no default)
-o OUTPUT --output OUTPUT     Output FITS file (no default)
--pix_list PIX_LIST           Comma separated list of pixels [Default: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
--binning BINT,BINC           Binning applied in the TIME,CHAN dimension [Default: 1,1]
-h --help                     show this help


MAKE_SPEC_FITS converts a specfile to a waterfall plot in the form of
a FITS file, one plane per pixel. Since the number of samples could be
slightly different per roach board (which covers 4 pxiels), there
could be a slight synchronization mis-match between some pixels.

"""

import sys
from docopt import docopt
from lmtslr.viewer.spec_file_viewer import SpecFileViewer
import lmtslr.utils.convert as acv

def main(argv):
    av = docopt(__doc__,options_first=True, version='0.3')

    nc_file    = av['--input']
    fits_file  = av['--output']
    pix_list   = acv.listi(av['--pix_list'],  16)
    binning    = acv.listi(av['--binning'],  2)

    SV = SpecFileViewer(nc_file)

    SV.write_fits(pix_list, fits_file, binning)
    
if __name__ == '__main__':
    main(sys.argv[1:])
