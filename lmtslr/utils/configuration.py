"""
Reads a default configuration file

This uses python's configobj library
"""

from configobj import ConfigObj, flatten_errors
from validate import Validator
import os

otf_config_spec_text = """
#lmtslr general items
[general]
# input path for scan files
path = string(max=500, default='/data_lmt')
# output filename
output = string(max=500, default='./output.nc')
# If reducing single Obsnum, the integer Observation number
obsnum = integer(1, 100000)

[spectra]
# spectral bank for processing
bank = integer(min=0, max=3, default=0)
# list of pixels to process
pix_list = int_list(min=1, max=32, default=list(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
# channels to eliminate
eliminate_list = int_list(min=1, default=list(4096,))
# spectral reduction type
stype = integer(min=0, max=3, default=1)
# If use_cal is True find and use the calibration scan
use_cal = boolean(default=False)
# if use_cal is False use this Tsys
tsys = float(min=10, default=200)
# use CAL within OTF scan
use_otf_cal = boolean(default=False)
# select spectral x axis
x_axis = option('VLSR', 'VSKY', 'VBARY', 'VSRC', 'FLSR', 'FSKY', 'FBARY', 'FSRC', default='VLSR')
# baseline order
b_order = integer(min=-1, max=4, default=0)
# list of lists for baselines
b_regions = string(default='[[-193, -93], [107,207]]')
# list of lists for line fit regions
l_regions = string(default='[[-93, 107]]')
# list to specify slice from spectrum for processing
slice = float_list(min=2, max=2, default=list(-200, 200))
"""

ps_config_spec_text = """
#lmtslr general items
[general]
# input path for scan files
path = string(max=500, default='/data_lmt')
# output filename
output = string(max=500, default='./output.nc')
# If reducing single Obsnum, the integer Observation number
obs_list = int_list()

[spectra]
# spectral bank for processing
bank = integer(min=0, max=3, default=0)
# list of pixels to process
pix_list = int_list(min=1, max=32, default=list(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
# channels to eliminate
eliminate_list = int_list(min=1, default=list(4096,))
# spectral reduction type
stype = integer(min=0, max=3, default=1)
# If use_cal is True find and use the calibration scan
use_cal = boolean(default=False)
# if use_cal is False use this Tsys
tsys = float(min=10, default=200)
# select spectral x axis
x_axis = option('VLSR', 'VSKY', 'VBARY', 'VSRC', 'FLSR', 'FSKY', 'FBARY', 'FSRC', default='VLSR')
# baseline order
b_order = integer(min=-1, max=4, default=0)
# list of lists for baselines
b_regions = string(default='[[-193, -93], [107,207]]')
# list of lists for line fit regions
l_regions = string(default='[[-93, 107]]')
# list to specify slice from spectrum for processing
slice = float_list(min=2, max=2, default=list(-200, 200))
"""

viewspec_config_spec_text = """
#lmtslr general items
[general]
# input path for scan files
input = string(max=500, default='./input.nc')

[spectra]
# show all pixels?
show_all_pixels = boolean(default=True)
# list of pixels to process
pix_list = int_list(min=1, max=32, default=list(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
# show a particular pixel
show_pixel = integer(default=None)
# rms cut off to use
rms_cut = float(default=10000)
# Y-range for plots
plot_range = int_list(min=2, max=2, default=list(-10, 10))
"""

viewcube_config_spec_text = """
#lmtslr general items
[general]
# input path for scan files
input = string(max=500, default='input.fits')

[spectra]
# velocity range in km/s for integrated intensity
v_range = float_list(min=2, max=2, default=list(-300, 200))
# scale factor for velocity (default 1/1000 to convert m/s to km/s
v_scale = float(default=0.001)
# offset in arcsec for spectrum
location = float_list(min=2, max=2, default=list(0, 0))
# scale factor for positional offset (default 1/3600 for arcsec to deg)
scale = float(default=2.78e-4)
# xy limits for final map 
limits = float_list(min=4, max=4, default=list(-100, 100, -100, 100))
# data range for tmax image
tmax_range = float_list(min=2, max=2, default=list(-1, 1))
# data range for tint image
tint_range = float_list(min=2, max=2, default=list(-1, 1))
# data to plot. one of TINT or TMAX
plot_type = option('TINT', 'TMAX', default='TINT')
# interpolation to use. 
interp = option('none', 'nearest', 'bilinear', 'bicubic', default='bilinear')
"""

grid_config_text = """
#lmtslr general items
[general]
# full path name to the gridding program
program_path = string(max=500, default='/usr/local/env/specenv/bin/spec_driver_fits')
# input NC filename
input = string(max=50000, default='./input.nc')
# output FITS filename
output = string(max=500, default='./output.fits')
# output weight FITS filename
weight = string(max=500, default='./weight.fits')

[cube]
# angular resolution to use for gridding (in arcsec)
resolution = float(min=0, default=14)
# cell size of output grid (in arcsec)
cell = float(min=0, default=7)
# list of pixels to process
pix_list = int_list(min=1, max=32, default=list(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
# rms cut off to use
rms_cut = float(default=10000)
# noise_sigma - apply rms weighting if > 0.0
noise_sigma = float(default=1.0)
# x extent of cube in arcsec (cube will go to +/- x_extent)
x_extent = float(min=1, default=300)
# y extent of cube in arcsec (cube will go to +/- y_extent)
y_extent = float(min=1, default=300)
# filter code for otf gridding, one of (0=box,1=jinc,2=gaussian)
otf_select = option(0, 1, 2, default=1)
# maximum radius of convulution (units lambda/D)
rmax = float(default=3)
# number of samples in convolution filter
n_samples = integer(default=256)
# otf_a parameter
otf_a = float(default=1.1)
# otf_b parameter
otf_b = float(default=4.75)
# otf_c parameter
otf_c = float(default=2.0)
"""


def default_config():
    cfg_spec = ConfigObj(config_spec_text.splitlines(), list_values=False)
    valid = Validator()
    cfg = ConfigObj(configspec=cfg_spec, stringify=True, list_values=True)
    test = cfg.validate(valid, copy=True)
    if test:
        return cfg

def validate_dictionary(cdic):
    """This function validates a dictionary against the config spec here"""
    cfg_spec = ConfigObj(config_spec_text.splitlines(), list_values=False)
    valid = Validator()
    cfg = ConfigObj(cdic, configspec=cfg_spec)
    rtn = cfg.validate(valid, preserve_errors=True)
    if type(rtn) is bool and rtn:
        return True
    else:
        res = flatten_errors(cfg, rtn)
        errortxt = ''
        for row in res:
            errortxt += 'In Section %s, key %s has error: %s' % (row[0], row[1], row[2])
            print(errortxt)
        return False

class Configuration:
    def __init__(self, filename, config_spec):
        """Initializes a config file if does not exist. If exists, uses
        it to validate the file, and setup default initial parameters"""
        self.cfg_spec = ConfigObj(config_spec.splitlines(), list_values=False)
        self.cfg_filename = filename
        valid = Validator()
        if not os.path.exists(self.cfg_filename):
            #no configuration file found
            print("File %s not found, so creating one from you from defaults" % self.cfg_filename)
            cfg = ConfigObj(configspec=self.cfg_spec, stringify=True, list_values=True)
            cfg.filename = self.cfg_filename
            test = cfg.validate(valid, copy=True)
            cfg.write()
        self.cfg = ConfigObj(self.cfg_filename, configspec=self.cfg_spec)
        rtn = self.cfg.validate(valid, preserve_errors=True)
        if type(rtn) is bool and rtn:
            print("Config file validated")
            self.tested = True
        else:
            self.tested = False
            res = flatten_errors(self.cfg, rtn)
            self.errortxt = ''
            for row in res:
                self.errortxt += 'In Section %s, key %s has error: %s' % (row[0], row[1], row[2])
            print(self.errortxt)
