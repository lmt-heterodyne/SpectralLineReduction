import argparse
from lmtslr.utils.configuration import Configuration, \
    otf_config_spec_text, viewspec_config_spec_text, \
    viewcube_config_spec_text, grid_config_text, \
    ps_config_spec_text

def comma_separated(instr):
    if isinstance(instr, str):
        return ','.join([a.strip() for a in instr.split(',')])

class HandleOptions:
    def __init__(self):
        self.attrs = set()

    def read_config_file(self, filename, cfg_spec_text):
        self.cfg = Configuration(filename, cfg_spec_text)
        if self.cfg.tested:
            # Validated
            for section, dic in self.cfg.cfg.items():
                for k, v in dic.items():
                    setattr(self, k, v)
                    self.attrs.add(k)
        if hasattr(self, 'path'):
            self.data_path = self.path
            self.attrs.add('data_path')
        if hasattr(self, 'output'):
            self.output_file_name = self.output
            self.attrs.add('output_file_name')
        if hasattr(self, 'input'):
            self.input_file_name = comma_separated(self.input)
            self.attrs.add('input_file_name')            
        if hasattr(self, 'b_regions'):
            print(self.b_regions)
            self.b_regions = eval(self.b_regions)
        if hasattr(self, 'l_regions'):
            self.l_regions = eval(self.l_regions)
            
    def print_all_options(self):
        print('program %s options' % (self.program))
        for attr in self.attrs:
            print("%s      =   %s" % (attr, getattr(self, attr)))
            
class HandleOTFProcessOptions(HandleOptions):
    def __init__(self):
        HandleOptions.__init__(self)

    def parse_options(self, args, program, arglevel=0, print_options=False):
        self.program = program
        self.arglevel = arglevel
        self.parser = argparse.ArgumentParser(prog=program,
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-c", "--config", dest="config", 
                            help="Name of configuration file to set parameters")
        self.parser.add_argument("-p", "--path", dest="path", default='/data_lmt', 
                            help="data path")
        self.parser.add_argument("-o", "--output", dest="output",
                            help="name of output SpecFile")
        self.parser.add_argument("-O", "--obsnum", dest="obsnum", type=int,
                            help="ObsNum of observation")
        #self.parser.add_argument("--obs_list", dest="obs_list", type=str,
        #                help="Comma separated list of ObsNums")
        self.parser.add_argument("-b", "--bank", dest="bank", type=int,
                            help="Spectral Bank for processing")
        self.parser.add_argument("--pix_list", dest="pix_list", type=str,
                            help="Comma separated list of pixels")
        self.parser.add_argument("--eliminate_list", dest="eliminate_list", type=str,
                            help="Comma separated list of pixels")
        self.parser.add_argument("--use_cal", dest="use_cal", action="store_true",
                            default=False, help="Use Calibration scan")
        self.parser.add_argument("--tsys", dest="tsys", type=float,
                            default=250.0,
                            help="If use_cal is False, value of Tsys to use")
        self.parser.add_argument("--use_otf_cal", dest="use_otf_cal",
                                 action="store_true", default=False, help="Use calibration within OTF scan")
        self.parser.add_argument("--stype", dest="stype", type=int,
                                 default=1, help="type of spectral line reduction; 0 - median; 1 - single ref spectra; 2 - bracketed ref")
        self.parser.add_argument("--x_axis", dest="x_axis",
                            default='VLSR',
                            help="select spectral x axis. options one of VLSR, VSKY, VBARY, VSRC, FLSR, FSKY, FBARY, FSRC")
        self.parser.add_argument("--b_order", dest="b_order", type=int,
                            default=0, help="set polynomial baseline order")
        self.parser.add_argument("--b_regions", dest="b_regions", type=str, default="[[],[]]",
                            help="enter list of lists for baseline regions")
        self.parser.add_argument("--l_regions", dest="l_regions", type=str, default="[[],[]]",
                            help="enter list of lists for line fit regions")
        self.parser.add_argument("--slice", dest="slice", type=str,
                            help="enter list to specify slice from spectrum for processing")

        args = self.parser.parse_args(args)
        if 'help' in args.__dict__:
            self.parser.print_help()
        for k, v in args.__dict__.items():
            if v is not None:
                setattr(self, k, v)
                self.attrs.add(k)
            if k == 'x_axis' and v is not None:
                if self.x_axis not in ['VLSR', 'VSKY', 'VBARY', 'VSRC', 'FLSR', 'FSKY', 'FBARY','FSRC']:
                    self.x_axis = None
            if k in ('pix_list', 'eliminate_list') and v is not None:
                setattr(self, k, list(map(int, v.split(','))))
            if k in ('b_regions', 'l_regions', 'slice') and v is not None:
                setattr(self, k, eval(v))
        if args.path:
            self.data_path = args.path
            self.attrs.add('data_path')
        if args.output:
            self.output_file_name = args.output
            self.attrs.add('output_file_name')
        if 'config' in args.__dict__ and args.config is not None:
            self.read_config_file(args.config, otf_config_spec_text)
        if print_options:
            self.print_all_options()

class HandlePSProcessOptions(HandleOptions):
    def __init__(self):
        HandleOptions.__init__(self)

    def parse_options(self, args, program, arglevel=0, print_options=False):
        self.program = program
        self.arglevel = arglevel
        self.parser = argparse.ArgumentParser(prog=program,
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-c", "--config", dest="config", 
                            help="Name of configuration file to set parameters")
        self.parser.add_argument("-p", "--path", dest="path",
                            help="data path")
        self.parser.add_argument("-o", "--output", dest="output",
                            help="name of output SpecFile")
        #self.parser.add_argument("-O", "--obsnum", dest="obsnum", type=int,
        #                    help="ObsNum of observation")
        self.parser.add_argument("--obs_list", dest="obs_list", type=str,
                        help="Comma separated list of ObsNums")
        self.parser.add_argument("-b", "--bank", dest="bank", type=int,
                            help="Spectral Bank for processing")
        self.parser.add_argument("--pix_list", dest="pix_list", type=str,
                            help="Comma separated list of pixels")
        self.parser.add_argument("--use_cal", dest="use_cal", action="store_true",
                            default=False, help="Use Calibration scan")
        self.parser.add_argument("--tsys", dest="tsys", type=float,
                            default=250.0,
                            help="If use_cal is False, value of Tsys to use")
        self.parser.add_argument("--stype", dest="stype", type=int,
                                 default=1, help="type of spectral line reduction; 0 - median; 1 - single ref spectra; 2 - bracketed ref")
        self.parser.add_argument("--x_axis", dest="x_axis",
                            default='VLSR',
                            help="select spectral x axis. options one of VLSR, VSKY, VBARY, VSRC, FLSR, FSKY, FBARY, FSRC")
        self.parser.add_argument("--b_order", dest="b_order", type=int,
                            default=0, help="set polynomial baseline order")
        self.parser.add_argument("--b_regions", dest="b_regions", type=str,
                            help="enter list of lists for baseline regions")
        self.parser.add_argument("--l_regions", dest="l_regions", type=str,
                            help="enter list of lists for line fit regions")
        self.parser.add_argument("--slice", dest="slice", type=str,
                            help="enter list to specify slice from spectrum for processing")

        args = self.parser.parse_args(args)
        if 'help' in args.__dict__:
            self.parser.print_help()
        for k, v in args.__dict__.items():
            if v is not None:
                setattr(self, k, v)
                self.attrs.add(k)
            if k == 'x_axis' and v is not None:
                if self.x_axis not in ['VLSR', 'VSKY', 'VBARY', 'VSRC', 'FLSR', 'FSKY', 'FBARY','FSRC']:
                    self.x_axis = None
            if k in ('pix_list', 'obs_list') and v is not None:
                setattr(self, k, list(map(int, v.split(','))))
            if k in ('b_regions', 'l_regions', 'slice') and v is not None:
                setattr(self, k, eval(v))
        if args.path:
            self.data_path = args.path
            self.attrs.add('data_path')
        if args.output:
            self.output_file_name = args.output
            self.attrs.add('output_file_name')
        if 'config' in args.__dict__ and args.config is not None:
            self.read_config_file(args.config, ps_config_spec_text)
        if print_options:
            self.print_all_options()
            
class HandleViewSpecFileOptions(HandleOptions):
    def __init__(self):
        HandleOptions.__init__(self)

    def parse_options(self, args, program, arglevel=0, print_options=False):
        self.program = program
        self.arglevel = arglevel
        self.parser = argparse.ArgumentParser(prog=program,
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-c", "--config", dest="config", 
                            help="Name of configuration file to set parameters")
        self.parser.add_argument("-i", "--input", dest="input",
                            help="set input filename")
        self.parser.add_argument("--pix_list", dest="pix_list", type=str,
                            help="Comma separated list of pixels")
        self.parser.add_argument("--show_all_pixels", dest="show_all_pixels",
                                 action="store_true",
                                 default=True, help="Show all pixels")
        self.parser.add_argument("-p", "--show_pixel", dest="show_pixel",
                                 type=int,
                                 help="Show a particular pixel")
        self.parser.add_argument("--rms_cut", dest="rms_cut",
                                 type=float,
                                 help="rms threshold for data")
        self.parser.add_argument("--plot_range", dest="plot_range", type=str,
                                 help="set plot range for plots")
        
        args = self.parser.parse_args(args)
        if 'help' in args.__dict__:
            self.parser.print_help()
        for k, v in args.__dict__.items():
            if v is not None:
                setattr(self, k, v)
                self.attrs.add(k)
            if k in ('pix_list', ) and v is not None:
                setattr(self, k, list(map(int, v.split(','))))
            if k in ('plot_range') and v is not None:
                setattr(self, k, list(map(float, v.split(','))))
        if hasattr(self, 'show_pixel') and self.show_pixel is not None:
            self.show_all_pixels = False
        if args.input:
            self.input_file_name = comma_separated(args.input)
            self.attrs.add('input_file_name')
        if 'config' in args.__dict__ and args.config is not None:
            self.read_config_file(args.config, viewspec_config_spec_text)
        if print_options:
            self.print_all_options()

class HandleViewCubeOptions(HandleOptions):
    def __init__(self):
        HandleOptions.__init__(self)

    def parse_options(self, args, program, arglevel=0, print_options=False):
        self.program = program
        self.arglevel = arglevel
        self.parser = argparse.ArgumentParser(prog=program,
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-c", "--config", dest="config", 
                            help="Name of configuration file to set parameters")
        self.parser.add_argument("-i", "--input", dest="input",
                            help="set input filename")
        self.parser.add_argument("--v_range", dest="v_range", type=str,
                            help="Comma separated pair of velocities in km/s for integrated intensity image")
        self.parser.add_argument("--v_scale", dest="v_scale", type=float,
                                 default=0.001,
                                 help="scale factor for velocity (use 1/1000 to convert m/s to km/s)")
        self.parser.add_argument("--location", dest="location", type=str,
                                 help="Comma separated pair of location positions dx, dy in arcsec for spectrum plot")
        self.parser.add_argument("--scale", dest="scale", type=float,
                                 default=1/3600.,
                                 help="scale factor for position offset (use 1/3600 to convert arcsec to degrees)")
        self.parser.add_argument("--limits", dest="limits", type=str,
                                 help="Comma separated xlo,xhi,ylo,yhi for final map")
        self.parser.add_argument("--tmax_range", dest="tmax_range", type=str,
                                 help="Comma separated data range data_low, data_high for Tmax image")
        self.parser.add_argument("--tint_range", dest="tint_range", type=str,
                                 help="Comma separated data range data_low, data_high for Tint image")        
        self.parser.add_argument("--plot_type", dest="plot_type", type=str,
                                 help="Plot type data to plot - valid options: TINT, TMAX")
        self.parser.add_argument("--interpolation", dest="interp", type=str,
                                 default='bilinear',
                                 help="interpolation scheme one of (none, nearest, bilinear, bicubic")
        
        args = self.parser.parse_args(args)
        if 'help' in args.__dict__:
            self.parser.print_help()
        for k, v in args.__dict__.items():
            if v is not None:
                setattr(self, k, v)
                self.attrs.add(k)
            if k in ('v_range', 'location', 'limits', 'tmax_range', 'tint_range') and v is not None:
                setattr(self, k, list(map(float, v.split(','))))
        if args.input:
            self.input_file_name = comma_separated(args.input)
            self.attrs.add('input_file_name')
        if 'config' in args.__dict__ and args.config is not None:
            self.read_config_file(args.config, viewcube_config_spec_text)
        if print_options:
            self.print_all_options()
    
class HandleGridOptions(HandleOptions):
    def __init__(self):
        HandleOptions.__init__(self)

    def parse_options(self, args, program, arglevel=0, print_options=False):
        self.program = program
        self.arglevel = arglevel
        self.parser = argparse.ArgumentParser(prog=program,
                                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-c", "--config", dest="config", 
                            help="Name of configuration file to set parameters")
        self.parser.add_argument("-p", "--program_path",  dest="program_path",
                                 default='/usr/local/env/specenv/bin/spec_driver_fits',
                                 help="full path name to grid program")
        self.parser.add_argument("-i", "--input", dest="input",
                            help="set input spec filename")
        self.parser.add_argument("-o", "--output", dest="output",
                            help="name of output fits file")        
        self.parser.add_argument("-w", "--weight", dest="weight",
                            help="name of output weight fits file")        
        self.parser.add_argument("--resolution", dest="resolution", type=float,
                                 default=14.0,
                                 help="resolution (arcsec)")
        self.parser.add_argument("--cell", dest="cell", type=float,
                                 default=7.0,
                                 help="cell size (arcsec)")
        self.parser.add_argument("--pix_list", dest="pix_list", type=str,
                            help="Comma separated list of pixels")        
        self.parser.add_argument("--rms_cut", dest="rms_cut",
                                 type=float,
                                 help="rms threshold for data")
        self.parser.add_argument("--noise_sigma", dest="noise_sigma",
                                 type=float, default=1.0,
                                 help="noise weighting - apply if > 0")                
        self.parser.add_argument("--x_extent", dest="x_extent", type=float,
                                 help="x extent of cube (arcsec) note: cube will go to +/- x_extent")
        self.parser.add_argument("--y_extent", dest="y_extent", type=float,
                                 help="y extent of cube (arcsec) note: cube will go to +/- y_extent")
        self.parser.add_argument("--otf_select", dest="otf_select", type=int,
                                 default=1,
                                 help="otf filter code one of (0=box, 1=jinc,2=gaussian)")
        self.parser.add_argument("--rmax", dest="rmax", type=float,
                                 default=3.0,
                                 help="maximum radius of convolution (units lambda/D)")
        self.parser.add_argument("--n_samples", dest="n_samples", type=int,
                                 default=256,
                                 help="number of samples in convolution filter")
        self.parser.add_argument("--otf_a", dest="otf_a", type=float,
                                 default=1.1,
                                 help="OTF A parameter")
        self.parser.add_argument("--otf_b", dest="otf_b", type=float,
                                 default=4.75,
                                 help="OTF B parameter")
        self.parser.add_argument("--otf_c", dest="otf_c", type=float,
                                 default=2.0,
                                 help="OTF C parameter")        
        
        args = self.parser.parse_args(args)
        if 'help' in args.__dict__:
            self.parser.print_help()
        for k, v in args.__dict__.items():
            if v is not None:
                setattr(self, k, v)
                self.attrs.add(k)
            if k in ('pix_list', ) and v is not None:
                setattr(self, k, list(map(int, v.split(','))))                
        if args.input:
            self.input_file_name = comma_separated(args.input)
            self.attrs.add('input_file_name')
        if args.output:
            self.output_file_name = args.output
            self.attrs.add('output_file_name')
        if args.weight:
            self.weight_file_name = args.weight
            self.attrs.add('weight_file_name')
        if 'config' in args.__dict__ and args.config is not None:
            self.read_config_file(args.config, grid_config_text)
        if hasattr(self, 'pix_list'):
            self.pix_list = '[' + ','.join(['%d' % d for d in self.pix_list]) + ']'
        if print_options:
            self.print_all_options()
        
    
