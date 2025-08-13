import os
import sys
import getopt


class HandleProcessOptions():
    """
    Class to handle parameters for the processing script.
    Args:
        none
    """
    def __init__(self):
        # initialization
        self.data_path = '/data_lmt/'
        self.obsnum = 0
        self.obs_list = []
        self.bank = 0
        self.pix_list = []
        self.eliminate_list = []
        self.use_cal = False
        self.tsys = 0.
        self.x_axis = 'VLSR'
        self.b_order = 0
        self.b_regions = [[],[]]
        self.l_regions = [[],[]]
        self.slice = []
        self.offx = 0
        self.offy = 0
        self.output_file_name = ''
        
    def read_config_file(self, filename):
        """
        Returns the options from the config file.
        Args:
            filename (str): name of config file
        Returns:
            s2 (list): list of options generated from config file
        """
        fd = open(filename,'r')
        f = fd.readlines()
        s = ''
        for line in f:
            x = line.split('#')
            s = s + '--' + x[0]
        s1 = s.replace('\n',' ')
        s2 = s1.split()
        return s2
        
    def parse_options(self, s, program, arglevel=0, print_options=False):
        """
        Sets options according to input arguments and returns the exit 
            code.
        Args:
            s (str): arguments to script
            program (script): program to return
            arglevel (int): argument level (default is 0)
            print_options (bool): option to print options
        Returns:
            result (int): exit code of setting options
        """
        self.program = program
        self.arglevel = arglevel
        try:
            opts, arg = getopt.getopt(s,
                                      "hc:p:O:i:o:",
                                      ["config=",
                                       "path=",
                                       "obsnum=",
                                       "obs_list=",
                                       "bank=",
                                       "pix_list=",
                                       "eliminate_list=",
                                       "use_cal=",
                                       "tsys=",
                                       "x_axis=",
                                       "b_order=",
                                       "b_regions=",
                                       "l_regions=",
                                       "slice=",
                                       "output="
                                      ]
            )
        except:
            result = -1
            print('error: argument error')
            return result
        
        for opt, arg in opts:
            if opt == '-h':
                self.print_help()
                return -1 
            elif opt in ("-c","--config"):
                ss = self.read_config_file(arg)
                self.parse_options(ss,self.program)
            elif opt in ("-p","--path"):
                self.data_path = arg
            elif opt in ("-O","--obsnum"):
                self.obsnum = eval(arg)
            elif opt in ("--obs_list"):
                self.obs_list = eval(arg)
            elif opt in ("--bank"):
                self.bank = eval(arg)
            elif opt in ("--pix_list"):
                self.pix_list = eval(arg)
            elif opt in ("--eliminate_list"):
                self.eliminate_list = eval(arg)
            elif opt in ("--use_cal"):
                self.decode_cal_string(arg)
            elif opt in ("--tsys"):
                self.tsys = eval(arg)
                
            elif opt in ("--x_axis"):
                self.x_axis = arg
            elif opt in ("--b_order"):
                self.b_order = eval(arg)
            elif opt in ("--b_regions"):
                self.b_regions = eval(arg)
            elif opt in ("--l_regions"):
                self.l_regions = eval(arg)
            elif opt in ("--slice"):
                self.slice = eval(arg)
            elif opt in ("--offx"):
                self.offx = eval(arg)
            elif opt in ("--offy"):
                self.offy = eval(arg)

            elif opt in ("-o", "--output"):
                self.output_file_name = arg
                
        if print_options:
            self.print_options()
        
        result = 0
        return result

    def decode_cal_string(self,s):
        """
        Sets the use_cal option according to input string.
        Args:
            s (str): string to set use_cal option
        Returns:
            none
        """
        if s[0] == 't' or s[0] == 'T':
             self.use_cal = True
        elif s[0] == 'f' or s[0] == 'F':
            self.use_cal = False
        else:
            self.use_cal = False
            print('unknown option for cal')
            
    def print_options(self):
        """
        Prints the current options.
        Args:
            none
        Returns:
            none
        """
        print('program %s options'%(self.program))
        print('data path        = ',self.data_path)
        print('obsnum           = ',self.obsnum)
        print('obsnum list      = ',self.obs_list)
        print('bank             = ',self.bank)
        print('pixel list       = ',self.pix_list)
        print('eliminate list   = ',self.eliminate_list)
        print('use cal          = ',self.use_cal)
        print('tsys             = ',self.tsys)
        print('baseline order   = ',self.b_order)
        print('baseline regions = ',self.b_regions)
        print('line regions     = ',self.l_regions)
        print('slice            = ',self.slice)
        print('offx             = ',self.offx)
        print('offy             = ',self.offy)
        print('output file      = ',self.output_file_name)
        print(' ')
        
    def print_help(self):
        """
        Prints explanations of options.
        Args:
            none
        Returns:
            none
        """
        print('program name = %s'%(self.program))
        print('--config [-c]    : name of configuration file to set \
            parameters')
        print('--output [-o]    : name of output SpecFile')
        print('DATA SPECIFICATION')
        print('--path [-p]      : set data path')
        print('--obsnum [-O]    : set obsnum')
        print('--obs_list       : enter list of obsnums')
        print('--bank           : spectral bank for processing')
        print('--pix_list       : enter list of pixels for processing')
        print('--eliminate_list : list of channels to be eliminated')
        print('CALIBRATION')
        print('--use_cal        : use calibration scan for cal')
        print('--tsys           : value for tsys if use_cal==False')
        print('SPECTRAL LINE REDUCTION')
        print('--x_axis         : select spectral x axis')
        print('--b_order        : set polynomial baseline order')
        print('--b_regions      : enter list of lists for baseline regions')
        print('--l_regions      : enter list of lists for line fit regions')
        print('--slice          : enter list to specify slice from spectrum for processing')
        print('--offx           : pointing offset in X')
        print('--offy           : pointing offset in Y')
        print(' ')


class HandleViewSpecFileOptions():
    """
    Class to handle parameters for viewing the spec file.
    Args:
        none
    """
    def __init__(self):
        # initialization
        self.input_file_name = ''
        self.show_all_pixels = True
        self.pix_list = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        self.show_pixel = 10
        self.rms_cut = 10000.
        self.plot_range = [-1.,1.]
        
    def read_config_file(self, filename):
        """
        Returns the options from the config file.
        Args:
            filename (str): name of config file
        Returns:
            s2 (list): list of options generated from config file
        """
        fd = open(filename,'r')
        f = fd.readlines()
        s = ''
        for line in f:
            x = line.split('#')
            s = s + '--' + x[0]
        s1 = s.replace('\n',' ')
        s2 = s1.split()
        return s2
        
    def parse_options(self, s, program, arglevel=0, print_options=False):
        """
        Sets options according to input arguments and returns the exit 
            code.
        Args:
            s (str): arguments to script
            program (script): program to return
            arglevel (int): argument level (default is 0)
            print_options (bool): option to print options
        Returns:
            result (int): exit code of setting options
        """
        self.program = program
        self.arglevel = arglevel
        try:
            opts, arg = getopt.getopt(s,
                                      "hc:i:p:",
                                      ["config=",
                                       "input=",
                                       "pix_list=",
                                       "show_pixel=",
                                       "rms_cut=",
                                       "plot_range="
                                      ]
            )
        except:
            result = -1
            print('error: argument error')
            return result
        
        for opt, arg in opts:
            if opt == '-h':
                self.print_help()
                exit(-1)
            elif opt in ("-c","--config"):
                ss = self.read_config_file(arg)
                self.parse_options(ss,self.program)
            elif opt in ("-i", "--input"):
                self.input_file_name = arg
            elif opt in ("--pix_list"):
                self.pix_list = eval(arg)
            elif opt in ("-p", "--show_pixel"):
                print('single')
                self.show_pixel = eval(arg)
                self.show_all_pixels = False
            elif opt in ("--rms_cut"):
                self.rms_cut = eval(arg)
            elif opt in ("--plot_range"):
                self.plot_range = eval(arg)

                
        if print_options:
            self.print_options()
        
        result = 0
        return result

    def print_options(self):
        """
        Prints the current options.
        Args:
            none
        Returns:
            none
        """
        print('program %s options'%(self.program))
        print('input file       = ',self.input_file_name)
        if self.show_all_pixels:  
            print('show all pixels  = ',self.show_all_pixels)
            print('pixel list       = ',self.pix_list)
        else:
            print('show_all_pixels  = ',self.show_all_pixels)
            print('show pixel       = ',self.show_pixel)
        print('rms cutoff       = ',self.rms_cut)
        print('plot range       = ',self.plot_range)

    def print_help(self):
        """
        Prints explanations of options.
        Args:
            none
        Returns:
            none
        """
        print('program name = %s'%(self.program))
        print('--config [-c]    : name of configuration file to set \
            parameters')
        print('--input [-i]     : set input file name')
        print('--pix_list       : enter list of pixels to be displayed')
        print('--show_pixel     : select specific pixel for display')
        print('--rms_cut        : set rms threshold for data')
        print('--plot_range     : list to set limits for data axis in plot')

        
class HandleGridOptions():
    """
    Class to handle parameters for running the gridding script.
    Args:
        none
    """
    def __init__(self):
        # initialization
        self.program_path = './spec_driver_fits'
        self.input_file_name = ''
        self.output_file_name = ''
        self.resolution = 14.
        self.cell = 7
        self.rms_cut = 10000.
        self.x_extent = 0.
        self.y_extent = 0.
        self.otf_select=1
        self.rmax = 3.0
        self.n_samples = 256
        self.otf_a = 1.1
        self.otf_b = 4.75
        self.otf_c = 2.0
        self.otf_filter_code = ['box', 'jinc', 'gauss']

    def read_config_file(self, filename):
        """
        Returns the options from the config file.
        Args:
            filename (str): name of config file
        Returns:
            s2 (list): list of options generated from config file
        """
        fd = open(filename,'r')
        f = fd.readlines()
        s = ''
        for line in f:
            x = line.split('#')
            s = s + '--' + x[0]
        s1 = s.replace('\n',' ')
        s2 = s1.split()
        return s2
        
    def parse_options(self, s, program, arglevel=0, print_options=False):
        """
        Sets options according to input arguments and returns the exit 
            code.
        Args:
            s (str): arguments to script
            program (script): program to return
            arglevel (int): argument level (default is 0)
            print_options (bool): option to print options
        Returns:
            result (int): exit code of setting options
        """
        self.program = program
        self.arglevel = arglevel
        try:
            opts, arg = getopt.getopt(s,
                                      "hc:i:o:p:",
                                      ["config=",
                                       "program_path=",
                                       "input=",
                                       "output=",
                                       "resolution=",
                                       "cell=",
                                       "pix_list=",
                                       "rms_cut=",
                                       "x_extent=",
                                       "y_extent=",
                                       "otf_select=",
                                       "rmax=",
                                       "n_samples=",
                                       "otf_a=",
                                       "otf_b=",
                                       "otf_c="
                                      ]
            )
        except:
            result = -1
            print('error: argument error')
            return result
        
        for opt, arg in opts:
            if opt == '-h':
                self.print_help()
                exit(-1)
            elif opt in ("-c","--config"):
                ss = self.read_config_file(arg)
                self.parse_options(ss,self.program)
            elif opt in ["-p", "--program_path"]:
                self.program_path = arg
            elif opt in ("-i", "--input"):
                self.input_file_name = arg
            elif opt in ("-o", "--output"):
                self.output_file_name = arg

            elif opt in ("--resolution"):
                self.resolution = eval(arg)
            elif opt in ("--cell"):
                self.cell = eval(arg)
            elif opt in ("--pix_list"):
                self.pix_list = arg
            elif opt in ("--rms_cut"):
                self.rms_cut = eval(arg)
            elif opt in ("--x_extent"):
                self.x_extent = eval(arg)
            elif opt in ("--y_extent"):
                self.y_extent = eval(arg)
            elif opt in ("--otf_select"):
                self.otf_select = eval(arg)
            elif opt in ("--rmax"):
                self.rmax = eval(arg)
            elif opt in ("--n_samples"):
                self.n_samples = eval(arg)
            elif opt in ("--otf_a"):
                self.otf_a = eval(arg)
            elif opt in ("--otf_b"):
                self.otf_b = eval(arg)
            elif opt in ("--otf_c"):
                self.otf_c = eval(arg)
                
        if print_options:
            self.print_options()
        
        result = 0
        return result

    def print_options(self):
        """
        Prints the current options.
        Args:
            none
        Returns:
            none
        """
        print('program %s options'%(self.program))
        print('program path     = ',self.program_path)
        print('input file name  = ',self.input_file_name)
        print('output file name = ',self.output_file_name)
        print('resolution       = ',self.resolution)
        print('cell size        = ',self.cell)
        print('pix list         = ',self.pix_list)
        print('rms cutoff       = ',self.rms_cut)
        print('x extent         = ',self.x_extent)
        print('y extent         = ',self.y_extent)
        print('otf filter select= ',self.otf_filter_code[self.otf_select])
        print('r max            = ',self.rmax)
        print('n conv samples   = ',self.n_samples)
        print('otf a            = ',self.otf_a)
        print('otf b            = ',self.otf_b)
        print('otf c            = ',self.otf_c)
        print(' ')

    def print_help(self):
        """
        Prints explanations of options.
        Args:
            none
        Returns:
            none
        """
        print('program %s options'%(self.program))
        print('--config [-c]      : configuration file')
        print('--program_path [-p]: full path name to grid program')
        print('--input [-i]       : input SpecFile name')
        print('--output [-o]      : output file name')
        print('--resolution       : resolution')
        print('--cell             : cell size')
        print('--pix_list         : list of pixels to process')
        print('--rms_cut          : rms threshold')
        print('--x_extent         : x extent of cube')
        print('--y_extent         : y extent of cube')
        print('--otf_select       : filter code (0=box,1=jinc,2=gaussian)')
        print('--rmax             : maximum radius of convolution \
            (units lambda/D)')
        print('--n_samples        : number of samples in convolution filter')
        print('--otf_a            : otf a parameter')
        print('--otf_b            : otf b parameter')
        print('--otf_c            : otf_c parameter')
        print(' ')

class HandleViewCubeOptions():
    """
    Class to handle parameters for viewing the FITS cube.
    Args:
        none
    """
    def __init__(self):
        # initialization
        self.input_file_name = ''
        self.v_range = []
        self.v_scale = 1000.
        self.location = []
        self.scale = 1.0/3600.
        self.limits = []
        self.tmax_range = []
        self.tint_range = []
        self.plot_type = 'TMAX'
        self.interp = 'bilinear'
        
    def read_config_file(self, filename):
        """
        Returns the options from the config file.
        Args:
            filename (str): name of config file
        Returns:
            s2 (list): list of options generated from config file
        """
        fd = open(filename,'r')
        f = fd.readlines()
        s = ''
        for line in f:
            x = line.split('#')
            s = s + '--' + x[0]
        s1 = s.replace('\n',' ')
        s2 = s1.split()
        return s2
        
    def parse_options(self, s, program, arglevel=0, print_options=False):
        """
        Sets options according to input arguments and returns the exit 
            code.
        Args:
            s (str): arguments to script
            program (script): program to return
            arglevel (int): argument level (default is 0)
            print_options (bool): option to print options
        Returns:
            result (int): exit code of setting options
        """
        self.program = program
        self.arglevel = arglevel
        try:
            opts, arg = getopt.getopt(s,
                                      "hc:i:p:",
                                      ["config=",
                                       "input=",
                                       "v_range=",
                                       "v_scale=",
                                       "location=",
                                       "scale=",
                                       "limits=",
                                       "tmax_range=",
                                       "tint_range=",
                                       "plot_type=",
                                       "interp="
                                      ]
            )
        except:
            result = -1
            print('error: argument error')
            return result
        
        for opt, arg in opts:
            if opt == '-h':
                self.print_help()
                exit(-1)
            elif opt in ("-c","--config"):
                ss = self.read_config_file(arg)
                self.parse_options(ss,self.program)
            elif opt in ("-i", "--input"):
                self.input_file_name = arg
            elif opt in ("--v_range"):
                self.v_range = eval(arg)
            elif opt in ("--v_scale"):
                self.v_scale = eval(arg)
            elif opt in ("--location"):
                self.location = eval(arg)
            elif opt in ("--scale"):
                self.scale = eval(arg)
            elif opt in ("--limits"):
                self.limits = eval(arg)
            elif opt in ("--tmax_range"):
                self.tmax_range = eval(arg)
            elif opt in ("--tint_range"):
                self.tint_range = eval(arg)
            elif opt in ("--plot_type"):
                self.decode_plot_type(arg)
            elif opt in ("--interp"):
                self.limits = arg
                
        if print_options:
            self.print_options()
        
        result = 0
        return result

    def decode_plot_type(self,arg):
        """
        Sets the plot_type option according to input string.
        Args:
            arg (str): string to set plot_type option
        Returns:
            none
        """
        if arg == 'TINT' or arg == 'tint':
            self.plot_type = 'TINT'
        elif arg == 'TMAX' or arg == 'tmax':
            self.plot_type = 'TMAX'
        else:
            print('%s is an invalid plot type'%(arg))
            self.plot_type = 'TMAX'
            
    def print_options(self):
        """
        Prints the current options.
        Args:
            none
        Returns:
            none
        """
        print('program %s options'%(self.program))
        print('input file name  = ',self.input_file_name)
        print('velocity range   = ',self.v_range)
        print('velocity scale   = ',self.v_scale)
        print('location         = ',self.location)
        print('scale            = ',self.scale)
        print('limits           = ',self.limits)
        print('tmax plot range  = ',self.tmax_range)
        print('tint plot range  = ',self.tint_range)
        print('plot type        = ',self.plot_type)
        print('interpolation    = ',self.interp)
        print('')
        
    def print_help(self):
        """
        Prints explanations of options.
        Args:
            none
        Returns:
            none
        """
        print('program %s options'%(self.program))
        print('--config [-c]    : name of configuration file to set parameters')
        print('--input          : input FITS file')
        print('--v_range        : [vlo,vhi] is velocity range for integrated \
            intensity (km/s)')
        print('--v_scale        : scale factor for velocity [default=1000 to \
            convert m/s to km/s]')
        print('--location       : [dx,dy] is location for spectrum plot \
            (offset in arcsec)')
        print('--scale          : scale factor for position offset \
            [default=1/3600 to convert arcsec to degrees]')
        print('--limits         : [xlo,xhi,ylo,yhi] limits for final map')
        print('--tmax_range     : [data_lo,data_hi] data range for tmax \
            image')
        print('--tint_range     : [data_lo,data_hi] data range for tint \
            image')
        print('--plot_type      : data to plot - valid options: TINT, TMAX')
        print('--interpolation  : valid options: none, nearest, bilinear, \
            bicubic. default=bilinear')
        print(' ')
     
