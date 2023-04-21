import os
import fnmatch
import glob

roach_pixels_all = [[i+j*4 for i in range(4)] for j in range(4)]

def lookup_roach_files(obsnum,
                       roach_list=['roach0', 'roach1', 'roach2', 'roach3', 'roach4', 'roach5', 'roach6', 'roach7'],
                       path=None,
                       debug=False):
    """
    Returns a tuple of the roach files which match a particular obsnum 
    and the number of those files.
    Args:
        obsnum (int): target obvservation number
        roach_list (list): list of the directories of roach files
        path (str): path to the where roach sub-directories are
        debug (boolean): if debug True, tends to print out more information
    Returns:
        (filenames (list), result (int)) : list of file names, number 
        of files found
    """
    debug=True
    if path == None:
        if 'DATA_LMT' in os.environ:
            path = os.environ['DATA_LMT'] + '/spectrometer'
        else:
            path = '/data_lmt/spectrometer/'
        
    if not os.path.isdir(path):
        print("Warning: path=%s does not exist" % path)
    nroach = len(roach_list)
    filenames = []
    result = 0
    for roach in roach_list:
        globs = os.path.join(path, roach, '%s_%d_*.nc' % (roach, obsnum))
        spec_filenames = glob.glob(globs)
                                   
        for filename in spec_filenames:
            if debug:
                print('found %s' % (filename))
            if  not 'allantest' in filename:
                if debug:
                    print('append %s' % (filename))
                filenames.append(filename)
                result = result + 1
    if filenames == []:
        if debug:
            print('lookup_roach_files: no files for obsnum', obsnum)
    else:
        if len(filenames) != len(roach_list):
            print('Warning: %d/%d missing roach files for obsnum=%d in %s' %
                  (len(roach_list)-len(filenames),len(roach_list),obsnum,path))
        
    return (filenames, result)


def find_roach_from_pixel(pixel_id):
    """
    Returns roach number on which target pixel is located.
    Args:
        pixel_id (int): target pixel number
    Returns:
        i (int): roach number on which target pixel is located
    """
    for i, lis in enumerate(roach_pixels_all):
        if pixel_id in lis:
            return [i]
    return []

def create_roach_list(pixel_list, bank=0):
    """
    Returns list of roach boards to be read given a list of pixels.
    Args:
        pixel_list (list): list of target pixels
        bank (integer):    bank=0 can loop over roach0..3
                           bank=1 can loop over roach4..7
    Returns:
        roach_list (list): list of roach boards to be read
    """
    rid = [0, 0, 0, 0, 0, 0, 0, 0]
    for pixel_id in pixel_list:
        r = find_roach_from_pixel(pixel_id)
        if r != []:
            rid[r[0]+4*bank] = 1
    roach_list = []
    for i in range(len(rid)):
        if rid[i] == 1:
            roach_list.append('roach%d' % (i))
    return roach_list
