import os
import fnmatch
import glob

roach_pixels_all = [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]]

def lookup_roach_files(obsnum,
                       roach_list=['roach0', 'roach1', 'roach2', 'roach3'],
                       path='/data_lmt/spectrometer/',
                       debug=False):
    """
    Returns a tuple of the roach files which match a particular obsnum 
    and the number of those files.
    Args:
        obsnum (int): target obvservation number
        roach_list (list): list of the directories of roach files
        path (str): path to the roach directories
        debug (boolean): if debug True, tends to print out more information
    Returns:
        (filenames (list), result (int)) : list of file names, number 
        of files found
    """
    nroach = len(roach_list)
    filenames = []
    result = 0
    for roach in roach_list:
        spec_filenames = glob.glob(os.path.join(path, roach, 
                                   '%s_%d_*.nc' % (roach, obsnum)))
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

def create_roach_list(pixel_list):
    """
    Returns list of roach boards to be read given a list of pixels.
    Args:
        pixel_list (list): list of target pixels
    Returns:
        roach_list (list): list of roach boards to be read
    """
    rid = [0, 0, 0, 0]
    for pixel_id in pixel_list:
        r = find_roach_from_pixel(pixel_id)
        if r != []:
            rid[r[0]] = 1
    roach_list = []
    for i in range(4):
        if rid[i] == 1:
            roach_list.append('roach%d' % (i))
    return roach_list
