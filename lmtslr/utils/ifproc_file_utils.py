import os
#import fnmatch
import glob

def lookup_ifproc_file(obsnum, path='/data_lmt/ifproc/'):
    """
    Returns the path to the NetCDF data file for a given obsnum.
    Args:
        obsnum (int): observation number of target observation
        path (str): path to the data directory (default is 
            '/data_lmt/ifproc/')
    Returns:
        filename (str): path to NetCDF data file of target observation
    """
    paths = [path]

    if 'ifproc' not in path:
        paths += ['/data_lmt/ifproc/']
    if 'lmtttpm' not in path:
        paths += ['/data_lmt/lmttpm/']
    if 'tel' not in path:
        paths += ['/data_lmt/tel/']

    print(paths)

    for path in paths:
        filenames = glob.glob(os.path.join(path, '*_%06d_*.nc' % obsnum))
        if len(filenames) > 0:
            print('found %s' % (filenames[0]))
            return filenames[0]
    return ''
    #filename = ''
    #for file in os.listdir(path):
    #    if fnmatch.fnmatch(file,'*_%06d_*.nc'%(obsnum)):
    #        print('found %s'%(file))
    #        filename = path+file
    #if filename == '':
    #print('lookup_ifproc_file: no file for obsnum ', obsnum)
    #if 'lmttpm' not in path:
    #    print('look in lmttpm')
    #    return lookup_ifproc_file(obsnum,path='/data_lmt/lmttpm/')
    #return(filename)
