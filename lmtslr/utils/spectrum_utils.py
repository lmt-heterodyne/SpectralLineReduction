import numpy

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return numpy.isnan(y), lambda z: z.nonzero()[0]

def interpolate_spectrum_with_nans(spec):
    """
    Given a 1-d spectrum, finds all NaNs and 
    replaces with interpolated values from adjacent channels
    """
    nans, x= nan_helper(spec)
    spec[nans]= numpy.interp(x(nans), x(~nans), spec[~nans])
    return spec
