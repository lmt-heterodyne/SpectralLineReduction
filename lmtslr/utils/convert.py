import sys

# convert strings to lists of ascii, int or float

# todo if nmax is given, a size check should be done.

def lista(sl, nmax=None):
    if nmax==1: return(sl)
    return sl.split(',')

def listi(sl, nmax=None):
    if nmax==1: return(int(sl))    
    return [int(s) for s in sl.split(',')]

def listf(sl, nmax=None):
    if nmax==1: return(float(sl))
    return [float(s) for s in sl.split(',')]

