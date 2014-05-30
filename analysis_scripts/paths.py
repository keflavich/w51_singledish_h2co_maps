import os
datapath = '/Users/adam/work/h2co/maps/W51/'
datapath_w51 = '/Users/adam/work/w51/'
figurepath = '/Users/adam/work/w51_singledish_maps/figures/'
datapath_spectra = os.path.join(datapath, 'spectralfits')
regionpath = '/Users/adam/work/w51_singledish_maps/regions/'

def dpath(x, datapath=datapath):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)

def dpath2(x, datapath=datapath_w51):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)

def rpath(x, datapath=regionpath):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)
