import os
datapath = '/Users/adam/work/h2co/maps/W51/'
datapath_w51 = '/Users/adam/work/w51/'
figurepath = '/Users/adam/work/w51_singledish_maps/figures/'
datapath_spectra = os.path.join(datapath, 'spectralfits')

def dpath(x, datapath=datapath):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)
