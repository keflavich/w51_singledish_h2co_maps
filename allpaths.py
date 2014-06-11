import os
root = os.path.expanduser('~/work/')
datapath_w51 = os.path.join(root, 'w51/')
datapath = os.path.join(datapath_w51, 'h2co_singledish')
source_root = os.path.join(root,'w51_singledish_maps')
analysis_figurepath = os.path.join(source_root, 'figures/')
figurepath = os.path.join(source_root, 'tex/figures/')
datapath_spectra = os.path.join(datapath, 'spectralfits')
regionpath = os.path.join(source_root, 'regions/')

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

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def afpath(x, figurepath=analysis_figurepath):
    return os.path.join(figurepath, x)
