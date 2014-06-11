import os
root = os.path.expanduser('~/work/')
datapath = os.path.join(root, 'h2co/maps/W51/')
datapath_w51 = os.path.join(root, 'w51/')
analysis_figurepath = os.path.join(root, 'w51_singledish_maps/figures/')
figurepath = os.path.join(root, 'w51_singledish_maps/tex/figures/')
datapath_spectra = os.path.join(datapath, 'spectralfits')
regionpath = os.path.join(root, 'w51_singledish_maps/regions/')

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

def afpath(x, figurepath=afigurepath):
    return os.path.join(figurepath, x)
