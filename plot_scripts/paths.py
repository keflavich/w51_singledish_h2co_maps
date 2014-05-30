import os
datapath = '/Users/adam/work/h2co/maps/W51/'
datapath_w51 = '/Users/adam/work/w51/'
figurepath = '/Users/adam/work/h2co/maps/paper/figures/'
figurepath = '/Users/adam/work/w51_singledish_maps/tex/figures/'
regionpath = '/Users/adam/work/w51_singledish_maps/regions/'

def dpath(x, datapath=datapath):
    return os.path.join(datapath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def rpath(x, regionpath=regionpath):
    return os.path.join(regionpath, x)
