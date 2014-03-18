from astropy import units as u
from astropy.io import fits

def beams():
    # smoothed by 20" during mapping
    aobeam = 54*u.arcsec
    gbbeam = 54*u.arcsec
    return aobeam,gbbeam

h2co11freq = 4.82966 * u.GHz
h2co22freq = 14.48848 * u.GHz

# etamb_111 = 0.51 # ALREADY ACCOUNTED FOR!
etamb_77 = 0.886

datapath = '/Users/adam/work/w51/'
datapath_cubes = '/Users/adam/work/w51/h2co_singledish/'
figpath = '/Users/adam/work/h2co/maps/paper/figures/'

tau11cubefn = 'W51_H2CO11_taucube_supersampled.fits'
tau22cubefn = "W51_H2CO22_pyproc_taucube_lores_supersampled.fits"

namedict = {'tau11cube':datapath_cubes+tau11cubefn,
            'tau22cube':datapath_cubes+tau22cubefn}

__cache__ = {}

def get_cached(name):
    if name not in __cache__:
        load_into_cache(name)
    return __cache__[name]

def load11taucube():
    return get_cached('tau11cube')

def load_into_cache(name):
    __cache__[name] = fits.open(namedict[name])
