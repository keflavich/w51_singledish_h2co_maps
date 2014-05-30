from astropy import units as u
from astropy.io import fits
from paths import datapath,datapath_w51,figurepath

def beams():
    # smoothed by 20" during mapping
    aobeam = 54*u.arcsec
    gbbeam = 54*u.arcsec
    return aobeam,gbbeam

h2co11freq = 4.82966 * u.GHz
h2co22freq = 14.48848 * u.GHz

TCMB = 2.7315

# etamb_111 = 0.51 # ALREADY ACCOUNTED FOR!
etamb_77 = 0.886

tau11cubefn = 'W51_H2CO11_taucube_supersampled.fits'
tau22cubefn = "W51_H2CO22_pyproc_taucube_lores_supersampled.fits"

namedict = {'tau11cube':datapath+tau11cubefn,
            'tau22cube':datapath+tau22cubefn}

def rrl(n,dn=1,amu=1.007825):    # compute Radio Recomb Line feqs in GHz
    # from Brown, Lockman & Knapp ARAA 1978 16 445
    nu = 3.289842e6*(1-5.48593e-4/amu)*(1/float(n)**2 - 1/float(n+dn)**2)
    return nu * u.GHz

zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

__cache__ = {}

def get_cached(name):
    if name not in __cache__:
        load_into_cache(name)
    return __cache__[name]

def load11taucube():
    return get_cached('tau11cube')

def load_into_cache(name):
    __cache__[name] = fits.open(namedict[name])
