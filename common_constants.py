from astropy import units as u
from astropy import constants
from astropy import coordinates
from astropy.io import fits
from paths import datapath,datapath_w51,figurepath
import numpy as np

def beams():
    # smoothed by 20" during mapping
    aobeam = 54*u.arcsec
    gbbeam = 54*u.arcsec
    return aobeam,gbbeam

aobeam, gbbeam = beams()
fwhm = np.sqrt(8*np.log(2))
aobeamarea = (2*np.pi*(aobeam/fwhm)**2)
gbbeamarea = (2*np.pi*(gbbeam/fwhm)**2)

h2co11freq = 4.82966 * u.GHz
h2co22freq = 14.48848 * u.GHz

TCMB = 2.7315

# etamb_111 = 0.51 # ALREADY ACCOUNTED FOR!
etamb_gbt = etamb_77 = 0.886

tau11cubefn = 'W51_H2CO11_taucube_supersampled.fits'
tau22cubefn = "W51_H2CO22_pyproc_taucube_lores_supersampled.fits"

namedict = {'tau11cube':datapath+tau11cubefn,
            'tau22cube':datapath+tau22cubefn}

def rrl(n,dn=1,amu=1.007825):    # compute Radio Recomb Line feqs in GHz
    # from Brown, Lockman & Knapp ARAA 1978 16 445
    nu = 3.289842e6*(1-5.48593e-4/amu)*(1/float(n)**2 - 1/float(n+dn)**2)
    return nu * u.GHz

zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

distance = 5.1*u.kpc

vrange1 = [40,62]*u.km/u.s
vrange2 = [62,75]*u.km/u.s

muh2 = 2.8
h2togm = constants.m_p * muh2

# The center for radial profiles
center = coordinates.Galactic(49.4904*u.deg, -0.3765*u.deg)

# CO constants
coprops12CO10 = {'Aul':7.203e-8*u.Hz, 'Be': 57.63596828e9*u.Hz, 'nu': 115.271202e9 *u.Hz, 'Ju': 1, 'h2toco':1e4}
coprops13CO10 = {'Aul':6.294e-8*u.Hz, 'Be': 55.1010138e9 *u.Hz, 'nu': 110.201370e9 *u.Hz, 'Ju': 1, 'h2toco':6e5}
coprops12CO21 = {'Aul':6.910e-7*u.Hz, 'Be': 57.63596828e9*u.Hz, 'nu': 230.5380e9   *u.Hz, 'Ju': 2, 'h2toco':1e4}
coprops13CO21 = {'Aul':6.038e-7*u.Hz, 'Be': 55.1010138e9 *u.Hz, 'nu': 220.3986765e9*u.Hz, 'Ju': 2, 'h2toco':6e5}
coprops12CO32 = {'Aul':2.497e-6*u.Hz, 'Be': 57.63596828e9*u.Hz, 'nu': 345.7959899e9*u.Hz, 'Ju': 3, 'h2toco':1e4}
coprops13CO32 = {'Aul':2.181e-6*u.Hz, 'Be': 55.1010138e9 *u.Hz, 'nu': 330.5879601e9*u.Hz, 'Ju': 3, 'h2toco':6e5}


def nj(tem,Ju,numlevs=50,Be=coprops13CO10['Be']):
    """
    Energy level population of upper state for CO
    """
    from astropy.constants import h, k_B
    J = np.arange(numlevs)
    ntotDn0 = np.sum( (2*J+1)*np.exp(-J*(J+1.0)*Be*h/(k_B*tem)) )
    if Ju == 0:
        return 1.0/ntotDn0
    expBe = np.exp(-h*Be*(Ju*(Ju+1))/ ( k_B * tem ) )
    nuDn0 = (2*Ju+1.0) * expBe
    #ntotDnu = ntotDn0 * (2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( k_B * tem ) ) 
    #ntotDnu = ntotDn0 * ((2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( k_B * tem ) ) - \
    #                    (2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( k_B * TCMB*u.K ) ) )
    return nuDn0/ntotDn0

# CO to H2 conversion factor
# Taken from co_tauthin.py (other repository)
def cotocol(tex=20,coprops=coprops13CO10):
    """
    Conversion ratio - ratio of particular CO molecule to H2
    """
    from astropy.constants import h, k_B, c
    from numpy import pi
    nu = coprops['nu']
    Ju = coprops['Ju']
    Be = coprops['Be']
    Aul = coprops['Aul']
    conversionratio = coprops['h2toco']
    if not isinstance(tex,np.ndarray):
        if hasattr(tex,'__len__') and len(tex) > 1:
            tex = np.array(tex)
        else:
            tex = np.array([tex])
    if not hasattr(tex,'unit') or tex.unit is None:
        tex = tex * u.K
    ntotDnu = 1.0 / np.array([ nj(T,Ju,Be=Be,numlevs=50) for T in tex ]) 
    expcmb = np.exp(h*nu/(k_B*TCMB*u.K))
    Nupper = ( 8 * pi * nu**2 ) / ( c**3 * Aul ) * k_B / h * (expcmb - 1) / ( expcmb - (np.exp(h*nu/(k_B*tex))) ) # - (exp(hplanck*nu/(k_B*TCMB*u.K))-1)**-1)
    return (conversionratio * ntotDnu * Nupper).to(u.cm**-2 / (u.K * u.km/u.s))


__cache__ = {}

def get_cached(name):
    if name not in __cache__:
        load_into_cache(name)
    return __cache__[name]

def load11taucube():
    return get_cached('tau11cube')

def load_into_cache(name):
    __cache__[name] = fits.open(namedict[name])
