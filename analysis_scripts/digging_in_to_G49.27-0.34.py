import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy import constants
#from astroquery.splatalogue import Splatalogue
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from pyspeckit_individual_spectra import load_rrlcubes

src = Simbad.query_region(coordinates.Galactic(49.2676*u.deg, -0.3372*u.deg),
                          radius=5*u.arcsec)
viz = Vizier.query_region(coordinates.Galactic(49.2676*u.deg, -0.3372*u.deg),
                          radius=5*u.arcsec)

bgps = viz['J/ApJS/188/123/table1']
snu = bgps['S'][0] * u.Jy
smin,smaj = bgps['minAxis'][0]*u.arcsec,bgps['majAxis'][0]*u.arcsec

# copied from pick_bgps_candidates.py in BoundHII
def mass_conversion_factor(TK=20*u.K):
    return 14.30 * (np.exp(13.01*u.K/TK) - 1) * u.M_sun / u.Jy

def sigma_virial(mass, radius, aparam=1):
    """
    Return the velocity dispersion for a given mass and radius assuming virialized
    mass in solar masses
    radius in pc

    The "a" parameter is a correction factor for non-sphericity.  If a spherical core
    has a power-law density profile, a is given by:
    a = (1-k/3) / (1-2k/5) = 5/3 * (3-k)/(5-2k)
    k=0 -> a=1
    k=1 -> a=1.111
    k=2 -> a=1.666
    """
    sigma = np.sqrt(3* (5-2*aparam) / (3-aparam) * 2 * constants.G * (mass) / (radius))
    return sigma.to(u.km/u.s)

dist = 5.1*u.kpc
mass_g49pt37 = (mass_conversion_factor() * snu * (dist/u.kpc)**2).to(u.M_sun)
area_g49pt37 = (smin*smaj*2*np.pi/u.rad**2 * dist**2).to(u.pc**2)
r_eff = (area_g49pt37/np.pi)**0.5
sigma_g49pt37 = sigma_virial(mass_g49pt37, r_eff)

print "Mass and virial line width of G49.37-0.34: ",mass_g49pt37, sigma_g49pt37
print "Density: ",mass_g49pt37 / (4/3. * np.pi * r_eff**3)
print "Density: ",(mass_g49pt37 / (4/3. * np.pi * r_eff**3)/constants.m_p).to(u.cm**-3)

hcubes = load_rrlcubes()
h77a = hcubes['H77$\\alpha$']
h110a = hcubes['H110$\\alpha$']

ap = (49.268682, -0.33721476, 25.2/3600)
h77l = h77a.get_apspec(ap, coordsys='galactic', wunit='degree')
h77l.specfit()
h110l = h110a.get_apspec(ap, coordsys='galactic', wunit='degree')
h110l.specfit()
print h77l.specfit.parinfo.SHIFT0, h110l.specfit.parinfo.SHIFT0
