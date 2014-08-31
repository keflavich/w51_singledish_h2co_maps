"""
Compute the luminosity from HiGal SED fits
"""

from common_constants import distance
from astropy import units as u
from astropy.io import fits
from astropy import wcs
import numpy as np
import paths

def integ_to_lum(fn):

    img = fits.getdata(fn)*u.erg/u.s/u.cm**2
    w = wcs.WCS(fits.getheader(fn))
    pixscale = np.abs(w.wcs.get_cdelt()[0])*u.deg
    pixarea_cm = ((pixscale * distance)**2).to(u.cm**2,
                                               u.dimensionless_angles())

    lum = pixarea_cm * img

    return lum.to(u.L_sun) * 4 * np.pi


fn = paths.dpath2('HIGAL_W51_mosaic_fit_070to500_integral.fits')
lum = integ_to_lum(fn)
ff = fits.open(fn)
ff[0].data = lum.value
ff[0].header['BUNIT'] = 'L_sun'
ff.writeto(fn.replace("_integral","_luminosity"), clobber=True)
