import numpy as np
import astrodendro
import FITS_tools
from astropy import units as u
from paths import datapath
from astropy.io import fits

oneonefn = datapath+'W51_H2CO11_cube_supersampled_continuum.fits'
twotwofn = datapath+'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits'
oneoned = (fits.getdata(oneonefn) * u.K).to(u.Jy,
                                            equivalencies=u.brightness_temperature(2*np.pi*(50*u.arcsec)**2/(8*np.log(2)), 4.829*u.GHz))
twotwod = fits.getdata(twotwofn)
oneoneh = fits.getheader(oneonefn)
twotwoh = fits.getheader(twotwofn)

dendone = astrodendro.dendrogram.Dendrogram.compute(oneoned, 0.1*u.Jy, 0.05*u.Jy, 6)

metadata = {}
#metadata['data_unit'] = u.K
metadata['spatial_scale'] = FITS_tools.header_tools.header_to_platescale(oneoneh) * u.deg
metadata['beam_major'] = 50 * u.arcsec
metadata['beam_minor'] = 50 * u.arcsec

cat = astrodendro.pp_catalog(dendone, metadata)
