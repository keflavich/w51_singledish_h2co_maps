from astropy.io import fits
import numpy as np
from sdpy import makecube
from set_headers import set_header_keywords
from paths import dpath

"""
H107, 108 are affected by the missing sampler on day 1
H111 has terrible baseline issues throughout.  I don't know why.
"""

fn = dpath('W51_h%ialpha_cube_supersampled_sub.fits')
halpha_6cm_cube = np.mean([fits.getdata(fn % i) for i in (110,112)],axis=0)
# use 111's header as an average of 110/112
halpha_6cm_hdr = set_header_keywords(fits.getheader(fn % 111), telescope='Arecibo')
halpha_6cm_hdr['RESTFREQ'] = (4.74418e9,"H110alpha rest frequency")
halpha_6cm_hdu = fits.PrimaryHDU(halpha_6cm_cube, halpha_6cm_hdr)
halpha_6cm_hdu.writeto(dpath('W51_Halpha_6cm_cube_supersampled_sub.fits'),clobber=True)
makecube.make_flats(dpath('W51_Halpha_6cm_cube_supersampled'),vrange=[35,85],noisevrange=[-15,30])

fn = dpath('W51_h%ialpha_cube_supersampled_continuum.fits')
halpha_6cm_cont = np.mean([fits.getdata(fn % i) for i in (110,112)],axis=0)
# use 111's header as an average of 110/112
halpha_6cm_hdr = set_header_keywords(fits.getheader(fn % 111), telescope='Arecibo')
halpha_6cm_hdr['RESTFREQ'] = (4.74418e9,"H110alpha rest frequency")
halpha_6cm_conthdu = fits.PrimaryHDU(halpha_6cm_cont, halpha_6cm_hdr)
halpha_6cm_conthdu.writeto(dpath('W51_Halpha_6cm_cube_supersampled_continuum.fits'),clobber=True)

fn = dpath('W51_he%ialpha_cube_supersampled_sub.fits')
healpha_6cm_cube = np.mean([fits.getdata(fn % i) for i in (110,112)],axis=0)
# use 111's header as an average of 110/112
healpha_6cm_hdr = set_header_keywords(fits.getheader(fn % 111), telescope='Arecibo')
healpha_6cm_hdr['NOTE'] = 'Average of He110 and He112alpha'
halpha_6cm_hdr['RESTFREQ'] = (4.74418e9,"He110alpha rest frequency")
healpha_6cm_hdu = fits.PrimaryHDU(healpha_6cm_cube, healpha_6cm_hdr)
healpha_6cm_hdu.writeto(dpath('W51_healpha_6cm_cube_supersampled_sub.fits'),clobber=True)
makecube.make_flats(dpath('W51_healpha_6cm_cube_supersampled'),vrange=[35,85],noisevrange=[-15,30])
