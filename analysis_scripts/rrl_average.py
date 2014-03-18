from astropy.io import fits
import numpy as np
from gbtpy import makecube

"""
H107, 108 are affected by the missing sampler on day 1
H111 has terrible baseline issues throughout.  I don't know why.
"""

fn = 'W51_h%ialpha_cube_supersampled_sub.fits'
halpha_6cm_cube = np.mean([fits.getdata(fn % i) for i in (110,112)],axis=0)
# use 111's header as an average of 110/112
halpha_6cm_hdr = fits.getheader(fn % 111)
halpha_6cm_hdu = fits.PrimaryHDU(halpha_6cm_cube, halpha_6cm_hdr)
halpha_6cm_hdu.writeto('W51_Halpha_6cm_cube_supersampled_sub.fits')
makecube.make_flats('W51_Halpha_6cm_cube_supersampled',vrange=[45,75],noisevrange=[-15,30])
