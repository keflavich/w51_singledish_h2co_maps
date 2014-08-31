from astropy.io import fits
import FITS_tools
from paths import datapath
import pyspeckit
import numpy as np
from paths import h2co11,h2co22

h2co11.xarr.convert_to_unit('km/s')
argvmax = h2co11.cube[h2co11.xarr > 20,:,:].argmax(axis=0)
velo11 = np.array(h2co11.xarr[h2co11.xarr > 20][argvmax])
velo11[velo11==h2co11.xarr[h2co11.xarr > 20].min()] = np.nan
velo11[(velo11 < 20) | (velo11 > 80)] = np.nan

h2co22.xarr.convert_to_unit('km/s')
argvmax = h2co22.cube[h2co22.xarr > 20,:,:].argmax(axis=0)
velo22 = np.array(h2co22.xarr[h2co22.xarr > 20][argvmax])
velo22[velo22==h2co22.xarr.min()] = np.nan
velo22[(velo22 < 20) | (velo22 > 80)] = np.nan

header11 = FITS_tools.strip_headers.flatten_header(h2co11.header)
hdu_v11 = fits.PrimaryHDU(data=velo11, header=header11)
hdu_v11.writeto(datapath+'H2CO11_central_velocity.fits',clobber=True)

header22 = FITS_tools.strip_headers.flatten_header(h2co22.header)
hdu_v22 = fits.PrimaryHDU(data=velo22, header=header22)
hdu_v22.writeto(datapath+'H2CO22_central_velocity.fits',clobber=True)
