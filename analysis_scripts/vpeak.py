from astropy.io import fits
import FITS_tools
from paths import datapath
import pyspeckit
import numpy as np

h2co11 = pyspeckit.Cube(datapath+'W51_H2CO11_taucube_supersampled.fits')
h2co22 = pyspeckit.Cube(datapath+'W51_H2CO22_pyproc_taucube_lores_supersampled.fits')

h2co11.xarr.convert_to_unit('km/s')
argvmax = h2co11.cube.argmax(axis=0)
velo11 = np.array(h2co11.xarr[argvmax])
velo11[velo11==h2co11.xarr.min()] = np.nan

h2co22.xarr.convert_to_unit('km/s')
argvmax = h2co22.cube.argmax(axis=0)
velo22 = np.array(h2co22.xarr[argvmax])
velo22[velo22==h2co22.xarr.min()] = np.nan

header11 = FITS_tools.strip_headers.flatten_header(h2co11.header)
hdu_v11 = fits.PrimaryHDU(data=velo11, header=header11)
hdu_v11.writeto(datapath+'H2CO11_central_velocity.fits',clobber=True)

header22 = FITS_tools.strip_headers.flatten_header(h2co22.header)
hdu_v22 = fits.PrimaryHDU(data=velo22, header=header22)
hdu_v22.writeto(datapath+'H2CO22_central_velocity.fits',clobber=True)
