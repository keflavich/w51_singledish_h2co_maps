from astropy.io import fits
import FITS_tools
import numpy as np

f1 = fits.open('grs-48-cube.fits')
f2 = fits.open('grs-50-cube.fits')

outheader = f2[0].header.copy()
outheader['NAXIS1'] = 650

r1 = FITS_tools.cube_regrid.regrid_cube_hdu(f1[0], outheader, order=1)
r2 = FITS_tools.cube_regrid.regrid_cube_hdu(f2[0], outheader, order=1)

stacked = fits.PrimaryHDU()
# assume ==0 means bad data
# cast out of boolean
stacked.data = (r1.data + r2.data)/np.add((r1.data != 0), (r2.data != 0), dtype='float')
stacked.header = outheader
stacked.writeto('grs_48and50_cube.fits', clobber=True)

outhead2 = fits.getheader('/Users/adam/work/w51/h2co_singledish/W51_H2CO11_cube_supersampled.fits')
stacked_reproj = FITS_tools.cube_regrid.regrid_cube_hdu(stacked, outhead2, order=1)
stacked_reproj.writeto('grs_48and50_cube_supersampledh2cogrid.fits', clobber=True)
