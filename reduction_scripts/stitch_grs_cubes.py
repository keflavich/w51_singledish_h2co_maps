from astropy.io import fits
from astropy.utils.data import download_file
import FITS_tools
import numpy as np
import os
import shutil
import paths
import warnings
from astropy import wcs

url = 'http://grunt.bu.edu/grs-stitch/source/grs-{0}-cube.fits'
files = {lon:download_file(url.format(lon), cache=True)
         if not os.path.exists(os.path.split(url.format(lon))[1])
         else os.path.abspath(os.path.split(url.format(lon))[1])
         for lon in (48,49,50)}
for lon,fn in files.items():
    outfn = os.path.split(url.format(lon))[1]
    try:
        os.link(fn, paths.dpath2(outfn))
    except OSError as ex:
        if ex.errno == 17 and ex.strerror == 'File exists':
            # All is good.  Continue.
            pass
        else:
            warnings.warn("Data is downloaded onto a different filesystem;"
                          " using symlinks instead of hardlinks")
            os.symlink(fn, paths.dpath2(outfn))


f1 = fits.open('grs-48-cube.fits')
f2 = fits.open('grs-49-cube.fits')
f3 = fits.open('grs-50-cube.fits')

lonmin = wcs.WCS(f1[0].header).sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world((f1[0].header['NAXIS1']),[0],0)[0][0]
lonmax = wcs.WCS(f3[0].header).sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([0],[0],0)[0][0]
dlon = lonmax-lonmin
npix = dlon / abs(wcs.WCS(f1[0].header).wcs.get_cdelt()[0])
refpix = 1
refval = lonmax

outheader = f3[0].header.copy()
outheader['NAXIS1'] = int(npix)
outheader['CRPIX1'] = refpix
outheader['CRVAL1'] = refval
outheader['EQUINOX'] = 2000
del outheader['CROTA1']
del outheader['CROTA2']
del outheader['CROTA3']

r1 = FITS_tools.cube_regrid.regrid_cube_hdu(f1[0], outheader, order=1)
r2 = FITS_tools.cube_regrid.regrid_cube_hdu(f2[0], outheader, order=1)
r3 = FITS_tools.cube_regrid.regrid_cube_hdu(f3[0], outheader, order=1)

stacked = fits.PrimaryHDU()
# assume ==0 means bad data
# cast out of boolean
addnan = lambda x,y: np.add(np.nan_to_num(x),np.nan_to_num(y),dtype='float')
add = lambda x,y: np.add(x,y,dtype='float')
stacked.data = reduce(addnan, (r1.data, r2.data, r3.data))/reduce(add, [(r1.data != 0), (r2.data != 0), (r3.data != 0)])
stacked.header = outheader
stacked.writeto(paths.dpath2('grs_48and50_cube.fits'), clobber=True)

outhead2 = fits.getheader(paths.dpath('W51_H2CO11_cube_supersampled.fits'))
stacked_reproj = FITS_tools.cube_regrid.regrid_cube_hdu(stacked, outhead2, order=1)
stacked_reproj.writeto(paths.dpath2('grs_48and50_cube_supersampledh2cogrid.fits'), clobber=True)
