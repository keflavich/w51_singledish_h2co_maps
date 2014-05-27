import pyspeckit
import numpy as np
import FITS_tools
from astropy.io import fits
from astropy import units as u

dpath = '/Users/adam/work/h2co/maps/W51/'

cube1 = pyspeckit.Cube(dpath+'W51_h77a_pyproc_cube_supersampled_sub.fits')
cube1.xarr.refX_units='GHz'
cube1.xarr.refX = 14.12861
cube1.xarr.convert_to_unit('km/s')
E = cube1.cube[cube1.xarr.as_unit('km/s') < 0].std(axis=0)
cube1.errorcube = np.repeat(np.reshape(E,(1,)+E.shape),cube1.shape[0],axis=0)

parcubefilename = dpath+'H77a_fit_parameters.fits'
cube1.fiteach(guesses=[1,60,3],
              absorption=False,
              integral=False,
              fittype='gaussian',
              multicore=4,
              signal_cut=3,
              parlimited=[(True,False), (True,True), (True,False)],
              parlimits=[(0.02,0), (30,75), (2,0)],
              start_from_point=[99,115])
cube1.write_fit(parcubefilename,clobber=True)

cube1.mapplot(estimator=1, vmin=30,vmax=75)

hdr = FITS_tools.strip_headers.flatten_header(cube1.header)
hdu = fits.PrimaryHDU(data=cube1.parcube[1,:,:], header=hdr)
hdr['BUNIT'] = 'km/s'
hdu.writeto('H77a_central_velocity.fits',clobber=True)

hdu.data=cube1.parcube[2,:,:]
hdu.writeto('H77a_velocity_width.fits',clobber=True)

hdr['BUNIT'] = 'K'
hdu.data=cube1.parcube[0,:,:]
hdu.writeto('H77a_amplitude.fits',clobber=True)

hdr['BUNIT'] = 'K km/s'
hdu.data = cube1.parcube[0] * np.sqrt(2*np.pi)*cube1.parcube[2,:,:]
hdu.writeto('H77a_integral.fits',clobber=True)

# Old version:
#cube2 = pyspeckit.Cube(dpath+'W51_h110alpha_cube_supersampled_sub.fits')
# New version, with sqrt(2) better noise:
# (includes H110, H112a)
cube2 = pyspeckit.Cube(dpath+'W51_Halpha_6cm_cube_supersampled_sub.fits')

cube2.xarr.refX_units='GHz'
cube2.xarr.refX = 4.87416
cube2.xarr.convert_to_unit('km/s')
E = cube2.cube[cube2.xarr.as_unit('km/s') < 0].std(axis=0)
cube2.errorcube = np.repeat(np.reshape(E,(1,)+E.shape),cube2.shape[0],axis=0)

parcubefilename = dpath+'H110a_fit_parameters.fits'
cube2.fiteach(guesses=[1,60,3],
              absorption=False,
              integral=False,
              fittype='gaussian',
              multicore=4,
              signal_cut=3,
              parlimited=[(True,False), (True,True), (True,False)],
              parlimits=[(0.02,0), (30,75), (2,0)],
              start_from_point=[99,115])
cube2.write_fit(parcubefilename,clobber=True)

cube2.parcube[cube2.parcube == 0] = np.nan
cube2.mapplot(estimator=1, vmin=30,vmax=75)
F = cube2.mapplot.FITSFigure
F.show_colorscale()
F.recenter(49.24,-0.33,width=0.8,height=0.5)
F.set_tick_labels_xformat('dd.d')
F.set_tick_labels_yformat('dd.d')
F.add_scalebar(((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value)
F.scalebar.set_label("10 pc")
F.scalebar.set_color((0.8,0.3,0.01,0.9))
F.scalebar.set_linewidth(3)

hdr = FITS_tools.strip_headers.flatten_header(cube2.header)
hdu = fits.PrimaryHDU(data=cube2.parcube[1,:,:], header=hdr)
hdr['BUNIT'] = 'km/s'
hdu.writeto('H110a_central_velocity.fits',clobber=True)

hdu.data=cube2.parcube[2,:,:]
hdu.writeto('H110a_velocity_width.fits',clobber=True)

hdr['BUNIT'] = 'K'
hdu.data=cube2.parcube[0,:,:]
hdu.writeto('H110a_amplitude.fits',clobber=True)

hdr['BUNIT'] = 'K km/s'
hdu.data = cube2.parcube[0] * np.sqrt(2*np.pi)*cube2.parcube[2,:,:]
hdu.writeto('H110a_integral.fits',clobber=True)
