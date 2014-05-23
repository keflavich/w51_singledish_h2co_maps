from spectral_cube import SpectralCube,BooleanArrayMask
from paths import datapath
from astropy import units as u
from astropy.convolution import convolve
from astropy.io import fits
import numpy as np
import os

h2co11 = SpectralCube.read(os.path.join(datapath,
                                        'W51_H2CO11_cube_supersampled_sub.fits'))
h2co22 = SpectralCube.read(os.path.join(datapath,
                                        'W51_H2CO22_pyproc_cube_lores_supersampled_sub.fits'))

vrange = [45,75]*u.km/u.s
h2co11slab = h2co11.spectral_slab(*vrange)
h2co22slab = h2co22.spectral_slab(*vrange)
noisevrange = [-50,0]*u.km/u.s
h2co11noiseslab = h2co11.spectral_slab(*noisevrange)
h2co22noiseslab = h2co22.spectral_slab(*noisevrange)

# TODO: replace _apply_numpy_function with either apply_numpy_function or std
h2co11noise = h2co11noiseslab._apply_numpy_function(np.std,axis=0)
h2co22noise = h2co22noiseslab._apply_numpy_function(np.std,axis=0)

# TODO: implement Cube (lazy?) Arithmetic: avoid filling data!
sn11 = SpectralCube(np.abs(h2co11.filled_data[:])/h2co11noise, wcs=h2co11.wcs)
sn22 = SpectralCube(np.abs(h2co22.filled_data[:])/h2co22noise, wcs=h2co22.wcs)
mask = ((sn11 > 2) & (sn22 > 2))
# TODO: implement MASK slabbing - will remove 2 lines here
sn11slab = sn11.spectral_slab(*vrange)
sn22slab = sn22.spectral_slab(*vrange)
maskslab = ((sn11slab > 2) & (sn22slab > 2))

# TODO: fix "velocity-based" mask wcs
maskslab._mask1._wcs.wcs.restfrq = 0
maskslab._mask2._wcs.wcs.restfrq = 0

dx = np.diff(h2co11.spectral_axis)[0]
h2co11integ = h2co11slab.with_mask(maskslab).sum(axis=0) * dx
h2co22integ = h2co22slab.with_mask(maskslab).sum(axis=0) * dx


# neighborly masking
filt = np.ones([3,3,3],dtype='bool')
filt[1,1,1] = 0
# TODO: provide access to boolean arrays!
boolean_maskslab = maskslab.include(data=sn11slab._data, wcs=maskslab._mask1._wcs)
nneighbors = convolve(boolean_maskslab, filt)
boolean_maskslab[(nneighbors<7)] = False
boolean_maskslab[(nneighbors>=10)] = True
# Grow it again
nneighbors = convolve(boolean_maskslab, filt)
boolean_maskslab[(nneighbors>=5)] = True
boolean_maskslab = BooleanArrayMask(boolean_maskslab, maskslab._mask1._wcs)

h2co11integ2 = h2co11slab.with_mask(boolean_maskslab).sum(axis=0) * dx
h2co22integ2 = h2co22slab.with_mask(boolean_maskslab).sum(axis=0) * dx

cont11 = fits.getdata(os.path.join(datapath,'W51_H2CO11_cube_supersampled_continuum.fits')) + 2.73
cont22 = fits.getdata(os.path.join(datapath,'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits')) + 2.73
cont11[cont11<2.73] = 2.73
cont22[cont22<2.73] = 2.73

hdu = fits.open(os.path.join(datapath,'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits'))[0]
hdu.data = (np.isfinite(h2co11integ2).astype('int'))
hdu.writeto(os.path.join(datapath,'mask_h2co_signal.fits'), clobber=True)
