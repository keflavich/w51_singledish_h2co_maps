from spectral_cube import SpectralCube,BooleanArrayMask
from paths import datapath
from astropy import units as u
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

h2co11noise = h2co11noiseslab.std(axis=0)
h2co22noise = h2co22noiseslab.std(axis=0)

sn11 = SpectralCube(h2co11.data/h2co11noise.data, wcs=h2co11.wcs)
sn22 = SpectralCube(h2co22.data/h2co22noise.data, wcs=h2co22.wcs)
mask = BooleanArrayMask( (sn11.data > 2) & (sn22.data > 2) )
maskslab = mask.spectral_slab(*vrange)

dx = np.diff(h2co11.spectral_axis)
h2co11integ = h2co11slab.with_mask(maskslab).sum(axis=0) * dx
h2co22integ = h2co22slab.with_mask(maskslab).sum(axis=0) * dx
