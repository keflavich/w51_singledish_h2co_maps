from common_constants import vrange1,vrange2
import spectral_cube
import paths

cube13 = spectral_cube.SpectralCube.read(paths.dpath2('grs_48and50_cube.fits'))

cube13_slab1 = cube13.spectral_slab(vrange1[0],vrange1[1])
cube13_slab2 = cube13.spectral_slab(vrange2[0],vrange2[1])
cube13_slab3 = cube13.spectral_slab(vrange1[0],vrange2[1])

cube13_slab1_masked_mom0 = cube13_slab1.with_mask(cube13_slab1 > 0.5).moment0()
cube13_slab2_masked_mom0 = cube13_slab2.with_mask(cube13_slab2 > 0.5).moment0()
cube13_slab3_masked_mom0 = cube13_slab3.with_mask(cube13_slab3 > 0.5).moment0()
