from spectral_cube import SpectralCube,BooleanArrayMask
from paths import (dpath, h2co11subfn, h2co22subfn, cont6cm, cont2cm, h213co11subfn)
from astropy import units as u
from astropy.convolution import convolve
from astropy.io import fits
import numpy as np
from build_mask import includemask,cubemask
from astropy import log

TCMB = 2.7315

h2co11 = SpectralCube.read(h2co11subfn).with_mask(cubemask)
h2co22 = SpectralCube.read(h2co22subfn).with_mask(cubemask)

h2co11_13 = SpectralCube.read(h213co11subfn)
#h2co11_13._wcs = h2co11._wcs # this is a hack; there appear to be differences that can't be there...
#h2co22_13 = SpectralCube.read(dpath('W51_H213CO22_pyproc_cube_lores_supersampled_sub.fits'))

cont11 = fits.getdata(cont6cm) + TCMB
cont22 = fits.getdata(cont2cm) + TCMB
cont11[cont11<TCMB] = TCMB
cont22[cont22<TCMB] = TCMB


for label,vrange in (("",[40,75]), ("lower",[40,62]), ("upper",[62,75])):
    log.info("Working on {0} velocity range: {1}".format(label,vrange))
    vrange = vrange*u.km/u.s
    h2co11slab = h2co11.spectral_slab(*vrange)
    h2co11_13slab = h2co11_13.spectral_slab(*vrange)
    h2co22slab = h2co22.spectral_slab(*vrange)

    # Compute noise from an appropriate subset of the data
    noisevrange = [-50,0]*u.km/u.s
    h2co11noiseslab = h2co11.spectral_slab(*noisevrange)
    h2co22noiseslab = h2co22.spectral_slab(*noisevrange)
    h2co11noise = h2co11noiseslab.apply_numpy_function(np.std,axis=0)
    h2co22noise = h2co22noiseslab.apply_numpy_function(np.std,axis=0)

    dx = np.diff(h2co11.spectral_axis)[0]
    h2co11integ = h2co11slab.sum(axis=0) * dx
    h2co22integ = h2co22slab.sum(axis=0) * dx
    h2co11_13integ = h2co11_13slab.sum(axis=0) * dx

    hdu = fits.open(cont2cm)[0]
    hdu.data = (np.isfinite(h2co11integ).astype('int'))
    hdu.writeto(dpath(label+'mask_h2co_signal.fits'), clobber=True)

    h2co11peak = h2co11slab.min(axis=0)
    h2co22peak = h2co22slab.min(axis=0)
    h2co11_13peak = h2co11_13slab.min(axis=0)

    depth11 = -np.log((h2co11peak.value+cont11) / cont11)
    depth22 = -np.log((h2co22peak.value+cont22) / cont22)

    hdu.header['BUNIT'] = ''
    hdu.data = depth11
    hdu.writeto(dpath(label+'peak_optdepth_11.fits'), clobber=True)
    hdu.data = depth22
    hdu.writeto(dpath(label+'peak_optdepth_22.fits'), clobber=True)
    hdu.data = (depth11/depth22)
    hdu.writeto(dpath(label+'peak_optdepth_ratio.fits'), clobber=True)

    hdu.data = h2co11_13integ.value
    hdu.writeto(dpath(label+"H213CO_SNmasked_integrated.fits"), clobber=True)

    # These are not particularly useful
    # Better to try and measure a T_ex per pixel
    #tex11 = 1.0
    #tex22 = 1.5
    #depth11_tex = -np.log((h2co11peak+cont11-tex11) / (cont11-tex11))
    #depth22_tex = -np.log((h2co22peak+cont22-tex22) / (cont22-tex22))

    #hdu.data = depth11_tex.value
    #hdu.writeto(dpath(label+'peak_optdepth_11_tex1.0.fits'), clobber=True)
    #hdu.data = depth22_tex.value
    #hdu.writeto(dpath(label+'peak_optdepth_22_tex1.5.fits'), clobber=True)
    #hdu.data = (depth11_tex/depth22_tex).value
    #hdu.writeto(dpath(label+'peak_optdepth_ratio_tex1.0_1.5.fits'), clobber=True)
