import os

import pylab as pl
import numpy as np
import matplotlib as mpl

from aplpy_figure_maker import FITSFigure
from astropy.io import fits
from paths import datapath, datapath_w51, figurepath
import paths

import spectral_cube
from spectral_cube import SpectralCube,BooleanArrayMask

chi2cube = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_bestchi2.fits"))
okmask = BooleanArrayMask(np.isfinite(chi2cube.filled_data[:]), wcs=chi2cube.wcs)
chi2cube = chi2cube.with_mask(okmask)

denscube = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_bestdens.fits"))
denscube = denscube.with_mask(okmask)
colcube  = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_bestcol.fits"))
colcube = colcube.with_mask(okmask)
goodmask = chi2cube < 1

denscol = spectral_cube.SpectralCube(denscube.filled_data[:] + colcube.filled_data[:], wcs=denscube.wcs,
                                     mask=okmask)
wtdmeandens = denscol.sum(axis=0) / colcube.sum(axis=0)

hdu = fits.PrimaryHDU(data=wtdmeandens.value,
                      header=denscube.wcs.dropaxis(2).to_header())
hdu.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density.fits"), clobber=True)

masked_wtdmeans = denscol.with_mask(goodmask).sum(axis=0) / colcube.with_mask(goodmask).sum(axis=0)
hdu = fits.PrimaryHDU(data=masked_wtdmeandens,
                      header=flatheader)
hdu.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits"), clobber=True)

for vrange in ([52,66], [65,75]):
    fig = pl.figure(6,figsize=(12,12))
    fig.clf()
    
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=0, vmax=30)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['ratio'])
    savefig(fn.replace("fits","png"), bbox_inches='tight')


