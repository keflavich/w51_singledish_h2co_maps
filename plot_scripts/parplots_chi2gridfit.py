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
goodmask = chi2cube < 0.1




denscol = spectral_cube.SpectralCube(denscube.filled_data[:] * colcube.filled_data[:], wcs=denscube.wcs,
                                     mask=okmask)
wtdmeandens = denscol.sum(axis=0) / colcube.sum(axis=0)

hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                      header=denscube.wcs.dropaxis(2).to_header())
hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density.fits"), clobber=True)

masked_wtdmeans = denscol.with_mask(goodmask).sum(axis=0) / colcube.with_mask(goodmask).sum(axis=0)
hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                      header=denscube.wcs.dropaxis(2).to_header())
hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits"), clobber=True)

cmhot = mpl.cm.cubehelix
zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

labels = {'dens':'log$_{10}$(n(H$_2$) cm$^{-3}$)',
          'lindens':'n(H$_2$) [cm$^{-3}]$',
          'velocity':'Velocity ($V_{LSR}$ km s$^{-1}$)',
          'ratio':r'Ratio $\tau_{obs} 1-1 / \tau_{obs} 2-2$',
          'width':'Line Width (km s$^{-1}$)',
          'column':'log$_{10}$(N(H$_2$) cm$^{-2}$)'}

def makefig(hdu, ii, name):
    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()
    
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=2, vmax=6.5)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['dens'])
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_log.pdf'.format(name)))

    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()
    
    hdu.data = 10**hdu.data
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=1e4, vmax=2e5)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['lindens'])
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_linear.pdf'.format(name)))


hdu3 = fits.open(paths.dpath('H2CO_ParameterFits_bestdens_max.fits'))[0]
hdu4 = fits.open(paths.dpath('H2CO_ParameterFits_meandens_max.fits'))[0]

makefig(hdu1, 0, name='weighted_mean_bestfit_density')
makefig(hdu2, 1, name='masked_weighted_mean_bestfit_density')
makefig(hdu3, 2, name='bestfit_max_density')
makefig(hdu4, 3, name='meanmatch_max_density')
