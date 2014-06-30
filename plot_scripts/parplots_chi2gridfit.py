import os

import pylab as pl
import numpy as np
import matplotlib as mpl

from aplpy_figure_maker import FITSFigure
from astropy.io import fits
from paths import datapath, datapath_w51, figurepath
import paths

from common_constants import vrange1,vrange2

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

denscube_lin = spectral_cube.SpectralCube(10**denscube.filled_data[:], wcs=denscube.wcs, mask=okmask)
colcube_lin = spectral_cube.SpectralCube(10**colcube.filled_data[:], wcs=denscube.wcs, mask=okmask)



denscol = spectral_cube.SpectralCube(denscube_lin.filled_data[:] * colcube_lin.filled_data[:], wcs=denscube.wcs,
                                     mask=okmask)
wtdmeandens = np.log10(denscol.sum(axis=0) / colcube_lin.sum(axis=0))

hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                       header=denscube.wcs.dropaxis(2).to_header())
hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density.fits"), clobber=True)

masked_wtdmeans = np.log10(denscol.with_mask(goodmask).sum(axis=0) / colcube_lin.with_mask(goodmask).sum(axis=0))
hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                       header=denscube.wcs.dropaxis(2).to_header())
hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits"), clobber=True)

cmhot = mpl.cm.gist_stern
#cmhot = mpl.cm.cubehelix
zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

labels = {'dens':'log$_{10}$(n(H$_2$) [cm$^{-3}$])',
          'lindens':'n(H$_2$) [$10^4$ cm$^{-3}]$',
          'velocity':'Velocity ($V_{LSR}$ km s$^{-1}$)',
          'ratio':r'Ratio $\tau_{obs} 1-1 / \tau_{obs} 2-2$',
          'width':'Line Width (km s$^{-1}$)',
          'column':'log$_{10}$(N(H$_2$) cm$^{-2}$)'}

def makefig(hdu, ii, name, vmax_lin=100):
    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()
    
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=2, vmax=6.5)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['dens'])
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_log.pdf'.format(name)), dpi=72)

    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()
    
    hdu.data = 10**hdu.data / 1e4
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=0, vmax=vmax_lin)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['lindens'])
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_linear.pdf'.format(name)), dpi=72)


hdu3 = fits.open(paths.dpath('H2CO_ParameterFits_bestdens_max.fits'))[0]
hdu4 = fits.open(paths.dpath('H2CO_ParameterFits_meandens_max.fits'))[0]

makefig(hdu1, 0, name='weighted_mean_bestfit_density')
makefig(hdu2, 1, name='masked_weighted_mean_bestfit_density')
makefig(hdu3, 2, name='bestfit_max_density')
makefig(hdu4, 3, name='meanmatch_max_density')

for ii,vrange in enumerate((vrange1,vrange2)):
    wtdmeandens = np.log10(denscol.spectral_slab(*vrange).sum(axis=0)/colcube_lin.spectral_slab(*vrange).sum(axis=0))
    hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                            header=denscube.wcs.dropaxis(2).to_header())
    hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density_v{0}to{1}.fits".format(vrange[0].value, vrange[1].value)), clobber=True)

    masked_wtdmeans = np.log10(denscol.with_mask(goodmask).spectral_slab(*vrange).sum(axis=0) / colcube_lin.with_mask(goodmask).spectral_slab(*vrange).sum(axis=0))
    hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                           header=denscube.wcs.dropaxis(2).to_header())
    hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked_v{0}to{1}.fits".format(vrange[0].value, vrange[1].value)), clobber=True)

    makefig(hdu1, 4+2*ii, name='weighted_mean_bestfit_density_v{0}to{1}'.format(vrange[0].value, vrange[1].value), vmax_lin=100)
    makefig(hdu2, 5+2*ii, name='masked_weighted_mean_bestfit_density_v{0}to{1}'.format(vrange[0].value, vrange[1].value), vmax_lin=100)
