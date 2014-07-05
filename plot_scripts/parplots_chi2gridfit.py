import os
import copy

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

import co_intmaps


def makefig(hdu, ii, name, vmax_lin=100):
    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()

    #cm = copy.copy(pl.matplotlib.cm.gray)
    #cm._segmentdata = {'blue': ((0.0, 0, 0, 0.5), (1.0, 1, 1, 0.5)),
    #                   'green': ((0.0, 0, 0, 0.5), (1.0, 1, 1, 0.5)),
    #                   'red': ((0.0, 0, 0, 0.5), (1.0, 1, 1, 0.5))}
    #F1a = FITSFigure(co_intmaps.cube13_slab3_masked_mom0.hdu, convention='calabretta',
    #                 figure=fig, im_zorder=-10, grid=False, color=False,)
    
    F1 = F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=2, vmax=6.5)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['dens'])
    
    F.grid.set_color('black')
    F.grid.set_alpha(0.5)
    F.grid.set_linestyle('dotted')
    F.grid.ax.gridlines.set_zorder(1)
    F.image.set_zorder(20)
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_log.pdf'.format(name)), dpi=72)

    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()
    
    hdu.data = 10**hdu.data / 1e4
    F2 = F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=0, vmax=vmax_lin)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['lindens'])
    F.grid.set_color('black')
    F.grid.set_linestyle('dotted')
    F.grid.set_alpha(0.5)
    F.grid.ax.gridlines.set_zorder(0)
    F.image.set_zorder(20)
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_{1}_linear.pdf'.format(name,meastype)), dpi=72)

    return F1,F2



chi2levels = {'best':0.1, 'mean':1.0}


for meastype in ('best','mean'):
    chi2cube = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_{0}chi2.fits".format(meastype)))
    okmask = BooleanArrayMask(np.isfinite(chi2cube.filled_data[:]), wcs=chi2cube.wcs)
    chi2cube = chi2cube.with_mask(okmask)

    # Try masking based on stddev: uncertainty of 1 order of magnitude isn't super interesting...
    stdcube = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_stddens.fits".format(meastype)))
    stdcube = stdcube.with_mask(okmask)

    denscube = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_{0}dens.fits".format(meastype)))
    denscube = denscube.with_mask(okmask)
    colcube  = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_{0}col.fits".format(meastype)))
    colcube = colcube.with_mask(okmask)
    goodmask = chi2cube < chi2levels[meastype]
    goodmask_std = stdcube < 0.5

    flathead = denscube.wcs.dropaxis(2).to_header()

    denscube_lin = spectral_cube.SpectralCube(10**denscube.filled_data[:], wcs=denscube.wcs, mask=okmask)
    colcube_lin = spectral_cube.SpectralCube(10**colcube.filled_data[:], wcs=denscube.wcs, mask=okmask)



    denscol = spectral_cube.SpectralCube(denscube_lin.filled_data[:] * colcube_lin.filled_data[:], wcs=denscube.wcs,
                                         mask=okmask)
    wtdmeandens = np.log10(denscol.sum(axis=0) / colcube_lin.sum(axis=0))

    hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                           header=flathead)
    hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{0}_density.fits".format(meastype)), clobber=True)

    masked_wtdmeans = np.log10(denscol.with_mask(goodmask).sum(axis=0) / colcube_lin.with_mask(goodmask).sum(axis=0))
    hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                           header=flathead)
    hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{0}_density_chi2masked.fits".format(meastype)), clobber=True)

    stdmasked_wtdmeans = np.log10(denscol.with_mask(goodmask_std).sum(axis=0) / colcube_lin.with_mask(goodmask_std).sum(axis=0))
    hdu5 = fits.PrimaryHDU(data=stdmasked_wtdmeans.value,
                           header=flathead)
    hdu5.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{0}_density_stdmasked.fits".format(meastype)), clobber=True)


    cmhot = mpl.cm.gist_stern
    #cmhot = mpl.cm.cubehelix
    zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

    labels = {'dens':'log$_{10}$(n(H$_2$) [cm$^{-3}$])',
              'lindens':'n(H$_2$) [$10^4$ cm$^{-3}]$',
              'velocity':'Velocity ($V_{LSR}$ km s$^{-1}$)',
              'ratio':r'Ratio $\tau_{obs} 1-1 / \tau_{obs} 2-2$',
              'width':'Line Width (km s$^{-1}$)',
              'column':'log$_{10}$(N(H$_2$) cm$^{-2}$)'}

    hdu3 = fits.open(paths.dpath('H2CO_ParameterFits_bestdens_max.fits'))[0]
    hdu4 = fits.open(paths.dpath('H2CO_ParameterFits_meandens_max.fits'))[0]

    makefig(hdu1, 0, name='weighted_mean_{0}fit_density'.format(meastype))
    makefig(hdu2, 1, name='masked_weighted_mean_{0}fit_density'.format(meastype))
    makefig(hdu3, 2, name='bestfit_max_density')
    makefig(hdu4, 3, name='meanmatch_max_density')
    makefig(hdu5, 4, name='stdmasked_weighted_mean_{0}fit_density'.format(meastype))

    for ii,vrange in enumerate((vrange1,vrange2)):
        wtdmeandens = np.log10(denscol.spectral_slab(*vrange).sum(axis=0)/colcube_lin.spectral_slab(*vrange).sum(axis=0))
        hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                                header=denscube.wcs.dropaxis(2).to_header())
        hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{2}_density_v{0}to{1}.fits".format(vrange[0].value, vrange[1].value, meastype)), clobber=True)

        masked_wtdmeans = np.log10(denscol.with_mask(goodmask).spectral_slab(*vrange).sum(axis=0) / colcube_lin.with_mask(goodmask).spectral_slab(*vrange).sum(axis=0))
        hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                               header=denscube.wcs.dropaxis(2).to_header())
        hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{2}_density_chi2masked_v{0}to{1}.fits".format(vrange[0].value, vrange[1].value, meastype)), clobber=True)

        stdmasked_wtdmeans = np.log10(denscol.with_mask(goodmask_std).spectral_slab(*vrange).sum(axis=0) / colcube_lin.with_mask(goodmask_std).spectral_slab(*vrange).sum(axis=0))
        hdu3 = fits.PrimaryHDU(data=stdmasked_wtdmeans.value,
                               header=flathead)
        hdu3.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{2}_density_stdmasked_v{0}to{1}.fits".format(vrange[0].value, vrange[1].value, meastype)), clobber=True)

        F1,F2 = makefig(hdu1, 5+2*ii, name='weighted_mean_{2}fit_density_v{0}to{1}'.format(vrange[0].value, vrange[1].value, meastype), vmax_lin=100)
        F1,F2 = makefig(hdu2, 6+2*ii, name='masked_weighted_mean_{2}fit_density_v{0}to{1}'.format(vrange[0].value, vrange[1].value, meastype), vmax_lin=100)
        F1,F2 = makefig(hdu3, 7+2*ii, name='stdmasked_weighted_mean_{2}fit_density_v{0}to{1}'.format(vrange[0].value, vrange[1].value, meastype), vmax_lin=100)
