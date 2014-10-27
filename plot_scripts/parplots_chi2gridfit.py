import os
import copy

import pylab as pl
import numpy as np
import matplotlib as mpl

from aplpy_figure_maker import FITSFigure
from astropy.io import fits
from astropy import units as u
from astropy import log
from paths import datapath, datapath_w51, figurepath, h2co11taufn
from common_constants import get_cached
import paths

from common_constants import vrange1,vrange2

import spectral_cube
from spectral_cube import SpectralCube,BooleanArrayMask

import co_intmaps

mpl.rc_file('pubfiguresrc')

tau11cube = get_cached('tau11cube')[0].data
tau22cube = get_cached('tau22cube')[0].data
noise11 = tau11cube[:50,:,:].std(axis=0)
noise22 = tau22cube[:50,:,:].std(axis=0)
sn11 = tau11cube/noise11
sn22 = tau11cube/noise22
sncube = SpectralCube.read(h2co11taufn)
sncube._data = sn11

def makefig(hdu, ii, name, vmin_lin=0, vmax_lin=100, vmin_log=2, vmax_log=6.5,
            labeltype='dens', linscale=1e4, weight_hdu=None):
    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()

    #cm = copy.copy(pl.matplotlib.cm.gray)
    #cm._segmentdata = {'blue': ((0.0, 0, 0, 0.5), (1.0, 1, 1, 0.5)),
    #                   'green': ((0.0, 0, 0, 0.5), (1.0, 1, 1, 0.5)),
    #                   'red': ((0.0, 0, 0, 0.5), (1.0, 1, 1, 0.5))}
    #F1a = FITSFigure(co_intmaps.cube13_slab3_masked_mom0.hdu, convention='calabretta',
    #                 figure=fig, im_zorder=-10, grid=False, color=False,)

    cmhot.set_bad((0.75,0.75,0.75,1))
    
    F1 = F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=vmin_log, vmax=vmax_log)
    F.recenter(**zoomargs)
    F.colorbar.set_axis_label_text(labels[labeltype])
    F.colorbar.set_axis_label_rotation(270)
    F.colorbar.set_axis_label_pad(30)

    if weight_hdu:
        color = (0.75,)*3 # should be same as background #888
        colors = [color + (1,)] + [color + (alpha,) for alpha in np.exp(-(np.logspace(0.20,2)-1.7)**2/(2.5**2*2.))]
        levels = [0]+np.logspace(0.20,2).tolist()
        # For more efficient plotting, leave out totally transparent layers
        colors,levels = zip(*[(c,l) for c,l in zip(colors,levels) if c[3]>1e-5])
        F.show_contour(weight_hdu, levels=levels,
                       colors=colors, #smooth=3,
                       filled=True,
                       #linewidths=[1.0]*5,
                       zorder=30, convention='calabretta')
            
    
    F.grid.set_color('black')
    F.grid.set_alpha(0.5)
    F.grid.set_linestyle('dotted')
    F.grid.ax.gridlines.set_zorder(1)
    F.image.set_zorder(20)
    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_log.pdf'.format(name)), dpi=72)

    fig = pl.figure(6+ii,figsize=(12,12))
    fig.clf()
    
    # Rescale density so that we have 1..100 10^4
    hdu.data = 10**hdu.data / linscale
    F2 = F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=vmin_lin, vmax=vmax_lin)
    F.recenter(**zoomargs)
    F.colorbar.set_axis_label_text(labels['lin'+labeltype])
    F.colorbar.set_axis_label_rotation(270)
    F.colorbar.set_axis_label_pad(30)
    F.grid.set_color('black')
    F.grid.set_linestyle('dotted')
    F.grid.set_alpha(0.5)
    F.grid.ax.gridlines.set_zorder(0)
    F.image.set_zorder(20)

    if weight_hdu:
        F.show_contour(weight_hdu, levels=levels,
                       colors=colors,
                       filled=True,
                       #linewidths=[1.0]*5,
                       zorder=31, convention='calabretta')
            

    F.save(paths.fpath('H2CO_ParameterFitPlot_{0}_{1}_linear.pdf'.format(name,meastype)), dpi=72)

    return F1,F2

def plot_denscol(denscube, colcube, center, radius, alpha=0.3, vmin=43, vmax=76,
                 **kwargs):
    region_mask = BooleanArrayMask((((denscube.spatial_coordinate_map[0] - center[1])**2 + 
                                     (denscube.spatial_coordinate_map[1] - center[0])**2)**0.5 < radius)[None,:,:]
                                   * np.ones(denscube.spectral_axis.shape,dtype='bool')[:,None,None],
                                   wcs=denscube.wcs)
    whmask = np.where(denscube.with_mask(region_mask).get_mask_array())[0]
    velocities = denscube.with_mask(region_mask).spectral_axis[whmask].to(u.km/u.s)
    pl.clf()
    pl.scatter(colcube.with_mask(region_mask).flattened().value,
               denscube.with_mask(region_mask).flattened().value,
               c=velocities.value, edgecolor='none', alpha=alpha,
               vmin=vmin, vmax=vmax, **kwargs)

    pl.xlabel("log $N$(o-H$_2$CO)/(km s$^{-1}$ pc $^{-1}$)")
    pl.ylabel("log $n$(H$_2$)")
    cb = pl.colorbar()
    cb.set_label('Velocity (km s$^{-1}$)')


chi2levels = {'best':0.1, 'mean':1.0}

for meastype in ('likewtd','mean','best',):

    # Try masking based on stddev: uncertainty of 1 order of magnitude isn't super interesting...
    # (this is done first because chi2 doesn't exist for likewtd)
    stdcube = SpectralCube.read(paths.dpath("H2CO_ParameterFits_stddens.fits".format(meastype)))

    if os.path.exists(paths.dpath("H2CO_ParameterFits_{0}chi2.fits".format(meastype))):
        chi2cube = SpectralCube.read(paths.dpath("H2CO_ParameterFits_{0}chi2.fits".format(meastype)))
        okmask = BooleanArrayMask(np.isfinite(chi2cube.filled_data[:]), wcs=chi2cube.wcs)
        chi2cube = chi2cube.with_mask(okmask)
        goodmask = chi2cube < chi2levels[meastype]
    else:
        goodmask = okmask = stdcube.mask

    stdcube = stdcube.with_mask(okmask)

    denscube = SpectralCube.read(paths.dpath("H2CO_ParameterFits_{0}dens.fits".format(meastype)))
    denscube = denscube.with_mask(okmask)
    colcube  = SpectralCube.read(paths.dpath("H2CO_ParameterFits_{0}col.fits".format(meastype)))
    colcube = colcube.with_mask(okmask)
    goodmask_std = stdcube < 0.5

    flathead = denscube.wcs.dropaxis(2).to_header()

    denscube_lin = SpectralCube(10**denscube.filled_data[:].value,
                                              wcs=denscube.wcs,
                                              mask=okmask)
    colcube_lin = SpectralCube(10**colcube.filled_data[:].value,
                                             wcs=denscube.wcs, mask=okmask)


    totalcol = colcube_lin.sum(axis=0)
    totalcolgood = colcube_lin.with_mask(goodmask).sum(axis=0)
    totalcolgoodstd = colcube_lin.with_mask(goodmask_std).sum(axis=0)

    denscol = SpectralCube(denscube_lin.filled_data[:] * colcube_lin.filled_data[:], wcs=denscube.wcs,
                                         mask=okmask)
    wtdmeandens = np.log10(denscol.sum(axis=0) / totalcol)

    mindens_std = denscube.with_mask(goodmask_std).min(axis=0)
    mindens_chi2 = denscube.with_mask(goodmask).min(axis=0)

    hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                           header=flathead)
    hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{0}_density.fits".format(meastype)), clobber=True)

    masked_wtdmeans = np.log10(denscol.with_mask(goodmask).sum(axis=0) / totalcolgood)
    hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                           header=flathead)
    hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{0}_density_chi2masked.fits".format(meastype)), clobber=True)

    stdmasked_wtdmeans = np.log10(denscol.with_mask(goodmask_std).sum(axis=0) / totalcolgoodstd)
    hdu5 = fits.PrimaryHDU(data=stdmasked_wtdmeans.value,
                           header=flathead)
    hdu5.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{0}_density_stdmasked.fits".format(meastype)), clobber=True)

    hdu6 = fits.PrimaryHDU(data=mindens_std.value,
                           header=flathead)
    hdu6.writeto(paths.dpath("H2CO_ParameterFits_min_{0}_density_stdmasked.fits".format(meastype)), clobber=True)

    abundcube = (denscube.filled_data[:]-colcube.filled_data[:]).value
    abundimg = np.nanmean(abundcube,axis=0)
    hdu7 = fits.PrimaryHDU(data=abundimg,
                           header=flathead)
    hdu7.writeto(paths.dpath("H2CO_ParameterFits_mean_{0}_abundance.fits".format(meastype)), clobber=True)

    hdu8 = fits.PrimaryHDU(data=np.log10(totalcolgood.value),
                          header=flathead)
    hdu8.writeto(paths.dpath("H2CO_ParameterFits_total_{0}_column.fits".format(meastype)), clobber=True)

    fig20 = pl.figure(num=20)
    fig20.clf()
    ax20 = fig20.gca()
    abund = abundimg[np.isfinite(abundimg)]
    h,l,p = ax20.hist(abund, bins=200, alpha=0.5)
    ax20.set_xlabel("log $X($o-H$_2$CO$)$", labelpad=20)
    ax20.set_ylabel("N(pixels)")
    import pyspeckit
    # Do a very bad thing: fit to binned data
    sp = pyspeckit.Spectrum(xarr=l[:-1]+(l[1]-l[0])/2, data=h)
    sp.specfit(fittype='gaussian', multifit=True, guesses=[30, -9, 1, 150,-8,0.5])
    totamp = (sp.specfit.parinfo.values[0] + sp.specfit.parinfo.values[3])
    w1 = (sp.specfit.parinfo.values[0])/totamp
    w2 = (sp.specfit.parinfo.values[3])/totamp
    m1 = (sp.specfit.parinfo.values[1])
    m2 = (sp.specfit.parinfo.values[4])
    c1 = (sp.specfit.parinfo.values[2])**2
    c2 = (sp.specfit.parinfo.values[5])**2
    #from sklearn import mixture
    # apparently sklearn can't do this?
    #fit = mixture.GMM(n_components=2, covariance_type='full', n_iter=10000).fit(abund)
    #m1, m2 = fit.means_[:,0]
    #w1, w2 = fit.weights_
    #c1, c2 = fit.covars_[:,0,0]
    x = np.linspace(abund.min(), abund.max(), 1000)
    ax20.plot(x, totamp*(w1*np.exp(-(x-m1)**2/(2*c1)) + (w2*np.exp(-(x-m2)**2/(2*c2)))), linewidth=2, color='r', alpha=0.75)
    ax20.plot(x, totamp*(w1*np.exp(-(x-m1)**2/(2*c1))), linewidth=2, color='b', alpha=0.75)
    ax20.plot(x, totamp*(w2*np.exp(-(x-m2)**2/(2*c2))), linewidth=2, color='b', alpha=0.75)
    log.info("Method: {0}".format(meastype))
    log.info("Fitted 2 gaussians: \n{0}\nWeights: {1},{2}".format(sp.specfit.parinfo,w1,w2))
    
    fig20.savefig(paths.fpath("H2CO_ParameterFits_mean_{0}_abundance.pdf".format(meastype)), bbox_inches='tight')

    pl.figure(0)
    plot_denscol(denscube.with_mask(goodmask_std).spectral_slab(43*u.km/u.s,62*u.km/u.s),
                 colcube.with_mask(goodmask_std).spectral_slab(43*u.km/u.s,62*u.km/u.s),
                 center=(49.5,-0.4)*u.deg, radius=6*u.arcmin, vmax=62,
                 s=20 if meastype=='mean' else 50)
    pl.savefig(paths.fpath('H2CO_ParameterFits_ColumnVsDensity_{0}_43to62kms_W51Main.pdf'.format(meastype)))
    plot_denscol(denscube.with_mask(goodmask_std).spectral_slab(62*u.km/u.s,76*u.km/u.s),
                 colcube.with_mask(goodmask_std).spectral_slab(62*u.km/u.s,76*u.km/u.s),
                 center=(49.5,-0.4)*u.deg, radius=6*u.arcmin, vmax=76, vmin=62,
                 s=20 if meastype=='mean' else 50)
    pl.savefig(paths.fpath('H2CO_ParameterFits_ColumnVsDensity_{0}_62to76kms_W51Main.pdf'.format(meastype)))

    cmhot = mpl.cm.gist_stern
    cmhot.set_bad((0.75,0.75,0.75,1))
    #cmhot = mpl.cm.hot
    #cmhot = mpl.cm.cubehelix
    zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

    labels = {'dens':'log$_{10}$(n(H$_2$) [cm$^{-3}$])',
              'lindens':'n(H$_2$) [$10^4$ cm$^{-3}]$',
              'velocity':'Velocity ($V_{LSR}$ km s$^{-1}$)',
              'ratio':r'Ratio $\tau_{obs} 1-1 / \tau_{obs} 2-2$',
              'width':'Line Width (km s$^{-1}$)',
              'column':'log$_{10}$(N(H$_2$) cm$^{-2}$)',
              'lincolumn':'N(H$_2$) cm$^{-2}$',
              'abundance':'log$_{10}$(X(o-H$_2$CO) cm$^{-2}$)',
              'linabundance':'X(o-H$_2$CO) cm$^{-2}$',
             }

    hdu3 = fits.open(paths.dpath('H2CO_ParameterFits_bestdens_max.fits'))[0]
    hdu4 = fits.open(paths.dpath('H2CO_ParameterFits_meandens_max.fits'))[0]

    weight = sncube.max(axis=0)
    makefig(hdu1, 0, name='weighted_mean_{0}fit_density'.format(meastype), weight_hdu=weight.hdu)
    makefig(hdu2, 1, name='masked_weighted_mean_{0}fit_density'.format(meastype), weight_hdu=weight.hdu)
    makefig(hdu3, 2, name='bestfit_max_density', weight_hdu=weight.hdu)
    makefig(hdu4, 3, name='meanmatch_max_density', weight_hdu=weight.hdu)
    makefig(hdu5, 4, name='stdmasked_weighted_mean_{0}fit_density'.format(meastype), weight_hdu=weight.hdu)
    makefig(hdu6, 5, name='stdmasked_min_{0}fit_density'.format(meastype), weight_hdu=weight.hdu)
    makefig(hdu7, 6, name='mean_{0}fit_abundance'.format(meastype),
            vmin_lin=10**-14, vmax_lin=10**-6, vmin_log=-14, vmax_log=-6,
            labeltype='abundance', linscale=1, weight_hdu=weight.hdu)
    makefig(hdu8, 7, name='total_{0}fit_column'.format(meastype),
            vmin_lin=10**12, vmax_lin=10**16, vmin_log=12, vmax_log=16,
            labeltype='column', linscale=1, weight_hdu=weight.hdu)

    nextfig = 8

    for ii,vrange in enumerate((vrange1,vrange2)):
        vrange = vrange[0]*1.01, vrange[1]*0.99
        log.info("Integrating over {0}".format(vrange))

        wtdmeandens = np.log10(denscol.spectral_slab(*vrange).sum(axis=0) /
                               colcube_lin.spectral_slab(*vrange).sum(axis=0))
        hdu1 = fits.PrimaryHDU(data=wtdmeandens.value,
                                header=denscube.wcs.dropaxis(2).to_header())
        hdu1.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{2}_density_v{0}to{1}.fits".format(vrange[0].value,
                                                                                                      vrange[1].value,
                                                                                                      meastype)),
                     clobber=True)

        masked_wtdmeans = np.log10(denscol.with_mask(goodmask).spectral_slab(*vrange).sum(axis=0) /
                                   colcube_lin.with_mask(goodmask).spectral_slab(*vrange).sum(axis=0))
        hdu2 = fits.PrimaryHDU(data=masked_wtdmeans.value,
                               header=denscube.wcs.dropaxis(2).to_header())
        hdu2.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{2}_density_chi2masked_v{0}to{1}.fits".format(vrange[0].value,
                                                                                                                 vrange[1].value,
                                                                                                                 meastype)),
                     clobber=True)

        stdmasked_wtdmeans = np.log10(denscol.with_mask(goodmask_std).spectral_slab(*vrange).sum(axis=0) /
                                      colcube_lin.with_mask(goodmask_std).spectral_slab(*vrange).sum(axis=0))
        hdu3 = fits.PrimaryHDU(data=stdmasked_wtdmeans.value,
                               header=flathead)
        hdu3.writeto(paths.dpath("H2CO_ParameterFits_weighted_mean_{2}_density_stdmasked_v{0}to{1}.fits".format(vrange[0].value,
                                                                                                                vrange[1].value,
                                                                                                                meastype)),
                     clobber=True)

        weight = sncube.spectral_slab(*vrange).max(axis=0)

        F1,F2 = makefig(hdu1, nextfig+0+2*ii,
                        name='weighted_mean_{2}fit_density_v{0}to{1}'.format(vrange[0].value,
                                                                             vrange[1].value,
                                                                             meastype),
                        vmax_lin=100,
                        weight_hdu=weight.hdu)
        F1,F2 = makefig(hdu2, nextfig+1+2*ii,
                        name='masked_weighted_mean_{2}fit_density_v{0}to{1}'.format(vrange[0].value,
                                                                                    vrange[1].value,
                                                                                    meastype),
                        vmax_lin=100,
                        weight_hdu=weight.hdu)
        F1,F2 = makefig(hdu3, nextfig+2+2*ii,
                        name='stdmasked_weighted_mean_{2}fit_density_v{0}to{1}'.format(vrange[0].value,
                                                                                       vrange[1].value,
                                                                                       meastype),
                        vmax_lin=100,
                        weight_hdu=weight.hdu)

        mindens_std = denscube.with_mask(goodmask_std).spectral_slab(*vrange).min(axis=0)
        hdu4 = fits.PrimaryHDU(data=mindens_std.value,
                               header=flathead)
        hdu4.writeto(paths.dpath("H2CO_ParameterFits_min_{2}_density_stdmasked_v{0}to{1}.fits".format(vrange[0].value,
                                                                                                      vrange[1].value,
                                                                                                      meastype)),
                     clobber=True)
        F1,F2 = makefig(hdu4, nextfig+3+2*ii,
                        name='stdmasked_min_{2}fit_density_v{0}to{1}'.format(vrange[0].value,
                                                                             vrange[1].value,
                                                                             meastype),
                        vmax_lin=100, weight_hdu=weight.hdu)
