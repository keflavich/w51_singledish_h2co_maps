import pylab as pl
import numpy as np
import astroML.plotting as ampl
from astropy.io import fits
from astropy import units as u

#from agpy.mad import MAD
#import itertools
import pyspeckit
#import mpl_toolkits
import aplpy
from common_constants import get_cached
from spectral_cube import SpectralCube, LazyMask
from paths import fpath,dpath,rpath,datapath_w51

pl.rcParams['font.size'] = 16

cube = dens = SpectralCube.read(dpath('W51_H2CO11_to_22_logdensity_supersampled.fits'))
header = fits.getheader(dpath('W51_H2CO11_cube_supersampled_continuum.fits'))

dens_peak = dens.max(axis=0)
data = dens.filled_data[:]
dens_mean = np.log10(np.nanmean(10**data,axis=0))
#dens_mean = dens._apply_along_axis(lambda x: 10**x

hdu_peak = fits.PrimaryHDU(data=dens_peak.value, header=header)
hdu_mean = fits.PrimaryHDU(data=dens_mean, header=header)

for fn in (1,2):
    pl.figure(fn)
    pl.clf()

F1 = aplpy.FITSFigure(hdu_peak, convention='calabretta', figure=pl.figure(1),
                      subplot=[0.1, 0.1, 0.5, 0.8])
F2 = aplpy.FITSFigure(hdu_mean, convention='calabretta', figure=pl.figure(2),
                      subplot=[0.1, 0.1, 0.5, 0.8])

F1.set_tick_labels_xformat('dd.d')
F2.set_tick_labels_xformat('dd.d')
F1.set_tick_labels_yformat('dd.d')
F2.set_tick_labels_yformat('dd.d')

nh2 = r'$\log(n(H_2))$ [cm$^{-3}$]'

for F in (F1,F2):
    F.show_colorscale()
    F.recenter(49.23,-0.28,width=1,height=0.5)
    F.add_colorbar()
    F.colorbar.set_axis_label_text(nh2)
    F.colorbar.set_axis_label_rotation(270)
    F.colorbar.set_axis_label_pad(25)
    F.recenter(49.23,-0.28,width=1,height=0.5)

ymin1,ymax1 = F1._ax1.bbox._bbox._points.T[1]
ymin2,ymax2 = F2._ax1.bbox._bbox._points.T[1]

pl.figure(1)
ax1 = pl.axes([0.68,ymin1,0.25,ymax1-ymin1])
counts,bins,patches = ampl.hist(dens_peak[dens_peak==dens_peak], bins=100,
                                log=True, histtype='step', linewidth=2,
                                alpha=0.8, color='k', ax=ax1)
sp = pyspeckit.Spectrum(xarr=(bins[1:]+bins[:-1])/2.,data=counts)
sp.specfit(guesses=[200,3,1])
p,m,w = sp.specfit.parinfo.values
g = p*np.exp(-(bins-m)**2/(2*(w)**2))
pl.plot(bins, g,'r--', label=r'$\mu=%0.1f, \sigma=%0.1f$' % (m,w))
pl.fill_between(bins,0,g,color='r',alpha=0.1,hatch='//',facecolor='none')
pl.annotate(r'$\mu=%0.1f$' % m   ,(4.5,200),)
pl.annotate(r'$\sigma=%0.1f$' % w,(4.5,100),)
#pl.legend(loc='upper right',fontsize=18, bbox_to_anchor=(1.0,1.25))
ax1.yaxis.set_ticks_position('right')
ax1.set_xticks([2,3,4,5])
ax1.set_ylim(0.8,400)
ax1.set_xlim(1.5,6)
ax1.set_xlabel(nh2)

pl.figure(2)
#ax2 = pl.axes([0.68,0.1,0.25,0.8])
ax2 = pl.axes([0.68,ymin2,0.25,ymax2-ymin2])
counts,bins,patches = ampl.hist(dens_mean[dens_mean==dens_mean], bins=100, log=True,
                                histtype='step', linewidth=2, alpha=0.8, color='k', ax=ax2)
sp = pyspeckit.Spectrum(xarr=(bins[1:]+bins[:-1])/2.,data=counts)
sp.specfit(guesses=[200,3,1])
p,m,w = sp.specfit.parinfo.values
g = p*np.exp(-(bins-m)**2/(2*(w)**2))
pl.plot(bins, g,'r--', label=r'$\mu=%0.1f, \sigma=%0.1f$' % (m,w))
pl.fill_between(bins,0,g,color='r',alpha=0.1,hatch='//',facecolor='none')
pl.annotate(r'$\mu=%0.1f$' % m   ,(4.5,200),)
pl.annotate(r'$\sigma=%0.1f$' % w,(4.5,100),)
#pl.legend(loc='upper right',fontsize=18, bbox_to_anchor=(1.0,1.25))
ax2.yaxis.set_ticks_position('right')
ax2.set_xticks([2,3,4,5])
ax2.set_ylim(0.8,400)
ax2.set_xlim(1.5,6)
ax2.set_xlabel(nh2)

pl.figure(1)
pl.savefig(fpath('density_peak_projection_withhist.pdf'),bbox_inches='tight')
pl.figure(2)
pl.savefig(fpath('density_mean_projection_withhist.pdf'),bbox_inches='tight')

dens = dens.filled_data[:]
dens[np.isnan(dens)] = -np.inf
peakpos = np.argmax(dens,axis=0)
peakvel = cube.spectral_axis[peakpos].to(u.km/u.s).value # cube.xarr[peakpos].copy() # to get around np read-only flags
peakvel[np.isnan(dens_peak)] = np.nan
dens[np.isinf(dens)+np.isnan(dens)] = np.inf
minvpos = np.argmin(dens,axis=0)
minvel = cube.spectral_axis[minvpos].copy().to(u.km/u.s).value # cube.xarr[minvpos].copy() # to get around np read-only flags
minvel[np.isnan(dens_peak)] = np.nan

hdu_peakvel = fits.PrimaryHDU(data=peakvel, header=header)

fn = 3
pl.figure(fn)
pl.clf()

F3 = aplpy.FITSFigure(hdu_peakvel,convention='calabretta',figure=pl.figure(fn))

F3.set_tick_labels_xformat('dd.d')
F3.set_tick_labels_yformat('dd.d')

nh2 = r'$\log(n(H_2))$ [cm$^{-3}$]'

F3.show_colorscale(vmin=45)
F3.recenter(49.23,-0.28,width=1,height=0.5)
F3.add_colorbar()
F3.colorbar.set_axis_label_text(r'$V_{LSR} [$km s$^{-1}]$')
F3.colorbar.set_axis_label_rotation(270)
F3.colorbar.set_axis_label_pad(30)
F3.recenter(49.23,-0.28,width=1,height=0.5)

F3.save(fpath('velocity_at_peak_density.pdf'))

tau11cube = get_cached('tau11cube')[0].data
tau22cube = get_cached('tau22cube')[0].data
noise11 = tau11cube[:50,:,:].std(axis=0)
noise22 = tau22cube[:50,:,:].std(axis=0)
sn11 = tau11cube/noise11
sn22 = tau11cube/noise22
mask_1 = (sn11 > 1) & (sn22 > 1)
mask_3 = (sn11 > 3) & (sn22 > 3)
mask_5 = (sn11 > 5) & (sn22 > 5)

dens.value[np.isinf(dens)] = np.nan
dens = dens.value

def makepv(cube):
    pv = np.nanmax(cube,axis=1)
    # mask out NaNs and 0's, 0's assumed to be non-signal
    pv = np.ma.masked_where(np.isnan(pv)+(pv==0),pv)
    return pv

pvdens = makepv(dens)
pvdens3 = makepv(dens * mask_3) # 0's ought to be not-max
pvdens5 = makepv(dens * mask_5) # 0's ought to be not-max
pv11 = makepv(tau11cube*mask_1)
pv22 = makepv(tau22cube*mask_1)
pvmin = np.nanmin(dens,axis=1)
pvmin = np.ma.masked_where(np.isnan(pvmin),pvmin)
#yax = cube.xarr
yax = cube.spectral_axis.to(u.km/u.s).value
xax = (np.arange(pvdens.shape[1])+1-header['CRPIX1'])*header['CD1_1']+header['CRVAL1']

fn = 4
pl.figure(fn)
for pv,qty in zip((pvdens,pv11,pv22,pvdens3,pvdens5),
                  ('density','tauoneone','tautwotwo',
                   'density3sigma','density5sigma')):
    pl.clf()
    ax = pl.gca()
    mappable = ax.pcolormesh(xax,yax,pv)
    ax.axis([49.7,48.8,40,80])
    cb = pl.colorbar(mappable)
    if 'dens' in qty:
        cb.set_label(nh2,rotation=270)
    else:
        cb.set_label("$\\tau$",rotation=270)
    cb.ax.yaxis.labelpad=30
    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel(r'$V_{LSR} [$km s$^{-1}]$')

    pl.savefig(fpath('pvdiagram_%s_max_along_latitude.pdf' % qty))

fn = 5
pl.figure(fn)
pl.clf()
ax = pl.gca()
mappable = ax.pcolormesh(xax,yax,pvmin)
ax.axis([49.7,48.8,40,80])
cb = pl.colorbar(mappable)
cb.set_label(nh2,rotation=270)
cb.ax.yaxis.labelpad=30
ax.set_xlabel("Galactic Longitude")
ax.set_ylabel(r'$V_{LSR} [$km s$^{-1}]$')


# slices
zax = (np.arange(dens.shape[1])+1-header['CRPIX2'])*header['CD2_2']+header['CRVAL2']
slice_pix = 70

for slice_pix in [50,60,70,80,90,100]:
    pv_slice = dens[:,slice_pix,:]
    pv_slice = np.ma.masked_where(np.isnan(pv_slice),pv_slice)
    fn = 6
    pl.figure(fn)
    pl.clf()
    ax = pl.gca()
    mappable = ax.pcolormesh(xax,yax,pv_slice)
    ax.axis([49.7,48.8,40,80])
    cb = pl.colorbar(mappable)
    cb.set_label(nh2,rotation=270)
    cb.ax.yaxis.labelpad=30
    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel(r'$V_{LSR} [$km s$^{-1}]$')
    ax.set_title(r"Slice $\ell=%0.2f$" % (zax[slice_pix]))

    pl.savefig(fpath('pvdiagram_slice_ell=%0.2f.pdf' % zax[slice_pix]),
               bbox_inches='tight')



oneone = SpectralCube.read(dpath('/W51_H2CO11_taucube_supersampled.fits'))
twotwo = SpectralCube.read(dpath('/W51_H2CO22_pyproc_taucube_lores_supersampled.fits'))
thirteen = SpectralCube.read(dpath('grs_48and50_cube_supersampledh2cogrid.fits',datapath_w51))
thirteen._mask = LazyMask(lambda x: 1, thirteen)
#sl13 = thirteen.spectral_slab(40*u.km/u.s, 80*u.km/u.s)
#sl13mom0 = sl13.moment0()
mask = (thirteen > 0.1) | (LazyMask(np.isnan, thirteen))
comasked = oneone.with_mask(mask)
slab_masked = comasked.spectral_slab(40*u.km/u.s, 80*u.km/u.s)
#spectral_cube.wcs_utils.check_equality(mask._wcs, comasked.wcs, verbose=True)
mom0 = slab_masked.moment0()
fig2 = pl.figure(2)
fig2.clf()
F = aplpy.FITSFigure(mom0.hdu, figure=fig2)
F.show_grayscale()
F.recenter(49.182984, -0.33612413, width=0.7, height=0.45)
F.show_regions(rpath('cyan_segments.reg'))


pl.draw()
pl.show()
