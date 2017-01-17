"""
Requires rgb_wcs branch of aplpy June 10, 2014
"""
from __future__ import print_function
from astropy.io import fits
import aplpy
#from aplpy_figure_maker import FITSFigure
from aplpy import FITSFigure
import PIL
import matplotlib.colors as mc
import numpy as np
import pylab as pl
from astropy import wcs
#from agpy.cubes import flatten_header
from FITS_tools.strip_headers import flatten_header
from astropy import units as u
import matplotlib
import pyregion
import re
from paths import datapath_w51 as datapath
from paths import figurepath as figpath
from paths import h2co11taufn
import paths
import os

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
pl.mpl.rc_file(os.path.join(paths.source_root, 'plot_scripts/pubfiguresrc'))


h2co11 = fits.open(h2co11taufn)
co32 = fits.open(datapath+'w51_bieging_13co32.fits')
#hdr = fits.getheader(datapath+'h2co_integ.fits')
hdr = flatten_header(h2co11[0].header)
cohdr = flatten_header(co32[0].header)

h2co_45to55 = h2co11[0].data[95:105,:,:].sum(axis=0)
h2co_55to60 = h2co11[0].data[105:110,:,:].sum(axis=0)
h2co_60to65 = h2co11[0].data[110:115,:,:].sum(axis=0)
h2co_65to75 = h2co11[0].data[115:125,:,:].sum(axis=0)

co_45to55 = co32[0].data[31:51,:,:].sum(axis=0)
co_55to60 = co32[0].data[51:61,:,:].sum(axis=0)
co_60to65 = co32[0].data[61:71,:,:].sum(axis=0)
co_65to75 = co32[0].data[71:91,:,:].sum(axis=0)

rgb = np.array([h2co_55to60,h2co_60to65,h2co_65to75])

L,H = np.nanmin(rgb),np.nanmax(rgb)
ncolors = 256
N = mc.BoundaryNorm(np.linspace(L,H,ncolors),ncolors)
rgbint = N(rgb[:,:,:]).data.swapaxes(0,1).swapaxes(1,2).astype('uint8')

alpha = ((np.max(rgb,axis=0) > 0.5) * 200).astype('uint8')
rgba = np.concatenate([rgbint,alpha[:,:,np.newaxis]],axis=2)

rgbI = PIL.Image.fromarray(rgba[::-1,:,:],mode='RGBA')
rgbI.save(figpath+'temp.png')

#hdu = fits.PrimaryHDU(data=rgb, header=h2co11[0].header)

#hdu = fits.PrimaryHDU(h2co_65to75, hdr)
WCS = wcs.WCS(hdr)

pl.close(1)
fig = pl.figure(1,figsize=(15,10))
fig.clf()
F = aplpy.FITSFigure(datapath+'v2.0_ds2_l050_13pca_map20_reproject.fits',
                     figure=fig, convention='calabretta', colorbar=False,
                     color=False)

H = fits.Header.fromtextfile(paths.pdpath('hdr4096.hdr'))
hwcs = wcs.WCS(H)
#F.show_rgb('/Volumes/128gbdisk/w51/pngs/W51_4096sq_WISE_bolo_mosaic_rotated_blackbg.png',wcs=hwcs)
F.show_rgb(paths.pdpath('W51_4096sq_WISE_bolo_mosaic_rotated_blackbg.png'),
           wcs=hwcs)

#F.show_regions(paths.rpath('large_scale_regions.reg'))
regions = pyregion.open(paths.rpath('large_scale_regions.reg'))
#text = re.compile("text={([^}]*)}")
for reg in regions:
    #t = text.search(reg.comment).groups()[0]
    t = reg.attr[1]['text']
    F.add_label(reg.coord_list[0], reg.coord_list[1], t, color='white',
                size=16, weight='bold')
F.set_tick_labels_xformat('dd.d')
F.set_tick_labels_yformat('dd.d')
F.recenter(49.27, -0.32, width=0.9, height=0.4)
F.save(paths.fpath('W51_wisecolor_largescale_labeled.pdf'), dpi=72)
F.show_rgb(paths.dpath("make_pretty_picture/W51_modified.png",paths.datapath_w51),
           wcs=hwcs)
F.save(paths.fpath('W51_wisecolor_modified_largescale_labeled.pdf'), dpi=72)

F.add_scalebar(((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value)
F.scalebar.set_label("10 pc")
F.scalebar.set_font_size(18)
F.scalebar.set_font_weight('bold')
F.scalebar.set_color('w')
F.scalebar.set_linewidth(3)
F.save(paths.fpath('W51_wisecolor_modified_largescale_labeled_scalebar.pdf'), dpi=150)
F.scalebar.hide()

for L in list(F._layers.keys()):
    if L in F._layers:
        F.remove_layer(L)

F.set_tick_labels_xformat('dd.d')
F.set_tick_labels_yformat('dd.d')
F.show_regions(paths.rpath('image_region_labels.reg'))
F.save(paths.fpath('W51_wisecolor_labeled_detail.pdf'), dpi=72)
F.recenter(49.436,-0.365,0.21)
F.save(paths.fpath('W51_wisecolor_zoom_W51Main.pdf'), dpi=72)
F.recenter(49.05,-0.33,0.20)
F.save(paths.fpath('W51_wisecolor_zoom_W51B.pdf'), dpi=72)


#F.show_grayscale()
#F.show_rgb('temp.png',wcs=WCS)
F.set_auto_refresh(False)
F.recenter(49.436,-0.365,0.2)

F.recenter(49.2743,-0.3439,width=1,height=0.45)
F.set_tick_labels_xformat('dd.d')
F.set_tick_labels_yformat('dd.d')


dens = fits.getdata(datapath+'h2co_singledish/W51_H2CO11_to_22_logdensity_supersampled.fits')
header = fits.getheader(datapath+'h2co_singledish/W51_H2CO11_cube_supersampled_continuum.fits')
dens_peak = np.nanmax(dens,axis=0)
hdu_peak = fits.PrimaryHDU(data=dens_peak, header=header)
F.show_contour(hdu_peak, levels=[2.5,3,3.5,4,4.5,5], colors=[(1,min([0.15+0.15*x,1]),0,min([0.15+0.15*x,1])) for x in range(1,6)], filled=True)

for ii,lev in enumerate([2.5,3,3.5,4,4.5]):
    F.add_label(48.85,-0.48+0.020*ii,
                "$%i\\times10^{%i}$" % (np.floor(10**(lev % 1)), np.floor(lev)), 
                color=(1,min([0.15+0.15*(ii+1),1]),0,1),
                weight='bold',
                size=18)

F.add_label(48.85,-0.48+0.020*(ii+1),
            "$n(H_2)$ cm$^{-3}$",
            color=(1,0.6,0,1),
            weight='bold',
            size=18)


scalebar_length = ((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value
try:
    F.add_scalebar(scalebar_length)
except:
    F.scalebar.show(scalebar_length)
F.scalebar.set_label("10 pc")
F.scalebar.set_font_size(18)
F.scalebar.set_font_weight('bold')
F.scalebar.set_color('w')
#F.scalebar.set_color((0.8,0.3,0.01,0.9))
#F.scalebar.set_linewidth(3)

print("Trying to save...")

F.refresh()
# PDF led to 
# python(29168,0x7fff7eb43310) malloc: *** error for object 0x117462790: pointer being freed was not allocated
# *** set a breakpoint in malloc_error_break to debug
# Abort trap: 6
#F.save(figpath+'w51_wisecolor_densityoverlay.pdf')
F.save(figpath+'w51_wisecolor_densityoverlay.png', dpi=72)
for L in list(F._layers.keys()):
    if L in F._layers:
        F.remove_layer(L)



#F.show_contour(fits.PrimaryHDU(h2co_45to55, hdr), levels=[0.15], colors=['b'])
#F.show_contour(fits.PrimaryHDU(h2co_55to60, hdr), levels=[0.40,5], filled=True, colors=[(0,0,1,0.3)])
#F.show_contour(fits.PrimaryHDU(h2co_60to65, hdr), levels=[0.75,5], filled=True, colors=[(0,1,0,0.3)])
#F.show_contour(fits.PrimaryHDU(h2co_65to75, hdr), levels=[0.50,5], filled=True, colors=[(1,0,0,0.3)])

F.show_contour(fits.PrimaryHDU(co_45to55, cohdr), levels=[15,55,85,500], filled=True, colors=[(0,0.5,0.5,0.3),(0,0.5,0.5,0.4),(0,0.5,0.5,0.5)])
F.refresh()
F.save(figpath+'w51_wisecolor_cooverlay_45to55.png')
for L in list(F._layers.keys()):
    F.remove_layer(L)

F.show_contour(fits.PrimaryHDU(co_55to60, cohdr), levels=[15,55,85,500], filled=True, colors=[(0,0,1,0.2),(0,0,1,0.3),(0,0,1,0.4)])
F.refresh()
F.save(figpath+'w51_wisecolor_cooverlay_55to60.png')
for L in list(F._layers.keys()):
    F.remove_layer(L)

F.show_contour(fits.PrimaryHDU(co_60to65, cohdr), levels=[15,55,85,500], filled=True, colors=[(0,0.5,0,0.3),(0,0.5,0,0.4),(0,0.5,0,0.6)])
F.refresh()
F.save(figpath+'w51_wisecolor_cooverlay_60to65.png')
for L in list(F._layers.keys()):
    F.remove_layer(L)

F.show_contour(fits.PrimaryHDU(co_65to75, cohdr), levels=[15,55,85,500], filled=True, colors=[(1,0,0,0.2),(1,0,0,0.3),(1,0,0,0.4)])
F.refresh()
F.save(figpath+'w51_wisecolor_cooverlay_65to75.png')
for L in list(F._layers.keys()):
    F.remove_layer(L)

F.refresh()
F.save(figpath+'w51_wisecolor_nooverlay.png')

F.show_rectangles([49.48640],[-0.37889],0.05688888888889344,0.05688888888889344)
F.refresh()
F.save(figpath+'w51_wisecolor_nooverlay_vlazoomregion.png')


for L in list(F._layers.keys()):
    F.remove_layer(L)
F.recenter(49.48,-0.38,0.05)
slength = ((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value
if slength != scalebar_length:
    try:
        F.scalebar.set_length(slength)
    except:
        F.scalebar._scalebar_settings.pop('linewidth')
        F.scalebar.set_length(slength)
        F.scalebar.set_linewidth(3)
F.save(figpath+'w51_wisecolor_zoomW51Main.png')
F.refresh()
F.show_regions(paths.rpath('w51main_spectral_apertures.reg'))
F.save(figpath+'w51_wisecolor_zoomW51Main_spectralapertures.png')
