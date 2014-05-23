import aplpy
import pyregion
import pylab as pl
import numpy as np
from astropy import units as u

datapath = '/Users/adam/work/w51/'
figpath = '/Users/adam/work/w51/figures/'

pl.figure(0)
pl.clf()

F = aplpy.FITSFigure(datapath+'v2.0_ds2_l050_13pca_map20.fits', convention='calabretta', figure=pl.figure(0))
F.show_grayscale(vmax=35,vmin=-1.0,stretch='arcsinh',invert=True)
F.recenter(49.39,-0.33,radius=15/60.)

endpoints_wcs = pyregion.open(datapath+'pvendpoints.reg')
for jj,color in enumerate(('green','red','blue','cyan','yellow','orange','purple')):
    coords = np.array([s.coord_list for s in endpoints_wcs if s.attr[1]['color'] == color])
    F.show_lines([coords.T], color=color)

F.add_scalebar(((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value)
F.scalebar.set_label("10 pc")
F.scalebar.set_color((0.8,0.3,0.01,0.9))
F.scalebar.set_linewidth(3)

F.save(figpath+'W51_PVDiagrams_PVmap.pdf')

############################
# New stuff added 12/12/2013
############################

from astropy.io import fits
from agpy.cubes import flatten_header
co32 = fits.open(datapath+'w51_bieging_13co32.fits')
cohdr = flatten_header(co32[0].header)
co_45to55 = co32[0].data[31:51,:,:].sum(axis=0)
co_55to60 = co32[0].data[51:61,:,:].sum(axis=0)
co_60to65 = co32[0].data[61:71,:,:].sum(axis=0)
co_65to75 = co32[0].data[71:91,:,:].sum(axis=0)

F.show_contour(fits.PrimaryHDU(co_45to55, cohdr), levels=[15,55,85,500], filled=True, colors=[(0,0.5,0.5,0.3),(0,0.5,0.5,0.4),(0,0.5,0.5,0.5)])
F.save(figpath+'w51_bgpsgrayscale_cooverlay_45to55.png')
for L in F._layers.keys(): 
    if 'contour' in L: F.remove_layer(L)

F.show_contour(fits.PrimaryHDU(co_55to60, cohdr), levels=[15,55,85,500], filled=True, colors=[(0,0,1,0.2),(0,0,1,0.3),(0,0,1,0.4)])
F.save(figpath+'w51_bgpsgrayscale_cooverlay_55to60.png')
for L in F._layers.keys(): 
    if 'contour' in L: F.remove_layer(L)

F.show_contour(fits.PrimaryHDU(co_60to65, cohdr), levels=[15,55,85,500], filled=True, colors=[(0,0.5,0,0.3),(0,0.5,0,0.4),(0,0.5,0,0.6)])
F.save(figpath+'w51_bgpsgrayscale_cooverlay_60to65.png')
for L in F._layers.keys(): 
    if 'contour' in L: F.remove_layer(L)

F.show_contour(fits.PrimaryHDU(co_65to75, cohdr), levels=[15,55,85,500], filled=True, colors=[(1,0,0,0.2),(1,0,0,0.3),(1,0,0,0.4)])
F.save(figpath+'w51_bgpsgrayscale_cooverlay_65to75.png')
for L in F._layers.keys(): 
    if 'contour' in L: F.remove_layer(L)

F.save(figpath+'w51_bgpsgrayscale_nooverlay.png')
