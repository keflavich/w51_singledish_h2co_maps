import os

import pylab as pl
import numpy as np
import matplotlib as mpl

import aplpy
from astropy.io import fits
from astroquery.vizier import Vizier
from astropy import coordinates
from paths import datapath, datapath_w51, figurepath

# Shortcut functions to get the full paths to files
def p1(x):
    return os.path.join(datapath, x)
def p2(x):
    return os.path.join(datapath_w51, x)
def p3(x):
    return os.path.join(figurepath, x)

def savefig(savename, **kwargs):
    """
    Save both a png and a pdf version of the image
    """
    pl.savefig(p3(savename.replace("pdf","png")), **kwargs)
    pl.savefig(p3(savename.replace("png","pdf")), **kwargs)

def setprops(F):
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')

zoomargs = {'x': 49.31, 'y':-0.35, 'width': 0.55, 'height':0.34}
# copied from projection_figures
zoomargs = {'x': 49.23, 'y': -0.28, 'width':1, 'height':0.5}

for fignum in (1,2,3):
    if pl.fignum_exists(fignum):
        pl.close(fignum)

Vizier.ROW_LIMIT = 999999
ysos = Vizier.get_catalogs('J/ApJ/706/83/ysos')['J/ApJ/706/83/ysos']
masers = Vizier.query_region('W51', radius='2 deg',
                             catalog='J/A+A/291/261/table1')['J/A+A/291/261/table1']
masers2 = Vizier.query_region('W51', radius='2 deg',
                              catalog='J/MNRAS/418/1689/table2')['J/MNRAS/418/1689/table2']
radec = [coordinates.SkyCoord(r,d,unit=('deg','deg'),frame='icrs')
         for r,d in zip(masers2['RAJ2000'],masers2['DEJ2000'])]
glon,glat = np.array(zip(*[(rd.galactic.l.degree,
                            rd.galactic.b.degree)
                           for rd in radec]))

# Added 5/23/2014
mask = fits.getdata(p1('mask_h2co_signal.fits'))

for suffix,extrastr in ((".fits",""), ("_prefiltered.fits", "filtered")):
    parcubefile = fits.open(p1('W51_taucube_fit_parcube_try10'+suffix))
    #parcubefile[0].header = fits.getheader('W51_H2CO11_taucube_integrated.fits')
    parcubefile[0].header['CTYPE3'] = 'linear'
    parcubefile[0].data *= mask

    cmhot = mpl.cm.gist_stern_r
    cmhot.set_bad('white')
    cmhot.set_under('white')
    cmhot.set_over('black')

    cmjet = mpl.cm.jet
    cmjet.set_bad('white')
    cmjet.set_under('white')
    cmjet.set_over('black')

    topbounds = [0.1,0.50,0.8,0.45]
    bottombounds = [0.1,0.05,0.8,0.464]
    # centerx, centery, radius, width, height

    pl.figure(1,figsize=(12,12))
    pl.clf()
    dens1 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[0],figure=pl.figure(1),subplot=topbounds)
    setprops(dens1)
    dens1.show_colorscale(vmin=2,vmax=6,cmap=cmhot)
    dens1.recenter(**zoomargs)
    dens1.add_colorbar()
    dens1.colorbar._colorbar_axes.set_ylabel('log$_{10}$(n(H$_2$) cm$^{-3}$)')
    dens1.hide_xaxis_label()
    dens1.hide_xtick_labels()
    velo1 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[4],figure=pl.figure(1),subplot=bottombounds)
    setprops(velo1)
    velo1.show_colorscale(vmin=52,vmax=66,cmap=cmjet)
    velo1.recenter(**zoomargs)
    velo1.add_colorbar()
    velo1.colorbar._colorbar_axes.set_ylabel('Velocity ($V_{LSR}$ km s$^{-1}$)')
    savefig('W51_H2CO_2parfittry10_v1_densityvelocity%s.png' % extrastr,bbox_inches='tight')

    dens1.remove_colorbar()
    dens1.hide_colorscale()
    col1 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[1],figure=pl.figure(1),subplot=topbounds)
    col1.show_colorscale(vmin=11,vmax=13.5,cmap=cmhot)
    col1.recenter(**zoomargs)
    col1.add_colorbar()
    col1.colorbar._colorbar_axes.set_ylabel('log$_{10}$(N(H$_2$) cm$^{-2}$)')
    col1.hide_xaxis_label()
    col1.hide_xtick_labels()
    setprops(col1)
    savefig('W51_H2CO_2parfittry10_v1_columnvelocity%s.png' % extrastr,bbox_inches='tight')

    pl.figure(2,figsize=(12,12))
    pl.clf()
    dens2 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[8],figure=pl.figure(2),subplot=topbounds)
    setprops(dens2)
    dens2.show_colorscale(vmin=2,vmax=6,cmap=cmhot)
    dens2.recenter(**zoomargs)
    dens2.add_colorbar()
    dens2.colorbar._colorbar_axes.set_ylabel('log$_{10}$(n(H$_2$) cm$^{-3}$)')
    dens2.hide_xaxis_label()
    dens2.hide_xtick_labels()
    velo2 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[12],figure=pl.figure(2),subplot=bottombounds)
    setprops(velo2)
    velo2.show_colorscale(vmin=65,vmax=73,cmap=cmjet)
    velo2.recenter(**zoomargs)
    velo2.add_colorbar()
    velo2.colorbar._colorbar_axes.set_ylabel('Velocity ($V_{LSR}$ km s$^{-1}$)')
    savefig('W51_H2CO_2parfittry10_v2_densityvelocity%s.png' % extrastr,bbox_inches='tight')

    dens2.remove_colorbar()
    dens2.hide_colorscale()
    col2 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[9],figure=pl.figure(2),subplot=topbounds)
    col2.show_colorscale(vmin=11,vmax=13.5,cmap=cmhot)
    col2.recenter(**zoomargs)
    col2.add_colorbar()
    col2.colorbar._colorbar_axes.set_ylabel('log$_{10}$(N(H$_2$) cm$^{-2}$)')
    col2.hide_xaxis_label()
    col2.hide_xtick_labels()
    savefig('W51_H2CO_2parfittry10_v2_columnvelocity%s.png' % extrastr,bbox_inches='tight')

    pl.figure(3)
    pl.clf()
    pl.rc('font',size=24)
    dens1 = aplpy.FITSFigure(parcubefile,convention='calabretta',slices=[0],figure=pl.figure(3))
    setprops(dens1)
    dens1.show_colorscale(vmin=2,vmax=6,cmap=cmhot)
    dens1.recenter(**zoomargs)
    dens1.add_colorbar()
    dens1.colorbar._colorbar_axes.set_ylabel('log$_{10}$(n(H$_2$) cm$^{-3}$)')
    dens1.set_tick_labels_format('dd.d','d.dd')
    #dens1.hide_xaxis_label()
    #dens1.hide_xtick_labels()
    dens1.refresh()
    savefig('W51_H2CO_2parfittry10_v1_justdensity%s.png' % extrastr,bbox_inches='tight')
    dens1.show_markers(np.array(ysos['_Glon'][ysos['Cl1']=='I']),
                       np.array(ysos['_Glat'][ysos['Cl1']=='I']),
                       edgecolor='k',marker='x')
    dens1.show_markers(np.array(masers['GLON']),
                       np.array(masers['GLAT']),
                       marker='+',edgecolor='k')
    dens1.show_markers(glon,glat,marker='+',edgecolor='k')
    dens1.refresh()
    savefig('W51_H2CO_2parfittry10_v1_justdensity%s_withYSOs.png' % extrastr,bbox_inches='tight')
