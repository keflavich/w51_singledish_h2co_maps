import os

import pylab as pl
import numpy as np
import matplotlib as mpl

from aplpy_figure_maker import FITSFigure
from astropy.io import fits
from astroquery.vizier import Vizier
from astropy import coordinates
from paths import datapath, datapath_w51, figurepath
from FITS_tools.strip_headers import flatten_header

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

zoomargs = {'x': 49.31, 'y':-0.35, 'width': 0.55, 'height':0.34}
# copied from projection_figures
zoomargs = {'x': 49.23, 'y': -0.28, 'width':1, 'height':0.5}
zoomargs = dict(x=49.27, y=-0.32, width=0.9, height=0.4)

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

for suffix,extrastr in ((".fits",""), ):#("_prefiltered.fits", "filtered")):
    parcubefile = fits.open(p1('W51_taucube_fit_parcube_try10'+suffix))
    #parcubefile[0].header = fits.getheader('W51_H2CO11_taucube_integrated.fits')
    parcubefile[0].header['CTYPE3'] = 'linear'
    parcubefile[0].data *= mask

    cmhot = mpl.cm.gist_stern_r
    cmhot.set_bad('white')
    cmhot.set_under('white')
    cmhot.set_over('black')

    cmjet = mpl.cm.jet
    cmjet = mpl.cm.RdYlBu
    cmjet.set_bad('white')
    cmjet.set_under('white')
    cmjet.set_over('black')

    topbounds = [0.1,0.50,0.8,0.45]
    bottombounds = [0.1,0.05,0.8,0.464]
    # centerx, centery, radius, width, height

    pl.figure(1,figsize=(12,12))
    pl.clf()
    dens1 = FITSFigure(parcubefile,convention='calabretta',slices=[0],figure=pl.figure(1),subplot=topbounds)
    dens1.show_colorscale(vmin=2,vmax=6,cmap=cmhot)
    dens1.recenter(**zoomargs)
    dens1.colorbar._colorbar_axes.set_ylabel('log$_{10}$(n(H$_2$) cm$^{-3}$)')
    dens1.hide_xaxis_label()
    dens1.hide_xtick_labels()
    velo1 = FITSFigure(parcubefile,convention='calabretta',slices=[4],figure=pl.figure(1),subplot=bottombounds)
    velo1.show_colorscale(vmin=52,vmax=66,cmap=cmjet)
    velo1.recenter(**zoomargs)
    velo1.colorbar._colorbar_axes.set_ylabel('Velocity ($V_{LSR}$ km s$^{-1}$)')
    savefig('W51_H2CO_2parfittry10_v1_densityvelocity%s.png' % extrastr,bbox_inches='tight')

    dens1.remove_colorbar()
    dens1.hide_colorscale()
    col1 = FITSFigure(parcubefile,convention='calabretta',slices=[1],figure=pl.figure(1),subplot=topbounds)
    col1.show_colorscale(vmin=11,vmax=13.5,cmap=cmhot)
    col1.recenter(**zoomargs)
    col1.colorbar._colorbar_axes.set_ylabel('log$_{10}$(N(H$_2$) cm$^{-2}$)')
    col1.hide_xaxis_label()
    col1.hide_xtick_labels()
    savefig('W51_H2CO_2parfittry10_v1_columnvelocity%s.png' % extrastr,bbox_inches='tight')

    pl.figure(2,figsize=(12,12))
    pl.clf()
    dens2 = FITSFigure(parcubefile,convention='calabretta',slices=[8],figure=pl.figure(2),subplot=topbounds)
    dens2.show_colorscale(vmin=2,vmax=6,cmap=cmhot)
    dens2.recenter(**zoomargs)
    dens2.colorbar._colorbar_axes.set_ylabel('log$_{10}$(n(H$_2$) cm$^{-3}$)')
    dens2.hide_xaxis_label()
    dens2.hide_xtick_labels()
    velo2 = FITSFigure(parcubefile,convention='calabretta',slices=[12],figure=pl.figure(2),subplot=bottombounds)
    velo2.show_colorscale(vmin=65,vmax=73,cmap=cmjet)
    velo2.recenter(**zoomargs)
    velo2.colorbar._colorbar_axes.set_ylabel('Velocity ($V_{LSR}$ km s$^{-1}$)')
    savefig('W51_H2CO_2parfittry10_v2_densityvelocity%s.png' % extrastr,bbox_inches='tight')

    dens2.remove_colorbar()
    dens2.hide_colorscale()
    col2 = FITSFigure(parcubefile,convention='calabretta',slices=[9],figure=pl.figure(2),subplot=topbounds)
    col2.show_colorscale(vmin=11,vmax=13.5,cmap=cmhot)
    col2.recenter(**zoomargs)
    col2.colorbar._colorbar_axes.set_ylabel('log$_{10}$(N(H$_2$) cm$^{-2}$)')
    col2.hide_xaxis_label()
    col2.hide_xtick_labels()
    savefig('W51_H2CO_2parfittry10_v2_columnvelocity%s.png' % extrastr,bbox_inches='tight')

    pl.figure(3)
    pl.clf()
    pl.rc('font',size=24)
    dens1 = FITSFigure(parcubefile,convention='calabretta',slices=[0],figure=pl.figure(3))
    dens1.show_colorscale(vmin=2,vmax=6,cmap=cmhot)
    dens1.recenter(**zoomargs)
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

    column1 = 10**parcubefile[0].data[1] * parcubefile[0].data[5]
    column2 = 10**parcubefile[0].data[9] * parcubefile[0].data[13]
    totalcolumn = np.log10(column1 + column2)
    totalcolumn[totalcolumn == 0] = np.nan
    colhdu = fits.PrimaryHDU(data=totalcolumn, header=flatten_header(parcubefile[0].header))

    pl.figure(4)
    pl.clf()
    tcol = FITSFigure(colhdu, figure=pl.figure(4))
    tcol.show_colorscale(cmap=cmhot)
    tcol.colorbar._colorbar_axes.set_ylabel('log$_{10}$(N(H$_2$) cm$^{-2}$)')
    savefig('W51_H2CO_2parfittry10_totalcolumn.png', bbox_inches='tight')


    labels = {'dens':'log$_{10}$(n(H$_2$) cm$^{-3}$)',
              'lindens':'n(H$_2$) [cm$^{-3}]$',
              'velocity':'Velocity ($V_{LSR}$ km s$^{-1}$)',
              'ratio':r'Ratio $\tau_{obs} 1-1 / \tau_{obs} 2-2$',
              'width':'Line Width (km s$^{-1}$)',
              'column':'log$_{10}$(N(H$_2$) cm$^{-2}$)'}
    cmaps = {'dens':cmhot,
             'velocity':cmjet,
             'width':cmjet,
             'column':cmhot,}

    for parnum, param in [(0,'dens1'),(1,'column1'),(4,'velocity1'),(5,'width1'),
                          (8,'dens2'),(9,'column2'),(12,'velocity2'),(13,'width2')]:
        fig = pl.figure(5,figsize=(12,12))
        pl.clf()
        F = FITSFigure(parcubefile,convention='calabretta',slices=[parnum],figure=fig)
        F.show_colorscale(cmap=cmaps[param[:-1]])
        F.recenter(**zoomargs)
        F.colorbar._colorbar_axes.set_ylabel(labels[param[:-1]])
        savefig('W51_H2CO_2parfittry10_{0}.png'.format(param),bbox_inches='tight')

rfiles = ['W51_H2CO_max_ratio.fits',
          'W51_H2CO_mid_ratio.fits',
          'W51_H2CO_min_ratio.fits',
          ]
for fn in rfiles:
    fig = pl.figure(6,figsize=(12,12))
    fig.clf()
    hdu = fits.open(p1(fn))[0]
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=0, vmax=30)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['ratio'])
    savefig(fn.replace("fits","png"), bbox_inches='tight')


for sigma in (0.1,0.5,1.0):
    densfiles = ['W51_H2CO_logdensity_textbg_max_ratio_sigma{0:0.1f}.fits'.format(sigma),
                 'W51_H2CO_logdensity_textbg_mid_ratio_sigma{0:0.1f}.fits'.format(sigma),
                 'W51_H2CO_logdensity_textbg_min_ratio_sigma{0:0.1f}.fits'.format(sigma),
                 ]

    cmhot = mpl.cm.cubehelix
    cmhot.set_bad('white')
    cmhot.set_under('black')
    cmhot.set_over('white')

    for fn in densfiles:
        fig = pl.figure(6,figsize=(12,12))
        fig.clf()
        hdu = fits.open(p1(fn))[0]
        hdu.data = 10**hdu.data
        F = FITSFigure(hdu,convention='calabretta',figure=fig)
        F.show_colorscale(cmap=cmhot, vmin=10**2.5, vmax=10**4.5)
        F.recenter(**zoomargs)
        F.colorbar._colorbar_axes.set_ylabel(labels['lindens'])
        savefig(fn.replace("fits","png"), bbox_inches='tight')

pl.show()
