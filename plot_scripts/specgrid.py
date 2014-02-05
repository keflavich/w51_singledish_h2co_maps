import pyspeckit
import itertools
from astropy import wcs
import pylab as pl
import numpy as np
from FITS_tools import strip_headers
do_some_plots=True
np.seterr(all='ignore')

datapath = '/Users/adam/work/h2co/maps/W51/'
figpath = '/Users/adam/work/h2co/figures/'

def spectral_grid(cube11=pyspeckit.Cube(datapath+'W51_H2CO11_taucube_supersampled.fits'),
                  #h213co11cube=pyspeckit.Cube(datapath+'.fits'),
                  cube22=pyspeckit.Cube(datapath+'W51_H2CO22_pyproc_taucube_lores_supersampled.fits'),
                  figure=pl.figure(1,figsize=(10,10)),
                  yrange=(-2.1,0.5),
                  ratio=False):

    for c in [cube11,cube22]: # h213cocube
        c.xarr.convert_to_unit('km/s')

    def some_plots(xx,yy,dolegend=False):
        #c13 = h213co11cube.get_spectrum(xx,yy)
        c12 = cube11.get_spectrum(xx,yy)
        c22 = cube22.get_spectrum(xx,yy)

        for c in (c12,c22):# c13,
            c.plotter.autorefresh=False

        c12.baseline(exclude=[-225,200],order=5)
        #c13.baseline(exclude=[-225,200],order=5)
        #c13.plotter(label='H$_{2}$$^{13}$CO 1-1',axis=pl.gca(),color='b',alpha=0.5)
        #c12.plotter(axis=c13.plotter.axis,clear=False,color='k',label="H$_{2}$CO 1-1")
        c12.plotter(axis=pl.gca(),color='k',label="H$_{2}$CO 1-1")
        #(c13*6).plotter(label='6$\\times$H$_{2}$$^{13}$CO',axis=pl.gca(),color='r',clear=False)
        c22.plotter(axis=c12.plotter.axis,clear=False,color='r',linewidth=2,alpha=0.8,
                    label='H$_2$CO 2-2')

        if ratio:
            r = c12.copy()
            r.data = c22.data/c12.data
            r.data[(r.data>1)] = np.nan
            r.data[(r.data<1/13.)] = np.nan
            r.plotter(axis=c12.plotter.axis,clear=False,color='r')
        if dolegend:
            pl.legend(loc='best')
        pl.gca().set_xlim(-100,100)
      

    pl.clf()

    pl.subplots_adjust(wspace=0,hspace=0,left=0.05,right=0.95,bottom=0.05,top=0.95)

    nx,ny = 6,6
    xc,yc = 52,48 # center
    xc,yc = 45,42 # bottom left
    w = wcs.WCS(strip_headers.flatten_header(cube11.header))
    for ii,(spy,spx) in enumerate(itertools.product(range(nx),range(ny))):
        sp = pl.subplot(nx,ny,ii+1)
        sp.zorder = nx-spx+ny*spy
        some_plots(spx*2+xc,(ny-spy-1)*2+yc,dolegend=False)#,dolegend=(ii==24))
        sp.set_ylim(*yrange)
        if spx > 0:
            sp.set_yticks([])
        if spy < ny-1:
            sp.set_xticks([])

        (l,b), = w.wcs_pix2world([[spx*2+xc,(ny-spy-1)*2+yc]],0)
        sp.annotate('%0.3f %+0.3f' % (l,b),(0.5,0.9),xycoords='axes fraction',fontsize=12,
                    horizontalalignment='center')

        pl.draw()

#spectral_grid()
#pl.savefig(figpath+'/spectralgrid_absorption.pdf')

spectral_grid(#cube11=pyspeckit.Cube('/Users/adam/work/gc/limabean/LimaBean_H2CO11_taucube.fits'),
              #cube22=pyspeckit.Cube('/Users/adam/work/gc/limabean/LimaBean_H2CO22_taucube_smoothtoCband.fits'),
              #h213co11cube=pyspeckit.Cube('/Users/adam/work/gc/limabean/LimaBean_H213CO_taucube.fits'),
              figure=pl.figure(2,figsize=(10,10)),
              yrange=(-0.05,0.2),
              ratio=False)
pl.savefig(figpath+'spectralgrid_optdepth.pdf')

