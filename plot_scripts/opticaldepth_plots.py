import pylab as pl

from aplpy_figure_maker import FITSFigure
from paths import dpath, fpath

files = [dpath('peak_optdepth_11.fits'),
         dpath('peak_optdepth_22.fits'),
         dpath('peak_optdepth_ratio.fits')]

for ii in xrange(1,4):
    pl.figure(ii)
    pl.clf()

figs = [FITSFigure(fn, figure=pl.figure(ii+1))
        for ii,fn in enumerate(files)]

figs[2].show_colorscale(cmap=pl.cm.gist_stern, vmin=0, vmax=15)

figs[0].save(fpath('peak_observed_opticaldepth_11.pdf'))
figs[1].save(fpath('peak_observed_opticaldepth_22.pdf'))
figs[2].save(fpath('peak_observed_opticaldepth_ratio.pdf'))

pl.show()
