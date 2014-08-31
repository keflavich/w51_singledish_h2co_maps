import pylab as pl

from aplpy_figure_maker import FITSFigure
from paths import dpath, fpath

for label in ("","lower","upper"):
    files = [dpath(label+'peak_optdepth_11.fits'),
             dpath(label+'peak_optdepth_22.fits'),
             dpath(label+'peak_optdepth_ratio.fits'),
             dpath(label+'peak_optdepth_11_tex1.0.fits'),
             dpath(label+'peak_optdepth_22_tex1.5.fits'),
             dpath(label+'peak_optdepth_ratio_tex1.0_1.5.fits')]

    for ii in xrange(1,len(files)+1):
        pl.figure(ii)
        pl.clf()

    figs = [FITSFigure(fn, figure=pl.figure(ii+1))
            for ii,fn in enumerate(files)]

    figs[2].show_colorscale(cmap=pl.cm.gist_stern, vmin=0, vmax=15)
    figs[5].show_colorscale(cmap=pl.cm.gist_stern, vmin=0, vmax=15)

    figs[0].save(fpath(label+'peak_observed_opticaldepth_11.pdf'), dpi=72)
    figs[1].save(fpath(label+'peak_observed_opticaldepth_22.pdf'), dpi=72)
    figs[2].save(fpath(label+'peak_observed_opticaldepth_ratio.pdf'), dpi=72)
    figs[5].save(fpath(label+'peak_observed_opticaldepth_ratio_tex.pdf'), dpi=72)

pl.show()
