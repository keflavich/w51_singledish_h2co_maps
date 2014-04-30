import aplpy_figure_maker
import pylab as pl
from common_constants import datapath_cubes,figpath

dpath = datapath_cubes

cm = pl.cm.jet
cm.set_bad('w')
cm.set_under('w')
cm.set_over('w')

pl.rcParams['font.size'] = 16

for rrl in (77,110):

    fig1 = pl.figure(1)
    fig1.clf()

    F = aplpy_figure_maker.FITSFigure(dpath+'H%ia_central_velocity.fits' % rrl,
                                      figure=fig1)
    F.colorbar.set_axis_label_text(r'$V_{LSR} [$km s$^{-1}]$')
    F.colorbar.set_axis_label_rotation(270)
    F.colorbar.set_axis_label_pad(30)
    F.show_colorscale(vmin=45,vmax=77,cmap=cm)
    F.save(figpath+'H%ia_central_velocity.png' % rrl)
