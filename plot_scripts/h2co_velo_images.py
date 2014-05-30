import aplpy_figure_maker
import pylab as pl
from paths import figurepath,datapath
from common_constants import zoomargs

dpath = datapath

cm = pl.cm.jet
cm.set_bad('w')
cm.set_under('w')
cm.set_over('w')

pl.rcParams['font.size'] = 16

for line in (11,22):

    fig1 = pl.figure(1)
    fig1.clf()

    F = aplpy_figure_maker.FITSFigure(dpath+'H2CO%i_central_velocity.fits' % line,
                                      figure=fig1)
    F.colorbar.set_axis_label_text(r'$V_{LSR} [$km s$^{-1}]$')
    F.colorbar.set_axis_label_rotation(270)
    F.colorbar.set_axis_label_pad(30)
    F.recenter(**zoomargs)
    F.show_colorscale(vmin=45,vmax=77,cmap=cm)
    F.save(figurepath+'H2CO%i_central_velocity.pdf' % line)
