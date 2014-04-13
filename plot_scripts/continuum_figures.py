import pylab as pl
import numpy as np
import aplpy
from paths import datapath

for fn in (1,2):
    pl.figure(fn)
    pl.clf()


oneonefn = datapath+'W51_H2CO11_cube_supersampled_continuum.fits'
twotwofn = datapath+'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits'
F1 = aplpy.FITSFigure(oneonefn,convention='calabretta',figure=pl.figure(1),)
F2 = aplpy.FITSFigure(twotwofn,convention='calabretta',figure=pl.figure(2),)

F1.set_tick_labels_xformat('dd.d')
F2.set_tick_labels_xformat('dd.d')
F1.set_tick_labels_yformat('dd.d')
F2.set_tick_labels_yformat('dd.d')

F1.show_grayscale(stretch='log',vmid=-1,invert=True,vmax=300)
F2.show_grayscale(stretch='log',vmid=-1,invert=True,vmax=25)
F1.recenter(49.23,-0.28,width=1,height=0.5)
F2.recenter(49.23,-0.28,width=1,height=0.5)
F1.add_colorbar()
F2.add_colorbar()

con11 = 10**np.linspace(-1,2.5,6)
con22 = 10**np.linspace(-1,1.5,6)
linewidths = np.linspace(2,0.5,6)

def loc1(x):
    return np.interp(x, F1.colorbar._colorbar._values, np.linspace(0,1,len(F1.colorbar._colorbar._values)))
def loc2(x):
    return np.interp(x, F2.colorbar._colorbar._values, np.linspace(0,1,len(F2.colorbar._colorbar._values)))

F1.colorbar.set_ticks(con11)
F2.colorbar.set_ticks(con22)
L1 = F1.colorbar._colorbar_axes.hlines(loc1(con11),0,1, color='r',linewidth=2,alpha=0.8)
L2 = F2.colorbar._colorbar_axes.hlines(loc2(con22),0,1, color='r',linewidth=2,alpha=0.8)

for LL in (L1,L2):
    LL.set_linewidths(linewidths)

F1.colorbar._colorbar_axes.figure.show()
F2.colorbar._colorbar_axes.figure.show()

F1.save('/Users/adam/work/h2co/maps/paper/figures/continuum_11.pdf',dpi=100)
F2.save('/Users/adam/work/h2co/maps/paper/figures/continuum_22.pdf',dpi=100)
F1.show_contour(oneonefn, levels=con11, colors=['r'], linewidths=linewidths)
F2.show_contour(twotwofn, levels=con22, colors=['r'], linewidths=linewidths)
F1.save('/Users/adam/work/h2co/maps/paper/figures/continuum_11_contours.pdf',dpi=100)
F2.save('/Users/adam/work/h2co/maps/paper/figures/continuum_22_contours.pdf',dpi=100)

if False:
    for fn in (1,2):
        pl.figure(fn)
        pl.draw()
        pl.show()
