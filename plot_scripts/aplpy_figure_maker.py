import aplpy
import matplotlib
from astropy import units as u

def FITSFigure(name, convention='calabretta', xcen=49.27, ycen=-0.32, width=0.9,
               height=0.4, colorbar=True, color=True, grid=True,
               cmap=matplotlib.cm.afmhot_r, im_zorder=20, 
               scale_zorder=40,
               scalebar=True, transparent_nan=True, **kwargs):
    F = aplpy.FITSFigure(name, convention=convention, **kwargs)

    refresh = F._parameters.auto_refresh
    F.set_auto_refresh(False)

    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')

    if transparent_nan:
        cmap.set_bad(color='w', alpha=0.0)

    if color:
        F.show_colorscale(cmap=cmap)
    else:
        F.show_grayscale()
    if colorbar:
        try:
            F.add_colorbar()
        except:
            pass

    if grid:
        F.add_grid()
        F.grid.set_color('black')
        F.grid.set_linestyle('dotted')
        F.grid.set_alpha(0.5)
        F.grid.ax.gridlines.set_zorder(0)

    F.image.set_zorder(im_zorder)

    F.recenter(xcen,ycen,width=width,height=height)

    if scalebar:
        F.add_scalebar(((10*u.pc)/(5.1*u.kpc)*u.radian).to(u.deg).value)
        F.scalebar.set_label("10 pc")
        F.scalebar.set_font_size(18)
        F.scalebar.set_font_weight('bold')
        F.scalebar.set_color((0.8,0.3,0.01,0.9))
        F.scalebar.set_linewidth(3)
        F.scalebar._scalebar.set_zorder(scale_zorder)

    F.set_auto_refresh(refresh)

    return F
