import aplpy
import matplotlib

def FITSFigure(name, convention='calabretta', xcen=49.27, ycen=-0.32, width=0.9,
               height=0.4, colorbar=True, color=True, grid=True,
               cmap=matplotlib.cm.afmhot_r, im_zorder=20, 
               transparent_nan=True, **kwargs):
    F = aplpy.FITSFigure(name, convention=convention, **kwargs)

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

    return F
