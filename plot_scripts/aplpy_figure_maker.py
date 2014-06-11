import aplpy
import matplotlib

def FITSFigure(name, convention='calabretta', xcen=49.27, ycen=-0.32, width=0.9,
               height=0.4, colorbar=True, color=True, **kwargs):
    F = aplpy.FITSFigure(name, convention=convention, **kwargs)

    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')

    if color:
        F.show_colorscale(cmap=matplotlib.cm.afmhot_r)
    else:
        F.show_grayscale()
    if colorbar:
        try:
            F.add_colorbar()
        except:
            pass

    F.recenter(xcen,ycen,width=width,height=height)

    return F
