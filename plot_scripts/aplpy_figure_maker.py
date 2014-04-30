import aplpy

def FITSFigure(name, convention='calabretta', xcen=49.23, ycen=-0.28, width=1,
               height=0.5, **kwargs):
    F = aplpy.FITSFigure(name, convention=convention, **kwargs)

    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')

    F.show_colorscale()
    F.recenter(xcen,ycen,width=width,height=height)
    F.add_colorbar()

    return F
