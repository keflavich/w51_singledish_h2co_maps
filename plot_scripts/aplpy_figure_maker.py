import aplpy

def FITSFigure(name, convention='calabretta', xcen=49.27, ycen=-0.32, width=0.9,
               height=0.4, **kwargs):
    F = aplpy.FITSFigure(name, convention=convention, **kwargs)

    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')

    F.show_colorscale()
    F.recenter(xcen,ycen,width=width,height=height)
    try:
        F.add_colorbar()
    except:
        pass

    return F
