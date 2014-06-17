import os

import pylab as pl
import numpy as np
import matplotlib as mpl

from aplpy_figure_maker import FITSFigure
from paths import datapath, datapath_w51, figurepath

import spectral_cube

denscube = spectral_cube.SpectralCube.read("H2CO_ParameterFits_bestdens.fits")
colcube  = spectral_cube.SpectralCube.read("H2CO_ParameterFits_bestcol.fits")
chi2cube = spectral_cube.SpectralCube.read("H2CO_ParameterFits_bestchi2.fits")
for vrange in ([52,66], [65,75]):
    fig = pl.figure(6,figsize=(12,12))
    fig.clf()
    
    F = FITSFigure(hdu,convention='calabretta',figure=fig)
    F.show_colorscale(cmap=cmhot, vmin=0, vmax=30)
    F.recenter(**zoomargs)
    F.colorbar._colorbar_axes.set_ylabel(labels['ratio'])
    savefig(fn.replace("fits","png"), bbox_inches='tight')


