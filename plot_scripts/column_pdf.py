import numpy as np
from astropy.io import fits
import pylab as pl
import paths
import matplotlib
import os
matplotlib.rc_file(paths.rcfilepath)

d = fits.getdata(paths.dpath2('HIGAL_W51_mosaic_fit_070to500_N.fits'))
pl.clf()
X = pl.hist(np.log10(d[d==d]), bins=50, histtype='stepfilled', linewidth=3,
            edgecolor='k', facecolor=(0,0,0,0.5), log=True)
ax = pl.gca()
ax.set_ylim(1,ax.get_ylim()[1])
pl.draw()
pl.show()
pl.xlabel("Column Density N(H$_2$)")
pl.ylabel("Number of pixels")
ax.set_xlim(20.95, 23.22)
pl.savefig(paths.fpath("HiGal_ColumnDensityHistogram.pdf"))
