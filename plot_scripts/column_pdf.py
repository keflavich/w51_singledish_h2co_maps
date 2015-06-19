import numpy as np
from astropy.io import fits
import pylab as pl
import paths
import matplotlib
import os
matplotlib.rc_file(paths.rcfilepath)

d = fits.getdata(paths.dpath2('HIGAL_W51_mosaic_fit_070to500_N.fits'))
f = fits.open(paths.dpath2('HIGAL_W51_mosaic_fit_070to500_N.fits'))
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

import pyregion
regs = pyregion.open(paths.rpath('main_and_b_squares.reg'))

reg1 = pyregion.ShapeList([regs[0]])
reg2 = pyregion.ShapeList([regs[1]])
mask1 = reg1.get_mask(f[0])
mask2 = reg2.get_mask(f[0])

X = pl.hist(np.log10(d[mask1]), bins=50, histtype='stepfilled', linewidth=3,
            edgecolor='g', facecolor=(0,1,0,0.5), log=True)
X = pl.hist(np.log10(d[mask2]), bins=50, histtype='stepfilled', linewidth=3,
            edgecolor='b', facecolor=(0,0,1,0.5), log=True)
pl.savefig(paths.fpath("HiGal_ColumnDensityHistogram_mainandb.pdf"))
