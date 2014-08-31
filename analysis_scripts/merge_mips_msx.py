import mpl_plot_templates
from astropy.io import fits
import numpy as np
import pylab as pl
import matplotlib

msx = fits.open('MSX_E_W51_mosaic_bgmatch.fits')
mips = fits.open('W51_mips1.fits')
msxd = msx[0].data
mipsd = mips[0].data

ok = np.isfinite(msxd) & np.isfinite(mipsd)

msxOK = msxd > 500
fitOK = msxOK * (mipsd<1400)*(mipsd>250)

pl.figure(1)
pl.clf()
mpl_plot_templates.adaptive_param_plot(mipsd.flat, msxd.flat, marker=',',
                                       bins=50, threshold=10)
mpl_plot_templates.adaptive_param_plot(mipsd[fitOK], msxd[fitOK], marker=',',
                                       bins=50, threshold=10, color='r')

pv = np.polyfit(mipsd[fitOK],msxd[fitOK],1)
print pv
x = np.linspace(0,2000)
pl.plot(x, np.polyval(pv, x), linewidth=3, alpha=0.5)

mixed = msxd
replace = msxd > 1000
replace = ((mipsd > 1500) & (msxd > 1000)) | (True - np.isfinite(mipsd))
mixed[True-replace] = mipsd[True-replace]*pv[0] + pv[1]

pl.figure(2)
pl.clf()
pl.imshow(mixed, norm=matplotlib.colors.LogNorm())

msx[0].data = mixed
msx.writeto('MSX_MIPS_merged.fits',clobber=True)

pl.show()
