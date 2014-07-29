import pylab as pl
import numpy as np
import astroML.plotting as ampl
from astropy.io import fits
from agpy.mad import MAD
from paths import dpath,fpath

pl.rcParams['font.size'] = 20

h2co11 = fits.getdata(dpath('W51_H2CO11_taucube_supersampled.fits'))
h2co22 = fits.getdata(dpath('W51_H2CO22_pyproc_taucube_lores_supersampled.fits'))

ratio = fits.getdata(dpath('W51_H2CO11_to_22_tau_ratio_supersampled_neighbors.fits'))

pl.close(1)
pl.figure(1, figsize=(10,10))
pl.clf()

ax1 = pl.subplot(3,1,1)
oneone = h2co11[h2co11==h2co11]
counts, bins, patches = ampl.hist(oneone, bins=100, log=True, histtype='step',
                                  linewidth=2, alpha=0.8, color='k')
ylim = ax1.get_ylim()
med, mad = np.median(oneone),MAD(oneone)
pl.plot(bins,counts.max()*np.exp(-(bins-med)**2/(2*mad**2)),'r--')
ax1.set_ylim(*ylim)
ax1.set_xlabel("$\\tau_{obs}($H$_2$CO 1-1$)$", labelpad=10)
ax1.set_ylabel("$N($voxels$)$")

ax2 = pl.subplot(3,1,2)
twotwo = h2co22[h2co22==h2co22]
counts, bins, patches = ampl.hist(twotwo, bins=100, log=True, histtype='step', linewidth=2, alpha=0.8, color='k')
ylim = ax2.get_ylim()
med, mad = np.median(twotwo),MAD(twotwo)
pl.plot(bins,counts.max()*np.exp(-(bins-med)**2/(2*mad**2)),'r--')
ax2.set_ylim(*ylim)
ax2.set_xlabel("$\\tau_{obs}($H$_2$CO 2-2$)$", labelpad=10)
ax2.set_ylabel("$N($voxels$)$")


ax3 = pl.subplot(3,1,3)
rr = ratio[ratio==ratio]
counts, bins, patches = ampl.hist(rr, bins=100, log=True, histtype='step', linewidth=2, alpha=0.8, color='k')
ax3.set_xlabel("Ratio (1-1)/(2-2)")
ax3.set_ylabel("$N($voxels$)$")
#ylim = ax3.get_ylim()
#med, mad = np.median(rr),MAD(rr)
#pl.plot(bins,counts.max()*np.exp(-(bins-med)**2/(2*mad**2)),'k--')
#ax3.set_ylim(*ylim)

pl.subplots_adjust(hspace=0.4)
pl.savefig(fpath('cube_histograms_tau_and_ratio.pdf'),
           bbox_inches='tight')
pl.show()
