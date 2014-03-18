import pylab as pl
import numpy as np
import astroML.plotting as ampl
from astropy.io import fits
import sys
sys.path.append("/Users/adam/work/h2co/lowdens/code/")
from smoothtau_models import select_data,generate_tau_functions
import itertools

pl.rcParams['font.size'] = 20

ratio = fits.getdata('W51_H2CO11_to_22_tau_ratio_supersampled_neighbors.fits')

colors = itertools.cycle(["#"+x for x in '348ABD, 7A68A6, A60628, 467821, CF4457, 188487, E24A33'.split(', ')])

abunds = (-8.5,-10)

dcs = {}

for abund in abunds:
    dcs[abund] = {}
    for sigma in (0.001,1.0):
        tau1,tau2,dens,col = select_data(abund)
        tau,vtau,vtau_ratio = generate_tau_functions(abundance=abund)

        tauratio = vtau_ratio(dens, line1=tau1, line2=tau2, sigma=sigma)

        ok = np.arange(tauratio.size) > np.argmax(tauratio)

        def ratio_to_dens(ratio):
            inds = np.argsort(tauratio[ok])
            return np.interp(ratio, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)

        if sigma == 0.001: sigma = 0

        r = ratio_to_dens(ratio)
        dcs[abund][sigma] = r[r==r]

pl.figure(3)
pl.clf()
ax = pl.axes([0.1,0.1,0.65,0.8])

for abund in abunds:
    for sigma in (0,1.0):
        rr = dcs[abund][sigma]
        counts, bins, patches = ampl.hist(rr, bins=100, log=True, histtype='step', 
                                          linewidth=2, alpha=1.0, color=colors.next(),
                                          label=r"$X=%s, \sigma=%i$" % (abund,sigma))
        pl.xlabel(r'$\log(n(H_2))$ [cm$^{-3}$]')

ax.set_ylim(1,2e3)
ax.set_xlim(1.5,6)
ax.set_ylabel("$N$(voxels)")

pl.legend(bbox_to_anchor=(1.0, 1.0),fontsize=18,loc='upper left')

pl.savefig('/Users/adam/work/h2co/maps/paper/figures/cube_histograms_density_ppv_multimodel.pdf',bbox_inches='tight')
