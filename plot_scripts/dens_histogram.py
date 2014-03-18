import pylab as pl
import numpy as np
import astroML.plotting as ampl
from astropy.io import fits
from agpy.mad import MAD
import itertools
import pyspeckit
import mpl_toolkits

pl.rcParams['font.size'] = 20


nh2 = r'$\log(n(H_2))$ [cm$^{-3}$]'

dens = fits.getdata('W51_H2CO11_to_22_logdensity_supersampled.fits')
cube = pyspeckit.Cube('W51_H2CO11_to_22_logdensity_supersampled.fits')
densOK = dens==dens

pl.figure(2)
pl.clf()
ax = pl.subplot(1,2,1)
densp = dens[densOK]
counts,bins,patches = ampl.hist(densp, bins=100, log=True, histtype='step', linewidth=2, alpha=0.8, color='k')
ylim = ax.get_ylim()

sp = pyspeckit.Spectrum(xarr=(bins[1:]+bins[:-1])/2.,data=counts)
sp.specfit(guesses= [660.23122694399035,
                     3.1516848752486522,
                     0.33836811902343894,
                     396.62714060001434,
                     2.5539176548294318,
                     0.32129608858734149,
                     199.13259679527025,
                     3.730112763513838,
                     0.4073913996012487])

def ntuples(lst, n):
    return zip(*[lst[i::n]+lst[:i:n] for i in range(n)])

dashes = itertools.cycle([ (16,8),  (4,4), (8,4,4,4),  ])
hatches = itertools.cycle(('//','--','||'))

med, mad = np.median(densp),MAD(densp)
composite = bins*0
for p,m,w in ntuples(sp.specfit.parinfo.values,3):
    g = p*np.exp(-(bins-m)**2/(2*(w)**2))
    pl.plot(bins, g,'r', dashes=dashes.next(), label=r'$\mu=%0.1f, \sigma=%0.1f$' % (m,w))
    #pl.fill_between(bins,0,g,color='r',alpha=0.1,hatch=hatches.next(),facecolor='none')
    composite += g

pl.plot(bins, composite, 'b', linewidth=3, alpha=0.5, dashes=(10,5))
pl.fill_between(bins,1,composite,color='b', hatch=r'\\', alpha=0.2, facecolor='none')
counts,bins,patches = ampl.hist(densp, bins=100, log=True, histtype='step', linewidth=2, alpha=0.8, color='k')
ax.set_ylim(0.9,ylim[1])
ax.set_xticklabels([])

pl.legend(loc='best',fontsize=14)
pl.ylabel("$N($voxels$)$")

divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
lowax = divider.append_axes('bottom', size='10%', pad=0)
sp.specfit.plotresiduals(fig=2, clear=False, axis=lowax, linewidth=2, yoffset=0.5)
lowax.set_title("")
#lowax.yaxis.set_ticks_position('right')
lowax.set_yticks([-50,0,50])
lowax.set_yticklabels([])
ylim2 = lowax.get_ylim()

pl.xlabel(r'$\log(n(H_2))$ [cm$^{-3}$]')

ax = pl.subplot(1,2,2)
densp = dens[densOK]
counts,bins,patches = ampl.hist(densp, bins=100, log=True, histtype='step', linewidth=2, alpha=0.8, color='k')
ylim = ax.get_ylim()

sp = pyspeckit.Spectrum(xarr=(bins[1:]+bins[:-1])/2.,data=counts)
sp.specfit(guesses= [1000,
                     3.1,
                     0.5,])

for p,m,w in ntuples(sp.specfit.parinfo.values,3):
    g = p*np.exp(-(bins-m)**2/(2*(w)**2))
    pl.plot(bins, g,'r', dashes=dashes.next(), label=r'$\mu=%0.1f, \sigma=%0.1f$' % (m,w))
    pl.fill_between(bins,0,g,color='r',alpha=0.1,hatch=hatches.next(),facecolor='none')

counts,bins,patches = ampl.hist(densp, bins=100, log=True, histtype='step', linewidth=2, alpha=0.8, color='k')
ax.set_ylim(0.9,ylim[1])
ax.set_xticklabels([])

pl.legend(loc='best',fontsize=14)
ax.set_yticks([])

divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
lowax = divider.append_axes('bottom', size='10%', pad=0)
sp.specfit.plotresiduals(fig=2, clear=False, axis=lowax, linewidth=2, yoffset=0.5)
lowax.set_title("")
#lowax.yaxis.set_ticks_position('right')
#lowax.set_yticks([])
lowax.yaxis.set_ticks_position('right')
lowax.set_yticks([-50,0,50])
lowax.set_ylim(*ylim2)
lowax.set_xticks(lowax.get_xticks()[1:])

pl.xlabel(r'$\log(n(H_2))$ [cm$^{-3}$]')

pl.subplots_adjust(wspace=0)

pl.savefig('/Users/adam/work/h2co/maps/paper/figures/cube_histograms_density_ppv.pdf',bbox_inches='tight')


pl.figure(6)
pl.clf()
ax = pl.gca()
n_of_v = np.array([(np.median(sl[sl==sl]),MAD(sl[sl==sl])) for sl in dens])

ax.plot(cube.xarr, n_of_v[:,0])
ax.fill_between(cube.xarr,n_of_v[:,0]-n_of_v[:,1],n_of_v[:,0]+n_of_v[:,1], alpha=0.2)
ax.set_xlim(42,75)
ax.set_xlabel(r'$V_{LSR} [$km s$^{-1}]$')
ax.set_ylabel(nh2)

pl.savefig('/Users/adam/work/h2co/maps/paper/figures/density_vs_velocity.pdf',bbox_inches='tight')








pl.show()
