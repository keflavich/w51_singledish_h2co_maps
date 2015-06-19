import numpy as np
import imf
import pylab as pl
from astropy.utils.console import ProgressBar
pl.matplotlib.rc_file('pubfiguresrc')

nstars_per_cluster = 100
maxmass = 3e4

cl_masses = np.logspace(2.5,np.log10(maxmass))
clusters = np.array([[imf.make_cluster(mass*(np.random.randn()/20.+1.),silent=True) for ii in range(nstars_per_cluster)] for mass in ProgressBar(cl_masses)])
luminosities = [imf.lum_of_cluster(cl) for section in ProgressBar(clusters) for cl in section]
#cl_masses_flat = np.array([x for mass in cl_masses for x in [mass for ii in range(nstars_per_cluster)]])
cl_masses_flat = np.array([sum(cl) for section in clusters for cl in section])
n_ostars = np.array([(cl>8).sum() for section in clusters for cl in section])

# guess, not measurement
nostars_w51 = n_ostars[(cl_masses_flat > 4e3) & (cl_masses_flat < 9e3)]


pl.figure(6).clf()
pl.semilogx(cl_masses_flat, luminosities, '.', alpha=0.5)
pl.axhline(np.log10(8.3e6), linewidth=10, alpha=0.1, zorder=-5)
pl.axhline(np.log10(2e7), linewidth=10, alpha=0.1, zorder=-5)
pl.xlim(10**2.5, (maxmass))
pl.ylim(5.5, 8.1)
pl.xlabel("Cluster mass ($M_\odot$)")
pl.ylabel("Cluster luminosity (log $L_\odot$)")
pl.savefig("cluster_mass_vs_luminosity.png")

pl.figure(9).clf()
pl.loglog(cl_masses_flat, n_ostars, '.', alpha=0.5)
pl.xlim(10**2.5, (maxmass))
pl.ylim(0, n_ostars.max())
pl.xlabel("Cluster mass ($M_\odot$)")
pl.ylabel("Number of OB-stars")
pl.savefig("cluster_mass_vs_n_ostars.png")


def lnprob(mass, luminosity=2e7, luminosity_error=5e6):
    if mass < 1:
        return -np.inf, 0
    if mass > 1e5:
        # use a best-fit relation; scatter is minimal
        lum = 10**(np.log10(mass) + 3.55)
        # n > 20 msun
        nostars = 10**(np.log10(mass) - 2.5)
    else:
        cluster = imf.make_cluster(mass, silent=True)
        lum = 10**imf.lum_of_cluster(cluster)
        nostars = (cluster > 20).sum()
    return -(luminosity-lum)**2/(2*luminosity_error**2), nostars

import emcee
ndim = 1
nwalkers = 48
p0 = np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim))*1000 + 5000
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=4)
pos,prob,state,blobs = sampler.run_mcmc(p0, 400)




pl.figure(7).clf()
pl.hist(nostars_w51, bins=30)

pl.figure(8)
pl.clf()
pl.hist(sampler.flatchain, bins=100)
pl.xlabel("Cluster Mass at $L=2 \pm 0.5\\times10^7$ L$_\odot$")

print("Cluster mass: {0} \pm {1}".format(sampler.flatchain.mean(),
                                         sampler.flatchain.std()))

pl.figure(9)
pl.clf()
flatblobs = np.array([x for y in sampler.blobs for x in y])
a,b,c = pl.hist(flatblobs, bins=np.arange(0,50,1))
x = np.linspace(0,50,500)
pl.plot(x, a.max()*np.exp(-(x-flatblobs.mean())**2 / (2*flatblobs.var())), '-', linewidth=2, alpha=0.5)
pl.xlabel("$N(M>20$M$_\odot)$ at $L=2 \pm 0.5\\times10^7$ L$_\odot$")

print("N(Ostars): {0} \pm {1}".format(flatblobs.mean(),
                                      flatblobs.std()))

