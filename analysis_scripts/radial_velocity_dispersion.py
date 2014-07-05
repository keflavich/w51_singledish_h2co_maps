import spectral_cube
import pyspeckit
import paths
from common_constants import distance, center
import numpy as np

sc = spectral_cube.SpectralCube.read(paths.dpath2('grs_48and50_cube_supersampledh2cogrid.fits'))
cube13 = pyspeckit.Cube(paths.dpath2('grs_48and50_cube_supersampledh2cogrid.fits'))

pars = []
spectra = []

import pylab as pl
for ii in range(1,25):
    pl.close(ii)

aprs = np.arange(50, 600, 30)

for apradius in aprs:
#apradius = 100
#if True:
    sp = cube13.get_apspec([center.l.deg, center.b.deg, apradius], coordsys='galactic')
    sp.plotter(xmin=30,xmax=90)
    guesses = [0.5,50,5,
               1,60,5,
               0.5,70,5]
    limits = [(0.1,3.5),(40,80),(1,15)]*3
    limits[1] = (45,55)
    limits[4] = (55,65)
    limits[7] = (65,75)
    limits[5] = (2,15)
    sp.specfit(fittype='gaussian', multifit=True, guesses=guesses, limits=limits, limited=[(True,True)]*9)
    sp.specfit.plot_fit(show_components=True)
    pars.append(sp.specfit.modelpars)
    spectra.append(sp)

[sp.specfit.plot_fit(show_components=True) for sp in spectra]

parsarr = np.array(pars)
pl.figure(0)
pl.clf()
pl.plot((8*np.log(2))**0.5*parsarr[:,2])
pl.plot((8*np.log(2))**0.5*parsarr[:,5])
pl.plot((8*np.log(2))**0.5*parsarr[:,8])
pl.plot((8*np.log(2))**0.5*(parsarr[:,5]**2+parsarr[:,8]**2)**0.5, '--')
geosum = (parsarr[:,2]**2+parsarr[:,5]**2+parsarr[:,8]**2)**0.5
pl.plot((8*np.log(2))**0.5*geosum, '--')
pl.show()

np.savetxt(paths.dpath('vwidth_vs_r.txt'),
           np.array([aprs, parsarr[:,2], parsarr[:,5], parsarr[:,8], geosum]).T,
           header=" ".join(["{0:>15s}".format(x) for x in ('rad_as','sigma_1','sigma_2','sigma_3','sigma_all')]),
           fmt="  %15i %15f %15f %15f %15f")
