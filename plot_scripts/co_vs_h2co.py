from astropy.io import fits
from FITS_tools.strip_headers import flatten_header
import numpy as np
from paths import datapath,datapath_w51
import pylab as pl
import matplotlib
matplotlib.rc_file('ggplotrc')

h2co11 = fits.getdata(datapath+'W51_H2CO11_taucube_supersampled.fits')
h2co22 = fits.getdata(datapath+'W51_H2CO22_pyproc_taucube_lores_supersampled.fits')

noise11 = h2co11[:50,:,:].std(axis=0)
noise22 = h2co22[:50,:,:].std(axis=0)

sn11 = h2co11/noise11
sn22 = h2co22/noise22

mask = (sn11 > 2) & (sn22 > 2)

cofiles = ['CO_12CO10_supersampledh2cogrid.fits',
           'CO_12CO21_supersampledh2cogrid.fits',
           'CO_12CO32_supersampledh2cogrid.fits',
           'CO_13CO10_supersampledh2cogrid.fits',
           'CO_13CO21_supersampledh2cogrid.fits',
           'CO_13CO32_supersampledh2cogrid.fits',
           'CO_13CO32_supersampledh2cogrid_smoothy.fits',
           'CO_C18O32_supersampledh2cogrid.fits',
           ]

twelveco10 = fits.getdata(datapath_w51+cofiles[0])
comask = twelveco10 > 3.5

for cf in cofiles:
    codata = fits.getdata(datapath_w51+cf)
    pl.figure(1)
    pl.clf()
    pl.plot(codata[comask], h2co11[comask], ',', alpha=0.5)

    pl.savefig(
