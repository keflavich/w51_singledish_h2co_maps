from astropy.io import fits
import numpy as np
from paths import datapath,datapath_w51,figurepath
import pylab as pl
import matplotlib
from mpl_plot_templates import adaptive_param_plot
import FITS_tools
from FITS_tools.spectral_regrid import spec_pix_to_world,spec_world_to_pix
from astropy import units as u
from astropy.wcs import WCS

matplotlib.rc_file('ggplotrc')

h2co11 = fits.getdata(datapath+'W51_H2CO11_taucube_supersampled.fits')
h2co22 = fits.getdata(datapath+'W51_H2CO22_pyproc_taucube_lores_supersampled.fits')

header = fits.getheader(datapath+'W51_H2CO22_pyproc_taucube_lores_supersampled.fits')
x_to_vel = lambda x: spec_pix_to_world(x, WCS(header), 2, u.m/u.s).to(u.km/u.s)
vel_to_x = lambda x: spec_world_to_pix(x, WCS(header), 2, u.m/u.s)

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
comask = twelveco10 > 5

vranges = [(40,50),(50,60),(60,65),(65,75)]

for cf in cofiles:
    codata = fits.getdata(datapath_w51+cf)
    for vr in vranges:
        slc = slice(vel_to_x(vr[0]*u.km/u.s),vel_to_x(vr[1]*u.km/u.s))
        cdat = codata[slc,:,:]
        hdat = h2co11[slc,:,:]
        cmsk = comask[slc,:,:]
        hmsk = mask[slc,:,:]
        pl.figure(1)
        pl.clf()
        adaptive_param_plot(cdat[cmsk], hdat[cmsk], bins=50, threshold=25,
                            marker='.', alpha=0.5, color='b')
        adaptive_param_plot(cdat[hmsk], hdat[hmsk], bins=50, threshold=25,
                            marker='.', alpha=0.5, color='r', cmap=pl.cm.bone)

        pl.savefig(figurepath+"co_vs_h2co_v{0}to{1}_".format(vr[0],vr[1])+cf.replace(".fits",".pdf"))
