import numpy as np
from astropy import coordinates
from astropy import wcs
from astropy.io import fits
from image_tools import radialprofile
from astropy import units as u
from astropy import constants
import pylab as pl
from common_constants import distance
import paths

pl.mpl.rc_file('/Users/adam/work/w51_singledish_maps/plot_scripts/ggplotrc')

center = coordinates.Galactic(49.4904*u.deg, -0.3765*u.deg)

column_image = fits.open('HIGAL_W51_mosaic_fit_160to500_N.fits')
w = wcs.WCS(column_image[0].header)

pixsize = np.abs(w.wcs.get_cdelt()[0]) * u.deg

x,y = w.wcs_world2pix([center.l.deg], [center.b.deg], 0)

# Background assumed ~4e21
result = radialprofile.azimuthalAverage((column_image[0].data-4e21), center=(x,y), return_nr=True)
nr, centers, rprof = result

distprof = centers*(pixsize*distance).to(u.pc,u.dimensionless_angles())
rprof = rprof*u.cm**-2

muh2 = 2.8 # AMU per particle
radOK = (distprof < 25*u.pc) & (distprof > (pixsize*distance).to(u.pc, u.dimensionless_angles()))
massprof = (rprof*constants.m_p*muh2*nr*(pixsize*distance)**2).to(u.M_sun, u.dimensionless_angles())
mean_rprof = np.cumsum(rprof) / np.cumsum(nr)

# Truncate profiles
distprof,massprof = distprof[radOK], massprof[radOK]

cumul_massprof = np.cumsum(massprof)
dens_profile = (cumul_massprof.to(constants.m_p)/(4/3.*np.pi*distprof.to(u.cm)**3)) / muh2

fig = pl.figure(1)
fig.clf()
ax1 = pl.subplot(2,1,1)
ax1.plot(distprof.value, massprof.value)
ax1.set_ylabel("Mass per shell")

ax2 = pl.subplot(2,1,2, sharex=ax1)
ax2.plot(distprof.value, cumul_massprof.value)
ax2.set_xlabel("Radius (pc)")
ax2.set_ylabel("Enclosed Mass")

pl.subplots_adjust(hspace=0)
# from http://matplotlib.org/examples/pylab_examples/subplots_demo.html
pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
ax2.set_xlim(distprof.value.min(), distprof.value.max())

pl.savefig(paths.fpath("HiGal_RadialProfile_Mass.pdf"))

fig = pl.figure(2)
fig.clf()
ax1 = pl.subplot(2,1,1)
ax1.loglog(distprof.value, (massprof.to(constants.m_p)/muh2/(4/3.*np.pi*distprof.to(u.cm)**3)).value)
ax1.set_ylabel("Shell-averaged $n(H_2)$")

ax2 = pl.subplot(2,1,2, sharex=ax1)
ax2.loglog(distprof.value, dens_profile.value)
ax2.plot(distprof.value, distprof.value**-2 / distprof.value[0]**-2 * dens_profile[0].value, '--')
ax2.plot(distprof.value, distprof.value**-1.5 / distprof.value[0]**-1.5 * dens_profile[0].value, ':')
ax2.set_xlabel("Radius (pc)")
ax2.set_ylabel("Density")
ax2.set_xlim(distprof.value.min(), distprof.value.max())

pl.subplots_adjust(hspace=0)
# from http://matplotlib.org/examples/pylab_examples/subplots_demo.html
pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
pl.savefig(paths.fpath("HiGal_RadialProfile_VolumeDensity.pdf"))


fig = pl.figure(3)
fig.clf()
ax1 = pl.subplot(2,1,1)
ax1.loglog(distprof.value, rprof[radOK].value)
ax1.set_ylabel("Shell-averaged N(H$_2$)")

ax2 = pl.subplot(2,1,2, sharex=ax1)
ax2.loglog(distprof.value, mean_rprof[radOK].value)
ax2.plot(distprof.value, distprof.value**-2 / distprof.value[0]**-2 * mean_rprof[0].value, '--')
ax2.plot(distprof.value, distprof.value**-1.5 / distprof.value[0]**-1.5 * mean_rprof[0].value, ':')
ax2.set_xlabel("Radius (pc)")
ax2.set_ylabel("Column Density")
ax2.set_xlim(distprof.value.min(), distprof.value.max())

pl.subplots_adjust(hspace=0)
# from http://matplotlib.org/examples/pylab_examples/subplots_demo.html
pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

pl.savefig(paths.fpath("HiGal_RadialProfile_ColumnDensity.pdf"))

h2coimg = fits.getdata(paths.dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits"))
h2conr, h2cocen, h2codensprof = radialprofile.azimuthalAverage(h2coimg, center=(x,y), return_nr=True, interpnan=True, mask=np.isfinite(h2coimg))



pl.show()
