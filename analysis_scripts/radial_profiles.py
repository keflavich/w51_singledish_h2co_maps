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
import FITS_tools

pl.mpl.rc_file('/Users/adam/work/w51_singledish_maps/plot_scripts/ggplotrc')

center = coordinates.Galactic(49.4904*u.deg, -0.3765*u.deg)

h2co_densfile = paths.dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits")
column_image = FITS_tools.hcongrid.hcongrid_hdu(fits.open(paths.dpath2('HIGAL_W51_mosaic_fit_160to500_N.fits'))[0],
                                                fits.getheader(h2co_densfile))
w = wcs.WCS(column_image.header)
dens_image = fits.open(h2co_densfile)[0]
w2 = wcs.WCS(dens_image.header)

pixsize = np.abs(w.wcs.get_cdelt()[0]) * u.deg
h2co_pixsize = np.abs(w2.wcs.get_cdelt()[0]) * u.deg

# Formaldehyde: mask out NaNs
x,y = w2.wcs_world2pix([center.l.deg], [center.b.deg], 0)
data = dens_image.data
mask = np.isfinite(data)
data[~mask] = 0
result = radialprofile.azimuthalAverage(data, center=(x,y), return_nr=True, mask=mask)
h2co_nr, h2co_centers, h2co_rprof = result


x,y = w.wcs_world2pix([center.l.deg], [center.b.deg], 0)

# Background assumed ~4e21
nr_all, centers_all, rprof_all = radialprofile.azimuthalAverage((column_image.data-4e21), center=(x,y),
                                                                return_nr=True)
# Use H2CO mask so that we measure profile in same apertures
result = radialprofile.azimuthalAverage((column_image.data-4e21), center=(x,y),
                                        return_nr=True, mask=mask)
nr, centers, rprof = result

distprof = centers*(pixsize*distance).to(u.pc,u.dimensionless_angles())
rprof = rprof*u.cm**-2
h2co_distprof = h2co_centers*(h2co_pixsize*distance).to(u.pc,u.dimensionless_angles())
h2co_rprof = h2co_rprof*u.cm**-3

muh2 = 2.8 # AMU per particle
radOK = (distprof < 25*u.pc) & (distprof > (pixsize*distance).to(u.pc, u.dimensionless_angles()))
massprof = (rprof*constants.m_p*muh2*nr*(pixsize*distance)**2).to(u.M_sun, u.dimensionless_angles())
mean_rprof = np.cumsum(rprof) / np.cumsum(nr)

# Truncate profiles
distprof,massprof = distprof[radOK], massprof[radOK]
dr = np.mean(np.diff(distprof))

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
shell_volume = (nr/nr_all.astype('float'))[radOK] * 4/3.*np.pi*((distprof+dr/2)**3-(distprof-dr/2)**3)
shell_density = (massprof.to(constants.m_p)/muh2/shell_volume.to(u.cm**3))
ax1.loglog(distprof.value, shell_density.value)
ax1.loglog(h2co_distprof.value, 10**h2co_rprof.value)
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
