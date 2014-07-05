import numpy as np
from astropy import coordinates
from astropy import wcs
from astropy.io import fits
from image_tools import radialprofile
from astropy import units as u
from astropy import constants
import pylab as pl
from common_constants import distance, center
import paths
import FITS_tools

pl.mpl.rc_file('/Users/adam/work/w51_singledish_maps/plot_scripts/pubfiguresrc')

h2co_densfile = paths.dpath("H2CO_ParameterFits_weighted_mean_mean_density_stdmasked.fits")
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
background = 4e21
bgsub = (column_image.data-background)
bgsub *= (bgsub > 0)
nr_all, centers_all, rprof_all = radialprofile.azimuthalAverage(bgsub, center=(x,y),
                                                                return_nr=True)
# Use H2CO mask so that we measure profile in same apertures
result = radialprofile.azimuthalAverage(bgsub, center=(x,y),
                                        return_nr=True, mask=mask)
nr, centers, rprof = result

distprof = centers*(pixsize*distance).to(u.pc,u.dimensionless_angles())
rprof_all = rprof_all*u.cm**-2
rprof = rprof*u.cm**-2
h2co_distprof = h2co_centers*(h2co_pixsize*distance).to(u.pc,u.dimensionless_angles())
h2co_rprof = h2co_rprof*u.cm**-3

muh2 = 2.8 # AMU per particle
radOK = (distprof < 25*u.pc) & (distprof > (pixsize*distance).to(u.pc, u.dimensionless_angles()))
massprof = (rprof*constants.m_p*muh2*nr*(pixsize*distance)**2).to(u.M_sun, u.dimensionless_angles())
massprof_all = (rprof_all*constants.m_p*muh2*nr_all*(pixsize*distance)**2).to(u.M_sun, u.dimensionless_angles())
mean_rprof = np.cumsum(rprof) / np.cumsum(nr)

# Truncate profiles
distprof,massprof = distprof[radOK], massprof[radOK]
massprof_all = massprof_all[radOK]
dr = np.mean(np.diff(distprof))

cumul_massprof = np.cumsum(massprof)
dens_profile = (cumul_massprof.to(constants.m_p)/(4/3.*np.pi*distprof.to(u.cm)**3)) / muh2

fig = pl.figure(1)
fig.clf()
ax1 = pl.subplot(2,1,1)
ax1.plot(distprof.value, massprof.value)
ax1.set_ylabel("Mass per shell")

ax2 = pl.subplot(2,1,2, sharex=ax1)
ax2.plot(distprof.value, cumul_massprof.value/1e3)
ax2.set_xlabel("Radius (pc)")
ax2.set_ylabel("Enclosed Mass ($10^3 M_\odot$)")

pl.subplots_adjust(hspace=0)
# from http://matplotlib.org/examples/pylab_examples/subplots_demo.html
pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
ax2.set_xlim(distprof.value.min(), distprof.value.max())

pl.savefig(paths.fpath("HiGal_RadialProfile_Mass.pdf"))


rad_as, sigma1, sigma2, sigma3, sigma_all = np.loadtxt(paths.dpath('vwidth_vs_r.txt')).T
rad_as *= u.arcsec
rad_pc = (rad_as*distance).to(u.pc, u.dimensionless_angles())
fwhm13co1 = sigma1*(8*np.log(2))**0.5*u.km/u.s
fwhm13co2 = sigma2*(8*np.log(2))**0.5*u.km/u.s
fwhm13co3 = sigma3*(8*np.log(2))**0.5*u.km/u.s
fwhm13co_all = sigma_all*(8*np.log(2))**0.5*u.km/u.s

fig1a = pl.figure(4)
fig1a.clf()
ax1 = fig1a.add_subplot(2,1,1)
ax1.plot(distprof.value, cumul_massprof.value/1e3)
ax1.set_ylabel("Enclosed Mass ($10^3 M_\odot$)")

ax2 = fig1a.add_subplot(2,1,2, sharex=ax1)
vesc = ((cumul_massprof * 2 * constants.G / distprof)**0.5).to(u.km/u.s)
ax2.plot(distprof.value, vesc.value, label='Escape Speed')
ax2.plot(rad_pc.value, fwhm13co_all.value, label='1D FWHM Velocity Dispersion')
ax2.plot(rad_pc.value, fwhm13co1.value)
ax2.plot(rad_pc.value, fwhm13co2.value)
ax2.plot(rad_pc.value, fwhm13co3.value)
ax2.set_xlabel("Radius (pc)")
ax2.set_ylabel("Speed (km/s)")
ax2.legend(loc='lower right', fontsize=18)
ax2.set_xlim(0,15)
pl.setp( ax1.get_xticklabels(), visible=False)
fig1a.subplots_adjust(hspace=0)

crossingpt = np.argmin(np.abs(vesc.value-np.interp(distprof.value, rad_pc.value, fwhm13co2.value)))
ylim1 = ax1.get_ylim()
ylim2 = ax2.get_ylim()
ax1.hlines([1e1,1e2], *ax1.get_xlim(), colors='k', linestyle=':', alpha=0.5, linewidth=1, zorder=0)
ax1.set_ylim(*ylim1)
ax2.set_ylim(*ylim2)

pl.savefig(paths.fpath("HiGal_RadialProfile_Mass_VelocityDispersion.pdf"))
ax1.vlines(distprof.value[crossingpt], *ylim1, linestyle='--', alpha=0.5, linewidth=1, zorder=0)
ax2.vlines(distprof.value[crossingpt], *ylim2, linestyle='--', alpha=0.5, linewidth=1, zorder=0)
ax1.set_ylim(*ylim1)
ax2.set_ylim(*ylim2)
pl.savefig(paths.fpath("HiGal_RadialProfile_Mass_VelocityDispersion_mark_levels.pdf"))
ax1.fill_between(ax1.get_xlim(), [10,10], [30,30], color='k', alpha=0.2)
pl.savefig(paths.fpath("HiGal_RadialProfile_Mass_VelocityDispersion_highlightMPC.pdf"))



fig = pl.figure(2)
fig.clf()
ax1 = pl.subplot(2,1,1)
shell_volume = (nr/nr_all.astype('float'))[radOK] * 4/3.*np.pi*((distprof+dr/2)**3-(distprof-dr/2)**3)
shell_density = (massprof.to(constants.m_p)/muh2/shell_volume.to(u.cm**3))
shell_volume_all = 4/3.*np.pi*((distprof+dr/2)**3-(distprof-dr/2)**3)
shell_density_all = (massprof_all.to(constants.m_p)/muh2/shell_volume_all.to(u.cm**3))
ax1.semilogy(distprof.value, shell_density.value, label='Masked')
ax1.semilogy(distprof.value, shell_density_all.value, label='Total')
h2co_line, = ax1.semilogy(h2co_distprof.value, 10**h2co_rprof.value, label='H$_2$CO')
ax1.set_ylabel("Shell-averaged $n(H_2)$")
ax1.legend(loc='lower left', fontsize=18)

# What is the cumulative density profile of formaldehyde if we assume the gas is volume-filling?
h2co_dens_filling = (10**h2co_rprof.value).cumsum() / h2co_nr.cumsum()

ax2 = pl.subplot(2,1,2, sharex=ax1)
ax2.semilogy(distprof.value, dens_profile.value)
matchloc = 30 # ~6 pc
ax2.plot(distprof.value, distprof.value**-2 / distprof[matchloc].value**-2 * dens_profile[matchloc].value, '--', label='$R^{-2}$')
ax2.plot(distprof.value, distprof.value**-2.5 / distprof[matchloc].value**-2.5 * dens_profile[matchloc].value, '--',
         dashes=(15,10), color='#009955', label='$R^{-2.5}$')
#ax2.plot(distprof.value, distprof.value**-1.5 / distprof[matchloc].value**-1.5 * dens_profile[matchloc].value, ':')
ax2.legend(loc='best')
ax2.set_xlabel("Radius (pc)")
ax2.set_ylabel("Volume-Averaged")
ax2.set_xlim(distprof.value.min(), distprof.value.max())
ax2.set_xlim(0, 15)
ax2.set_ylim(dens_profile.value[np.argmin(np.abs(distprof.value-15))], dens_profile.value[0], )

pl.subplots_adjust(hspace=0)
# from http://matplotlib.org/examples/pylab_examples/subplots_demo.html
pl.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
pl.savefig(paths.fpath("HiGal_RadialProfile_VolumeDensity.pdf"))
ax2.semilogy(h2co_distprof.value, h2co_dens_filling, color=h2co_line.get_color())
pl.savefig(paths.fpath("HiGal_RadialProfile_VolumeDensity_h2cofilling.pdf"))



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
