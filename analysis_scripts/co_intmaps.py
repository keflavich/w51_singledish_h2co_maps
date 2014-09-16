from common_constants import (vrange1, vrange2, cotocol, h2togm, distance,
                             h2co22freq, gbbeamarea,)
import spectral_cube
from spectral_cube import SpectralCube,BooleanArrayMask
from build_mask import cubemask,includemask
import numpy as np
import paths
from astropy import units as u
from astropy import constants
from astropy import coordinates
from astropy.table import Table,Column
from astropy.io import fits
from astropy import log
import copy
import os
import warnings

grsfilepath = paths.dpath2('grs_48and50_cube.fits')
grsSSfilepath = paths.dpath2('grs_48and50_cube_supersampledh2cogrid.fits')
if not (os.path.exists(grsfilepath) and os.path.exists(grsSSfilepath)):
    warnings.warn("{0} not found.  Running reduction/stitch_grs_cubes.py")
    execfile(os.path.join(paths.source_root,'reduction/stitch_grs_cubes.py'))

cube13 = spectral_cube.SpectralCube.read(grsfilepath)

cube13_slab1 = cube13.spectral_slab(vrange1[0],vrange1[1])
cube13_slab2 = cube13.spectral_slab(vrange2[0],vrange2[1])
cube13_slab3 = cube13.spectral_slab(vrange1[0],vrange2[1])

snmask1 = cube13_slab1 > 0.5
snmask2 =  cube13_slab2 > 0.5
snmask3 =  cube13_slab3 > 0.5

cube13_slab1_masked_mom0 = cube13_slab1.with_mask(snmask1).moment0()
cube13_slab2_masked_mom0 = cube13_slab2.with_mask(snmask2).moment0()
cube13_slab3_masked_mom0 = cube13_slab3.with_mask(snmask3).moment0()

cube13ss = spectral_cube.SpectralCube.read(grsSSfilepath)

# Try masking based on stddev: uncertainty of 1 order of magnitude isn't super interesting...
stdcube = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_stddens.fits"))
okmask = BooleanArrayMask(np.isfinite(stdcube.filled_data[:]), wcs=stdcube.wcs)
stdcube = stdcube.with_mask(okmask)
goodmask_std = stdcube < 0.5
denscube_mean = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_likewtddens.fits")).with_mask(okmask)
denscube_min = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_mindens.fits"))
denscube_max = spectral_cube.SpectralCube.read(paths.dpath("H2CO_ParameterFits_maxdens.fits"))
high5e4dens = denscube_mean > np.log10(5e4)
high1e4dens = denscube_mean > np.log10(1e4)
high5e3dens = denscube_mean > np.log10(5e3)
high1e3dens = denscube_mean > np.log10(1e3)

# WCSes are offset: HACK!
cube13ss._wcs = stdcube.wcs

cube13ss_slab1 = cube13ss.spectral_slab(vrange1[0],vrange1[1])
cube13ss_slab2 = cube13ss.spectral_slab(vrange2[0],vrange2[1])
cube13ss_slab3 = cube13ss.spectral_slab(vrange1[0],vrange2[1])

total_co = cube13ss.sum()
total_co_slab = cube13ss_slab3.sum()

m1 = cube13ss_slab1 > 0.5
m2 = cube13ss_slab2 > 0.5
m3 = cube13ss_slab3 > 0.5

cube13ss_slab1_masked_mom0 = cube13ss_slab1.with_mask(m1).moment0()
cube13ss_slab2_masked_mom0 = cube13ss_slab2.with_mask(m2).moment0()
cube13ss_slab3_masked_mom0 = cube13ss_slab3.with_mask(m3).moment0()

cube13ss_slab1_mom0 = cube13ss_slab1.moment0()
cube13ss_slab2_mom0 = cube13ss_slab2.moment0()
cube13ss_slab3_mom0 = cube13ss_slab3.moment0()

pixsize = np.abs(np.prod((cube13ss_slab1_mom0.wcs.wcs.get_cdelt() * cube13ss_slab1_mom0.wcs.wcs.get_pc().diagonal())))**0.5 * u.deg
pixarea = ((distance * pixsize)**2).to(u.cm**2, u.dimensionless_angles())
log.info("Mass in slab 1 ({v1} to {v2} km/s): {m:0.5g} Msun".format(v1=vrange1[0], v2=vrange1[1],
                                                            m=(cube13ss_slab1_mom0*h2togm*pixarea*cotocol()).sum().to(u.M_sun).value))
log.info("Mass in slab 2 ({v1} to {v2} km/s): {m:0.5g} Msun".format(v1=vrange2[0], v2=vrange2[1],
                                                            m=(cube13ss_slab2_mom0*h2togm*pixarea*cotocol()).sum().to(u.M_sun).value))

high5e4dens_co = cube13ss.with_mask(high5e4dens)
high5e4dens_co_slab3 = high5e4dens_co.spectral_slab(vrange1[0],vrange2[1])
high5e4dens_co_slab3_mom0 = high5e4dens_co_slab3.moment0()

high1e4dens_co = cube13ss.with_mask(high1e4dens)
high1e4dens_co_slab3 = high1e4dens_co.spectral_slab(vrange1[0],vrange2[1])
high1e4dens_co_slab3_mom0 = high1e4dens_co_slab3.moment0()

high5e3dens_co = cube13ss.with_mask(high5e3dens)
high5e3dens_co_slab3 = high5e3dens_co.spectral_slab(vrange1[0],vrange2[1])
high5e3dens_co_slab3_mom0 = high5e3dens_co_slab3.moment0()

high1e3dens_co = cube13ss.with_mask(high1e3dens)
high1e3dens_co_slab3 = high1e3dens_co.spectral_slab(vrange1[0],vrange2[1])
high1e3dens_co_slab3_mom0 = high1e3dens_co_slab3.moment0()

dgmf1e4 = fits.PrimaryHDU(data=(high1e4dens_co_slab3_mom0 /
                                cube13ss_slab3_mom0).value,
                          header=high1e4dens_co_slab3_mom0.hdu.header)
dgmf5e4 = fits.PrimaryHDU(data=(high5e4dens_co_slab3_mom0 /
                                cube13ss_slab3_mom0).value,
                          header=high5e4dens_co_slab3_mom0.hdu.header)
dgmf5e3 = fits.PrimaryHDU(data=(high5e3dens_co_slab3_mom0 /
                                cube13ss_slab3_mom0).value,
                          header=high5e3dens_co_slab3_mom0.hdu.header)
dgmf1e3 = fits.PrimaryHDU(data=(high1e3dens_co_slab3_mom0 /
                                cube13ss_slab3_mom0).value,
                          header=high1e3dens_co_slab3_mom0.hdu.header)

# Purely detection-based thresholding
h2co11 = SpectralCube.read(paths.dpath('W51_H2CO11_cube_supersampled_sub.fits'))
h2co22 = SpectralCube.read(paths.dpath('W51_H2CO22_pyproc_cube_lores_supersampled_sub.fits'))

noisevrange = [-50,0]*u.km/u.s
h2co11noiseslab = h2co11.spectral_slab(*noisevrange)
h2co22noiseslab = h2co22.spectral_slab(*noisevrange)
h2co11noise = h2co11noiseslab.apply_numpy_function(np.std,axis=0)
h2co22noise = h2co22noiseslab.apply_numpy_function(np.std,axis=0)

# TODO: implement Cube (lazy?) Arithmetic: avoid filling data!
sn11 = SpectralCube(np.abs(h2co11.filled_data[:])/h2co11noise, wcs=h2co11.wcs)
sn22 = SpectralCube(np.abs(h2co22.filled_data[:])/h2co22noise, wcs=h2co22.wcs)
sn11_slab3 = sn11.spectral_slab(vrange1[0], vrange2[1])
sn22_slab3 = sn22.spectral_slab(vrange1[0], vrange2[1])

h2co11_total = np.array([cube13ss_slab3.with_mask(sn11_slab3>threshold).sum().value
                         for threshold in np.arange(1,6)])
h2co22_total = np.array([cube13ss_slab3.with_mask(sn22_slab3>threshold).sum().value
                         for threshold in np.arange(1,6)])
h2co11_fraction = h2co11_total / total_co_slab.value
h2co22_fraction = h2co22_total / total_co_slab.value
both_fraction = cube13ss_slab3.with_mask((sn11_slab3 > 2) &
                                         (sn22_slab3 > 2)).sum().value / total_co_slab.value

log.info("H2CO Fractions: 1-1: {0}".format(h2co11_fraction))
log.info("H2CO Fractions: 2-2: {0}".format(h2co22_fraction))
log.info("Both: {0}".format(both_fraction))


w51a_corners = denscube_mean.wcs.dropaxis(denscube_mean.wcs.wcs.spec).wcs_world2pix([49.4,49.6],[-0.5,-0.3],0)
w51b_corners = denscube_mean.wcs.dropaxis(denscube_mean.wcs.wcs.spec).wcs_world2pix([48.0,49.4],[-0.5, 0.0],0)
region_slices = {'All':(slice(None),slice(None),slice(None)),
                 'W51 Main':(
                             slice(None),
                             slice(w51a_corners[1][0],w51a_corners[1][1]), # y
                             slice(w51a_corners[0][1],w51a_corners[0][0]), # x
                 ),
                 'W51 B':(slice(None),
                          slice(w51b_corners[1][0],w51b_corners[1][1]), # y
                          slice(w51b_corners[0][1],w51b_corners[0][0]), # x
                          ),
                }

tables = {}

for region in region_slices:
    tbl = Table()
    tbl.meta['Region'] = region
    slices = region_slices[region]

    total_co_reg = cube13ss[slices].sum()
    total_co_slab_reg = cube13ss_slab3[slices].sum()

    tbl.add_column(Column(data=np.linspace(3.0, 6), name="Density Threshold", dtype='float'))
    for stat in ('mean','min','max'):
        denscube = eval('denscube_{0}'.format(stat))[slices].spectral_slab(vrange1[0], vrange2[1])
        tbl.add_column(Column(data=[cube13ss_slab3[slices].with_mask(denscube>threshold).sum().value
                                    for threshold in tbl['Density Threshold']],
                              dtype='float',
                              name='{0} Flux Over Threshold'.format(stat)))
        tbl.add_column(Column(data=tbl['{0} Flux Over Threshold'.format(stat)] / total_co_reg.value,
                              dtype='float', name='{0} F_total'.format(stat)))
        tbl.add_column(Column(data=tbl['{0} Flux Over Threshold'.format(stat)] / total_co_slab_reg.value,
                              dtype='float', name='{0} F_total(slab)'.format(stat)))

    pkfrac = cube13ss.with_mask(cubemask).spectral_slab(vrange1[0], vrange2[1])[slices].sum().value / total_co_slab_reg.value
    tbl.meta['Peak Fraction'] = pkfrac

    tables[region] = tbl


if __name__ == "__main__":
    from aplpy_figure_maker import FITSFigure
    import aplpy
    import pylab as pl
    fig = pl.figure(1, figsize=(12,12))
    fig.clf()
    F1 = FITSFigure(cube13ss_slab3_masked_mom0.hdu, subplot=[0.05,0.5,0.4,0.4], figure=fig)
    F1.show_grayscale()
    F1._ax1.set_title("$^{13}$CO masked with H$_2$CO")
    F2 = FITSFigure(cube13_slab3_masked_mom0.hdu, subplot=[0.50,0.5,0.4,0.4], figure=fig)
    F2.show_grayscale()
    F2._ax1.set_title("$^{13}$CO masked by S/N")
    F3 = FITSFigure(high5e4dens_co_slab3_mom0.hdu, subplot=[0.05,0.05,0.4,0.4], figure=fig)
    F3.show_grayscale()
    F3._ax1.set_title("$^{13}$CO masked by $n>5\\times10^4$")
    F4 = FITSFigure(high1e4dens_co_slab3_mom0.hdu, subplot=[0.50,0.05,0.4,0.4], figure=fig)
    F4.show_grayscale()
    F4._ax1.set_title("$^{13}$CO masked by $n>1\\times10^4$")

    pl.matplotlib.rc_file('pubfiguresrc')

    for region in tables:
        tbl = tables[region]
        pk_frac = tbl.meta['Peak Fraction']
        fig = pl.figure(2)
        fig.clf()
        pl.semilogx(10**tbl['Density Threshold'], tbl['mean F_total(slab)'])
        #pl.semilogx(10**tbl['Density Threshold'], tbl['min F_total(slab)'])
        #pl.semilogx(10**tbl['Density Threshold'], tbl['max F_total(slab)'])
        pl.fill_between(10**tbl['Density Threshold'], np.nan_to_num(tbl['min F_total(slab)']),
                        tbl['max F_total(slab)'], color='b', alpha=0.25)
        pl.vlines(10**4, 0,  pk_frac*1.05, color='k', linestyle='--')
        pl.hlines(pk_frac, 1e3, 1e6, color='k', linestyle='--')
        pl.xlabel("Density Threshold (cm$^{-3}$)")
        pl.ylabel("Fraction of cloud above density")
        pl.title(region)
        pl.ylim(0, pk_frac*1.05)
        pl.savefig(paths.fpath("FractionOfMassAboveDensity_{0}.pdf".format(region.replace(" ",""))))

    def show_con(dgmf, fignum=3, label="$f(n>10^4$ cm$^{-3})$"):
        figN = pl.figure(fignum)
        figN.clf()
        gray = copy.copy(pl.cm.gray)
        gray.set_bad('black')
        gray.set_under('black')
        FF = FITSFigure(cube13ss_slab3_masked_mom0.hdu, figure=figN, colorbar=False,
                        cmap=gray, transparent_nan=False)
        #FF = FITSFigure(dgmf1e4, figure=fig4)
        #FF.show_grayscale()
        FF.show_contour(dgmf, levels=[0.1,0.3,0.5,0.7,0.9], zorder=1000, smooth=1)
        FF.add_colorbar()
        cax = FF.colorbar._colorbar_axes
        cax.collections[0].set_visible(False)
        FF.colorbar = pl.colorbar(FF._layers['contour_set_1'], cax=cax)
        FF.colorbar.set_label("Dense Fraction", rotation=270, labelpad=30)
        FF._ax1.set_title(label)
        #FF.scalebar._scalebar.set_zorder(40)
        #FF.refresh()
        return FF

    F5 = show_con(dgmf5e4, 3, "$f(n>5\\times10^4$ cm$^{-3})$")
    F5.save(paths.fpath('DGMF_5e4_Contours_on_13CO.pdf'))
    F6 = show_con(dgmf1e4, 4, "$f(n>1\\times10^4$ cm$^{-3})$")
    F6.save(paths.fpath('DGMF_1e4_Contours_on_13CO.pdf'))
    F7 = show_con(dgmf5e3, 5, "$f(n>5\\times10^3$ cm$^{-3})$")
    F8 = show_con(dgmf1e3, 6, "$f(n>1\\times10^3$ cm$^{-3})$")

    # Make a "star forming mass map"
    # use some simple assumptions for CO
    co_to_mass = ((8e20*u.cm**-2 * (15*u.arcsec * 5.1*u.kpc).to(u.pc,
                                                               u.dimensionless_angles())**2
                  * constants.m_p * 2.8).to(u.M_sun)).value
    co_to_mass_surf = ((8e20*u.cm**-2 * constants.m_p *
                        2.8).to(u.M_sun/u.pc**2)).value
    totalmassmap = (cube13ss_slab3_masked_mom0.to(u.K*u.km/u.s).value *
                    co_to_mass_surf)
    sfmassmap = (dgmf1e4.data *
                 cube13ss_slab3_masked_mom0.to(u.K*u.km/u.s).value *
                 co_to_mass_surf)
    badsfmass = ~np.isfinite(sfmassmap)
    sfmassmap[badsfmass] = 0
    sfmasshdu = fits.PrimaryHDU(sfmassmap, header=dgmf1e4.header)

    # Make an SFR map from the radio continuum data
    cont2cm = fits.getdata(paths.cont2cm)
    hdr2cm = fits.getheader(paths.cont2cm)
    # Use Murphy 2011 (2011ApJ...737...67M) calibration: eqn 11
    te = 7.5e3*u.K
    nu = h2co22freq
    sfrperergshz = (4.6e-28 *u.Msun/u.yr * (te/(1e4*u.K))**-0.45 *
                    (nu.to(u.GHz).value)**0.1) / u.erg / u.s**-1 / u.Hz**-1
    pixel_area = (hdr2cm['CDELT2']*u.deg*distance).to(u.cm, u.dimensionless_angles())**2
    sfr2cm = ((cont2cm*u.K).to(u.Jy, u.brightness_temperature(gbbeamarea,
                                                              h2co22freq)) *
              distance**2 * 
              sfrperergshz).to(u.M_sun/u.yr)
    sfr2cmd = (sfr2cm / pixel_area).to(u.M_sun/u.yr/u.kpc**2)
    ok2cm = cont2cm > 0.1



    from astroquery.vizier import Vizier
    catalog_list = Vizier.find_catalogs('Kang W51')
    ysos = Vizier(row_limit=1e6).query_region(coordinates.SkyCoord.from_name('W51'),
                                              radius=1*u.deg,
                                              catalog=catalog_list.keys())
    cl1 = ((ysos[0]['Cl1'] == 'I') |
           (ysos[0]['Cl1'] == 'F') |
           (ysos[0]['Cl2'] == 'I') |
           (ysos[0]['Cl2'] == 'F'))
    cl1coords = coordinates.SkyCoord(ysos[0][cl1]['_RAJ2000'],
                                     ysos[0][cl1]['_DEJ2000'],
                                     frame='fk5').galactic
    fig7 = pl.figure(7)
    fig7.clf()
    gray = copy.copy(pl.cm.gray_r)
    gray.set_bad('white')
    gray.set_under('white')
    FMM = FITSFigure(sfmasshdu, figure=fig7, cmap=gray, stretch='log', vmax=99.95)
    FMM.colorbar._colorbar.set_label("Star-Forming Gas Surface Density\n[$M_{\odot}$ pc$^{-2}$",
                                     rotation=270, labelpad=50)
    FMM.show_markers(cl1coords.l, cl1coords.b, marker='x', edgecolor='r', zorder=1100)
    FMM.save(paths.fpath('StarFormingMassMap_Kang2009ClassIYSOs.pdf'))

    fig8 = pl.figure(8)
    fig8.clf()
    ax8 = fig8.gca()
    ax8.loglog(totalmassmap[~badsfmass & ok2cm], sfr2cmd[~badsfmass & ok2cm], 'r.', alpha=0.1)
    ax8.loglog(sfmassmap[~badsfmass & ok2cm], sfr2cmd[~badsfmass & ok2cm], 'k.', alpha=0.1)
    ax8.loglog([1,1e4],np.array([1,1e4])*1.2e-8*1e6,'b--', alpha=0.5)
    ax8.loglog([700,2900],[2,900],'g:', alpha=0.7)
    ax8.set_xlabel("Star Forming Gas Surface Density ($M_{\\odot}$ pc$^{-2}$)")
    ax8.set_ylabel("2 cm continuum SFR ($M_{\\odot}$ yr$^{-1}$ kpc$^{-2}$)")
    ax8.set_ylim(0.5,1e3)
    fig8.savefig(paths.fpath('ksrelation_2cm_densegas.pdf'))
