from astropy.wcs import WCS
from astropy import units as u
from astropy.io import fits
import wcsaxes
from astropy import coordinates
import numpy as np
import pylab as pl
pl.rcParams['font.size'] = 20
import paths
import os
import matplotlib
from agpy import asinh_norm
from matplotlib.patches import Polygon
import pvextractor
import aplpy
import warnings

cubefiles = ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
             'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits',
             'grs_48and50_cube_supersampledh2cogrid.fits',
             'H2CO_ParameterFits_bestdens.fits',
             'w51_bieging_12co32.fits',
             'w51_12co10_carpenter_rightaxes.fits',
             )

if 'colorpvs' in locals():
    pvs = colorpvs['cyan']

else:
    from pvdiagrams import get_pvs
    import pyregion
    endpoints_wcs = pyregion.open(paths.rpath('pvendpoints.reg'))
    colorpvs = {}
    for color in ['cyan', 'purple']:
        pvs = {}
        coords = np.array([s.coord_list for s in endpoints_wcs if
                           s.attr[1]['color'] == color])
        endpoints = coordinates.Galactic(coords[:,0],coords[:,1], unit=(u.deg,u.deg))
        for cubefn in cubefiles:
            pvs[cubefn] = get_pvs(cubefn, endpoints)
        colorpvs[color] = pvs

def pvplots(pvs, color='cyan', coname='13CO', h2coname='h2co11_22', extranumber=0, width=0.85,
            cubefiles=cubefiles):
    fig11 = pl.figure(11+extranumber)
    fig11.clf()
    cm = pl.cm.gray_r
    #aspects = []

    hdus = []
    for ii,cfn in enumerate(cubefiles[:3]):
        pvhdu, npts, velo, dv = pvs[cfn]
        hdus.append(pvhdu)

        wcs = WCS(pvhdu.header)
        #no works:
        # wcs.wcs.cunit[wcs.wcs.spec] = 'km/s'
        # wcs.wcs.crval[wcs.wcs.spec] = wcs.wcs.crval[wcs.wcs.spec] / 1e3
        # print wcs.wcs.cdelt
        # wcs.wcs.cdelt[wcs.wcs.spec] = wcs.wcs.cdelt[wcs.wcs.spec] / 1e3
        # print wcs.wcs.cdelt

        ax = wcsaxes.WCSAxesSubplot(fig11, 3,1, ii+1, wcs=wcs)
        fig11.add_axes(ax)

        #ax.imshow(pvhdu.data, cmap=cm, norm=asinh_norm.AsinhNorm())
        vmin = np.percentile(pvhdu.data[np.isfinite(pvhdu.data)], 10)
        vmax = np.percentile(pvhdu.data[np.isfinite(pvhdu.data)], 99.95)
        ax.imshow(pvhdu.data, cmap=cm, vmin=vmin, vmax=vmax, norm=asinh_norm.AsinhNorm())
        ax.set_ylim(*wcs.wcs_world2pix([0,0],[35e3,80e3],0)[1])
        #ax.set_ylim(35e3, 80e3)
        #ax.set_ylim(35, 80)
        #ax.set_yticklabels(ax.get_yticks()/1e3)

        ax.coords[1].set_format_unit(u.km / u.s)
        if ii<2:
            ax.coords[0].ticklabels.set_visible(False)
            #ax.coords[1].set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x: x/1e3)) # notimplemented =(
            #ax.set_xticklabels(["" for x in ax.get_xticklabels()])

        #pl.draw; pl.show()

        #ax.coords[1].set_ticks(spacing=10*u.km/u.s)

        #for k in ax.coords[1].ticklabels.text:
        #    text = [str(int(tt)/1000) for tt in ax.coords[1].ticklabels.text[k]]
        #    ax.coords[1].ticklabels.text[k] = text
        #    ax.coords[1].ticklabels.set_text(text)
        #    ax.yaxis.set_ticklabels(['{0:d}'.format(int(yy)) for yy in ax.yaxis.get_ticklocs()])

        #aspects.append(np.diff(ax.get_ylim()) / np.diff(ax.get_xlim()))

    fig11.subplots_adjust(hspace=0)
    #ax.set_xlabel("Offset (degrees)")
    ax.coords[0].set_axislabel("Offset (degrees)")
    #fig11.axes[0].set_xlabel("Offset (pc)")
    #fig11.axes[0].coords[0].set_ticks(spacing=(10*u.pc/(5.1*u.kpc)).to(u.deg, u.dimensionless_angles()).value)
    #fig11.axes[0].coords[0].set_ticks(spacing=(10*u.pc/(5.1*u.kpc)).to(u.deg, u.dimensionless_angles()))
    #for ax in fig11.axes:
    #    ax.yaxis.set_ticklabels(['{0:d}'.format(int(yy)) for yy in ax.yaxis.get_ticklocs()])
    #fig11.axes[1].set_ylabel("Radio $V_{LSR}$ (m/s)")
    fig11.axes[1].coords[1].set_axislabel("Radio $V_{LSR}$ (km/s)")

    #import ipdb; ipdb.set_trace()

    # the -0.5 is just to force the colors to the high-end
    con = ax.contour(hdus[0].data, levels=[-0.5,0.025,0.05,0.1,0.2,0.4])
    fig11.savefig(os.path.join(paths.figurepath, color+"_filaments_H2CO_"+coname+"_pvslice.pdf"))
    for cc in ax.collections:
        cc.set_visible(False)

    #mask = hdus[2].data > 0.2
    # Ratio figure
    #rfig = hdus[0].data / hdus[1].data
    #rfig[~mask] = np.nan
    #rhdu = pvhdu.copy()
    #rhdu.data = rfig
    # Density Figure
    nhdu, nnpts, nvelo, ndv = pvs[cubefiles[3]]


    auto_levels = aplpy.image_util.percentile_function(hdus[0].data)
    vmin = auto_levels(0.25)
    vmax = auto_levels(99.75)
    levels = np.linspace(vmin, vmax, 5)
    c1 = ax.contour(hdus[0].data, levels=levels)
    #fig11.savefig(os.path.join(paths.figurepath, color+"_filaments_H2CO_"+coname+"_density_pvslice.pdf"))
    fig11.savefig(os.path.join(paths.figurepath, color+"_wide_PV_"+h2coname+"_"+coname+"_wcsaxes.pdf"))

    for cc in c1.collections:
        cc.set_visible(False)
    c = ax.contourf(nhdu.data, colors=[(1,0,0.1,0.7), (1,0,0.0,0.6),
                                       (1,0,0,0.5), (1.0,0.33,0,0.5,),
                                       (0.75,0.6,0,0.5), (0.6,0.9,0,0.3),
                                       (0.3,0.95,0,0.3)][::-1],
                levels=[1,2,3,4,5,6]) # density
                #levels=[0,1,2,5,10,15,35]) # ratio
    fig11.savefig(os.path.join(paths.figurepath, color+"_filaments_"+h2coname+"_"+coname+"_density_pvslice.pdf"))
    ((x0,y0),(x1,y1)) = ax.bbox._bbox.get_points()
    A = matplotlib.axes.Axes(fig11, [x1,y0,0.025,y1-y0])
    pl.colorbar(c, cax=A)
    fig11.add_axes(A)
    fig11.savefig(os.path.join(paths.figurepath, color+"_filaments_"+h2coname+"_"+coname+"_density_pvslice_colorbar.pdf"))



    #fig.axes[-1].set_aspect(aspects[0]/2)

    fig = pl.figure(12+extranumber)
    fig.clf()

    ffs = []
    for ii,cfn in enumerate(cubefiles[:3]):
        pvhdu, npts, velo, dv = pvs[cfn]
        hdus.append(pvhdu)

        F = aplpy.FITSFigure(pvhdu, subplot=[0.1, 0.05+(0.9-0.3*(ii+1)), 0.88,
                                             0.3], figure=fig)
        F.show_grayscale(invert=True)

        F.recenter(width/2, 60000, width=width, height=35000)

        ffs.append(F)
    levels = ffs[2].show_contour(hdus[0], return_levels=True)
    print levels

    # change labels AFTER overplotting contours
    for F in ffs:
        F._wcs.wcs.cunit[1] = 'km/s' # u.km/u.s
        F._wcs.wcs.cdelt[1] = F._wcs.wcs.cdelt[1]/1000.
        F._wcs.wcs.crval[1] = F._wcs.wcs.crval[1]/1000.

    ffs[0].hide_xtick_labels()
    ffs[1].hide_xtick_labels()
    ffs[0].axis_labels.hide_x()
    ffs[1].axis_labels.hide_x()
    ffs[0].axis_labels.hide_y()
    ffs[1].axis_labels.set_ytext('$V_{LSR}$ km $\mathrm{s}^{-1}$')
    ffs[2].axis_labels.hide_y()

    F.save(os.path.join(paths.figurepath, color+"_wide_PV_"+h2coname+"_"+coname+".pdf"))

    return ffs


if __name__ == "__main__":

    ffsca = pvplots(colorpvs['cyan'], coname='12CO21',
                   color='cyan', width=0.85, cubefiles=
                    ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
                     'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits',
                     'w51_bieging_12co32.fits',
                     'H2CO_ParameterFits_bestdens.fits',
                    ),
                   )

    ffscb = pvplots(colorpvs['cyan'], coname='12CO21',
                    h2coname='h2co11_13CO10',
                    color='cyan', width=0.85, cubefiles=
                    ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
                     'grs_48and50_cube_supersampledh2cogrid.fits',
                     'w51_bieging_12co32.fits',
                     'H2CO_ParameterFits_bestdens.fits',
                    ),
                   )

    ffscc = pvplots(colorpvs['cyan'], coname='12CO10',
                   color='cyan', width=0.85, cubefiles=
                    ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
                     'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits',
                     'w51_12co10_carpenter_rightaxes.fits',
                     'H2CO_ParameterFits_bestdens.fits',
                    ),
                   )

    ffscd = pvplots(colorpvs['cyan'], coname='13CO10',
                    h2coname='h2co11_13CO10',
                    color='cyan', width=0.85, cubefiles=
                    ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
                     'grs_48and50_cube_supersampledh2cogrid.fits',
                     'grs_48and50_cube_supersampledh2cogrid.fits',
                     'H2CO_ParameterFits_bestdens.fits',
                    ),
                   )

    if False:
        ffsc = pvplots(colorpvs['cyan'], color='cyan', width=0.85)
        ffsp = pvplots(colorpvs['purple'], color='purple', extranumber=-2, width=0.45)




        bgps = fits.open(os.path.join(paths.datapath_w51, 'v2.0_ds2_l050_13pca_map20.fits'))[0]
        bgps_wcs = WCS(bgps.header)
        column = fits.open(os.path.join(paths.datapath_w51, 'HIGAL_W51_mosaic_fit_160to500_N.fits'))[0]
        column_wcs = WCS(column.header)

        fig13 = pl.figure(13)
        fig13.clf()
        F = aplpy.FITSFigure(column, convention='calabretta', figure=fig13)
        #F.show_grayscale(invert=True, vmax=25, stretch='log', vmin=0, vmid=-0.2)
        F.show_grayscale(invert=True, vmax=1e23, stretch='log', vmin=2e21, vmid=1e21)
        F.recenter(49.22, -0.33612413, width=0.85, height=0.55)
        F.show_regions(os.path.join(paths.datapath_w51, 'cyan_segments.reg'))
        F.show_regions(os.path.join(paths.datapath_w51, 'bluered_segments.reg'))

        pve_path = pvextractor.geometry.path.Path

        warnings.warn("This next step may cause an error on some matplotlib backends,"
                      "e.g. MacOSX.  Try another backend, e.g. Qt4Agg")

        F.remove_layer('region_set_1')
        F.remove_layer('region_set_2')

        for color in ['cyan', 'purple']:
            pvs = {}
            coords = np.array([s.coord_list for s in endpoints_wcs if
                               s.attr[1]['color'] == color])
            endpoints = coordinates.Galactic(coords[:,0],coords[:,1], unit=(u.deg,u.deg))
            polygons = pve_path(endpoints, width=25*u.arcsec).sample_polygons(1,
                                                                              wcs=column_wcs)
            for poly in polygons:
                F._ax1.add_patch(Polygon(zip(poly.x, poly.y), ec=color, fc='none',
                                           transform=F._ax1.transData, clip_on=True,
                                           clip_box=F._ax1.bbox, zorder=10))
            #    F._ax1.draw_artist(Polygon(zip(poly.x, poly.y), ec='green', fc='none',
            #                               transform=F._ax1.transData, clip_on=True,
            #                               clip_box=F._ax1.bbox, zorder=10))
        F.set_tick_labels_xformat('dd.d')
        F.set_tick_labels_yformat('dd.d')
        F.refresh()
        F.save(os.path.join(paths.figurepath, 'filament_extraction_region_on_HiGal.pdf'))

        pl.show()
