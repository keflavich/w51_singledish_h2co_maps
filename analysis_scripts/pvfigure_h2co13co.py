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

cubefiles = ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
             'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits',
             'grs_48and50_cube_supersampledh2cogrid.fits')

if 'colorpvs' in locals():
    pvs = colorpvs['cyan']

else:
    from pvdiagrams import get_pvs
    import pyregion
    pvs = {}
    endpoints_wcs = pyregion.open('pvendpoints.reg')
    color = 'cyan'
    coords = np.array([s.coord_list for s in endpoints_wcs if s.attr[1]['color'] == color])
    endpoints = coordinates.Galactic(coords[:,0],coords[:,1], unit=(u.deg,u.deg))
    for cubefn in cubefiles:
        pvs[cubefn] = get_pvs(cubefn, endpoints)

def pvplots(pvs, color='cyan', extranumber=0, width=0.85):
    fig11 = pl.figure(11+extranumber)
    fig11.clf()
    cm = pl.cm.gray_r
    #aspects = []

    hdus = []
    for ii,cfn in enumerate(cubefiles):
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
        ax.set_ylim(*wcs.wcs_world2pix([0,0],[40e3,80e3],0)[1])
        #ax.set_yticklabels(ax.get_yticks()/1e3)

        if ii<2:
            ax.coords[0].ticklabels.set_visible(False)
            #ax.coords[1].set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x: x/1e3)) # notimplemented =(
            for k in ax.coords[1].ticklabels.text:
                text = [str(int(tt)/1000) for tt in ax.coords[1].ticklabels.text[k]]
                ax.coords[1].ticklabels.text[k] = text
                ax.coords[1].ticklabels.set_text(text)
            ax.set_xticklabels(["" for x in ax.get_xticklabels()])

        #aspects.append(np.diff(ax.get_ylim()) / np.diff(ax.get_xlim()))

    fig11.subplots_adjust(hspace=0)
    ax.set_xlabel("Offset (degrees)")
    fig11.axes[1].set_ylabel("Radio $V_{LSR}$ (m/s)")
    fig11.axes[0].set_xlabel("Offset (pc)")
    fig11.axes[0].coords[0].set_ticks(spacing=(10*u.pc/(5.1*u.kpc)).to(u.deg, u.dimensionless_angles()).value)

    # the -0.5 is just to force the colors to the high-end
    con = ax.contour(hdus[0].data, levels=[-0.5,0.025,0.05,0.1,0.2,0.4])
    fig11.savefig(os.path.join(paths.figurepath, "filaments_H2CO_13CO_pvslice.pdf"))
    for cc in ax.collections:
        cc.set_visible(False)

    mask = hdus[2].data > 0.2
    rfig = hdus[0].data / hdus[1].data
    rfig[~mask] = np.nan
    rhdu = pvhdu.copy()
    rhdu.data = rfig
    c = ax.contourf(rfig, colors=[(1,0,0.5,0.3), (1,0,0,0.5), (1.0,0.33,0,0.5,),
                              (0.75,0.6,0,0.5), (0.6,0.9,0,0.3), (0.3,0.95,0,0.3)],
                levels=[0,1,2,5,10,15,35])
    fig11.savefig(os.path.join(paths.figurepath, color+"_filaments_H2CO_13CO_density_pvslice.pdf"))
    A = matplotlib.axes.Axes(fig11, [0.9,0.1,0.025,0.25])
    pl.colorbar(c, cax=A)
    fig11.add_axes(A)
    fig11.savefig(os.path.join(paths.figurepath, color+"_filaments_H2CO_13CO_density_pvslice_colorbar.pdf"))



    #fig.axes[-1].set_aspect(aspects[0]/2)

    fig = pl.figure(12+extranumber)
    fig.clf()

    ffs = []
    for ii,cfn in enumerate(cubefiles):
        pvhdu, npts, velo, dv = pvs[cfn]
        hdus.append(pvhdu)

        F = aplpy.FITSFigure(pvhdu, subplot=[0.1, 0.05+(0.9-0.3*(ii+1)), 0.88, 0.3], figure=fig)
        F.show_grayscale(invert=True)

        F.recenter(width/2, 65000, width=width, height=30000)

        ffs.append(F)
    ffs[2].show_contour(hdus[0])

    F.save(os.path.join(paths.figurepath, color+"_wide_PV_h2co11_22_13co.pdf"))


pvplots(pvs, color='cyan', width=0.85)
pvplots(colorpvs['purple'], color='purple', extranumber=-2, width=0.45)




bgps = fits.open(os.path.join(paths.datapath_w51, 'v2.0_ds2_l050_13pca_map20.fits'))[0]
bgps_wcs = WCS(bgps.header)
column = fits.open(os.path.join(paths.datapath_w51, 'higalsedfit_70to500_l048_beta1.75N.fits'))[0]
column_wcs = WCS(column.header)

fig13 = pl.figure(13)
fig13.clf()
F = aplpy.FITSFigure(column, convention='calabretta', figure=fig13)
#F.show_grayscale(invert=True, vmax=25, stretch='log', vmin=0, vmid=-0.2)
F.show_grayscale(invert=True, vmax=1e23, stretch='log', vmin=2e21, vmid=1e21)
F.recenter(49.22, -0.33612413, width=0.85, height=0.55)
F.show_regions(os.path.join(paths.datapath_w51, 'cyan_segments.reg'))
F.show_regions(os.path.join(paths.datapath_w51, 'bluered_segments.reg'))
polygons = pvextractor.geometry.path.Path(endpoints, width=25*u.arcsec).sample_polygons(1, wcs=column_wcs)
for poly in polygons:
    F._ax1.draw_artist(Polygon(zip(poly.x, poly.y), ec='green', fc='none',
                               transform=F._ax1.transData, clip_on=True,
                               clip_box=F._ax1.bbox))
F.set_tick_labels_xformat('dd.d')
F.set_tick_labels_yformat('dd.d')
F.save(os.path.join(paths.figurepath, 'filament_extraction_region_on_HiGal.pdf'))

pl.show()
