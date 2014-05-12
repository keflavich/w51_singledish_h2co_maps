from astropy.wcs import WCS
from astropy import units as u
import wcsaxes
from astropy import coordinates
import numpy as np
import pylab as pl
import paths
import os

cubefiles = ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
             'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits',
             'grs_48and50_cube_supersampledh2cogrid.fits')

if 'colorpvs' in locals():
    pvs = colorpvs['cyan']

else:
    from pvdiagrams import get_pvs
    import pyregion
    pvs = []
    endpoints_wcs = pyregion.open('pvendpoints.reg')
    color = 'cyan'
    coords = np.array([s.coord_list for s in endpoints_wcs if s.attr[1]['color'] == color])
    endpoints = coordinates.Galactic(coords[:,0],coords[:,1], unit=(u.deg,u.deg))
    for cubefn in cubefiles:
        pvs[cubefn] = get_pvs(cubefn, endpoints)


fig = pl.figure(11)
fig.clf()
cm = pl.cm.gray_r
aspects = []

for ii,cfn in enumerate(cubefiles):
    pvhdu, npts, velo, dv = pvs[cfn]

    wcs = WCS(pvhdu.header)
    #no works:
    # wcs.wcs.cunit[wcs.wcs.spec] = 'km/s'
    # wcs.wcs.crval[wcs.wcs.spec] = wcs.wcs.crval[wcs.wcs.spec] / 1e3
    # print wcs.wcs.cdelt
    # wcs.wcs.cdelt[wcs.wcs.spec] = wcs.wcs.cdelt[wcs.wcs.spec] / 1e3
    # print wcs.wcs.cdelt

    ax = wcsaxes.WCSAxesSubplot(fig, 3,1, ii+1, wcs=wcs)
    fig.add_axes(ax)

    ax.imshow(pvhdu.data, cmap=cm)
    ax.set_ylim(*wcs.wcs_world2pix([0,0],[40e3,80e3],0)[1])

    aspects.append(np.diff(ax.get_ylim()) / np.diff(ax.get_xlim()))

fig.axes[-1].set_aspect(aspects[0]/2)

import aplpy
fig = pl.figure(12)
fig.clf()

ffs = []
hdus = []
for ii,cfn in enumerate(cubefiles):
    pvhdu, npts, velo, dv = pvs[cfn]
    hdus.append(pvhdu)

    F = aplpy.FITSFigure(pvhdu, subplot=[0.1, 0.05+(0.9-0.3*(ii+1)), 0.9, 0.3], figure=fig)
    F.show_grayscale(invert=True)

    F.recenter(0.325, 65000, width=0.65, height=30000)

    ffs.append(F)
ffs[2].show_contour(hdus[0])

F.save(os.path.join(paths.figurepath, "wide_PV_h2co11_22_13co.pdf"))
pl.show()
