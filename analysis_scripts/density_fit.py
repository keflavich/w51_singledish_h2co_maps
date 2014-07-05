from paths import datapath,dpath,rpath
from astropy.convolution import convolve
from astropy.io import fits
from astropy import wcs
import pyregion
import numpy as np

from astropy.utils.console import ProgressBar
import itertools
from pyspeckit import parallel_map

from h2co_modeling.grid_fitter import grid_2p_getmatch

from load_pyspeckit_cubes import (both, T, F, cont11, cont22, h2co11filename,
                                  cube1, cube2, texgrid1,  taugrid1,  texgrid2,
                                  taugrid2,  hdr, fit_a_pixel)
from common_constants import TCMB

### copied from tau_ratio_cube
# These are not used for fitting, only for masking!
# The absorption lines themselves (from cube1,cube2 above) are used for the
# fitting
h2co11 = fits.getdata(dpath('W51_H2CO11_taucube_supersampled.fits'))
h2co22 = fits.getdata(dpath('W51_H2CO22_pyproc_taucube_lores_supersampled.fits'))

noise11 = h2co11[:50,:,:].std(axis=0)
noise22 = h2co22[:50,:,:].std(axis=0)

sn11 = h2co11/noise11
sn22 = h2co22/noise22

includemask = (sn11 > 2) & (sn22 > 2)
ratio = h2co11/h2co22
ratio[True-includemask] = np.nan

filt = np.ones([3,3,3],dtype='bool')
filt[1,1,1] = 0
nneighbors = convolve(np.isfinite(ratio), filt)
ratio[(nneighbors<7) + (True-np.isfinite(nneighbors))] = np.nan
includemask[(nneighbors<7) + (True-np.isfinite(nneighbors))] = False
# End masking
###



data_iterator = [np.array(cube1.cube[includemask]),
                 np.array(cube1.errorcube[includemask]),
                 np.repeat(cont11[None,:],cube1.shape[0],axis=0).reshape(cube1.cube.shape)[includemask],
                 np.array(cube2.cube[includemask]),
                 np.array(cube2.errorcube[includemask]),
                 np.repeat(cont22[None,:],cube2.shape[0],axis=0).reshape(cube2.cube.shape)[includemask],]

args = zip(*data_iterator)[0]
print(args,fit_a_pixel(args))

pb = ProgressBar(len(data_iterator[0]))

#results = map(fit_a_pixel, zip(*(data_iterator + [itertools.cycle((pb,))])),)

results = parallel_map.parallel_map(fit_a_pixel,
                                    zip(*(data_iterator + [itertools.cycle((pb,))])),
                                    numcores=8)

cubenames = ['bestdens','bestcol','besttemp','bestopr','bestchi2',
             'mindens','mincol','mintemp','minopr','minchi2',
             'maxdens','maxcol','maxtemp','maxopr','maxchi2',
             'meandens','meancol','meantemp','meanopr','meanchi2',
             'stddens','stdcol','stdtemp','stdopr','stdchi2']
cubetargets = [ (0,0), (0,1), (0,2), (0,3), (0,4), 
                (1,0), (1,1), (1,2), (1,3), (1,4),
                (2,0), (2,1), (2,2), (2,3), (2,4),
                (3,0), (3,1), (3,2), (3,3), (3,4),
                (4,0), (4,1), (4,2), (4,3), (4,4),]

cubes = {c: np.empty(h2co11.shape)+np.nan for c in cubenames}

header = wcs.WCS(cube1.header).to_header()
flatheader = wcs.WCS(cube1.header).sub([wcs.WCSSUB_CELESTIAL]).to_header()

for cc,(i1,i2) in zip(cubenames,cubetargets):
    cubes[cc][includemask] = [r[i1][i2] for r in results]
    hdu = fits.PrimaryHDU(data=cubes[cc],
                          header=header)
    hdu.writeto(dpath("H2CO_ParameterFits_{0}.fits".format(cc)), clobber=True)
    hdu = fits.PrimaryHDU(data=np.nanmax(cubes[cc],axis=0),
                          header=flatheader)
    hdu.writeto(dpath("H2CO_ParameterFits_{0}_max.fits".format(cc)), clobber=True)
    hdu = fits.PrimaryHDU(data=np.nanmin(cubes[cc],axis=0),
                          header=flatheader)
    hdu.writeto(dpath("H2CO_ParameterFits_{0}_min.fits".format(cc)), clobber=True)

wtdmeandens = np.log10((10**np.nansum(cubes['bestdens'] *
                                      10**cubes['bestcol'],axis=0)) /
                       10**np.nansum(cubes['bestcol'],axis=0))
hdu = fits.PrimaryHDU(data=wtdmeandens,
                      header=flatheader)
hdu.writeto(dpath("H2CO_ParameterFits_weighted_mean_density.fits"), clobber=True)

masked_wtdmeandens = np.log10((np.nansum(10**cubes['bestdens'] *
                                         10**cubes['bestcol'] *
                                         (cubes['bestchi2'] < 1), axis=0) /
                               np.nansum(10**cubes['bestcol'] *
                                         (cubes['bestchi2'] < 1), axis=0)))
hdu = fits.PrimaryHDU(data=masked_wtdmeandens,
                      header=flatheader)
hdu.writeto(dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits"), clobber=True)

#chi2 = np.empty(h2co11.shape) + np.nan
#chi2[includemask] = [r[5] for r in results]
#
#hdu = fits.PrimaryHDU(data=chi2,
#                      header=header)
#hdu.writeto(dpath("H2CO_ParameterFits_{0}.fits".format('chi2')), clobber=True)
