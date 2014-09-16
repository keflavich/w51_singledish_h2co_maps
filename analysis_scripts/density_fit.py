from paths import datapath,dpath,rpath
from astropy.convolution import convolve
from astropy.io import fits
from astropy import wcs
import pyregion
import numpy as np
from astropy import log

from astropy.utils.console import ProgressBar
import itertools
from pyspeckit import parallel_map

from h2co_modeling.grid_fitter import grid_2p_getmatch

from load_pyspeckit_cubes import (both, T, F, cont11, cont22, 
                                  cube1, cube2, texgrid1,  taugrid1,  texgrid2,
                                  taugrid2,  hdr, fit_a_pixel)
from common_constants import TCMB

from build_mask import includemask


###
# Create a new tau ratio cube
ratioF = fits.open(datapath+'W51_H2CO11_to_22_tau_ratio_supersampled.fits')
tau_ratio_fgbg11 = -np.log((cont11+cube1.cube)/(cont11))
tau_ratio_fgbg22 = -np.log((cont22+cube2.cube)/(cont22))
newratio = np.array(tau_ratio_fgbg11/tau_ratio_fgbg22)
newratio[~includemask] = np.nan
ratioF[0].data = newratio
ratioF.writeto(datapath+'W51_H2CO11_to_22_tau_ratio_supersampled_fgbg.fits', clobber=True)
###




data_iterator = [np.array(cube1.cube[includemask]),
                 np.array(cube1.errorcube[includemask]),
                 np.repeat(cont11[None,:],cube1.shape[0],axis=0).reshape(cube1.cube.shape)[includemask],
                 np.array(cube2.cube[includemask]),
                 np.array(cube2.errorcube[includemask]),
                 np.repeat(cont22[None,:],cube2.shape[0],axis=0).reshape(cube2.cube.shape)[includemask],]

for d in data_iterator:
    if np.any(np.isnan(d)):
        raise ValueError("There are NaNs in the data.  This is strictly impossible.")

args = zip(*data_iterator)[0]
print(args,fit_a_pixel(args))

pb = ProgressBar(len(data_iterator[0]))

#results = map(fit_a_pixel, zip(*(data_iterator + [itertools.cycle((pb,))])),)

# This line contains all the meat!  It can take a long time
results = parallel_map.parallel_map(fit_a_pixel,
                                    zip(*(data_iterator + [itertools.cycle((pb,))])),
                                    numcores=8)

cubenames = ['bestdens','bestcol','besttemp','bestopr','bestchi2',
             'mindens','mincol','mintemp','minopr','minchi2',
             'maxdens','maxcol','maxtemp','maxopr','maxchi2',
             'meandens','meancol','meantemp','meanopr','meanchi2',
             'stddens','stdcol','stdtemp','stdopr','stdchi2',
             'likewtddens','likewtdcol','likewtdtemp','likewtdopr',
             'likestddens','likestdcol','likestdtemp','likestdopr',
            ]
cubetargets = [ (0,0), (0,1), (0,2), (0,3), (0,4), 
                (1,0), (1,1), (1,2), (1,3), (1,4),
                (2,0), (2,1), (2,2), (2,3), (2,4),
                (3,0), (3,1), (3,2), (3,3), (3,4),
                (4,0), (4,1), (4,2), (4,3), (4,4),
                (6,0), (6,1), (6,2), (6,3), 
                (7,0), (7,1), (7,2), (7,3), 
              ]

cubes = {c: np.empty(cube1.cube.shape)+np.nan for c in cubenames}

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
