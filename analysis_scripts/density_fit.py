from load_pyspeckit_cubes import (both, T, F, cont11, cont22, h2co11filename,
                                  cube1, cube2, texgrid1,  taugrid1,  texgrid2,
                                  taugrid2,  hdr)
from paths import datapath,dpath,rpath
from astropy.convolution import convolve
from astropy.io import fits
from astropy import wcs
import numpy as np

from astropy.utils.console import ProgressBar
import itertools
from pyspeckit import parallel_map

from h2co_modeling.grid_fitter import grid_2p_getmatch

### copied from tau_ratio_cube
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
###

# copied from pyspeckit.spectrum.models.formaldehyde
winds,zinds,yinds,xinds = np.indices(taugrid1.shape)
densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
temparr    = (zinds+hdr['CRPIX3']-1)*hdr['CDELT3']+hdr['CRVAL3'] # temperature
oprarr     = (winds+hdr['CRPIX4']-1)*hdr['CDELT4']+hdr['CRVAL4'] # log ortho/para ratio

def fit_a_pixel(args):
    if len(args) == 6:
        tline1, etline1, cont1, tline2, etline2, cont2 = args
    elif len(args) == 7:
        tline1, etline1, cont1, tline2, etline2, cont2, pb = args
        pb.update()

    pargrid1 = (cont1*np.exp(-taugrid1) + (1-np.exp(-taugrid1))*texgrid1)
    pargrid2 = (cont2*np.exp(-taugrid2) + (1-np.exp(-taugrid2))*texgrid2)

    match, indbest, chi2 = grid_2p_getmatch(tline1+cont1, etline1, pargrid1,
                                                tline2+cont2, etline2, pargrid2)

    chi2best = chi2.flat[indbest]
    chi2_all = chi2[match]
    dens_best = densityarr.flat[indbest]
    dens_all = densityarr[match]
    col_best = columnarr.flat[indbest]
    col_all = columnarr[match]
    temp_best = temparr.flat[indbest]
    temp_all = temparr[match]
    opr_best = oprarr.flat[indbest]
    opr_all = oprarr[match]

    best = (dens_best, col_best, temp_best, opr_best, chi2best)
    mins = map(min, (dens_all, col_all, temp_all, opr_all, chi2_all))
    maxs = map(max, (dens_all, col_all, temp_all, opr_all, chi2_all))
    mean = map(np.mean, (dens_all, col_all, temp_all, opr_all, chi2_all))
    std = map(np.std, (dens_all, col_all, temp_all, opr_all, chi2_all))

    return best,mins,maxs,mean,std,chi2best


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
cubetargets = [ (0,0), (0,1), (0,2), (0,3),
                (1,0), (1,1), (1,2), (1,3),
                (2,0), (2,1), (2,2), (2,3),
                (3,0), (3,1), (3,2), (3,3),
                (4,0), (4,1), (4,2), (4,3)]

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

wtdmeandens = ((np.nansum(cubes['bestdens'] * cubes['bestcol'],axis=0)) 
               / np.nansum(cubes['bestcol'],axis=0))
hdu = fits.PrimaryHDU(data=wtdmeandens,
                      header=flatheader)
hdu.writeto(dpath("H2CO_ParameterFits_weighted_mean_density.fits"), clobber=True)

masked_wtdmeandens = (np.nansum(cubes['bestdens'] * cubes['bestcol'] * (cubes['bestchi2'] < 1), axis=0)
                      / np.nansum(cubes['bestcol'] * (cubes['bestchi2'] < 1), axis=0))
hdu = fits.PrimaryHDU(data=masked_wtdmeandens,
                      header=flatheader)
hdu.writeto(dpath("H2CO_ParameterFits_weighted_mean_density_chi2masked.fits"), clobber=True)

#chi2 = np.empty(h2co11.shape) + np.nan
#chi2[includemask] = [r[5] for r in results]
#
#hdu = fits.PrimaryHDU(data=chi2,
#                      header=header)
#hdu.writeto(dpath("H2CO_ParameterFits_{0}.fits".format('chi2')), clobber=True)
