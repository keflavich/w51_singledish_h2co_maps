import pyspeckit
from pyspeckit.spectrum import models
from pyspeckit.wrappers import fith2co
from astropy.io import fits
import numpy as np
import pyregion
import paths
from paths import (datapath, dpath, rpath, mpath, h2co11subfn, h2co22subfn,
                   cont2cm, cont6cm)
from common_constants import TCMB, etamb_gbt
from h2co_modeling.grid_fitter import grid_2p_getmatch
from astropy import log
import warnings
import os

# These are old versions preserved for posterity; the files have been
# restructured and the Arecibo data are pre-corrected for etamb (they have to
# be correct on-the-fly because of a slight elevation dependence in the main
# beam correction)
#cube1 = pyspeckit.Cube('W51_H2CO11_cube_sub.fits')        / 0.51 # eta_mb = 0.51 for arecibo @ c-band according to outergal paper
#cube2 = pyspeckit.Cube('W51_H2CO22_pyproc_cube_sess22_sub.fits') / 0.886 # from both outergal and pilot

# this is no longer a relevant todo item
# the abundance is degenerate with the line-of-sight length scale / velocity
# dispersion.
# those can vary.
#print "TO DO: fix (hold in place) abundance"

plot=False

#etamb already accounted for
# July 6: Really?  Is it?  I'm pretty damned sure it's not.  Jackass.
# etamb is dealt with in the reduction script for H2CO 1-1 in
# reduce_map.pro:accum_map.
# For 2-2, it is never touched: it is only ever mentioned in make_taucube
# The problem is that I made the transition from tau->sub without checking this...
# even though it's obviously indicated in the post-import lines above
cube1 = pyspeckit.Cube(h2co11subfn)
cube2 = pyspeckit.Cube(h2co22subfn)
cube2.data /= etamb_gbt # etamb scaling *crucial* for 2-2!
cube1.xarr.refX_units='GHz'
cube1.xarr.refX = 4.829659400
cube1.xarr.velocity_convention = 'radio'
E1 = cube1.cube[cube1.xarr.as_unit('km/s') < 0].std(axis=0)
cube1.errorcube = np.repeat(np.reshape(E1,(1,)+E1.shape),cube1.shape[0],axis=0)
cube2.xarr.refX_units='GHz'
cube2.xarr.refX = 14.48847881
cube2.xarr.velocity_convention = 'radio'
E2 = cube2.cube[cube2.xarr.as_unit('km/s') < 0].std(axis=0)
cube2.errorcube = np.repeat(np.reshape(E2,(1,)+E2.shape),cube2.shape[0],axis=0)
both = pyspeckit.CubeStack([cube1,cube2])
if plot:
    both.mapplot()
else:
    both.mapplot.makeplane()
#both.units = 'Optical Depth $\\tau$'
#both.header['BUNIT'] = 'Optical Depth $\\tau$'
# need continua now
cont11filename = cont6cm
cont22filename = cont2cm
cont11 = fits.getdata(cont11filename) + TCMB
cont22 = fits.getdata(cont22filename) / etamb_gbt + TCMB
cont11[cont11<TCMB] = TCMB
cont22[cont22<TCMB] = TCMB

contfrontregions = pyregion.open(rpath('continuum_in_the_front.reg'))
header = fits.getheader(cont11filename)
contfrontmask = contfrontregions.get_mask(fits.PrimaryHDU(data=cont11,header=header))
cont11[contfrontmask] = TCMB
cont22[contfrontmask] = TCMB


#path_to_data = "/Users/adam/work/h2co/radex/troscompt_grid_March2012"

texgrid1 = fits.getdata(paths.model11tex)
taugrid1 = fits.getdata(paths.model11tau)
texgrid2 = fits.getdata(paths.model22tex)
taugrid2 = fits.getdata(paths.model22tau)
hdr    = fits.getheader(paths.model11tex)
# # this deserves a lot of explanation:
# # models.formaldehyde.formaldehyde_radex is the MODEL that we are going to fit
# # models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks,
# # and annotations
# # all of the parameters after the first are passed to the model function
T,F = True,False
formaldehyde_radex_fitter = models.model.SpectralModel(
        models.formaldehyde.formaldehyde_radex_orthopara_temp, 8,
        parnames=['density','column','orthopara','temperature','center','width','tbackground1','tbackground2'], 
        parvalues=[4,12,0.0,15.0,0,1,2.73,2.73],
        parlimited=[(True,True), (True,True), (True,True), (True,True), (False,False), (True,False), (True,False), (True,False)], 
        parlimits=[(1,8), (11,16), (-3,np.log10(3.0)), (5,55), (0,0), (0,0), (2.73,0), (2.73,0)],
        fixed=[F,F,T,T,F,F,T,T],
        parsteps=[0.01,0.01,0,0,0,0,0,0],
        fitunits='Hz',
        texgrid=((4,5,texgrid1),(14,15,texgrid2)),
        taugrid=((4,5,taugrid1),(14,15,taugrid2)),
        hdr=hdr,
        shortvarnames=("n","N",'OP','T',"v","\\sigma",'T_{bg,1}','T_{bg,2}'),
        )


#both.Registry.add_fitter('formaldehyde_radex',formaldehyde_radex_fitter,4,multisingle='multi')
both.Registry.add_fitter('formaldehyde_radex',formaldehyde_radex_fitter,8,multisingle='multi')

both.plot_special = fith2co.plotter_override
both.plot_special_kwargs = {'vrange':[30,90],'fignum':3,'reset_xlimits':True}


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
    #spec = (1.0-np.exp(-np.array(tau_nu_cumul)))*(Tex-Tbackground)

    match, indbest, chi2 = grid_2p_getmatch(tline1+cont1, etline1, pargrid1,
                                            tline2+cont2, etline2, pargrid2)

    likelihood = np.exp(-chi2/2.)
    def lwt(x, mask=match):
        """ Likelihood-weighted mean """
        return (x*likelihood)[mask].sum() / likelihood[mask].sum()
    def lstd(x, mask=match):
        """ Likelihood-weighted stddev """
        return (((x-lwt(x,mask))**2*likelihood)[mask].sum() / likelihood[mask].sum() * (mask.sum()/(mask.sum()-1)))**0.5

    chi2best = chi2.flat[indbest]
    dens_best = densityarr.flat[indbest]
    col_best = columnarr.flat[indbest]
    temp_best = temparr.flat[indbest]
    opr_best = oprarr.flat[indbest]

    best = (dens_best, col_best, temp_best, opr_best, chi2best)

    if np.count_nonzero(match) == 0:
        log.warn(str(args))
        warnings.warn("Found no matches.  Returning NaNs")
        return best, np.nan, np.nan, np.nan, np.nan, np.nan


    chi2_all = chi2[match]
    dens_all = densityarr[match]
    col_all = columnarr[match]
    temp_all = temparr[match]
    opr_all = oprarr[match]

    mins = map(min, (dens_all, col_all, temp_all, opr_all, chi2_all))
    maxs = map(max, (dens_all, col_all, temp_all, opr_all, chi2_all))
    mean = map(np.mean, (dens_all, col_all, temp_all, opr_all, chi2_all))
    std = map(np.std, (dens_all, col_all, temp_all, opr_all, chi2_all))
    lweighted = map(lwt, (densityarr, columnarr, temparr, oprarr))
    lerror = map(lstd, (densityarr, columnarr, temparr, oprarr))

    return best,mins,maxs,mean,std,chi2best,lweighted,lerror
