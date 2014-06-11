import pyspeckit
from pyspeckit.spectrum import models
from pyspeckit.wrappers import fith2co
from astropy.io import fits
import numpy as np
import pyregion
from paths import datapath,dpath,rpath
from common_constants import TCMB
#cube1 = pyspeckit.Cube('W51_H2CO11_cube_sub.fits')        / 0.51 # eta_mb = 0.51 for arecibo @ c-band according to outergal paper
#cube2 = pyspeckit.Cube('W51_H2CO22_pyproc_cube_sess22_sub.fits') / 0.886 # from both outergal and pilot

# this is no longer a relevant todo item
# the abundance is degenerate with the line-of-sight length scale / velocity
# dispersion.
# those can vary.
#print "TO DO: fix (hold in place) abundance"

plot=False

#etamb already accounted for
h2co11filename = '/Users/adam/work/h2co/maps/W51/W51_H2CO11_cube_supersampled_sub.fits'
cube1 = pyspeckit.Cube(h2co11filename)
cube2 = pyspeckit.Cube('/Users/adam/work/h2co/maps/W51/W51_H2CO22_pyproc_cube_lores_supersampled_sub.fits')
cube1.xarr.refX_units='GHz'
cube1.xarr.refX = 4.829659400
E1 = cube1.cube[cube1.xarr.as_unit('km/s') < 0].std(axis=0)
cube1.errorcube = np.repeat(np.reshape(E1,(1,)+E1.shape),cube1.shape[0],axis=0)
cube2.xarr.refX_units='GHz'
cube2.xarr.refX = 14.48847881
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
cont11filename = '/Users/adam/work/h2co/maps/W51/W51_H2CO11_cube_supersampled_continuum.fits'
cont11 = fits.getdata(cont11filename) + TCMB
cont22 = fits.getdata('/Users/adam/work/h2co/maps/W51/W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits') + TCMB
cont11[cont11<TCMB] = TCMB
cont22[cont22<TCMB] = TCMB

contfrontregions = pyregion.open(rpath('continuum_in_the_front.reg'))
header = fits.getheader(datapath+'W51_H2CO11_cube_supersampled_continuum.fits')
contfrontmask = contfrontregions.get_mask(fits.PrimaryHDU(data=cont11,header=header))
cont11[contfrontmask] = TCMB
cont22[contfrontmask] = TCMB


path_to_data = "/Users/adam/work/h2co/radex/troscompt_grid_March2012"

texgrid1 = fits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tex1.fits')
taugrid1 = fits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau1.fits')
texgrid2 = fits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tex2.fits')
taugrid2 = fits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau2.fits')
hdr    = fits.getheader(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau2.fits')
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

