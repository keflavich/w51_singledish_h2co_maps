"""
Fit the continuum and absorption together appropriately (line depth depends on
tau and tex)

Step 1: Just fit tau and tex for each line independently (no RADEX needed)
Step 2: Combine 1-1 and 2-2

***THIS WILL NOT WORK***
Tau and Tex are very degenerate in all interesting cases.  Only the joint
restriction provided by the LVG models provides useful constraints.

It is necessary to fit *both* lines using *both* tex and tau

Therefore, the most efficient and sensible route is to fit pure absorption
profiles with no degenerate parameters, then try to fit RADEX models to those.
"""
import pyspeckit
import numpy as np
from astropy.io import fits
# shortcuts to minimize tedium
T,F = True,False

dpath = '/Users/adam/work/h2co/maps/W51/'
pfx11 = 'W51_H2CO11_cube_supersampled'
pfx22 = 'W51_H2CO22_pyproc_cube_lores_supersampled'

cube1 = pyspeckit.Cube(dpath+pfx11+'.fits')
C1 = fits.getdata(dpath+pfx11+"_continuum.fits")
cube1.Registry.add_fitter('formaldehyde_vtau', pyspeckit.models.formaldehyde.formaldehyde_vtau.background_fitter, 5, multisingle='multi', override=True)

# these should be set from the header...
cube1.xarr.refX_units='GHz'
cube1.xarr.refX = 4.829659400

# add CMB
cube1.cube += 2.73
C1 += 2.73

cube1.plot_spectrum(104,93)
cube1.specfit(fittype='formaldehyde_vtau',
              guesses=[2.73,1,0.1,60,1], # Tbg,Tex,tau,v,w
              absorption=True,
              limited=[(T,F), (T,T), (T,T), (T,T), (T,F)],
              limits=[(0,0),(0.6,1.6), (0,1.5), (40,75), (0,0)])

# compute error from low velocity data
E1 = cube1.cube[cube1.xarr.as_unit('km/s') < 0].std(axis=0)
#cube1.errorcube = np.repeat(np.reshape(E,(1,)+E.shape),cube1.shape[0],axis=0)

parcubefilename = dpath+'H2CO11_h2cofit_parameters.fits'
# parnames=['Tex','tau','center','width'],
cube1.fiteach(guesses=[2.73,1,0.1,70,1],#,0,1,1,70,1],
              errmap=E1,
              continuum_map=C1,
              absorption=True,
              integral=False,
              fittype='formaldehyde_vtau',
              multicore=8,
              signal_cut=5,
              parlimited=[(True,False),(True,True), (True,True), (True,True), (True,False)],
              parlimits=[(0,0),(0,100), (0,1.5), (65,75), (0,0)],#+[(0,0), (0,0), (0,1.5), (65,75), (0,0)],
              start_from_point=[150,100],
              verbose_level=1)
cube1.write_fit(parcubefilename,clobber=True)


cube2 = pyspeckit.Cube(dpath+pfx22+'.fits')
cube2.Registry.add_fitter('formaldehyde_vtau', pyspeckit.models.formaldehyde.formaldehyde_vtau_background_fitter, 5, multisingle='multi')

# these should be set from the header...
cube2.xarr.refX_units='GHz'
cube2.xarr.refX = 14.48847881

# compute error from low velocity data
E2 = cube2.cube[cube2.xarr.as_unit('km/s') < 0].std(axis=0)
#E2 = fits.getdata(datapath+pfx22+
#cube2.errorcube = np.repeat(np.reshape(E,(1,)+E.shape),cube2.shape[0],axis=0)

parcubefilename = dpath+'H2CO22_h2cofit_parameters.fits'
# parnames=['Tbg', 'Tex','tau','center','width'],
cube2.fiteach(guesses=[2.73,1,0.1,60,1],
              errmap=E2,
              absorption=True,
              integral=False,
              fittype='formaldehyde_vtau',
              multicore=1,
              signal_cut=5,
              parlimited=[(True,True), (True,True), (True,True), (True,False)]*2,
              parlimits=[(0,100), (1.05,1.95) (0,1.5), (40,65), (0,0)],
              start_from_point=[56,50])
cube2.write_fit(parcubefilename,clobber=True)
