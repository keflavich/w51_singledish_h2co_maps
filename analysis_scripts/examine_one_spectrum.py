execfile('/Users/adam/work/w51_singledish_maps/analysis_scripts/pyspeckit_individual_spectra.py')
import paths

#if not 'spectra' in locals():
#    spectra = w51main()
#    s17=spectra[-1]
from load_pyspeckit_cubes import formaldehyde_radex_fitter

s17 = pyspeckit.Spectrum(paths.dpath('spectralfits/spectralfits_w51main_aperture_southdensespot.fits'))
s17.Registry.add_fitter('formaldehyde_radex',formaldehyde_radex_fitter,8,multisingle='multi')

dofit(s17, s17.header['CONT11'], s17.header['CONT22'], 2)
plotitem(s17, 17, dolegend=True)

verbose = False
fixed=np.array([F,F,T,T,F,F,T,T]*3)
limits, limited = modelpars()
vguesses=[67.3,58.3]
c11=s17.header['CONT11']
c22=s17.header['CONT22']
gg = [4,13,-3,20,vguesses[0],1, c11, c22]*2
gg[12] = vguesses[1]

gg[0] = 4.415
gg[1] = 12.7407
gg[5] = 2.1752
gg[4] = 67.3357
gg[8] = 5.6
gg[9] = 13.4
gg[10] = np.log10(3)
gg[13] = 3.5
fixed[:8] = True
fixed[10] = False
fixed[11] = False

# Add a 3rd component
gg += [4,13,-3,20,58.3,2.5, c11, c22]
limited = limited + limited[:8]
limited[12] = (True,True)
limited[20] = (True,True)
limited[13] = (True,True)
limited[21] = (True,True)
limited[8] = (True,True)
limited[16] = (True,True)
limited[9] = (True,True)
limited[17] = (True,True)
limits = limits + limits[:8]
limits[12] = (56,60)
limits[20] = (56,60)
limits[13] = (2,4)
limits[21] = (2,4)
limits[8] = (2,7)
limits[16] = (2,5.5)
limits[9] = (12.25,14)
limits[17] = (12.25,14)

s17.specfit(fittype='formaldehyde_radex',guesses=gg,
         fixed=fixed, multifit=True,quiet=True,verbose=verbose,
         limits=limits,
         limited=limited,
         use_window_limits=False, fit_plotted_area=False)
plotitem(s17, 17, dolegend=True)
s17.specfit.add_sliders()
