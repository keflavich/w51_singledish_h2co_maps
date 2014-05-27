"""
Script to do the parameter fitting for the H2CO cubes

It has gone through at least 11 iterations, some of which tried fitting the
observed optical depth cubes (bad) and others which do the full line modeling.
"""
from load_pyspeckit_cubes import both,T,F,cont11,cont22,h2co11filename
from astropy.io import fits
import os
import numpy as np
import types

#sp = both.get_spectrum(14,31)
#sp.specfit(fittype='formaldehyde_radex',guesses=[4,13,-20,1],multifit=True,quiet=False,verbose=True,negamp=True)
x,y = 88,77
both.plot_spectrum(x,y)
both.specfit(fittype='formaldehyde_radex',guesses=[4,13,-3,20,65,1,
                                                   cont11[y,x],
                                                   cont22[y,x]]*2,
             fixed=[F,F,T,T,F,F,T,T]*2,
             multifit=True,quiet=False,verbose=True,
             use_window_limits=False,
             fit_plotted_area=False)
both.plot_spectrum(x,y, errstyle='fill', residfignum=5)
#both.specfit.add_sliders()

doextra=False

# Added 5/23/2014
mask = fits.getdata('mask_h2co_signal.fits')
both.maskmap = mask

twopar = True
if True:

    #fith2co.plotter_override(sp,vrange=[-50,10],fignum=5,reset_xlimits=True)
    #sp.specfit(fittype='formaldehyde_radex_tau',guesses=[4,13,-20,2],quiet=False,verbose=True,multifit=True)
    #fith2co.plotter_override(sp,vrange=[-50,10],fignum=5,reset_xlimits=True)

    # try2 uses obsfreq instead of restfreq in 2-2 reduction.  Weird artifacts in sess22 though.  Darn.
    # try3 uses sess22 without the failed second receiver
    # try4 uses different ranges for 1-1 polyfit to try to recover those wide lines in the west
    # try5 restricts velocity components to one side or the other
    # try6 includes errors
    # try7 has a lower signal cut
    # try8 includes tbackground
    # try9 includes tbackground with corrected velocities
    # try10 includes tbackground and tries two independent pars.  Scary.
    # try11 includes tbackground and tries two independent pars.  It is pre-masked. 5/23/2014
    parcubefilename = outcube = 'W51_taucube_fit_parcube_try11.fits'
    # parcube is messed?
    refit=not os.path.exists(outcube)
    #refit=True
    if refit:
        if twopar:
            guesses = np.empty((16,) + cont11.shape)
            guesses[:6,:,:] = np.array([4.5,13.1,-3,20,60,1.0,]).reshape((6,1,1,))
            guesses[6,:,:] = cont11
            guesses[7,:,:] = cont22
            guesses[8:14,:,:] = np.array([4.5,13.1,-3,20,70,1.0,]).reshape((6,1,1,))
            guesses[14,:,:] = cont11
            guesses[15,:,:] = cont22

            # correct the continuum in the right-side filament
            guesses[np.array([6,7,14,15]),89:114,196:240] = 2.73

            both.fiteach(guesses=guesses,
                         absorption=True,
                         integral=False,
                         fittype='formaldehyde_radex',
                         multicore=4,
                         signal_cut=4,
                         fixed=[F,F,T,T,F,F,T,T]*2,
                         parlimited=[(T,T), (T,T), (T,T), (T,T), (F,F), (T,F), (T,F), (T,F)]*2,
                         parlimits=[(1,8), (11,16), (-3,np.log10(3.0)), (5,55), (0,0), (0,0), (2.73,0), (2.73,0)]*2,
                         start_from_point=[x,y])
        else:
            guesses = np.empty((8,) + cont11.shape)
            guesses[:6,:,:] = np.array([4.5,13.1,-3,20,60,1.0,]).reshape((6,1,1,))
            guesses[6,:,:] = cont11
            guesses[7,:,:] = cont22
            both.fiteach(guesses=guesses,
                         absorption=True,
                         integral=False,
                         fittype='formaldehyde_radex',
                         multicore=8,
                         signal_cut=3,
                         fixed=[F,F,T,T,F,F,T,T],
                         parlimited=[(T,T), (T,T), (T,T), (T,T), (F,F), (T,F), (T,F), (T,F)],
                         parlimits=[(1,8), (11,16), (-3,np.log10(3.0)), (5,55), (0,0), (0,0), (2.73,0), (2.73,0)],
                         start_from_point=[x,y])
        both.write_fit(parcubefilename,clobber=True)
        cubeheader = fits.getheader(h2co11filename)
        tempf = fits.open(parcubefilename) # error in pyspeckit header tracking
        tempf[0].header = cubeheader
        tempf.writeto(parcubefilename,clobber=True)
        # filter out bad fits for each component independently
        #tempf[0].data[0:4, (tempf[0].data[3,:,:] < 1 )]  = 0  # narrow lines unacceptable
        #tempf[0].data[4:8, (tempf[0].data[7,:,:] < 1 )]  = 0  # narrow lines unacceptable
        tempf.writeto(parcubefilename.replace(".fits","_prefiltered.fits"),clobber=True)
        #tempf[0].data[0, (tempf[0].data[8,:,:]  == 0)|(tempf[0].data[8,:,:] > 1 )]  = 0  # density
        #tempf[0].data[4, (tempf[0].data[12,:,:] == 0)|(tempf[0].data[12,:,:] > 1)] = 0
        #tempf[0].data[1, (tempf[0].data[9,:,:]  == 0)|(tempf[0].data[9,:,:] > 1 )]  = 0  # column
        #tempf[0].data[5, (tempf[0].data[13,:,:] == 0)|(tempf[0].data[13,:,:] > 1)] = 0
        #tempf[0].data[2, (tempf[0].data[10,:,:] == 0)|(tempf[0].data[10,:,:] > 5)] = 0  # velocity
        #tempf[0].data[6, (tempf[0].data[14,:,:] == 0)|(tempf[0].data[14,:,:] > 5)] = 0
        #tempf[0].data[3, (tempf[0].data[11,:,:] == 0)|(tempf[0].data[11,:,:] > 2)] = 0  # width
        #tempf[0].data[7, (tempf[0].data[15,:,:] == 0)|(tempf[0].data[15,:,:] > 2)] = 0
        #tempf.writeto(parcubefilename.replace(".fits","_filtered.fits"),clobber=True)
    else:
        both.load_model_fit(parcubefilename, 8,
                            fittype='formaldehyde_radex',
                            _temp_fit_loc=(x,y))
        #both.parcube[4:,both.parcube[-1,:,:]==0] = 0
        # alternative
        # both.load_model_fit('W51_scaled_parcube.fits',4,fittype='formaldehyde_radex_tau',_temp_fit_loc=(15,30))
        #
    both.mapplot(estimator=0, vmin=3, vmax=6)

    # Individual spectra stuff - not set up for W51 yet
    # sp1 = pyspeckit.Spectrum('G173.47+2.44_h2co.fits',wcstype='D') / 0.51
    # sp1.baseline(exclude=[-25,-15],exclude_units='km/s')
    # sp2 = pyspeckit.Spectrum('G173.47+2.44_h2co_Tastar.fits') / 0.886
    # sp2.crop(-400,400,units='km/s')
    # spectra = pyspeckit.Spectra([sp1,sp2])
    # spectra.Registry.add_fitter('formaldehyde_radex_tau',
    #                 formaldehyde_radex_fitter,4,multisingle='multi')
    # spectra.specfit(fittype='formaldehyde_radex',guesses=[4,13,-20,2],quiet=False,verbose=True,multifit=True,negamp=True)
    # fith2co.plotter_override(spectra,vrange=[-50,10],fignum=6,reset_xlimits=True)
    # sp1.specfit.Registry.add_fitter('formaldehyde_radex',
    #                 formaldehyde_radex_fitter,4,multisingle='multi')
    # sp2.specfit.Registry.add_fitter('formaldehyde_radex',
    #                 formaldehyde_radex_fitter,4,multisingle='multi')
    # sp1.specfit(fittype='formaldehyde_radex',guesses=[4,13,-20,2],quiet=False,verbose=True,multifit=True,negamp=True)
    # sp2.specfit(fittype='formaldehyde_radex',guesses=[4,13,-20,2],quiet=False,verbose=True,multifit=True,negamp=True)
    # sp1.specfit.seterrspec()
    # sp2.specfit.seterrspec()
    # sp1.error=sp1.specfit.errspec
    # sp2.error=sp2.specfit.errspec
    # spectra = pyspeckit.Spectra([sp1,sp2])
    # spectra.Registry.add_fitter('formaldehyde_radex',
    #                 formaldehyde_radex_fitter,4,multisingle='multi')
    # spectra.specfit(fittype='formaldehyde_radex',guesses=[3.3,12.5,-16.2,2,5.5,13.5,-17,0.8],quiet=False,verbose=True,multifit=True,negamp=True)
    # fith2co.plotter_override(spectra,vrange=[-50,10],fignum=6,reset_xlimits=True, residfignum=8,show_components=True)

    if doextra:
        sp = both.get_spectrum(56,50)
        sp.specfit(fittype='formaldehyde_radex_tau',guesses=[3.3,12.5,-16.2,2,5.5,13.5,-17,0.8],quiet=False,verbose=True,multifit=True)
        fith2co.plotter_override(sp,vrange=[-50,10],fignum=10,reset_xlimits=True, residfignum=11, show_components=True)

        #both.fiteach(guesses=[3.6,13.3,-16.66,1.82,5.7,13.3,-17.5,1.3],absorption=True,integral=False,fittype='formaldehyde_radex',multicore=8)

        for ii in xrange(44,50):
            both.plot_spectrum(54,ii)
            pylab.savefig("W51_bestfit_spec54_%i.png" % ii)

        both.plot_spectrum(53,49)
        pylab.savefig("W51_bestfit_spec53_49_IRS2.png")
        both.plot_spectrum(53,46)
        pylab.savefig("W51_bestfit_spec53_49_W51e2.png")

        both.mapplot(estimator=0,vmin=2)

        # Analyze the low-velocity cloud
        err = cube1.cube[cube1.xarr < 0,:,:].std(axis=0)
        wtdavg1 =  (cube1.cube / err**2).sum(axis=1).sum(axis=1) / (1/err**2).sum()
        err2 = cube2.cube[cube2.xarr < 0,:,:].std(axis=0)
        wtdavg2 =  (cube2.cube / err2**2).sum(axis=1).sum(axis=1) / (1/err2**2).sum()
        sp1w = cube1.get_spectrum(55,55)
        sp1w.data = wtdavg1
        sp2w = cube2.get_spectrum(55,55)
        sp2w.data = wtdavg2

        sp1w.crop(-10,20)
        sp2w.crop(-10,20)
        both1dW = pyspeckit.Spectra([sp1w,sp2w])
        both1dW.plot_special = types.MethodType(fith2co.plotter_override, both1dW, both1dW.__class__)
        both1dW.plot_special_kwargs = {'vrange':[-10,20],'fignum':14,'reset_xlimits':True}
        both1dW.Registry.add_fitter('formaldehyde_radex_tau',formaldehyde_radex_tau_fitter,4,multisingle='multi')
        both1dW.specfit(fittype='formaldehyde_radex_tau',guesses=[4,12,5,1],verbose=True,quiet=False)
        both1dW.plot_special(vrange=[-10,20])

        # two IRDCs on the west side that have no fits; they are apparently surprisingly low column and density
        irdc1 = both.get_apspec((48.921,-0.28267,0.0247*3600),coordsys='galactic')
        irdc2 = both.get_apspec((48.99173,-0.304066,0.0247*3600),coordsys='galactic')
        for ii,irdc in enumerate((irdc1,irdc2)):
            irdc.error[irdc.xarr.as_unit('GHz') > 10] = irdc[irdc.xarr.as_unit('GHz') > 10].stats()['std']
            irdc.error[irdc.xarr.as_unit('GHz') < 10] = irdc[irdc.xarr.as_unit('GHz') < 10].stats()['std']
            irdc.plot_special = types.MethodType(fith2co.plotter_override, irdc, irdc.__class__)
            irdc.plot_special_kwargs = {'vrange':[30,90],'fignum':145+ii,'reset_xlimits':True}
            irdc.Registry.add_fitter('formaldehyde_radex_tau',formaldehyde_radex_tau_fitter,4,multisingle='multi')
            irdc.specfit(fittype='formaldehyde_radex_tau',guesses=[4.5,13.1,47,2.0,4,13,67,2],verbose=True,quiet=False)
            irdc.plot_special(vrange=[30,90],fignum=145+ii,reset_xlimits=True)
        print "Something is clearly wrong with how the model is being calculated... the residuals are right, but the models just aren't plotting right."



        # tooling around with pymc
        from pymodelfit import FunctionModel1DAuto

        class DoubleFormaldehydeModel(FunctionModel1DAuto):
            def f(self,x,density0=4,column0=13,xoff_v0=0.0,width0=1.0,density1=4,column1=13,xoff_v1=0.0,width1=1.0):
                return pyspeckit.spectrum.models.formaldehyde.formaldehyde_radex(x, density=density0, column=column0, xoff_v=xoff_v0,width=width0,texgrid=((4,5,texgrid1),(14,15,texgrid2)),taugrid=((4,5,taugrid1),(14,15,taugrid2)),hdr=hdr) + pyspeckit.spectrum.models.formaldehyde.formaldehyde_radex(x, density=density1, column=column1, xoff_v=xoff_v1,width=width1,texgrid=((4,5,texgrid1),(14,15,texgrid2)),taugrid=((4,5,taugrid1),(14,15,taugrid2)),hdr=hdr)

        class FormaldehydeModel(FunctionModel1DAuto):
            def f(self,x,density=4,column=13,xoff_v=0.0,width=1.0):
                return pyspeckit.spectrum.models.formaldehyde.formaldehyde_radex(x, density=density, column=column, xoff_v=xoff_v,width=width,texgrid=((4,5,texgrid1),(14,15,texgrid2)),taugrid=((4,5,taugrid1),(14,15,taugrid2)),hdr=hdr)

        class FormaldehydeTauModel(FunctionModel1DAuto):
            def f(self,x,density=4,column=13,xoff_v=0.0,width=1.0):
                return formaldehyde_radex_tau(x, density=density, column=column,
                        xoff_v=xoff_v, width=width,
                        taugrid=((4, 5, taugrid1), (14, 15, taugrid2)),
                        hdr=hdr)

        class DoubleFormaldehydeModelFF(FunctionModel1DAuto):
            def f(self,x,density0=4,column0=13,xoff_v0=0.0,width0=1.0,density1=4,column1=13,xoff_v1=0.0,width1=1.0,ff2=1.0):
                return (pyspeckit.spectrum.models.formaldehyde.formaldehyde_radex(x, density=density0, column=column0, xoff_v=xoff_v0,width=width0,texgrid=((4,5,texgrid1),(14,15,texgrid2)),taugrid=((4,5,taugrid1),(14,15,taugrid2)),hdr=hdr) + 
                        ff2*pyspeckit.spectrum.models.formaldehyde.formaldehyde_radex(x, density=density1, column=column1, xoff_v=xoff_v1,width=width1,texgrid=((4,5,texgrid1),(14,15,texgrid2)),taugrid=((4,5,taugrid1),(14,15,taugrid2)),hdr=hdr))

        compFMtau = FormaldehydeTauModel()
        compFMCtau = compFMtau.getMCMC(both1dW.xarr, 
                both1dW.data,
                priors={'density':(1,6),'width':(0.1,4),'xoff_v':(0,10),'column':(10,16)})
        compFMCtau.sample(10000,burn=100)
        compFMCs = compFMC.stats()

        pylab.figure(1)
        pylab.clf()
        agpy.pymc_plotting.hist2d(compFMCtau,'density','column',bins=50,doerrellipse=False,varslice=(1000,None,None))
        pylab.ylabel('Column per km/s per pc')
        pylab.xlabel('Density (H$_2$ cm$^{-3}$)')

        spec67_64 = both.get_spectrum(67,64)
        spec67_64_mcmc = compFMtau.getMCMC(spec67_64.xarr, 
                spec67_64.data,
                priors={'density':(1,6),'width':(0.1,4),'xoff_v':(50,70),'column':(11,16)})
        spec67_64_mcmc.sample(20000,burn=100)
        spec67_64_mcmc_stats = spec67_64_mcmc.stats()
        spec67_64.specfit.modelpars = [spec67_64_mcmc_stats[x]['quantiles'][50] for x in ['density', 'column',  'xoff_v', 'width',]]
        spec67_64.specfit.modelerrs = [spec67_64_mcmc_stats[x]['standard deviation'] for x in ['density', 'column',  'xoff_v', 'width',]]
        spec67_64.specfit.parinfo = spec67_64.specfit.fitter._make_parinfo(values=spec67_64.specfit.modelpars,errors=spec67_64.specfit.modelerrs)[0]
        spec67_64.specfit.npeaks = 1
        spec67_64.units = 'Optical Depth $\\tau$'
        spec67_64.plot_special = types.MethodType(fith2co.plotter_override, spec67_64, spec67_64.__class__)
        spec67_64.plot_special(vrange=(40,85),fignum=148)
        spec67_64.plotter.savefig('spec67_64_bestfit_MCMC.png')

        pylab.figure(2)
        pylab.clf()
        agpy.pymc_plotting.hist2d(spec67_64_mcmc,'density','column',bins=20,doerrellipse=False,varslice=(1000,None,None),fignum=2)
        pylab.ylabel('Column per km/s per pc')
        pylab.xlabel('Density (H$_2$ cm$^{-3}$)')
        pylab.savefig('MCMC_DensColplot_67_64.png')

        spec65_60 = both.get_spectrum(65,60)
        spec65_60_mcmc = compFMtau.getMCMC(spec65_60.xarr, 
                spec65_60.data,
                priors={'density':(1,6),'width':(0.1,4),'xoff_v':(50,70),'column':(11,16)})
        spec65_60_mcmc.sample(20000,burn=100)
        spec65_60_mcmc_stats = spec65_60_mcmc.stats()
        spec65_60.specfit.modelpars = [spec65_60_mcmc_stats[x]['quantiles'][50] for x in ['density', 'column',  'xoff_v', 'width',]]
        spec65_60.specfit.modelerrs = [spec65_60_mcmc_stats[x]['standard deviation'] for x in ['density', 'column',  'xoff_v', 'width',]]
        spec65_60.specfit.parinfo = spec65_60.specfit.fitter._make_parinfo(values=spec65_60.specfit.modelpars,errors=spec65_60.specfit.modelerrs)[0]
        spec65_60.specfit.npeaks = 1
        spec65_60.units = 'Optical Depth $\\tau$'
        spec65_60.plot_special = types.MethodType(fith2co.plotter_override, spec65_60, spec65_60.__class__)
        spec65_60.plot_special(vrange=(40,85),fignum=147)
        spec65_60.plotter.savefig('spec65_60_bestfit_MCMC.png')

        pylab.figure(3)
        pylab.clf()
        agpy.pymc_plotting.hist2d(spec65_60_mcmc,'density','column',bins=20,doerrellipse=False,varslice=(1000,None,None),fignum=3)
        pylab.ylabel('Column per km/s per pc')
        pylab.xlabel('Density (H$_2$ cm$^{-3}$)')
        pylab.savefig('MCMC_DensColplot_65_60.png')


        """
        compFM = DoubleFormaldehydeModel()
        nsamples = 500

        compFMC = compFM.getMCMC(spectra2.xarr,spectra2.data,priors={'density0':(1,6),'width0':(0.1,20),'xoff_v0':(-20,-10),'column0':(10,16),'density1':(4,8),'xoff_v1':(-20,-10),'width1':(0.1,20),'column1':(10,16)})
        compFMC.sample(nsamples,burn=100)
        compFMCs = compFMC.stats()
        for p in compFM.params: setattr(compFM,p,compFMCs[p]['mean'])

        compFMC2 = compFM.getMCMC(spectra2.xarr,spectra2.data,priors={'density0':(2.5,3.9),'width0':(0.1,2.5),'xoff_v0':(-19,-15),'column0':(12.6,13.6),'density1':(3,8),'xoff_v1':(-19,-15),'width1':(0.1,2.5),'column1':(10,16)})
        compFMC2.sample(nsamples,burn=100)
        compFMC2s = compFMC2.stats()
        for p in compFM.params: setattr(compFM,p,compFMC2s[p]['mean'])

        # this one looks OK
        compFMC4 = compFM.getMCMC(spectra2.xarr,spectra2.data,priors={'density0':(1.5,2.5),
                    'width0':(1.0,1.5), 'xoff_v0':(-17,-16), 'column0':(13.5,14.1),
                    'density1':(4,6), 'xoff_v1':(-18,-16), 'width1':(0.4,0.8), 'column1':(13,15)})
        compFMC4.sample(nsamples,burn=100)
        compFMC4s = compFMC4.stats()
        for p in compFM.params: setattr(compFM,p,compFMC4s[p]['mean'])
        # In [437]: compFM.pardict
        # Out[437]: 
        # {'column0': 13.887521642187115,
        #  'column1': 13.963400410810015,
        #  'density0': 2.1024458697974047,
        #  'density1': 4.5207616120955221,
        #  'width0': 1.3137357237158149,
        #  'width1': 0.60555978959148793,
        #  'xoff_v0': -16.450329019841522,
        #  'xoff_v1': -17.131943593967165}


        # this one looks absolutely terrible (it's not even a fit)
        # maybe the weights are wrong because of differing bin sizes and errors?
        compFMFF = DoubleFormaldehydeModelFF()
        compFMC3 = compFMFF.getMCMC(spectra2.xarr,spectra2.data,priors={'density0':(2.0,3.9),
            'width0':(0.1,2.5), 'xoff_v0':(-19,-15), 'column0':(13,15),
            'density1':(4.0,8), 'xoff_v1':(-19,-15), 'width1':(0.1,2.5), 'column1':(13.25,15),
            'ff2':(0.1,1.0)})
        compFMC3.sample(nsamples,burn=100)
        compFMC3s = compFMC3.stats()
        for p in compFMFF.params: setattr(compFMFF,p,compFMC3s[p]['mean'])

        spectra.specfit.modelpars = list(compFMFF.parvals[:-1])
        spectra.specfit.fitter.mpp = spectra.specfit.modelpars
        spectra.specfit.fitter.mpperr = [compFMC3s[p]['standard deviation'] for p in compFMFF.params]
        for par in spectra.specfit.parinfo:
            parname = par['parname'].lower() if not 'CENTER' in par['parname'] else par['parname'].replace('CENTER','xoff_v')
            par['value'] = compFMC3s[parname]['mean']
            par['error'] = compFMC3s[parname]['standard deviation']
        spectra.specfit.fitter.parinfo = spectra.specfit.parinfo

        """
