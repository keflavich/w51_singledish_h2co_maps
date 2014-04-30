from __future__ import print_function
from load_pyspeckit_cubes import both,T,F,cont11,cont22,cont11filename
import pyspeckit
from astropy.io import fits
import pyregion
import types
from pyspeckit.wrappers import fith2co
import numpy as np
import pylab as pl

regfiledict = {"/Users/adam/work/w51/dense_filament_spectral_apertures.reg":
               'spectralfits/spectralfits_70kmscloud',
               "/Users/adam/work/w51/filament_leftside_spectral_apertures.reg":
               'spectralfits/spectralfits_70kmscloudLeft',
               "/Users/adam/work/w51/w51main_spectral_apertures.reg":
               'spectralfits/spectralfits_w51main',
               "/Users/adam/work/w51/maus_spectral_apertures.reg":
               'spectralfits/spectralfits_maus',
               "/Users/adam/work/w51/middlechunk_spectral_apertures.reg":
               'spectralfits/spectralfits_63kmscloud'}


def modelpars():
    
    limits = [(1, 8),
              (11, 16),
              (-3, 0.47712125471966244),
              (5, 55),
              (50, 80), # velo
              (0, 0),
              (2.73, 0),
              (2.73, 0),
              (1, 8),
              (11, 16),
              (-3, 0.47712125471966244),
              (5, 55),
              (50, 80),
              (0, 0),
              (2.73, 0),
              (2.73, 0)]

    limited = [(True, True),
               (True, True),
               (True, True),
               (True, True),
               (True, True), # velo
               (True, False),
               (True, False),
               (True, False),
               (True, True),
               (True, True),
               (True, True),
               (True, True),
               (True, True),
               (True, False),
               (True, False),
               (True, False)]

    return limits,limited

def plotitem(sp, ii=0, **kwargs):
    sp.plot_special = types.MethodType(fith2co.plotter_override, sp, sp.__class__)
    plotkwargs = {'vrange':[40,80],'fignum':ii+1,'reset_xlimits':True,
                  'annotate': False, 'errstyle':'fill'}
    plotkwargs.update(kwargs)
    sp.plot_special_kwargs = plotkwargs
    sp.plot_special(**plotkwargs)
    # plot cleanup junk
    sp.plotter.figure.subplots_adjust(hspace=0)
    # move h2co2-2 title down
    sp.axisdict['twotwo'].title.set_y(0.9)
    #sp.plotter.figure.subplots_adjust(right=0.8)
    #sp.specfit.annotate(bbox_to_anchor=(1.33,1.0))
    #sp.specfit.annotate(loc='lower right')
    #sp.specfit.annotate(bbox_to_anchor=(1.33,1.0))
    # chop off the top ylabel
    sp.axisdict['twotwo'].set_yticks(sp.axisdict['twotwo'].get_yticks()[:-1])
    # remove top plot xticks
    sp.axisdict['oneone'].set_xticks([])
    sp.axisdict['oneone'].set_xticklabels([])

def do_indiv_fits(regfilename, outpfx, ncomp=1):
    regions = pyregion.open(regfilename)

    spectra = [both.get_apspec(r.coord_list, coordsys=r.coord_format,
                               wunit='degree')
               for r in regions]

    cont11hdu = fits.open(cont11filename)[0]

    limits, limited = modelpars()

    for cont in ('front','back',):
        print("Continuum in the %s" % cont)
        for ii,(sp,r) in enumerate(zip(spectra,regions)):
            fig = pl.figure(ii, figsize=(12,8))
            fig.clf()
            name = r.attr[1]['text']

            if cont == 'back':
                M = pyregion.ShapeList([r]).get_mask(cont11hdu)
                c11 = cont11[M].mean()
                c22 = cont22[M].mean()
                sp.header['CONT11'] = c11
                sp.header['CONT22'] = c22
            elif cont == 'front':
                c11 = 2.7315
                c22 = 2.7315

            gg = [4,13,-3,20,68,1, c11, c22]*2
            gg[8] = 5
            gg[12] = 60
            fixed = np.array([F,F,T,T,F,F,T,T]*2)

            if ncomp == 1:
                sp.specfit(fittype='formaldehyde_radex',guesses=gg[:8],
                           fixed=fixed[:8], multifit=True,quiet=True,verbose=False,
                           limits=limits[:8],
                           limited=limited[:8],
                           use_window_limits=False, fit_plotted_area=False)
            elif ncomp == 2:
                sp.specfit(fittype='formaldehyde_radex',guesses=gg,
                           fixed=fixed, multifit=True,quiet=True,verbose=False,
                           limits=limits,
                           limited=limited,
                           use_window_limits=False, fit_plotted_area=False)
            else:
                raise ValueError('ncomp must be 1 or 2')

            sp.specname = name
            print("ap%i %s: ni=%i, X^2%0.1f, X^2/n=%0.1f, n=%0.1f" % (ii,
                                                                      sp.specname,
                                                                      sp.specfit.fitter.mp.niter,
                                                                      sp.specfit.chi2,
                                                                      sp.specfit.chi2/sp.specfit.dof,
                                                                      sp.specfit.parinfo.DENSITY0.value))

            #if sp.specfit.chi2 / sp.specfit.dof > 2:
            #    gg[4] = sp.specfit.parinfo.CENTER0.value
            #    fixed[4] = True
            #    sp.specfit(fittype='formaldehyde_radex',guesses=gg,
            #               fixed=fixed, multifit=True,quiet=False,verbose=True,
            #               limits=limits,
            #               limited=limited,
            #               use_window_limits=False, fit_plotted_area=False)

            plotitem(sp, ii)

            sp.plotter.savefig(outpfx+'_aperture_%s_%s.pdf' % (name,cont))
            sp.write(outpfx+"_aperture_%s.fits" % name)

    return spectra

def filaments_right():
    return do_indiv_fits("/Users/adam/work/w51/dense_filament_spectral_apertures.reg",
                         'spectralfits/spectralfits_70kmscloud')

def filaments_left():
    return do_indiv_fits("/Users/adam/work/w51/filament_leftside_spectral_apertures.reg",
                         'spectralfits/spectralfits_70kmscloudLeft')

def w51main():
    return do_indiv_fits("/Users/adam/work/w51/w51main_spectral_apertures.reg",
                         'spectralfits/spectralfits_w51main',
                         ncomp=2)

def maus():
    return do_indiv_fits("/Users/adam/work/w51/maus_spectral_apertures.reg",
                         'spectralfits/spectralfits_maus')

def middlechunk():
    spectra = do_indiv_fits("/Users/adam/work/w51/middlechunk_spectral_apertures.reg",
                            'spectralfits/spectralfits_63kmscloud')
    return spectra

def do_all_h2co():
    return (filaments_right() +
            filaments_left() +
            middlechunk() +
            maus() +
            w51main()
            )

def fit_twocomp_foregroundbackground(sp):
    """
    Fit a two-component model where one component is in the foreground and one
    in the background relative to the continuum
    """
    
    limits,limited = modelpars()
    c11,c22 = sp.header['CONT11'],sp.header['CONT22']
    gg = [4,13,-3,20,69,1, 2.7315, 2.7315]*2
    gg[8] = 5
    gg[12] = 62
    gg[14] = c11
    gg[15] = c22
    fixed = np.array([F,F,T,T,F,F,T,T]*2)

    sp.specfit(fittype='formaldehyde_radex', guesses=gg, fixed=fixed,
               multifit=True, limits=limits, limited=limited,
               use_window_limits=False, fit_plotted_area=False)

    plotitem(sp, 0, show_components=True)

    print("ni=%i, X^2%0.1f, X^2/n=%0.1f, n=%0.1f" %
          (sp.specfit.fitter.mp.niter, sp.specfit.chi2,
           sp.specfit.chi2/sp.specfit.dof, sp.specfit.parinfo.DENSITY0.value))
    print(sp.specfit.parinfo)

    return sp

def do_ap2_middlechunk(spectra_middle=None):
    if spectra_middle is None:
        spectra_middle = middlechunk()

    sp = fit_twocomp_foregroundbackground(spectra_middle[2])

    sp.plotter.savefig('spectralfits/middlechunk_63kms_aperture2_twocomponentfit.pdf')

    """
    Left region selection unchanged.  xminpix, xmaxpix: 0,402
    ni=16, X^22047.3, X^2/n=5.2, n=4.2
    Param #0     DENSITY0 =      4.20893 +/-        0.100445   Range:     [1,8]
    Param #1      COLUMN0 =       12.372 +/-       0.0287282   Range:   [11,16]
    Param #2   ORTHOPARA0 =           -3 (fixed)  Range:[-3,0.477121]
    Param #3 TEMPERATURE0 =           20 (fixed)  Range:    [5,55]
    Param #4      CENTER0 =      68.2007 +/-        0.097412   Range:   [50,80]
    Param #5       WIDTH0 =      1.36213 +/-        0.104359   Range:   [0,inf)
    Param #6 TBACKGROUND0 =       2.7315 (fixed)  Range:[2.73,inf)
    Param #7 TBACKGROUND1 =       2.7315 (fixed)  Range:[2.73,inf)
    Param #8     DENSITY1 =      4.69134 +/-      0.00987699   Range:     [1,8]
    Param #9      COLUMN1 =      12.8271 +/-      0.00661199   Range:   [11,16]
    Param #10   ORTHOPARA1 =           -3 (fixed)  Range:[-3,0.477121]
    Param #11 TEMPERATURE1 =           20 (fixed)  Range:    [5,55]
    Param #12      CENTER1 =      61.9618 +/-       0.0123433   Range:   [50,80]
    Param #13       WIDTH1 =      1.73291 +/-        0.013299   Range:   [0,inf)
    Param #14 TBACKGROUND2 =      12.6773 (fixed)  Range:[2.73,inf)
    Param #15 TBACKGROUND3 =      3.56966 (fixed)  Range:[2.73,inf)
    """

class LoadCOCubes:

    def __call__(self):

        datapath = '/Users/adam/work/w51/'
        if not hasattr(self,'co1332'):
            self.co1210 = pyspeckit.Cube(datapath+'w51_12co10_carpenter_rightaxes.fits')
            self.co1332 = pyspeckit.Cube(datapath+'13co_final_cube_c.fits')
            self.co1232 = pyspeckit.Cube(datapath+'12co_final_cube_c.fits')
            self.co1321 = pyspeckit.Cube(datapath+'w51_bieging_13co32.fits')
            self.co1221 = pyspeckit.Cube(datapath+'w51_bieging_12co32.fits')
            self.co1310 = pyspeckit.Cube(datapath+'grs-50-cube_supersampledh2cogrid.fits')

        cocubes = {'$^{12}$CO 1-0': self.co1210,
                   '$^{13}$CO 1-0': self.co1310,
                   '$^{13}$CO 2-1': self.co1321,
                   '$^{13}$CO 3-2': self.co1332,
                   '$^{12}$CO 2-1': self.co1221,
                   '$^{12}$CO 3-2': self.co1232,}
        
        for c in cocubes.values():
            if c.xarr.frame is None:
                c.xarr.frame = 'LSRK'

        return cocubes

load_cocubes = LoadCOCubes()

class LoadRRLCubes:
    def __call__(self):

        datapath = '/Users/adam/work/h2co/maps/W51/'
        if not hasattr(self,'h77a'):
            self.h77a = pyspeckit.Cube(datapath+'W51_h77a_pyproc_cube_supersampled_sub.fits')
            self.h110a = pyspeckit.Cube(datapath+'W51_h110alpha_cube_supersampled_sub.fits')

        hcubes = {'H110$\\alpha$': self.h110a,
                  'H77$\\alpha$': self.h77a}

        return hcubes

load_rrlcubes = LoadRRLCubes()

class LoadHICubes:
    def __call__(self):

        datapath = '/Volumes/128gbdisk/w51/'
        if not hasattr(self,'h_one'):
            self.h_one = pyspeckit.Cube(datapath+'MOS_049_Tb_reprojsupersampled.fits')

        hicubes = {'HI': self.h_one}

        for c in hicubes.values():
            if c.xarr.frame is None:
                c.xarr.frame = 'LSRK'


        return hicubes

load_HIcubes = LoadHICubes()

def do_cospectra():
    cocubes = load_cocubes()

    order = ['$^{13}$CO 1-0',
             '$^{13}$CO 2-1',
             '$^{13}$CO 3-2',
             '$^{12}$CO 3-2',
             '$^{12}$CO 1-0',
             '$^{12}$CO 2-1',
             ]

    pl.close(1)
    fig = pl.figure(1, figsize=(12,8))

    colors = ['k','r','g','b','m','#11AAAA']

    for regfn in regfiledict:
        regions = pyregion.open(regfn)

        for r in regions:
            spectra = [(c,cocubes[c].get_apspec(r.coord_list,
                                                coordsys=r.coord_format,
                                                wunit='degree'))
                       for c in order]

            fig.clf()

            for ii,(coline,sp) in enumerate(spectra):
                sp.xarr.convert_to_unit('km/s')
                sp.specname = r.attr[1]['text']
                sp.unit = 'K'
                if '12' in coline:
                    linestyle='steps--'
                else:
                    linestyle='steps-'
                sp.plotter(figure=fig, xmin=40, xmax=80, color=colors[ii], clear=False,
                           label=coline, linestyle=linestyle)

            sp.plotter.axis.set_ylabel("$T_A^*$ (K)")
            sp.plotter.savefig(regfiledict[regfn]+"_CO_aperture_%s.pdf" % sp.specname)
            sp.write(regfiledict[regfn]+"_CO_aperture_%s.fits" % sp.specname)

def do_rrlspectra():
    hcubes = load_rrlcubes()
    order = ['H77$\\alpha$', 'H110$\\alpha$']

    pl.close(1)
    fig = pl.figure(1, figsize=(12,8))

    colors = ['k','r','g','b','m','#11AAAA']

    for regfn in regfiledict:
        regions = pyregion.open(regfn)

        for r in regions:
            spectra = [(c,hcubes[c].get_apspec(r.coord_list,
                                               coordsys=r.coord_format,
                                               wunit='degree'))
                       for c in order]

            fig.clf()

            for ii,(hline,sp) in enumerate(spectra):
                sp.xarr.convert_to_unit('km/s')
                sp.specname = r.attr[1]['text']
                sp.unit = 'K'
                sp.specfit(fittype='gaussian', guesses=[1,60,5], multifit=True)
                pars = sp.specfit.modelpars
                pars[2] /= 5
                sp.plotter(figure=fig, xmin=40, xmax=80, color=colors[ii], clear=False,
                           label=hline)
                fakespec = pyspeckit.Spectrum(xarr=sp.xarr,
                                              data=sp.specfit.get_model_frompars(sp.xarr,
                                                                                 pars),
                                              header=sp.header)
                fakespec.plotter(axis=sp.plotter.axis, xmin=40, xmax=80,
                                 color=colors[ii], clear=False,
                                 linestyle='steps--')

            sp.plotter.axis.set_ylabel("$T_A^*$ (K)")
            sp.plotter.savefig(regfiledict[regfn]+"_RRL_aperture_%s.pdf" % sp.specname)
            sp.write(regfiledict[regfn]+"_RRL_aperture_%s.fits" % sp.specname)

def do_HIspectra():
    hcubes = load_HIcubes()
    order = ['HI']

    pl.close(1)
    fig = pl.figure(1, figsize=(12,8))

    colors = ['k','r','g','b','m','#11AAAA']

    for regfn in regfiledict:
        regions = pyregion.open(regfn)

        for r in regions:
            spectra = [(c,hcubes[c].get_apspec(r.coord_list,
                                               coordsys=r.coord_format,
                                               wunit='degree'))
                       for c in order]

            fig.clf()

            for ii,(hline,sp) in enumerate(spectra):
                sp.xarr.convert_to_unit('km/s')
                sp.specname = r.attr[1]['text']
                sp.unit = 'K'
                sp.plotter(figure=fig, xmin=40, xmax=80, color=colors[ii],
                           clear=False, label=hline)

            sp.plotter.axis.set_ylabel("$T_A^*$ (K)")
            sp.plotter.savefig(regfiledict[regfn]+"_HI_aperture_%s.pdf" % sp.specname)
            sp.write(regfiledict[regfn]+"_HI_aperture_%s.fits" % sp.specname)


def do_all3spectra():
    hcubes = load_rrlcubes()
    Horder = ['H77$\\alpha$', 'H110$\\alpha$']
    cocubes = load_cocubes()
    COorder = ['$^{13}$CO 1-0',
               '$^{13}$CO 2-1',
               '$^{13}$CO 3-2',
               '$^{12}$CO 3-2',
               '$^{12}$CO 1-0',
               '$^{12}$CO 2-1',
               ]
    for c in both.cubelist:
        c.xarr.convert_to_unit('km/s')
    h2cocubes = {'H$_2$CO 1-1':both.cubelist[0],
                 'H$_2$CO 2-2':both.cubelist[1]}
    h2coorder = ['H$_2$CO 2-2',
                 'H$_2$CO 1-1']
    pl.close(1)
    fig = pl.figure(1, figsize=(12,8))

    colors = ['k','r','g','b','m','#11AAAA']

    for regfn in regfiledict:
        regions = pyregion.open(regfn)

        for r in regions:
            fig.clf()
            subplots = [pl.subplot(3,1,i) for i in range(1,4)]

            for ax,cubes,order in zip(subplots,
                                      [h2cocubes,cocubes,hcubes],
                                      [h2coorder,COorder,Horder]):


                spectra = [(c,cubes[c].get_apspec(r.coord_list,
                                                  coordsys=r.coord_format,
                                                  wunit='degree'))
                           for c in order]

                for ii,(line,sp) in enumerate(spectra):
                    sp.xarr.convert_to_unit('km/s')
                    sp.specname = r.attr[1]['text']
                    sp.unit = 'K'
                    if 'alpha' in line:
                        sp.specfit(fittype='gaussian', guesses=[1,60,5], multifit=True)
                        pars = sp.specfit.modelpars
                        pars[2] /= 5
                        fakespec = pyspeckit.Spectrum(xarr=sp.xarr,
                                                      data=sp.specfit.get_model_frompars(sp.xarr,
                                                                                         pars),
                                                      header=sp.header)
                        fakespec.plotter(axis=ax, xmin=40, xmax=80,
                                         color=colors[ii], clear=False,
                                         linestyle='steps--')
                    linestyle = 'steps--' if '12' in line else 'steps-'
                    sp.plotter(figure=fig, xmin=40, xmax=80, color=colors[ii], clear=False,
                               label=line, linestyle=linestyle, axis=ax)

            fig.subplots_adjust(hspace=0)
            subplots[0].set_xticks([])
            subplots[1].set_xticks([])
            subplots[1].set_title("")
            subplots[2].set_title("")
            sp.plotter.axis.set_ylabel("$T_A^*$ (K)")
            sp.plotter.savefig(regfiledict[regfn]+"_allthree_aperture_%s.pdf" % sp.specname)

def do_all():
    do_all_h2co()
    do_cospectra()
    do_all3spectra()
    do_rrlspectra()
    do_HIspectra()
