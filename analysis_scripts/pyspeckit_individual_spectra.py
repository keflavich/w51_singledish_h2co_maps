from __future__ import print_function
import glob
from load_pyspeckit_cubes import both,T,F,cont11,cont22,cont11filename
import matplotlib as mpl
import pyspeckit
from astropy.io import fits
import astropy.table
import pyregion
import types
from pyspeckit.wrappers import fith2co
import numpy as np
import pylab as pl
import os
from paths import datapath,datapath_w51,figurepath,datapath_spectra

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
    
    limits = [(1, 8), # dens
              (11, 16), # col
              (-3, 0.47712125471966244), # opr
              (5, 55), # temp
              (50, 80), # velo
              (0, 0), # width
              (2.73, 0), # cont11
              (2.73, 0), # cont22
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

def plotitem(sp, ii=0, errstyle='fill', vrange=[40,80], resid=False,
             dolegend=False, refresh=False, show_components=True,
             residkwargs={}, **kwargs):
    sp.plot_special = types.MethodType(fith2co.plotter_override, sp, sp.__class__)
    plotkwargs = {'vrange':vrange,'fignum':ii+1,'reset_xlimits':True,
                  'annotate': False, 'errstyle':errstyle,
                  'show_components': show_components}
    plotkwargs.update(kwargs)

    if resid:
        plotkwargs['resid_overlay'] = True
        plotkwargs['residkwargs'] = {'zeroline': True,
                                     }
        plotkwargs['residkwargs'].update(residkwargs)

    sp.plot_special_kwargs = plotkwargs
    sp.plot_special(**plotkwargs)

    # plot cleanup junk
    sp.plotter.figure.subplots_adjust(hspace=0)
    # move h2co2-2 title down CANCELED by set_position below
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

    # move both titles to bottom right
    for label, ax in sp.axisdict.items():
        ax.title.set_position((0.85,0.07))

    if dolegend:
        for ax in sp.axisdict.values():
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        sp.specfit.annotate(bbox_to_anchor=(1.01,0.9), loc='upper left',
                            chi2='allthree', fontsize=14)
            
    if refresh:
        sp.plotter.refresh()

def initialize_table():
    parnames = ['DENSITY', 'COLUMN', 'ORTHOPARA',
                'TEMPERATURE', 'CENTER', 'WIDTH', 'TBACKGROUND0',
                'TBACKGROUND1',]
    colnames = (['name', 'component_number', 'frontback', 'frontbackbest',
                 'chi2', 'dof', 'opt_chi2', 'opt_red_chi2', 'ra', 'dec',
                 'radius'] +
                [x for y in zip(parnames, ['e'+p for p in parnames]) for x in y])
    dtypes = ['S20',np.int,'S5','S5'] + [np.float]*(7+len(parnames)*2)
    table = astropy.table.Table(names=colnames, dtypes=dtypes)
    return table

def add_parinfo_to_table(table, parinfo, chi2, dof, opt_chi2, opt_red_chi2, ra,
                         dec, radius, frontback="", frontbackbest="", name=""):
    new_row = [name]
    new_row.append(0)
    new_row += [frontback, frontbackbest, chi2, dof, opt_chi2, opt_red_chi2,
                ra, dec, radius]
    for pp in parinfo[:8]:
        new_row.append(pp.value)
        new_row.append(np.nan if pp.fixed else pp.error)
    table.add_row(new_row)

    if len(parinfo)>8:
        new_row = [name]
        new_row.append(1)
        new_row += [frontback, frontbackbest, chi2, dof, opt_chi2,
                    opt_red_chi2, ra, dec, radius]
        for pp in parinfo[8:]:
            new_row.append(pp.value)
            new_row.append(np.nan if pp.fixed else pp.error)
        table.add_row(new_row)

def do_indiv_fits(regfilename, outpfx, ncomp=1, dobaseline=False, table=None,
                  tableprefix="", **kwargs):
    regions = pyregion.open(regfilename)

    spectra = [both.get_apspec(r.coord_list, coordsys=r.coord_format,
                               wunit='degree')
               for r in regions]

    cont11hdu = fits.open(cont11filename)[0]

    pl.ioff()

    for ii,(sp,r) in enumerate(zip(spectra,regions)):

        if hasattr(ncomp,'__len__'):
            nc = ncomp[ii]
        else:
            nc = ncomp

        chi2, parinfo = {},{}
        sp.header['REGION'] = "{shape}({ra},{dec},{radius})".format(shape=r.name,
                                                                    ra=r.coord_list[0],
                                                                    dec=r.coord_list[1],
                                                                    radius=r.coord_list[2])

        ## this is pretty hacky =(
        #if dobaseline:
        #    spdict = fith2co.BigSpectrum_to_H2COdict(sp)
        #    for n,s in spdict.items():
        #        s.baseline(exclude=[2,8,43,75])
        #    sp = pyspeckit.Spectra(spdict.values())
        #    spectra[ii] = sp



        for cont in ('front','back',):
            print("Continuum in the %s" % cont)
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

            dofit(sp, c11, c22, nc, **kwargs)

            sp.specname = name
            sp.header['OBJECT'] = name
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

            sp.plotter.autorefresh=False
            plotitem(sp, ii)

            sp.plotter.savefig(outpfx+'_aperture_%s_%s.pdf' % (name,cont))

            plotitem(sp, ii, dolegend=True)

            # seriously, something aint'n't right here
            pl.figure(sp.plotter.figure.number)
            pl.savefig(outpfx+'_aperture_%s_%s_legend.pdf' % (name,cont),
                       bbox_extra_artists=[sp.specfit.fitleg])

            # for writing to file, select the best-fit
            chi2[cont] = sp.specfit.chi2
            parinfo[cont] = sp.specfit.parinfo
            if cont == 'back': # second one...
                if sp.specfit.chi2 > chi2['front']:
                    print("Back chi^2 > front chi^2.  Replacing parameters with Continuum in the Front")
                    best = 'front'
                else:
                    best = 'back'
                # set the previous, unset to match this...
                # since we're guaranteed to be in the 2nd of 2 in a loop here,
                # no danger
                table[-1]['frontbackbest'] = best
                if nc == 2:
                    # or previous *two* if twocomp
                    table[-2]['frontbackbest'] = best
                # at this stage, 'frontbackbest' should be universally assigned....
                if np.any(table['frontbackbest'] == ''):
                    import ipdb; ipdb.set_trace()
            else:
                best = ''

            if table is not None:
                add_parinfo_to_table(table, sp.specfit.parinfo, sp.specfit.chi2,
                                     sp.specfit.dof, sp.specfit.optimal_chi2(reduced=False),
                                     sp.specfit.optimal_chi2(reduced=True),
                                     ra=r.coord_list[0],
                                     dec=r.coord_list[1],
                                     radius=r.coord_list[2],
                                     frontback=cont,
                                     frontbackbest=best,
                                     name=tableprefix+sp.specname)

            if best != '':
                sp.specfit.parinfo = parinfo[best]
                sp.specfit.chi2 = chi2[best]

            sp.write(outpfx+"_aperture_%s.fits" % name)

    table.write(os.path.join(datapath_spectra,
                             tableprefix+"spectralfit_table.ipac"),
                format='ascii.ipac')

    pl.ion()

    return spectra

def get_last_chi2_fromheader(header):
    """
    This is a fragile way to recover the best-fit chi2
    """
    for line in header:
        if 'Chi^2' in line:
            chi2 = line.split()[4]
    return chi2

def reload_fits():
#    spectra_front = []
#    spectra_back = []
    spectra = []
    for fn in glob.glob("spectralfits/*fits"):
        if not ('CO' in fn or 'HI' in fn or 'RRL' in fn):
            sp = pyspeckit.Spectrum(fn)
            spectra.append(sp)

    return spectra
#            if 'front' in fn:
#                spectra_front.append(fn)
#            else:
#                spectra_back.append(fn)
#
#    spectra = []
#    for spb,spf in zip(spectra_front,spectra_back):
#        if spb.specname != spf.specname:
#            raise ValueError("Mismatch between spectral names")
#        chi2_b = get_last_chi2_fromheader(spb.header)
#        chi2_f = get_last_chi2_fromheader(spf.header)
#        spectra.append(spb if chi2_b < chi2_f else chi2_f)
#
#    return spectra,spectra_front, spectra_back


limits, limited = modelpars()

def dofit(sp, c11, c22, ncomp, fixed=np.array([F,F,T,T,F,F,T,T]*2),
          limits=limits, limited=limited, verbose=False,
          c11b=None, c22b=None,
          vguesses=[68,60]):
    gg = [4,13,-3,20,vguesses[0],1, c11, c22]*2
    gg[8] = 5
    gg[12] = vguesses[1]

    # ALWAYS fix the continuum
    fixed[6::8] = True
    fixed[7::8] = True

    if c11b is not None:
        gg[14] = c11b
    if c22b is not None:
        gg[15] = c22b

    if ncomp == 1:
        sp.specfit(fittype='formaldehyde_radex',guesses=gg[:8],
                   fixed=fixed[:8], multifit=True,quiet=True,verbose=verbose,
                   limits=limits[:8],
                   limited=limited[:8],
                   use_window_limits=False, fit_plotted_area=False)
    elif ncomp == 2:
        sp.specfit(fittype='formaldehyde_radex',guesses=gg,
                   fixed=fixed, multifit=True,quiet=True,verbose=verbose,
                   limits=limits,
                   limited=limited,
                   use_window_limits=False, fit_plotted_area=False)
    else:
        raise ValueError('ncomp must be 1 or 2')


def filaments_right(table=None):
    if table is None:
        table = initialize_table()
    ncomp = [1,1,1,2,2,1,1]
    return do_indiv_fits("/Users/adam/work/w51/dense_filament_spectral_apertures.reg",
                         'spectralfits/spectralfits_70kmscloud', ncomp=ncomp,
                         tableprefix="filamentsright_",
                         table=table)

def filaments_left(table=None):
    if table is None:
        table = initialize_table()
    ncomp = [2,2,2,2,1,2,1]
    outpfx = 'spectralfits/spectralfits_70kmscloudLeft'
    spectra = do_indiv_fits("/Users/adam/work/w51/filament_leftside_spectral_apertures.reg",
                            outpfx, ncomp=ncomp,
                            tableprefix="filamentsleft_",
                            table=table)
    dofit(spectra[0], spectra[0].header['CONT11'], spectra[0].header['CONT22'], vguesses=[64, 68], ncomp=2)
    plotitem(spectra[0], 0, dolegend=True)
    pl.figure(spectra[0].plotter.figure.number)
    pl.savefig(outpfx+'_aperture_%s_%s_legend.pdf' % (spectra[0].specname,'back'),
               bbox_extra_artists=[spectra[0].specfit.fitleg])

    return spectra
    

def w51main(table=None):
    if table is None:
        table = initialize_table()
    return do_indiv_fits("/Users/adam/work/w51/w51main_spectral_apertures.reg",
                         'spectralfits/spectralfits_w51main',
                         ncomp=2,
                         tableprefix="w51main_",
                         table=table)

def maus(table=None):
    if table is None:
        table = initialize_table()

    limits, limited = modelpars()
    limits[12] = (40, 80)
    limits[13] = (1,3)
    limits[5] = (1,3)
    limited[5] = (True,True)
    limited[13] = (True,True)

    ncomp = [2,2,2,2,2,2]
    # front2nd means "continuum in front of 2nd component?"
    # front2nd[3] is sketchy, 5 is uncertain
    front2nd = [True, False, True, False, True, True]

    outpfx = 'spectralfits/spectralfits_maus'

    spectra = do_indiv_fits("/Users/adam/work/w51/maus_spectral_apertures.reg",
                            outpfx,
                            ncomp=ncomp,
                            limits=limits,
                            limited=limited,
                            vguesses=[64,50],
                            tableprefix="maus_",
                            table=table)


    for ii,sp in enumerate(spectra):

        if front2nd[ii]:
            dofit(sp, sp.header['CONT11'], sp.header['CONT22'],
                  vguesses=[68,50],
                  limits=limits,
                  limited=limited,
                  c11b=2.7315, c22b=2.7315,
                  ncomp=ncomp[ii])

            plotitem(sp, ii, dolegend=True)

            pl.figure(sp.plotter.figure.number)
            pl.savefig(outpfx+'_aperture_%s_%s_legend.pdf' % (sp.specname,'both'),
                       bbox_extra_artists=[sp.specfit.fitleg])


def middlechunk(table=None):
    if table is None:
        table = initialize_table()
    ncomp = [2,2,1,1,1,1,1]
    ncomp = [2,2,2,2,1,1,2] # in 3,4,7, second comp is at 70 kms
    outpfx = 'spectralfits/spectralfits_63kmscloud'
    spectra = do_indiv_fits("/Users/adam/work/w51/middlechunk_spectral_apertures.reg",
                            outpfx,
                            ncomp=ncomp,
                            tableprefix="middlechunk_",
                            table=table)

    sp = spectra[0]
    dofit(sp, sp.header['CONT11'], sp.header['CONT22'],
          vguesses=[59,63], ncomp=2)

    plotitem(sp, 0, dolegend=True)

    pl.figure(sp.plotter.figure.number)
    pl.savefig(outpfx+'_aperture_%s_%s_legend.pdf' % (sp.specname,'back'),
               bbox_extra_artists=[sp.specfit.fitleg])

    for ii,sp in enumerate(spectra):
        if ii == 0:
            continue
        dofit(sp, sp.header['CONT11'], sp.header['CONT22'],
              vguesses=[60,68],
              c11b=2.7315, c22b=2.7315,
              ncomp=ncomp[ii])

        plotitem(sp, ii, dolegend=True)

        pl.figure(sp.plotter.figure.number)
        pl.savefig(outpfx+'_aperture_%s_%s_legend.pdf' % (sp.specname,'both'),
                   bbox_extra_artists=[sp.specfit.fitleg])

    return spectra

def do_all_h2co():
    table = initialize_table()
    result = (filaments_right(table=table) + filaments_left(table=table) +
              middlechunk(table=table) + maus(table=table) +
              w51main(table=table))
    table.write(os.path.join(datapath_spectra, "spectralfit_table.ipac"),
                format='ascii.ipac')
    table2 = add_tex_tau_to_table(table)
    table2.write(os.path.join(datapath_spectra, "spectralfit_table_withtextau.ipac"),
                 format='ascii.ipac')
    return result,table

def split_table(table):
    comp0 = table[table['component_number']==0]
    comp0best = comp0[comp0['frontback'] == comp0['frontbackbest']]
    comp1 = table[table['component_number']==1]
    comp1best = comp1[comp1['frontback'] == comp1['frontbackbest']]
    return comp0,comp1,comp0best,comp1best

def split_table_bycolumn(table):
    best = table[table['frontback'] == table['frontbackbest']]
    deepest = [t1 if ((t2['component_number'] == 1 and t2['COLUMN'] < t1['COLUMN'])
                      or t2['component_number'] == 0) else t2
               for t1,t2 in zip(best[:-1],best[1:])
               if (t2['component_number']==1 or
                   (t1['component_number']==0 and t2['component_number']==0))
               ]
    return astropy.table.Table(rows=deepest, names=table.colnames)

def split_table_byvelo(table):
    best = table[table['frontback'] == table['frontbackbest']]
    highv = [t1 if ((abs(t2['CENTER']-68) > abs(t1['CENTER']-68)) or
                    t2['component_number'] == 0) else t2
             for t1,t2 in zip(best[:-1],best[1:])
             if (t2['component_number']==1 or
                 (t1['component_number']==0 and t2['component_number']==0))
             ]
    return astropy.table.Table(rows=highv, names=table.colnames)

def table_to_reg(table, regfilename, system='galactic'):
    with open(regfilename, 'w') as outf:
        outf.write(system+"\n")
        for row in table:
            rowdict = dict(zip(row.colnames, row.data))
            #rowdict['radius'] = "{0}\"".format(rowdict['radius']*3600)
            outf.write("circle({ra},{dec},{radius}) # text={{{DENSITY:0.2f}}}\n".format(**rowdict))


def add_tex_tau_to_table(table):
    import pyradex
    R = pyradex.Radex(species='o-h2co_troscompt', h2column=1e21, abundance=10**-8.5)

    newcols = {'tau1':[],
               'tau2':[],
               'tex1':[],
               'tex2':[]}
    for t in table:
        R.temperature = t['TEMPERATURE']
        R.column = 10**t['COLUMN']
        orthofrac = 10**t['ORTHOPARA']/(1+10**t['ORTHOPARA'])
        R.density = {'oH2': 10**t['DENSITY']*orthofrac,
                     'pH2': 10**t['DENSITY']*(1-orthofrac)}
        R.run_radex()
        newcols['tex1'].append(R.tex[0].value)
        newcols['tex2'].append(R.tex[2].value)
        newcols['tau1'].append(R.tau[0])
        newcols['tau2'].append(R.tau[2])

    for k in newcols:
        table.add_column(astropy.table.Column(name=k, data=newcols[k]))

    return table

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

    print("ni=%i, X^2=%0.1f, X^2/n=%0.1f, n=%0.1f" %
          (sp.specfit.fitter.mp.niter, sp.specfit.chi2,
           sp.specfit.chi2/sp.specfit.dof, sp.specfit.parinfo.DENSITY0.value))
    print(sp.specfit.parinfo)

    return sp

def paperfigure_filament_demonstrate_frontback(fixedTO=True):
    """
    fixedTO: fixed temerature/orthopara
    """
    with mpl.rc_context(fname='/Users/adam/.matplotlib/pubfiguresrc'):
        sp = pyspeckit.Spectrum('spectralfits/spectralfits_70kmscloud_aperture_ap6.fits')
        from load_pyspeckit_cubes import formaldehyde_radex_fitter
        sp.Registry.add_fitter('formaldehyde_radex',formaldehyde_radex_fitter,8,multisingle='multi')


        c11 = sp.header["CONT11"]
        c22 = sp.header["CONT22"]
        gg = [4,13,-3,20,68,1, c11, c22]*2
        gg[8] = 5
        gg[12] = 60
        if fixedTO:
            fixed = np.array([F,F,T,T,F,F,T,T]*2)
        else:
            fixed = np.array([F,F,F,F,F,F,T,T]*2)
        limits, limited = modelpars()

        # CHEATER baselining:
        sp.data[sp.xarr.as_unit('GHz') < 5] -= 0.05

        sp.specfit(fittype='formaldehyde_radex',guesses=gg[:8],
                   fixed=fixed[:8], multifit=True,quiet=True,verbose=False,
                   limits=limits[:8],
                   limited=limited[:8],
                   use_window_limits=False, fit_plotted_area=False, plot=False,
                   clear=False)

        print("Continuum in the back.")
        print(sp.specfit.parinfo)
        print(sp.specfit.chi2, sp.specfit.optimal_chi2(reduced=False), sp.specfit.optimal_chi2(reduced=False, threshold=0))
        print(sp.specfit.chi2/sp.specfit.dof, sp.specfit.optimal_chi2(reduced=True), sp.specfit.optimal_chi2(reduced=True, threshold=0))
        #parinfo_back = sp.specfit.parinfo
        #model_back = sp.specfit.model
        plotitem(sp,0,clear=True, vrange=[50,90])
        plotitem(sp,1,clear=True, vrange=[50,90], dolegend=True)
        plotitem(sp,3,clear=True, vrange=[50,90], resid=True,
                 resid_yoffsets={'oneone':0.15, 'twotwo': 0.04},
                 residkwargs={'color':'r'})

        gg[6] = 2.7315
        gg[7] = 2.7315
        sp.specfit(fittype='formaldehyde_radex',guesses=gg[:8],
                   fixed=fixed[:8], multifit=True,quiet=True,verbose=False,
                   limits=limits[:8],
                   limited=limited[:8],
                   use_window_limits=False, fit_plotted_area=False, plot=False,
                   clear=False)
        #parinfo_front = sp.specfit.parinfo
        #model_front = sp.specfit.model

        plotitem(sp, 3, clear=False, plot_fit_kwargs=dict(composite_fit_color='g'),
                 errstyle='none', vrange=[50,90], resid=True,
                 resid_yoffsets={'oneone':0.15, 'twotwo': 0.04},
                 residkwargs={'color':'g'})
        plotitem(sp, 0, clear=False, plot_fit_kwargs=dict(composite_fit_color='g'),
                 errstyle='none', vrange=[50,90])
        plotitem(sp,2,clear=True, vrange=[50,90],dolegend=True)
        #sp.specfit.plot_fit(composite_fit_color='g')
        #plotitem(sp,0,clear=False)
        print("Continuum in the front.")
        print(sp.specfit.parinfo)
        print(sp.specfit.chi2, sp.specfit.optimal_chi2(reduced=False), sp.specfit.optimal_chi2(reduced=False, threshold=0))
        print(sp.specfit.chi2/sp.specfit.dof, sp.specfit.optimal_chi2(reduced=True), sp.specfit.optimal_chi2(reduced=True, threshold=0))

        pl.figure(1)
        if fixedTO:
            pl.savefig('spectralfits/spectralfits_70kmscloud_aperture_ap6_modelcomparison.pdf')
        else:
            pl.savefig('spectralfits/spectralfits_70kmscloud_aperture_ap6_modelcomparison_freed.pdf')

        pl.figure(4)
        pl.subplot(211).set_ylim(-0.45,0.25)
        pl.subplot(212).set_ylim(-0.08,0.06)
        if fixedTO:
            pl.savefig('spectralfits/spectralfits_70kmscloud_aperture_ap6_modelcomparison_withresiduals.pdf')
        else:
            pl.savefig('spectralfits/spectralfits_70kmscloud_aperture_ap6_modelcomparison_withresiduals_freed.pdf')

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
                # fit twice: first time, sets the errors
                sp.specfit(fittype='gaussian', guesses=[1,60,5], multifit=True)
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
                rrlid = 'H77a' if '77' in hline else '110'
                sp.write(regfiledict[regfn]+"_RRL%s_aperture_%s.fits" % (rrlid,sp.specname))

            sp.plotter.axis.set_ylabel("$T_A^*$ (K)")
            sp.plotter.savefig(regfiledict[regfn]+"_RRL_aperture_%s.pdf" % sp.specname)

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
