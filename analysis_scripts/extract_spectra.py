from paths import datapath,rpath
from paths import figurepath
import os
import pyspeckit
import pyregion
import pylab as pl

scube11 = pyspeckit.Cube(datapath+'W51_H2CO11_cube_supersampled_sub.fits')
taucube11 = pyspeckit.Cube(datapath+'W51_H2CO11_taucube_supersampled.fits')
scube11.xarr.convert_to_unit('km/s')
#taucube11_13 = pyspeckit.Cube(datapath+'W51_H213CO_taucube.fits')
#scube11_13 = pyspeckit.Cube(datapath+'W51_H213CO_cube_sub.fits')
#scube11_13.xarr.convert_to_unit('km/s')
scube22 = pyspeckit.Cube(datapath+'W51_H2CO22_pyproc_cube_lores_supersampled_sub.fits')
taucube22 = pyspeckit.Cube(datapath+'W51_H2CO22_pyproc_taucube_lores_supersampled.fits')
scube22.xarr.convert_to_unit('km/s')

regfn = rpath('w51main_spectral_apertures.reg')
regions = pyregion.open(regfn)

spectra11 = {}
spectra22 = {}

# Loading
for ii, reg in enumerate(regions):
    sp11 = scube11.get_apspec(reg.coord_list, coordsys=reg.coord_format, wunit='degree')
    sp22 = scube22.get_apspec(reg.coord_list, coordsys=reg.coord_format, wunit='degree')
    taup11 = taucube11.get_apspec(reg.coord_list, coordsys=reg.coord_format, wunit='degree')
    taup22 = taucube22.get_apspec(reg.coord_list, coordsys=reg.coord_format, wunit='degree')

    sp11.specname = reg.attr[1]['text']
    sp22.specname = reg.attr[1]['text']
    taup11.specname = reg.attr[1]['text']
    taup22.specname = reg.attr[1]['text']
    spectra11[ii] = taup11
    spectra22[ii] = taup22
    #sp11.write(outpath+'H2CO_11_%s.fits' % sp11.specname.replace(" ","_"))
    #sp22.write(outpath+'H2CO_22_%s.fits' % sp22.specname.replace(" ","_"))

twotwooffset = -0.025

# Plotting
for ii in spectra11:
    sp11,sp22 = [s[ii] for s in (spectra11,spectra22)]

    fig = pl.figure(ii)
    pl.clf()
    sp11.plotter(figure=fig, clear=True,zorder=50,xmin=-75,xmax=150, linewidth=2,
                 alpha=0.6)
    sp22.plotter(axis=sp11.plotter.axis, clear=False, color='r',
                 offset=twotwooffset, use_window_limits=True,zorder=25,
                 linewidth=2)
    sp11.plotter.axis.hlines([twotwooffset,0.0],
                             -100,100,color='purple',linestyle='--',alpha=0.5,zorder=5,
                             linewidth=2)
    sp11.plotter.axis.set_ylim(min([sp22.data.min()+twotwooffset]),
                               max([sp11.data.max()]))
    sp11.plotter.axis.set_xlim(40,90)

    sp11.plotter.refresh()

    ratio = sp11.copy()
    ratio.data = sp11.data/sp22.data
    ratio.units=r'$\tau$ ratio'

    inset = pl.axes([0.65,0.65,0.25,0.25])
    xmin,xmax = 40,90
    ratio.plotter(axis=inset, xmin=xmin, xmax=xmax, ymin=0,ymax=20, 
                  linewidth=2, alpha=0.6)
    ratio.plotter.axis.set_xlim(40,90)

    if 'off' not in sp11.specname.lower():
        sp11.plotter.savefig(figurepath+sp11.specname.replace(" ","_")+"_11_22_spectra.pdf")

