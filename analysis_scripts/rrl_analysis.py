import aplpy
from astropy.io import fits
import pylab as pl
import itertools
import multiplot
import mpl_plot_templates
import numpy as np
from agpy import fit_a_line
import astropy.units as u
import common_constants
aobeam,gbbeam = common_constants.beams()
from common_constants import h2co11freq,h2co22freq,etamb_77

pl.mpl.rcParams['axes.color_cycle'] = ["#"+x for x in "348ABD, 7A68A6, A60628, 467821, CF4457, 188487, E24A33".split(", ")]

# Have to convert to Jy to compare fairly
fwhm = np.sqrt(8*np.log(2))
ktojy111 = (1*u.K).to(u.Jy,u.brightness_temperature((2*np.pi*(aobeam/fwhm)**2), h2co11freq))
ktojy77 = (1*u.K).to(u.Jy,u.brightness_temperature((2*np.pi*(gbbeam/fwhm)**2), h2co22freq))

aofn = 'W51_Halpha_6cm_cube_supersampled.fits'
h112i = fits.getdata(aofn.replace("cube","integrated")) * ktojy111
h112c = fits.getdata('W51_h112alpha_cube_supersampled_continuum.fits') * ktojy111
h112h = fits.getheader('W51_h112alpha_cube_supersampled_continuum.fits')

h77iname = 'W51_h77a_pyproc_integrated_supersampled.fits'
h77i = fits.getdata(h77iname) * ktojy77 / etamb_77
h77cname = 'W51_h77a_pyproc_cube_supersampled_continuum.fits'
h77cname = 'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits'
h77c = fits.getdata(h77cname) * ktojy77 / etamb_77
h77h = fits.getheader(h77cname)

h77lc = h77i/h77c * (h77i.value > 0.01)
h112lc = h112i/h112c * (h112i.value > 0.01)
rrlmask = (h77i.value > 0.01) * (h112i.value > 0.01)
h77th112 = h77i/h112i * rrlmask
c2cmtc6cm = h77c/h112c * rrlmask

figs = []

titles = {'h77lc':r'H77$\alpha$ Line/Continuum',
          'h112lc':r'H112$\alpha$ Line/Continuum',
          'h77th112':r'H77$\alpha$ / H112$\alpha$',
          'c2cmtc6cm':r'$S_{15 GHz}/S_{5 GHz}$'}

for ii,img in enumerate('h77lc,h112lc,h77th112,c2cmtc6cm'.split(',')):
    name = img
    img = locals()[name].value
    hdu = fits.PrimaryHDU(data=img, header=h112h)
    pl.figure(ii)
    pl.clf()
    F = aplpy.FITSFigure(hdu, figure=pl.figure(ii))
    F.show_colorscale(vmin=0,vmax=0.1 if 'lc' in name else 2.5)
    F.recenter(49.235,-0.303,width=0.952,height=0.433)
    F.add_colorbar()
    F._ax1.set_title(titles[name])
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    #F.set_tick_xspacing(0.3)
    figs.append(F)
    F.save('/Users/adam/work/h2co/maps/paper/figures/ratiomap_cont_%s.pdf' % name)
    F.save('/Users/adam/work/h2co/maps/paper/figures/ratiomap_cont_%s.png' % name)

pl.figure(6)
pl.clf()
F = aplpy.FITSFigure(h77cname, figure=pl.figure(6), convention='calabretta')
F.show_grayscale(vmin=-0.1,vmid=-0.4,vmax=8,invert=True,stretch='log')
F.recenter(49.235,-0.303,width=0.952,height=0.433)
F.show_contour(hdu, levels=[0.3,0.5,0.7], colors=['b',(0.2,1,0.4,0.8),'r'])
F.add_colorbar()

pl.figure(7)
pl.clf()
F = aplpy.FITSFigure(h77iname, figure=pl.figure(7), convention='calabretta')
F.show_grayscale(vmin=-0.1/23,vmid=-0.4/23,vmax=8./23,invert=True,stretch='log')
F.recenter(49.235,-0.303,width=0.952,height=0.433)
F.show_contour(hdu, levels=[0.3,0.5,0.7], colors=['b',(0.2,1,0.4,0.8),'r'])
F.add_colorbar()

keys = [r'H77$\alpha$',
       r'$S_{15 GHz}$',
       r'H112$\alpha$',
       r'$S_{5 GHz}$']

data = {r'H77$\alpha$':  h77i,
        r'$S_{15 GHz}$': h77c,
        r'H112$\alpha$': h112i,
        r'$S_{5 GHz}$':  h112c}

combs = itertools.combinations(keys,2)

xpos = {r'H77$\alpha$':0,
        r'$S_{15 GHz}$':1,
        r'H112$\alpha$':2,}

ypos = {r'$S_{15 GHz}$':0,
        r'H112$\alpha$':1,
        r'$S_{5 GHz}$':2}

pl.figure(5)
pl.clf()
mp = multiplot.multipanel(dims=(3,3),diagonal=False,figID=5)
for ii,(xd,yd) in zip(mp.grid.keys(),combs):
    print mp.axis_number(ypos[yd],xpos[xd],), xpos[xd], ypos[yd]
    ax = mp.grid[mp.axis_number(ypos[yd],xpos[xd])]
    #ax.plot(data[xd][rrlmask],data[yd][rrlmask],'.')
    mpl_plot_templates.adaptive_param_plot(data[xd][rrlmask],
                                           data[yd][rrlmask],
                                           bins=30,
                                           threshold=5,
                                           fill=False,
                                           alpha=0.8,
                                           axis=ax,
                                           cmap=pl.mpl.cm.spectral)
    axlims = ax.axis()
    factor = fit_a_line.total_least_squares(data[xd][rrlmask],data[yd][rrlmask],intercept=False)
    ax.plot(np.linspace(0,20),factor*np.linspace(0,20),'k--',linewidth=2,alpha=0.5,label="$y=%0.2fx$" % factor)
    ax.axis(axlims) # reset plot limits

    ax.legend(loc='best',fontsize=16)

    # tweaks to prevent tick overlaps
    if xpos[xd] > 0:
        ax.set_yticklabels([])
        ax.set_ylabel("")
    else:
        ax.set_ylabel(yd)
    ax.set_xlabel(xd)
    if ypos[yd] > 0:
        ax.set_yticks(ax.get_yticks()[:-1])
        #if ypos[yd] == 2:
        #    ax.set_yticks(ax.get_yticks()[:-1])
    if ypos[yd] < 2:
        ax.set_xticklabels([])
    if xpos[xd] < 2:
        ax.set_xticks(ax.get_xticks()[:-1])
    if len(ax.get_yticks()) > 6:
        ax.set_yticks(ax.get_yticks()[::2])
    if len(ax.get_xticks()) > 6:
        ax.set_xticks(ax.get_xticks()[::2])
#for ii,(xd,yd) in enumerate(combs):
#    print ii
#    ax = pl.subplot(3,3,ii+1)
#    ax.plot(xd[rrlmask],yd[rrlmask],'.')
ax = mp.grid[mp.axis_number(0,2)]
mpl_plot_templates.adaptive_param_plot(h77th112[rrlmask],
                                       c2cmtc6cm[rrlmask],
                                       bins=np.linspace(0,2.5,25),
                                       threshold=8,
                                       fill=False,
                                       alpha=0.8,
                                       axis=ax,
                                       ncontours=10,
                                       cmap=pl.mpl.cm.spectral)
ax.plot(np.linspace(0,3),1.0*np.linspace(0,3),'k--',alpha=0.5)
ax.plot(np.linspace(0,3),0.6*np.linspace(0,3),'k:',alpha=0.5)
ax.set_xlabel(r"H77$\alpha$/H112$\alpha$")
ax.set_ylabel(r"$S_{15 GHz} / S_{5 GHz}$")
ax.set_ylim(0,2.5)
ax.set_xlim(0,2.5)
ax.set_xticks(ax.get_xticks()[1:])
ax.set_yticks(ax.get_yticks()[1:])

mp.grid[mp.axis_number(1,2)].set_visible(False)
mp.grid[mp.axis_number(0,1)].set_visible(False)
mp.grid[mp.axis_number(2,1)].set_xticks(mp.grid[mp.axis_number(2,1)].get_xticks()[:-1])

pl.savefig('/Users/adam/work/h2co/maps/paper/figures/continuum_rrl_comparisongrid.pdf')
pl.show()
