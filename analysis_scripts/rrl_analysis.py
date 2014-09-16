import os
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
from common_constants import h2co11freq,h2co22freq,etamb_77,rrl,aobeamarea,gbbeamarea
from paths import datapath,fpath
import matplotlib

matplotlib.rc_file('pubfiguresrc')
pl.mpl.rcParams['axes.color_cycle'] = ["#"+x for x in "348ABD, 7A68A6, A60628, 467821, CF4457, 188487, E24A33".split(", ")]

# Have to convert to Jy to compare fairly
ktojy111 = (1*u.K).to(u.Jy,u.brightness_temperature(aobeamarea, h2co11freq))
ktojy77 = (1*u.K).to(u.Jy,u.brightness_temperature(gbbeamarea, h2co22freq))

aofn = os.path.join(datapath,'W51_Halpha_6cm_cube_supersampled.fits')
h112i = fits.getdata(os.path.join(datapath,'H110a_integral.fits')) * ktojy111 # ,aofn.replace("cube","integrated"))) * ktojy111
h112a = fits.getdata(os.path.join(datapath,'H110a_amplitude.fits')) * ktojy111
h112a.value[h112a == 0] = np.nan
h112c = fits.getdata(os.path.join(datapath,'W51_Halpha_6cm_cube_supersampled_continuum.fits')) * ktojy111
h112c.value[h112c == 0] = np.nan
h112h = fits.getheader(os.path.join(datapath,'W51_Halpha_6cm_cube_supersampled_continuum.fits'))

h77iname = os.path.join(datapath,'H77a_integral.fits') #'W51_h77a_pyproc_integrated_supersampled.fits')
h77aname = os.path.join(datapath,'H77a_amplitude.fits') #'W51_h77a_pyproc_integrated_supersampled.fits')
h77i = fits.getdata(h77iname) * ktojy77 / etamb_77
h77a = fits.getdata(h77aname) * ktojy77 / etamb_77
h77a.value[h77a == 0] = np.nan
#h77cname = os.path.join(datapath,'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits')
h77cname = os.path.join(datapath,'W51_h77a_pyproc_cube_supersampled_continuum.fits')
h77c = fits.getdata(h77cname) * ktojy77 / etamb_77
h77c.value[h77c == 0] = np.nan
h77h = fits.getheader(h77cname)

he77a_name = os.path.join(datapath,'W51_he77a_pyproc_integrated_supersampled.fits')
he77a = fits.getdata(he77a_name)
he112a_name = os.path.join(datapath,'W51_healpha_6cm_integrated_supersampled.fits')
he112a = fits.getdata(he112a_name)

amplitude_threshold111 = 0.2 * ktojy111
amplitude_threshold77 = 0.03 * ktojy77
h77lc = h77a/h77c * (h77a > amplitude_threshold77)
h112lc = h112a/h112c * (h112a > amplitude_threshold111)
rrlmask = (h77a > amplitude_threshold77) * (h112a > amplitude_threshold111)
h77th112 = h77a/h112a * rrlmask
c2cmtc6cm = (h77c/h112c * rrlmask)
c2cmtc6cm.value[c2cmtc6cm <= 0] = np.nan

# Electron Temperatures
# Compute using Wilson 2009, eqn 14.58
vwidth77 = fits.getdata(os.path.join(datapath,'H77a_velocity_width.fits')) * u.km/u.s
vwidth110 = fits.getdata(os.path.join(datapath,'H110a_velocity_width.fits')) * u.km/u.s
anut = 1.0 # Gaunt correction factor, eqn 10.35, pg 251
hefrac = 0.1 # Wilson 2009, p369
h77te = (6.985e3/anut * (rrl(77).to(u.GHz).value)**1.1 * (1./(1.+hefrac)) * (vwidth77*fwhm).to(u.km/u.s).value**-1 * h77lc**-1)**0.87
h112te = (6.985e3/anut * (rrl(112).to(u.GHz).value)**1.1 * (1./(1.+hefrac)) * (vwidth110*fwhm).to(u.km/u.s).value**-1 * h112lc**-1)**0.87
h77tename = os.path.join(datapath,'H77a_electrontemperature.fits')
h112tename = os.path.join(datapath,'H112a_electrontemperature.fits')
hdu = fits.PrimaryHDU(data=h77te.value, header=h77h)
hdu.writeto(h77tename, clobber=True)
hdu = fits.PrimaryHDU(data=h112te.value, header=h112h)
hdu.writeto(h112tename, clobber=True)

figs = []

vmaxes = {'77lc': 0.4,
          '112lc': 0.2,
          'te': 2.0e4,
          'th': 3.5,
          'tc': 2.5}

vmins = {'lc': 0,
         'te': 5000,
         'th': 0,
         'tc': -0.1}

titles = {'h77lc':r'H77$\alpha$ Line/Continuum',
          'h112lc':r'H112$\alpha$ Line/Continuum',
          'h77te':r'H77$\alpha$ Electron Temperature',
          'h112te':r'H112$\alpha$ Electron Temperature',
          'h77th112':r'H77$\alpha$ / H112$\alpha$',
          'c2cmtc6cm':r'$S_{15 GHz}/S_{5 GHz}$'}

cmap = pl.cm.jet

for ii,img in enumerate('h77te,h112te,h77lc,h112lc,h77th112,c2cmtc6cm'.split(',')):
    name = img
    img = locals()[name].value
    hdu = fits.PrimaryHDU(data=img, header=h112h)
    pl.figure(ii)
    pl.clf()
    F = aplpy.FITSFigure(hdu, figure=pl.figure(ii))
    vmax = [vmaxes[k] for k in vmaxes if k in name][0]
    vmin = [vmins[k] for k in vmins if k in name][0]
    F.show_colorscale(vmin=vmin,vmax=vmax,cmap=cmap)
    F.recenter(49.235,-0.303,width=0.952,height=0.433)
    F.add_colorbar()
    F._ax1.set_title(titles[name])
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    #F.set_tick_xspacing(0.3)
    figs.append(F)
    F.save(fpath('ratiomap_cont_%s.pdf' % name))
    F.save(fpath('ratiomap_cont_%s.png' % name))
    F.show_contour(he77a_name, convention='calabretta', levels=[0.0125, 0.025, 0.05,0.1,0.15,0.2],
                   colors='w', linewidth=0.5)
    F.save(fpath('ratiomap_cont_%s_HeContours.pdf' % name))

pl.figure(ii+1)
pl.clf()
F = aplpy.FITSFigure(h77cname, figure=pl.figure(ii+1), convention='calabretta')
F.show_grayscale(vmin=-0.1,vmid=-0.4,vmax=550,invert=True,stretch='log')
F.recenter(49.235,-0.303,width=0.952,height=0.433)
F.show_contour(hdu, levels=[0.3,0.5,0.7], colors=['b',(0.2,1,0.4,0.8),'r'])
F.add_colorbar()

pl.figure(ii+2)
pl.clf()
F = aplpy.FITSFigure(h77iname, figure=pl.figure(ii+2), convention='calabretta')
F.show_grayscale(vmin=-0.1,vmid=-0.4,vmax=180.,invert=True,stretch='log')
F.recenter(49.235,-0.303,width=0.952,height=0.433)
F.show_contour(hdu, levels=[0.3,0.5,0.7], colors=['b',(0.2,1,0.4,0.8),'r'])
F.add_colorbar()

pl.figure(ii+3)
pl.clf()
ax = pl.gca()
pl.plot([000,4e4],[000,4e4],'k--')
h77te_med = np.median(h77te.value[np.isfinite(h77te.value)])
h112te_med = np.median(h112te.value[np.isfinite(h112te.value)])
h77te_mean = np.mean(h77te.value[np.isfinite(h77te.value)])
h112te_mean = np.mean(h112te.value[np.isfinite(h112te.value)])
mpl_plot_templates.adaptive_param_plot(h77te.value, h112te.value, bins=50,
                                       ncontours=10, threshold=5, fill=True,
                                       alpha=0.5, cmap=pl.mpl.cm.spectral,
                                       axis=ax)
ax.plot(h77te_med, h112te_med, 'kx', markersize=15, markeredgewidth=3)
ax.plot(h77te_mean, h112te_mean, 'k+', markersize=15, markeredgewidth=3)
ax.set_xlabel(r"$T_e(\mathrm{H}77\alpha)$")
ax.set_ylabel(r"$T_e(\mathrm{H}112\alpha)$")
ax.axis([0,2e4,0,2e4])
pl.savefig(fpath('electron_temperature_77vs111.pdf'))


keys = [
       r'$S_{5 GHz}$',
       r'$S_{15 GHz}$',
       r'H112$\alpha$',
       r'H77$\alpha$',
        ]

data = {r'H77$\alpha$':  h77a,
        r'$S_{15 GHz}$': h77c,
        r'H112$\alpha$': h112a,
        r'$S_{5 GHz}$':  h112c}

combs = list(itertools.combinations(keys,2))

xpos = {r'$S_{5 GHz}$':0,
        r'$S_{15 GHz}$':1,
        r'H112$\alpha$':2,
        }

ypos = {r'$S_{15 GHz}$':0,
        r'H112$\alpha$':1,
        r'H77$\alpha$':2,
        }

pl.figure(ii+4)
pl.clf()
mp = multiplot.multipanel(dims=(3,3),diagonal=False,figID=ii+4)
for ii,(xd,yd) in zip(mp.grid.keys(),combs):
    if not(yd in ypos and xd in xpos):
        continue
    print mp.axis_number(ypos[yd],xpos[xd],), xpos[xd], ypos[yd]
    ax = mp.grid[mp.axis_number(ypos[yd],xpos[xd])]
    #ax.plot(data[xd][rrlmask],data[yd][rrlmask],'.')
    mpl_plot_templates.adaptive_param_plot(data[xd][rrlmask].to(u.Jy).value,
                                           data[yd][rrlmask].to(u.Jy).value,
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
mpl_plot_templates.adaptive_param_plot(h77th112[rrlmask].value,
                                       c2cmtc6cm[rrlmask].value,
                                       bins=np.linspace(0,vmaxes['th'],25),
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
ax.set_ylim(0,vmaxes['tc'])
ax.set_xlim(0,vmaxes['th'])
ax.set_xticks(ax.get_xticks()[1:])
ax.set_yticks(ax.get_yticks()[1:])

mp.grid[mp.axis_number(1,2)].set_visible(False)
mp.grid[mp.axis_number(0,1)].set_visible(False)
mp.grid[mp.axis_number(2,1)].set_xticks(mp.grid[mp.axis_number(2,1)].get_xticks()[:-1])

pl.savefig(fpath('continuum_rrl_comparisongrid.pdf'))
pl.show()
