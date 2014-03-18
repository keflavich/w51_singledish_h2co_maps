import pylab as pl
import numpy as np
import FITS_tools
from astropy.io import fits
from astropy import units as u
from agpy import smooth
import mpl_plot_templates
import aplpy

oneonefn = 'W51_H2CO11_cube_supersampled_continuum.fits'
dcfn = 'w51.iuq.fits'

ac = fits.getdata(oneonefn)
dc = fits.getdata(dcfn,ext=1) / 1e3 # mK -> K
ac_hdr = fits.getheader(oneonefn)
dc_hdr = FITS_tools.strip_headers.flatten_header(fits.getheader(dcfn,ext=1))

r_dc = 9.5 * u.arcmin
r_ao = 50*u.arcsec
r_sm = ((r_dc**2-r_ao**2)**0.5).to(u.deg)
pixsize = abs(ac_hdr['CD1_1'])
fwhm = np.sqrt(8*np.log(2))

sm_ac = smooth(ac, r_sm.value / pixsize / fwhm, interpolate_nan=True)

okimg = np.isfinite(ac).astype('float')

resampled_ac = FITS_tools.hcongrid.hcongrid(sm_ac, ac_hdr, dc_hdr)
resampled_ok = FITS_tools.hcongrid.hcongrid(okimg, ac_hdr, dc_hdr)

xok = (np.isfinite(resampled_ac) * (resampled_ac != 0)).max(axis=0)
yok = (np.isfinite(resampled_ac) * (resampled_ac != 0)).max(axis=1)

croprange_x = np.argmax(xok),np.argmax(xok)+xok.sum()
croprange_y = np.argmax(yok),np.argmax(yok)+yok.sum()
slice_x = slice(*croprange_x)
slice_y = slice(*croprange_y)

cropped_dc = dc[slice_y,slice_x]
cropped_ac = resampled_ac[slice_y,slice_x]
cropped_rsok = resampled_ok[slice_y,slice_x]

crop_hdr = dc_hdr.copy()
crop_hdr['CRPIX1'] -= croprange_x[0]
crop_hdr['CRPIX2'] -= croprange_y[0]

ac_hdu = fits.PrimaryHDU(cropped_ac, crop_hdr)
dc_hdu = fits.PrimaryHDU(cropped_dc, crop_hdr)

fig = pl.figure(1)
sp1 = pl.subplot(2,2,1)
sp1b = sp1.bbox._bbox._points[0].tolist() + (sp1.bbox._bbox._points[1]-sp1.bbox._bbox._points[0]).tolist()
sp2 = pl.subplot(2,2,2)
sp2b = sp2.bbox._bbox._points[0].tolist() + (sp2.bbox._bbox._points[1]-sp2.bbox._bbox._points[0]).tolist()
pl.clf()
F1 = aplpy.FITSFigure(dc_hdu, convention='calabretta', subplot=sp1b, figure=fig)
F1._ax1.set_title("Urumqi 25m")
F2 = aplpy.FITSFigure(ac_hdu, convention='calabretta', subplot=sp2b, figure=fig)
F2._ax1.set_title("Arecibo")
for F in (F1,F2):
    F.show_grayscale(vmin=0.0, vmax=20, stretch='log', vmid=-0.5, invert=True)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.set_tick_xspacing(0.3)
    F.show_contour(fits.PrimaryHDU(okimg,ac_hdr), convention='calabretta', levels=[0.5], colors=['r'])
    F.add_colorbar()
F1.remove_colorbar() # hack to preserve figure size
F2.axis_labels.hide_y()
F2.tick_labels.hide_y()
F2.colorbar.set_ticks([0.1,0.5,1,2.5,5,10,20])


pl.subplot(2,1,2)
ok = np.isfinite(cropped_ac) # dc is all "ok"
pl.plot(np.linspace(0,20),np.linspace(0,20),'k--',linewidth=2,alpha=0.5)
pl.plot(cropped_dc[ok],cropped_ac[ok],',')
pl.plot(cropped_dc[np.round(cropped_rsok).astype('bool')],cropped_ac[np.round(cropped_rsok).astype('bool')],'.',color='r')
#mpl_plot_templates.adaptive_param_plot(cropped_dc[ok],cropped_ac[ok],bins=30,threshold=10,fill=True)
pl.xlabel(r'$T_B(K)$ Urumqi')
pl.ylabel(r'$T_B(K)$ Arecibo')
pl.axis([0,20,0,20])
pl.savefig('/Users/adam/work/h2co/maps/paper/figures/comparison_to_urumqi_6cm.pdf')

pl.show()
