import pylab as pl
import numpy as np
import FITS_tools
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy import convolution
import mpl_plot_templates
import aplpy
import paths

pl.matplotlib.rc_file(paths.rcfilepath)
pl.rcParams['font.size'] = 18

oneonefn = paths.cont6cm
dcfn = paths.dpath('w51.iuq.fits')

ac = fits.getdata(oneonefn)
# July 29: I didn't do anything, I SWEAR!
# This used to be ext=0, and now it's not, and I DON'T KNOW WHY
# (ext=1 yields an ERROR now?!)
# Re-downloading fresh fixes it - I think using WCSUTILS corrupts the file
dc = fits.getdata(dcfn,ext=1) / 1e3 # mK -> K
ac_hdr = fits.getheader(oneonefn)
ac_wcs = wcs.WCS(ac_hdr)
dc_hdr = FITS_tools.strip_headers.flatten_header(fits.getheader(dcfn,ext=1))

r_dc = 9.5 * u.arcmin
r_ao = 50*u.arcsec
r_sm = ((r_dc**2-r_ao**2)**0.5).to(u.deg)
pixsize = np.abs(np.prod((ac_wcs.wcs.get_cdelt() * ac_wcs.wcs.get_pc().diagonal())))**0.5
fwhm = np.sqrt(8*np.log(2))

kernel = convolution.Gaussian2DKernel(r_sm.value / pixsize / fwhm,
                                      x_size=ac.shape[1], y_size=ac.shape[0])
kernel.normalize()
sm_ac = convolution.convolve_fft(ac, kernel, interpolate_nan=True,
                                 ignore_edge_zeros=True)
# from AG_fft_tools import smooth
# sm_ac = smooth(ac, r_sm.value / pixsize / fwhm, interpolate_nan=True)

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

pl.close(1)
fig = pl.figure(1, figsize=(12,9))
sp1 = pl.subplot(2,2,1)
sp1b = sp1.bbox._bbox._points[0].tolist() + (sp1.bbox._bbox._points[1]-sp1.bbox._bbox._points[0]).tolist()
sp2 = pl.subplot(2,2,2)
sp2b = sp2.bbox._bbox._points[0].tolist() + (sp2.bbox._bbox._points[1]-sp2.bbox._bbox._points[0]).tolist()
pl.clf()
F1 = aplpy.FITSFigure(dc_hdu, convention='calabretta', subplot=sp1b, figure=fig)
F1.set_auto_refresh(False)
F1._ax1.set_title("Urumqi 25m")
F2 = aplpy.FITSFigure(ac_hdu, convention='calabretta', subplot=sp2b, figure=fig)
F2.set_auto_refresh(False)
F2._ax1.set_title("Arecibo")
for F in (F1,F2):
    F.show_grayscale(vmin=0.0, vmax=25, stretch='log', vmid=-0.5, invert=True)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.set_tick_xspacing(0.3)
    F.show_contour(fits.PrimaryHDU(okimg,ac_hdr), convention='calabretta', levels=[0.5], colors=['r'])
    F.add_colorbar()
F1.remove_colorbar() # hack to preserve figure size
F2.axis_labels.hide_y()
F2.tick_labels.hide_y()
F2.colorbar.set_ticks([0.1,0.5,1,2.5,5,10,20])
F2.colorbar._colorbar.set_label("$T_{MB}$ (K)",rotation=270, labelpad=20)


pl.subplot(2,1,2)
ok = np.isfinite(cropped_ac) # dc is all "ok"
pl.plot(np.linspace(0,25),np.linspace(0,25),'k--',linewidth=2,alpha=0.5, label=r'$y=x$')
pl.plot(np.linspace(0,25),np.linspace(0,25)*1.2,'k:' ,
        linewidth=2,alpha=0.5,label=r'$y=1.2 x$')
pl.plot(np.linspace(0,25),np.linspace(0,25)*1.4,'k-.' ,
        linewidth=2,alpha=0.5,label=r'$y=1.4 x$')
#pl.plot(cropped_dc[ok],cropped_ac[ok],',')
pl.plot(cropped_dc[np.round(cropped_rsok).astype('bool')],cropped_ac[np.round(cropped_rsok).astype('bool')],'.',color='r')
#mpl_plot_templates.adaptive_param_plot(cropped_dc[ok],cropped_ac[ok],bins=30,threshold=10,fill=True)
pl.xlabel(r'$T_{MB}(K)$ Urumqi', labelpad=15)
pl.ylabel(r'$T_{MB}(K)$ Arecibo')
pl.axis([0,20,0,25])
pl.legend(loc='upper left',fontsize=18)

F1.refresh()
F2.refresh()
pl.savefig(paths.fpath('comparison_to_urumqi_6cm.pdf'))

# This is just for checking purposes; this figure is not saved
diffmap = fits.PrimaryHDU(cropped_dc-cropped_ac, crop_hdr)
fig2 = pl.figure(2)
fig2.clf()
F3 = aplpy.FITSFigure(diffmap, convention='calabretta', figure=fig2)
F3.show_grayscale(invert=True)
F3.show_contour(fits.PrimaryHDU(cropped_ac, crop_hdr), convention='calabretta')
F3.add_colorbar()

pl.show()
