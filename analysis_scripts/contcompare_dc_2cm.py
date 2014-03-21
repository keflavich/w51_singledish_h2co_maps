import pylab as pl
import numpy as np
import FITS_tools
from astropy.io import fits
from astropy import units as u
from agpy import smooth
import mpl_plot_templates
import aplpy
import image_registration

twotwofn = 'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits'
dcfn = 'langston_14ghz_gpa-bk1.fit' # has already been corrected
im1,im2 = FITS_tools.match_images.match_fits(twotwofn,dcfn,use_montage=False)
etamb_twotwo = 0.886

downsample_factor = 5

ac = im1 / etamb_twotwo #fits.getdata(twotwofn)
# unknown whether Langston images need etamb correction!  I couldn't tell from
# examining the data OR from reading the website
# T_MB per beam?  Probably.
dc = im2 #fits.getdata(dcfn).squeeze()
ac_hdr = fits.getheader(twotwofn)
dc_hdr = FITS_tools.strip_headers.flatten_header(fits.getheader(dcfn))

r_big = 9.5 * u.arcmin
r_dc = 6.6*u.arcmin
r_gb = 50*u.arcsec
r_sm = ((r_big**2-r_gb**2)**0.5).to(u.deg)
r_sm_dc = ((r_big**2-r_dc**2)**0.5).to(u.deg)
pixsize = abs(ac_hdr['CD1_1'])
fwhm = np.sqrt(8*np.log(2))
ac_kernsize = (r_sm.value/pixsize)/fwhm
dc_kernsize = (r_sm_dc.value/pixsize)/fwhm

sm_ac = smooth(ac, ac_kernsize, interpolate_nan=True, downsample=True,downsample_factor=downsample_factor)
sm_dc = dc # alternate option
# this is slow; should use astropy version
sm_dc = smooth(dc, dc_kernsize, downsample=True,downsample_factor=downsample_factor) # (8.9**2-6.6**2)**0.5/2.35 == 2.54

# correction not needed
# sm_ac_hdu = fits.PrimaryHDU(sm_ac, ac_hdr)
# 
# result = image_registration.FITS_tools.register_fits(sm_ac_hdu, dcfn, return_cropped_images=True, verbose=True, upsample_factor=10)
# sm_ac,dc = result[-2:]

okimg = np.isfinite(ac).astype('float')

#resampled_ac = FITS_tools.hcongrid.hcongrid(sm_ac, ac_hdr, dc_hdr)
#resampled_ok = FITS_tools.hcongrid.hcongrid(okimg, ac_hdr, dc_hdr)
resampled_ac = sm_ac
resampled_ok = okimg[::downsample_factor,::downsample_factor]

#xok = (np.isfinite(resampled_ac) * (resampled_ac != 0)).max(axis=0)
#yok = (np.isfinite(resampled_ac) * (resampled_ac != 0)).max(axis=1)

#croprange_x = np.argmax(xok),np.argmax(xok)+xok.sum()
#croprange_y = np.argmax(yok),np.argmax(yok)+yok.sum()
#slice_x = slice(*croprange_x)
#slice_y = slice(*croprange_y)

cropped_dc = sm_dc
cropped_ac = resampled_ac
cropped_rsok = np.round(resampled_ok).astype('bool')

crop_hdr = ac_hdr.copy()
crop_hdr['CD1_1'] *= downsample_factor
crop_hdr['CD2_2'] *= downsample_factor
crop_hdr['CRPIX1'] = (crop_hdr['CRPIX1']-1)/float(downsample_factor)+1
crop_hdr['CRPIX2'] = (crop_hdr['CRPIX2']-1)/float(downsample_factor)+1
#crop_hdr['CRPIX1'] -= croprange_x[0]
#crop_hdr['CRPIX2'] -= croprange_y[0]

ac_hdu = fits.PrimaryHDU(cropped_ac, crop_hdr)
dc_hdu = fits.PrimaryHDU(cropped_dc, crop_hdr)

gpa_gt = (cropped_dc < cropped_ac).astype('float')
gpa_gt_hdu = fits.PrimaryHDU(gpa_gt, crop_hdr)

fig = pl.figure(1)
sp1 = pl.subplot(2,2,1)
sp1b = sp1.bbox._bbox._points[0].tolist() + (sp1.bbox._bbox._points[1]-sp1.bbox._bbox._points[0]).tolist()
sp2 = pl.subplot(2,2,2)
sp2b = sp2.bbox._bbox._points[0].tolist() + (sp2.bbox._bbox._points[1]-sp2.bbox._bbox._points[0]).tolist()
pl.clf()
F1 = aplpy.FITSFigure(dc_hdu, convention='calabretta', subplot=sp1b, figure=fig)
F1._ax1.set_title("NRAO 300ft")
F2 = aplpy.FITSFigure(ac_hdu, convention='calabretta', subplot=sp2b, figure=fig)
F2._ax1.set_title("Green Bank")
for F in (F1,F2):
    F.show_grayscale(vmin=0.0, vmax=1.5, stretch='log', vmid=-0.1, invert=True)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.set_tick_xspacing(0.3)
    F.show_contour(fits.PrimaryHDU(okimg,ac_hdr), convention='calabretta', levels=[0.5], colors=['r'])
    F.show_contour(gpa_gt_hdu, convention='calabretta', levels=[0.5], colors=['g'])
    F.add_colorbar()
F1.remove_colorbar() # hack to preserve figure size
F2.axis_labels.hide_y()
F2.tick_labels.hide_y()
F2.colorbar.set_ticks([0.05,0.1,0.25,0.5,1])

pl.subplot(2,1,2)
ok = np.isfinite(cropped_ac) # dc is all "ok"
pl.plot(np.linspace(0,1.5),np.linspace(0,1.5)*1.0,'k--',linewidth=2,alpha=0.5,label=r'$y=x$')
pl.plot(np.linspace(0,1.5),np.linspace(0,1.5)*0.8,'k:' ,linewidth=2,alpha=0.5,label=r'$y=0.8 x$')
pl.plot(np.linspace(0,1.5),np.linspace(0,1.5)*0.6,'k-.',linewidth=2,alpha=0.5,label=r'$y=0.6 x$')
pl.plot(cropped_dc[ok],cropped_ac[ok],',')
#pl.plot(cropped_dc[gpa_gt.astype('bool')],cropped_ac[gpa_gt.astype('bool')],',',color='g',alpha=0.5)
pl.plot(cropped_dc[cropped_rsok],cropped_ac[cropped_rsok],'.',color='r')
#pl.plot(cropped_dc[cropped_rsok*gpa_gt.astype('bool')],cropped_ac[cropped_rsok*gpa_gt.astype('bool')],'.',color='g',alpha=0.5)
#mpl_plot_templates.adaptive_param_plot(cropped_dc[ok],cropped_ac[ok],bins=30,threshold=10,fill=True)
pl.xlabel(r'$T_B(K)$ NRAO 300ft')
pl.ylabel(r'$T_B(K)$ Green Bank')
pl.axis([0,1.5,0,1.5])
pl.legend(loc='upper left',fontsize=18)
pl.savefig('/Users/adam/work/h2co/maps/paper/figures/comparison_to_gpa.pdf')

pl.show()