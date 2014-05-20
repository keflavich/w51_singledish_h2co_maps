import glob
from astropy import units as u
from astropy.io import fits
from higal_sedfitter import smooth,fit,higal_beams
from higal_sedfitter.fit import PixelFitter

pixelfitter = PixelFitter(bfixed=True)

for fn in glob.glob("HIGAL048*fits"):
    smooth.add_beam_information_to_higal_header(fn, name_to_um=higal_beams.num_to_um)

target_header = fits.getheader('HIGAL0482p012_500_RM.fits')
smooth.smooth_images_toresolution(45*u.arcsec, skip_existing=False,
                                  globs=['HIGAL048*'], regrid=True,
                                  target_header=target_header, clobber=True)
fit.fit_modified_blackbody_tofiles('HIGAL048*_{0}_RM_smregrid45.fits',
                                   pixelfitter=pixelfitter,
                                   name_to_um=higal_beams.num_to_um,
                                   out_prefix='HIGAL_L048_', integral=True)

    
for fn in glob.glob("HIGAL050*fits"):
    smooth.add_beam_information_to_higal_header(fn,
                                                name_to_um=higal_beams.num_to_um)

target_header = fits.getheader('HIGAL0504p012_500_RM.fits')
smooth.smooth_images_toresolution(45*u.arcsec, skip_existing=False,
                                  globs=['HIGAL050*'], regrid=True,
                                  target_header=target_header,
                                  clobber=True)
fit.fit_modified_blackbody_tofiles('HIGAL050*_{0}_RM_smregrid45.fits',
                                   pixelfitter=pixelfitter,
                                   name_to_um=higal_beams.num_to_um,
                                   out_prefix='HIGAL_L050_', integral=True)

import montage
header = fits.getheader('HIGAL0504p012_500_RM.fits')
header['CRPIX1'] = 1100
header['CRVAL1'] = 50.0
header['CRPIX2'] = 550
header['CRVAL2'] = 0.0
header['NAXIS1'] = 2200
header['NAXIS2'] = 1100
header.totextfile('header.hdr',clobber=True)
montage.wrapper(['HIGAL*T.fits'], outfile='HIGAL_W51_TemperatureFit_Beta1.75.fits', header='header.hdr', copy=True, combine='mean')
montage.wrapper(['HIGAL*N.fits'], outfile='HIGAL_W51_ColumnFit_Beta1.75.fits', header='header.hdr', copy=True, combine='mean')
montage.wrapper(['HIGAL*integral.fits'], outfile='HIGAL_W51_LuminosityFit_Beta1.75.fits', header='header.hdr', copy=True, combine='mean')
