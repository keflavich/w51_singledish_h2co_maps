import astropy.io.fits as pyfits
import astropy.io.fits as fits
import numpy as np
import os

from gbtpy import makecube

def fix_TDIM_in_header(fn):
    """
    Apparently TDIM that's stored in headers started being checked against the
    data in a recent version of astropy.io.fits, and that leads to files being
    unreadable if they have TDIM parameters.  I apparently discovered this
    during the September 2013 observing run on NGC 1333 and made the
    appropriate changes in the Arecibo reduction files but never re-reduced the
    W51 data
    """
    f = pyfits.open(fn)

    fixed=False
    for HDU in f:
        for kw in HDU.header.keys():
            if any([x in kw for x in ('TDIM',)]):
                fixed = True
                del HDU.header[kw]

    if fixed:
        outf = fn.replace(".fits","_fixed.fits")
        f.writeto(outf,output_verify='fix',clobber='fixed' in outf)
        print "Fixed file ",fn," -> ",outf
    else:
        outf = fn
        print "Did not fix file ",fn

    return outf

prefix = '/Users/adam/observations/arecibo/20120910/'

linefreq = 4.8296594e9

vmin = 50
vmax = 80
velocityrange = [-50,150] # match GBT exactly

cubename_supersampled = '/Users/adam/work/h2co/maps/w51/W51_H213CO11_cube_supersampled'
velocityrange = [-50,150] # match GBT exactly
cd3 = 1.0 # match GBT exactly
naxis3 = int((velocityrange[1]-velocityrange[0]) / cd3) + 1
crval3 = 50.0
vels = crval3+cd3*(np.arange(naxis3)+1-naxis3/2-1)
makecube.generate_header(49.209553,-0.277137,naxis1=308,naxis2=205,pixsize=15,naxis3=int(naxis3),cd3=cd3,crval3=crval3,clobber=True,
                         restfreq=4.8296594e9)
makecube.make_blank_images(cubename_supersampled,clobber=True)

files = ['W51_h213coW_spectra_0910.fits',
         'W51_h213coW_spectra_0911.fits',
         'W51_h213coW_spectra_0912.fits',
         'W51_h213coW_spectra_0915.fits',
         ]

for fn in files:
    if not os.path.exists(prefix+fn):
        print "File not found: ",prefix+fn
        continue
    try:
        fits.getdata(prefix+fn)
    except Exception as ex:
        print "Failed to load ",prefix+fn
        print ex
        continue

    fullfn = fix_TDIM_in_header(prefix+fn)
    makecube.add_file_to_cube(fullfn,
                              cubename_supersampled+'.fits',
                              add_with_kernel=True,
                              kernel_fwhm=20./3600.,
                              nhits=cubename_supersampled+'_nhits.fits',wcstype='D',
                              diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
                              velocityrange=velocityrange,excludefitrange=[vmin,vmax],linefreq=linefreq)

flat_vrange = [45,75]

makecube.runscript(cubename_supersampled)
makecube.make_flats(cubename_supersampled,vrange=flat_vrange,noisevrange=[-50,-1])
makecube.make_taucube(cubename_supersampled, continuum=cubename_supersampled+"_continuum.fits") # etamb accounted for already , etamb=0.51)
makecube.make_flats(cubename_supersampled.replace("cube","taucube"),vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')




