import numpy as np
import os

from fix_tdim import fix_TDIM_in_header
from sdpy import makecube
from makecube_pars import (glon, glat, naxis1, naxis2, naxis3, vels, crval3,
                           cd3, vmin, vmax, velocityrange)

prefix = '/Users/adam/observations/arecibo/20120910/'

linefreq = 4.59309e9

cubename_supersampled = '/Users/adam/work/h2co/maps/w51/W51_H213CO11_cube_supersampled'
makecube.generate_header(glon, glat, naxis1=naxis1, naxis2=naxis2,
                         pixsize=15, naxis3=int(naxis3), cd3=cd3,
                         crval3=crval3, clobber=True, 
                         restfreq=linefreq)
makecube.make_blank_images(cubename_supersampled,clobber=True)

files = ['W51_h213coW_spectra_0910.fits',
         'W51_h213coW_spectra_0911.fits',
         'W51_h213coW_spectra_0912.fits',
         'W51_h213coW_spectra_0915.fits',
         ]

for fn in files:
    fullfn = fix_TDIM_in_header(prefix+fn)
    makecube.add_file_to_cube(fullfn,
                              cubename_supersampled+'.fits',
                              add_with_kernel=True,
                              kernel_fwhm=20./3600.,
                              nhits=cubename_supersampled+'_nhits.fits',
                              diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
                              velocityrange=velocityrange,excludefitrange=[vmin,vmax],linefreq=linefreq)

flat_vrange = [45,75]

makecube.runscript(cubename_supersampled)
makecube.make_flats(cubename_supersampled,vrange=flat_vrange,noisevrange=[-50,-1])
makecube.make_taucube(cubename_supersampled,
                      continuum=cubename_supersampled+"_continuum.fits",
                      linefreq=linefreq) # etamb accounted for already , etamb=0.51)
makecube.make_flats(cubename_supersampled.replace("cube","taucube"),
                    vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')




