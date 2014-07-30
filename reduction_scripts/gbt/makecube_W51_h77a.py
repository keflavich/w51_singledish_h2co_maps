import pyfits,pywcs,coords
import numpy as np
import matplotlib
import pylab as pl
#import aplpy
import os
from agpy import asinh_norm
#from agpy import readcol,asinh_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import log

import sys
#sys.path.append('/Users/adam/repos/casaradio/branches/python/ginsburg/')
from sdpy import makecube

filelist = [
# session 10 is really messed up
#"Session10_32to93_B17_F1.fits",
#"Session10_32to93_B21_F1.fits", # bad in H2CO?
#"Session10_32to93_D33_F2.fits",
#"Session10_32to93_D37_F2.fits", # bad in H2CO?
"Session11_126to157_B17_F1.fits",
"Session11_126to157_B21_F1.fits",
"Session11_126to157_D33_F2.fits",
"Session11_126to157_D37_F2.fits",
"Session11_158to189_B17_F1.fits",
"Session11_158to189_B21_F1.fits",
"Session11_158to189_D33_F2.fits",
"Session11_158to189_D37_F2.fits",
"Session11_21to52_B17_F1.fits",
"Session11_21to52_B21_F1.fits",
"Session11_21to52_D33_F2.fits",
"Session11_21to52_D37_F2.fits",
"Session11_53to78_B17_F1.fits",
"Session11_53to78_B21_F1.fits",
"Session11_53to78_D33_F2.fits",
#"Session11_53to78_D37_F2.fits",  # really bad "jump" artifacts
"Session11_89to120_B17_F1.fits",
"Session11_89to120_B21_F1.fits",
"Session11_89to120_D33_F2.fits",
"Session11_89to120_D37_F2.fits",
"Session14_12to43_B17_F1.fits",
"Session14_12to43_B21_F1.fits",
#"Session14_12to43_D33_F2.fits", # jump artifacts
"Session14_12to43_D37_F2.fits",
"Session14_130to161_B17_F1.fits",
"Session14_130to161_B21_F1.fits",
"Session14_130to161_D33_F2.fits",
"Session14_130to161_D37_F2.fits",
"Session14_167to198_B17_F1.fits",
"Session14_167to198_B21_F1.fits",
#"Session14_167to198_D33_F2.fits",  # Moderate, but still unacceptable, "jump" artifacts
"Session14_167to198_D37_F2.fits",
"Session14_44to75_B17_F1.fits",
"Session14_44to75_B21_F1.fits",
"Session14_44to75_D33_F2.fits",
#"Session14_44to75_D37_F2.fits", # jump artifacts
"Session14_81to96_B17_F1.fits",
"Session14_81to96_B21_F1.fits",
"Session14_81to96_D33_F2.fits",
"Session14_81to96_D37_F2.fits",
"Session14_98to129_B17_F1.fits",
"Session14_98to129_B21_F1.fits",
"Session14_98to129_D33_F2.fits",
"Session14_98to129_D37_F2.fits",
"Session16_10to41_B17_F1.fits",
"Session16_10to41_B21_F1.fits",
"Session16_10to41_D33_F2.fits",
"Session16_10to41_D37_F2.fits",
"Session16_42to73_B17_F1.fits",
"Session16_42to73_B21_F1.fits",
"Session16_42to73_D33_F2.fits",
"Session16_42to73_D37_F2.fits",
"Session16_74to105_B17_F1.fits",
"Session16_74to105_B21_F1.fits",
"Session16_74to105_D33_F2.fits",
"Session16_74to105_D37_F2.fits",
"Session17_10to41_B17_F1.fits",
"Session17_10to41_B21_F1.fits",
"Session17_10to41_D33_F2.fits",
"Session17_10to41_D37_F2.fits",
"Session17_42to73_B17_F1.fits",
"Session17_42to73_B21_F1.fits",
"Session17_42to73_D33_F2.fits",
"Session17_42to73_D37_F2.fits",
"Session17_111to186_B17_F1.fits",
"Session17_111to186_B21_F1.fits",
"Session17_111to186_D33_F2.fits",
"Session17_111to186_D37_F2.fits",
]

filelist2=[
    "Session22_10to41_B17_F1.fits",
    "Session22_10to41_B21_F1.fits",
    "Session22_10to41_D33_F2.fits",
    "Session22_10to41_D37_F2.fits",
    "Session22_183to203_B17_F1.fits",
    "Session22_183to203_B21_F1.fits",
    "Session22_183to203_D33_F2.fits",
    "Session22_183to203_D37_F2.fits",
    "Session22_42to73_B17_F1.fits",
    "Session22_42to73_B21_F1.fits",
    "Session22_42to73_D33_F2.fits",
    "Session22_42to73_D37_F2.fits",
    "Session22_74to105_B17_F1.fits",
    "Session22_74to105_B21_F1.fits",
    "Session22_74to105_D33_F2.fits",
    "Session22_74to105_D37_F2.fits",
    ]

cubename_supersampled = '/Users/adam/work/h2co/maps/w51/W51_h77a_pyproc_cube_supersampled'
velocityrange=[-50,150]
cd3 = 1.0
#cd3 = 1.0 # Arecibo is limited to 0.64 because one of the receivers went bad at hi-res mode once
#naxis3 = int(np.ceil((velocityrange[1]-velocityrange[0])/cd3))+4 # +4 is BAD!  don't do that.
naxis3 = int((velocityrange[1]-velocityrange[0]) / cd3) + 1 # +1 is good: include -50
crval3 = 50.0
# dumb debug stuff
vels = crval3+cd3*(np.arange(naxis3)+1-naxis3/2-1)
# this will probably cause an error but I must insist..
#if velocityrange[0]<=vels.min() or velocityrange[1]>=vels.max():
#    raise ValueError("Add more points.  Something's going to be out of range for starlink")
makecube.generate_header(49.209553, -0.277137, naxis1=308, naxis2=205,
                         pixsize=15, naxis3=naxis3, cd3=cd3, crval3=crval3,
                         clobber=True, restfreq=14.12861587343394e9)
makecube.make_blank_images(cubename_supersampled,clobber=True)


for fn in filelist+filelist2:
    log.info("Adding file %s" % fn)
    fullfn = '/Users/adam/observations/gbt/W51map/'+fn
    d = pyfits.getdata(fullfn)
    pl.clf()
    pl.imshow(d['DATA'],norm=asinh_norm.AsinhNorm())
    pl.savefig(fullfn.replace(".fits","_data.png"))
    pl.clf()
    dsub = d['DATA']-np.median(d['DATA'],axis=1)[:,None]
    pl.imshow(dsub)
    pl.savefig(fullfn.replace(".fits","_data_subbed.png"))
    makecube.add_file_to_cube(fullfn,
                              cubename_supersampled+".fits",
                              add_with_kernel=True,
                              kernel_fwhm=20./3600.,
                              nhits=cubename_supersampled+"_nhits.fits", 
                              wcstype='V',
                              diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
                              velocityrange=velocityrange,excludefitrange=[20,90],noisecut=50,
                              continuum_prefix=cubename_supersampled+fn.replace(".fits",''),
                              negative_mean_cut=-1)
                              # noisecut was probably a bad idea: it cut SIGNAL  
                              # more aggressive noisecut

makecube.runscript(cubename_supersampled)
makecube.make_flats(cubename_supersampled,vrange=[20,90],noisevrange=[-15,20])
