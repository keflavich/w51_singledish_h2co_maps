import pyfits,pywcs,coords
import numpy as np
import matplotlib
import pylab as pl
#import aplpy
import os
from agpy import asinh_norm
#from agpy import readcol,asinh_norm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
#sys.path.append('/Users/adam/repos/casaradio/branches/python/ginsburg/')
from gbtpy import makecube
cubename = '/Users/adam/work/h2co/maps/w51/W51_h77a_pyproc_cube'
velocityrange=[-50,150]
cd3 = 1
naxis3 = int(np.ceil((velocityrange[1]-velocityrange[0])/cd3))+4
crval3 = (velocityrange[1]+velocityrange[0])/2.
# dumb debug stuff
vels = crval3+cd3*(np.arange(naxis3)+1-naxis3/2)
if velocityrange[0]<=vels.min() or velocityrange[1]>=vels.max():
    raise ValueError("Add more points.  Something's going to be out of range for starlink")
makecube.generate_header(49.209553,-0.277137,naxis1=192,naxis2=128,pixsize=24,naxis3=naxis3,cd3=cd3,crval3=crval3,clobber=True)
makecube.make_blank_images(cubename,clobber=True)

# probably these samplers:
# 2: ["B17","B21","D33","D37"],

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
"Session11_53to78_D37_F2.fits",
"Session11_89to120_B17_F1.fits",
"Session11_89to120_B21_F1.fits",
"Session11_89to120_D33_F2.fits",
"Session11_89to120_D37_F2.fits",
"Session14_12to43_B17_F1.fits",
"Session14_12to43_B21_F1.fits",
"Session14_12to43_D33_F2.fits",
"Session14_12to43_D37_F2.fits",
"Session14_130to161_B17_F1.fits",
"Session14_130to161_B21_F1.fits",
"Session14_130to161_D33_F2.fits",
"Session14_130to161_D37_F2.fits",
"Session14_167to198_B17_F1.fits",
"Session14_167to198_B21_F1.fits",
"Session14_167to198_D33_F2.fits",
"Session14_167to198_D37_F2.fits",
"Session14_44to75_B17_F1.fits",
"Session14_44to75_B21_F1.fits",
"Session14_44to75_D33_F2.fits",
"Session14_44to75_D37_F2.fits",
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
        
pl.figure(99)
for fn in filelist:
    print "Adding file %s" % fn
    fullfn = '/Users/adam/observations/gbt/W51map/'+fn
    d = pyfits.getdata(fullfn)
    pl.clf()
    pl.imshow(d['DATA'],norm=asinh_norm.AsinhNorm())
    pl.savefig(fullfn.replace(".fits","_data.png"),bbox_inches='tight')


    dsub = d['DATA']-np.median(d['DATA'],axis=1)[:,None]

    pl.clf()
    ax = pl.gca()
    dstd = dsub.std()
    dmean = dsub.mean()
    im = ax.imshow(dsub,vmin=dmean-5*dstd,vmax=dmean+5*dstd,norm=asinh_norm.AsinhNorm())
    divider = make_axes_locatable(ax)

    right = divider.append_axes("right", size="15%", pad=0.05)
    vright = divider.append_axes("right", size="15%", pad=0.05)
    meany = dsub.mean(axis=1)
    erry = dsub.std(axis=1)
    right.plot(meany,np.arange(meany.size))
    right.set_ylim(0,meany.size-1)
    right.set_yticks([])
    right.set_xticks([meany.min(),np.median(meany),meany.max()])
    pl.setp( right.xaxis.get_majorticklabels(), rotation=70 )
    right.set_title("$\mu$")
    vright.plot(erry,np.arange(erry.size))
    vright.set_ylim(0,erry.size-1)
    vright.set_yticks([])
    vright.set_xticks([erry.min(),np.median(erry),erry.max()])
    vright.set_xlabel("$\sigma$")
    vright.xaxis.set_ticks_position('top')
    pl.setp( vright.xaxis.get_majorticklabels(), rotation=70 )
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)

    top = divider.append_axes("top", size="15%", pad=0.05)
    vtop = divider.append_axes("top", size="15%", pad=0.05)
    meanx = dsub.mean(axis=0)
    errx = dsub.std(axis=0)
    top.plot(np.arange(meanx.size),meanx)
    top.set_xlim(0,meanx.size-1)
    top.set_xticks([])
    top.set_yticks([meanx.min(),np.median(meanx),meanx.max()])
    pl.setp(top.yaxis.get_majorticklabels(), rotation=20 )
    top.set_title("$\mu$")
    vtop.plot(np.arange(errx.size),errx,)
    vtop.set_xlim(0,errx.size-1)
    vtop.set_xticks([])
    vtop.set_yticks([errx.min(),np.median(errx),errx.max()])
    vtop.set_ylabel("$\sigma$")
    pl.setp(vtop.yaxis.get_majorticklabels(), rotation=-20 )

    pl.savefig(fullfn.replace(".fits","_data_subbed.png"),bbox_inches='tight')
    makecube.add_file_to_cube(fullfn,
        cubename+".fits",
        nhits=cubename+"_nhits.fits", wcstype='V',
        diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
        negative_mean_cut=-1,
        velocityrange=velocityrange,excludefitrange=[-90,-15,30,90],noisecut=1.0)

makecube.runscript(cubename)
makecube.make_flats(cubename,vrange=[30,90],noisevrange=[-15,30])
# silly makecube.make_taucube(cubename,cubename+"_continuum.fits",etamb=0.886)

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

import shutil
cubename2=cubename+"_sess22"
shutil.copy(cubename+".fits",cubename2+".fits")
shutil.copy(cubename+"_nhits.fits",cubename2+"_nhits.fits")

for fn in filelist2:
    print "Adding file %s" % fn
    fullfn = '/Volumes/disk4/gbt/w51_mosaic/'+fn
    d = pyfits.getdata(fullfn)
    pl.clf()
    pl.imshow(d['DATA'],norm=asinh_norm.AsinhNorm())
    pl.savefig(fullfn.replace(".fits","_data.png"),bbox_inches='tight')

    dsub = d['DATA']-np.median(d['DATA'],axis=1)[:,None]

    pl.clf()
    ax = pl.gca()
    dstd = dsub.std()
    dmean = dsub.mean()
    im = ax.imshow(dsub,vmin=dmean-5*dstd,vmax=dmean+5*dstd,norm=asinh_norm.AsinhNorm())
    divider = make_axes_locatable(ax)

    right = divider.append_axes("right", size="15%", pad=0.05)
    vright = divider.append_axes("right", size="15%", pad=0.05)
    meany = dsub.mean(axis=1)
    erry = dsub.std(axis=1)
    right.plot(meany,np.arange(meany.size))
    right.set_ylim(0,meany.size-1)
    right.set_yticks([])
    right.set_xticks([meany.min(),np.median(meany),meany.max()])
    pl.setp( right.xaxis.get_majorticklabels(), rotation=70 )
    right.set_title("$\mu$")
    vright.plot(erry,np.arange(erry.size))
    vright.set_ylim(0,erry.size-1)
    vright.set_yticks([])
    vright.set_xticks([erry.min(),np.median(erry),erry.max()])
    vright.set_xlabel("$\sigma$")
    vright.xaxis.set_ticks_position('top')
    pl.setp( vright.xaxis.get_majorticklabels(), rotation=70 )
    cax = divider.append_axes("right", size="5%", pad=0.05)
    pl.colorbar(im, cax=cax)

    top = divider.append_axes("top", size="15%", pad=0.05)
    vtop = divider.append_axes("top", size="15%", pad=0.05)
    meanx = dsub.mean(axis=0)
    errx = dsub.std(axis=0)
    top.plot(np.arange(meanx.size),meanx)
    top.set_xlim(0,meanx.size-1)
    top.set_xticks([])
    top.set_yticks([meanx.min(),np.median(meanx),meanx.max()])
    pl.setp(top.yaxis.get_majorticklabels(), rotation=20 )
    top.set_title("$\mu$")
    vtop.plot(np.arange(errx.size),errx,)
    vtop.set_xlim(0,errx.size-1)
    vtop.set_xticks([])
    vtop.set_yticks([errx.min(),np.median(errx),errx.max()])
    vtop.set_ylabel("$\sigma$")
    pl.setp(vtop.yaxis.get_majorticklabels(), rotation=-20 )

    pl.savefig(fullfn.replace(".fits","_data_subbed.png"),bbox_inches='tight')

    makecube.add_file_to_cube(fullfn,
        cubename2+".fits",nhits=cubename2+"_nhits.fits",wcstype='V',
        diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
        negative_mean_cut=-1,
        velocityrange=velocityrange,excludefitrange=[-90,-15,30,90],noisecut=1.0)

makecube.runscript(cubename2)
makecube.make_flats(cubename2,vrange=[30,90],noisevrange=[-15,20])
#makecube.make_taucube(cubename2,cubename2+"_continuum.fits",etamb=0.886)



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
    print "Adding file %s" % fn
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
                              nhits=cubename_supersampled+"_nhits.fits", wcstype='V',
                              diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
                              velocityrange=velocityrange,excludefitrange=[-15,30],noisecut=1.0)
                              # more aggressive noisecut

makecube.runscript(cubename_supersampled)
makecube.make_flats(cubename_supersampled,vrange=[30,90],noisevrange=[-15,30])
