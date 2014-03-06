import pylab as pl
import numpy as np
from astropy.io import fits as pyfits
from agpy import asinh_norm
from sdpy import makecube

linefreq = 14.488479e9

filelist = [
    "Session10_32to93_A13_F1.fits",
    # gigantic jump = bad "Session10_32to93_A9_F1.fits",
    # blech what the heck things on the left "Session10_32to93_C25_F2.fits",
    # negative continuum "Session10_32to93_C29_F2.fits",
    "Session11_126to157_A13_F1.fits",
    "Session11_126to157_A9_F1.fits",
    "Session11_126to157_C25_F2.fits",
    "Session11_126to157_C29_F2.fits",
    "Session11_158to189_A13_F1.fits",
    "Session11_158to189_A9_F1.fits",
    "Session11_158to189_C25_F2.fits",
    #"Session11_158to189_C29_F2.fits", # Blurred (triples)!   Could these be issues in the FITS packing?
    "Session11_21to52_A13_F1.fits",
    "Session11_21to52_A9_F1.fits",
    "Session11_21to52_C25_F2.fits", # small blur jump (not worth discarding?)
    "Session11_21to52_C29_F2.fits", # Major TSYS dip?
    "Session11_53to78_A13_F1.fits",
    "Session11_53to78_A9_F1.fits",
    "Session11_53to78_C25_F2.fits",
    "Session11_53to78_C29_F2.fits",
    "Session11_89to120_A13_F1.fits",
    "Session11_89to120_A9_F1.fits",
    #"Session11_89to120_C25_F2.fits",  Blurred in the top half
    "Session11_89to120_C29_F2.fits",
    "Session14_12to43_A13_F1.fits",
    "Session14_12to43_A9_F1.fits",
    "Session14_12to43_C25_F2.fits",
    "Session14_12to43_C29_F2.fits",
    "Session14_130to161_A13_F1.fits",
    "Session14_130to161_A9_F1.fits",
    "Session14_130to161_C25_F2.fits",
    "Session14_130to161_C29_F2.fits",
    "Session14_167to198_A13_F1.fits",
    "Session14_167to198_A9_F1.fits",
    "Session14_167to198_C25_F2.fits",
    "Session14_167to198_C29_F2.fits",
    "Session14_44to75_A13_F1.fits",
    "Session14_44to75_A9_F1.fits",
    "Session14_44to75_C25_F2.fits",
    #"Session14_44to75_C29_F2.fits", Blurred (triples)
    "Session14_81to96_A13_F1.fits",
    "Session14_81to96_A9_F1.fits",
    "Session14_81to96_C25_F2.fits",
    "Session14_81to96_C29_F2.fits",
    "Session14_98to129_A13_F1.fits",
    "Session14_98to129_A9_F1.fits",
    "Session14_98to129_C25_F2.fits",
    "Session14_98to129_C29_F2.fits",
    "Session16_10to41_A13_F1.fits",
    "Session16_10to41_A9_F1.fits",
    "Session16_10to41_C25_F2.fits",
    "Session16_10to41_C29_F2.fits",
    "Session16_42to73_A13_F1.fits",
    "Session16_42to73_A9_F1.fits",
    "Session16_42to73_C25_F2.fits",
    "Session16_42to73_C29_F2.fits",
    "Session16_74to105_A13_F1.fits",
    "Session16_74to105_A9_F1.fits",
    "Session16_74to105_C25_F2.fits",
    "Session16_74to105_C29_F2.fits",
    "Session17_10to41_A13_F1.fits",
    "Session17_10to41_A9_F1.fits",
    "Session17_10to41_C25_F2.fits",
    "Session17_10to41_C29_F2.fits",
    "Session17_42to73_A13_F1.fits",
    "Session17_42to73_A9_F1.fits",
    "Session17_42to73_C25_F2.fits",
    "Session17_42to73_C29_F2.fits",
    "Session17_110to186_A13_F1.fits",
    "Session17_110to186_A9_F1.fits",
    #"Session17_110to186_C25_F2.fits", Blurred
    "Session17_110to186_C29_F2.fits",
]
filelist2=[
    "Session22_10to41_A13_F1.fits",
    "Session22_10to41_A9_F1.fits",
    # neg stripes "Session22_10to41_C25_F2.fits",
    # neg stripes "Session22_10to41_C29_F2.fits",
    "Session22_183to203_A13_F1.fits",
    "Session22_183to203_A9_F1.fits",
    # neg stripes "Session22_183to203_C25_F2.fits",
    # neg stripes "Session22_183to203_C29_F2.fits",
    "Session22_42to73_A13_F1.fits",
    "Session22_42to73_A9_F1.fits",
    # neg stripes "Session22_42to73_C25_F2.fits",
    # neg stripes "Session22_42to73_C29_F2.fits",
    "Session22_74to105_A13_F1.fits",
    "Session22_74to105_A9_F1.fits",
    # neg stripes "Session22_74to105_C25_F2.fits",
    # neg stripes "Session22_74to105_C29_F2.fits",
]

if False: # old cube, don't care any more
    cubename_lores = '/Users/adam/work/h2co/maps/w51/W51_H2CO22_pyproc_cube_lores'
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
    #    raise ValueError("Add more points.  Something's going to be out of range for stupid star stupid link")
    makecube.generate_header(49.209553, -0.277137, naxis1=192, naxis2=128,
                             pixsize=24, naxis3=naxis3, cd3=cd3, crval3=crval3,
                             clobber=True, restfreq=linefreq)
    makecube.make_blank_images(cubename_lores,clobber=True)


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
                                  cubename_lores+".fits",
                                  nhits=cubename_lores+"_nhits.fits", wcstype='V',
                                  diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
                                  velocityrange=velocityrange,excludefitrange=[40,75],noisecut=np.inf)
                                  # more aggressive noisecut

    makecube.runscript(cubename_lores)
    makecube.make_flats(cubename_lores,vrange=[45,75],noisevrange=[-15,30])
    makecube.make_taucube(cubename_lores, cubename_lores+"_continuum.fits",
                          etamb=0.886, linefreq=linefreq*u.Hz, tex=2.0)


cubename_lores_supersampled = '/Users/adam/work/h2co/maps/w51/W51_H2CO22_pyproc_cube_lores_supersampled'
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
#    raise ValueError("Add more points.  Something's going to be out of range for stupid star stupid link")
makecube.generate_header(49.209553, -0.277137, naxis1=308, naxis2=205,
                         pixsize=15, naxis3=naxis3, cd3=cd3, crval3=crval3,
                         clobber=True, restfreq=14.488479e9)
makecube.make_blank_images(cubename_lores_supersampled,clobber=True)


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
                              cubename_lores_supersampled+".fits",
                              add_with_kernel=True,
                              kernel_fwhm=20./3600.,
                              nhits=cubename_lores_supersampled+"_nhits.fits", wcstype='V',
                              diagnostic_plot_name=fullfn.replace('.fits','_data_scrubbed.png'),
                              velocityrange=velocityrange,excludefitrange=[40,75],noisecut=np.inf,
                              continuum_prefix=cubename_lores_supersampled+fn.replace(".fits",''),
                              negative_mean_cut=-1)

makecube.runscript(cubename_lores_supersampled)
makecube.make_flats(cubename_lores_supersampled,vrange=[45,75],noisevrange=[-15,30])
makecube.make_taucube(cubename_lores_supersampled,cubename_lores_supersampled+"_continuum.fits",
                      etamb=0.886, linefreq=linefreq*u.Hz, tex=2.0)
