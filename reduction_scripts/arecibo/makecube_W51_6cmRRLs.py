# from astropy.io import fits
import os
import re
from agpy import readcol
from astropy import log
from astropy import units as u

from sdpy import makecube

bsgs = readcol('arecibo_bsg_freqref.txt',skipafter=1,asRecArray=True)

for line in bsgs:
    bsg = line['bsg']
    linefreq = line['restfreq']*u.MHz
    linename_ha = linename = line['linename']
    if linefreq == 0:
        continue
    # debug
    #if linename.lower() != 'h112ab':
    #    continue

    vmin = 30
    vmax = 90
    velocityrange = [-150,300]

    #makecube.generate_header(49.523158,-0.34987466,naxis1=96,naxis2=96,pixsize=20,naxis3=1600,cd3=0.5,clobber=True,restfreq=4.8296594e9)
    #makecube.generate_header(49.353568,-0.2982199,naxis1=144,naxis2=96,pixsize=20,naxis3=1600,cd3=0.5,clobber=True,restfreq=4.8296594e9)
    # reduced to CD3 = 1.0, naxis3 = 350 because of size and because the arecibo spectra looked artificially smoothed
    cd3 = 1.0
    naxis3 = int((velocityrange[1]-velocityrange[0]) / cd3)
    # Updated June 18, 2014 to use smaller pixels for more useful comparison to other data
    makecube.generate_header(49.209553, -0.277137, naxis1=308, naxis2=205,
                             pixsize=15, naxis3=int(naxis3), cd3=cd3,
                             crval3=50.0, clobber=True, author='Adam Ginsburg')
    cubename = '/Users/adam/work/h2co/maps/w51/W51_%slpha_cube' % linename
    makecube.make_blank_images(cubename,clobber=True)

    for date in ('0910','0911','0912','0915',):
        fn = '/Users/adam/observations/arecibo/2012{date}/W51_{line}_spectra_{date}.fits'.format(line=linename,date=date)
        if not os.path.exists(fn):
            fn = '/Users/adam/observations/arecibo/2012{date}/W51_{line}_spectra_{date}.fits'.format(line=linename.replace("h","H"),date=date)
        if not os.path.exists(fn):
            fn = '/Users/adam/observations/arecibo/20120910/W51_{line}_spectra_{date}.fits'.format(line=linename.replace("h","H"),date=date)
        log.info(" ".join([str(x) for x in (fn, velocityrange, linename, linefreq)]))
        if os.path.exists(fn):
            makecube.add_file_to_cube(fn,
                                      cubename+'.fits',
                                      nhits=cubename+'_nhits.fits',
                                      velocityrange=velocityrange,
                                      excludefitrange=[vmin,vmax], 
                                      progressbar=True,
                                      add_with_kernel=True,
                                      kernel_fwhm=20./3600.,
                                      chmod=True,
                                      linefreq=linefreq)
        else:
            raise IOError("Did not find file {0}".format(fn))
    #    f = pyfits.open(fn)
    #
    #    fixed=False
    #    for kw in ('TDIM%i' % ii for ii in xrange(5)):
    #        if kw in f[0].header:
    #            fixed = True
    #            del f[0].header[kw]
    #    if fixed:
    #        f[0].writeto(fn.replace(".fits","_fixed.fits"),output_verify='fix')
    #        fn = fn.replace(".fits","_fixed.fits")


    flat_vrange = [30,90]

    if os.path.exists(cubename+".fits"):
        makecube.runscript(cubename)
        makecube.make_flats(cubename,vrange=flat_vrange,noisevrange=[-50,-1])
        makecube.make_flats(cubename,vrange=[-100,-40],noisevrange=[-50,-1],
                            out_suffix='_Helium')
        #makecube.make_flats(cubename,vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')
    else:
        print "Failure at ",line
        
    # 1/10/2014: add "superresolution" maps
    vmin = 30
    vmax = 90
    velocityrange = [-50,150]

    #makecube.generate_header(49.523158,-0.34987466,naxis1=96,naxis2=96,pixsize=20,naxis3=1600,cd3=0.5,clobber=True,restfreq=4.8296594e9)
    #makecube.generate_header(49.353568,-0.2982199,naxis1=144,naxis2=96,pixsize=20,naxis3=1600,cd3=0.5,clobber=True,restfreq=4.8296594e9)
    # reduced to CD3 = 1.0, naxis3 = 350 because of size and because the arecibo spectra looked artificially smoothed
    cd3 = 1.0
    crval3 = 50.0
    naxis3 = int((velocityrange[1]-velocityrange[0]) / cd3)
    makecube.generate_header(49.209553,-0.277137,naxis1=308,naxis2=205,pixsize=15,naxis3=int(naxis3),cd3=cd3,crval3=crval3,clobber=True,
                             restfreq=linefreq, author='Adam Ginsburg')
    cubename = '/Users/adam/work/h2co/maps/w51/W51_%slpha_cube_supersampled' % linename
    makecube.make_blank_images(cubename,clobber=True)

    for date in ('0910','0911','0912','0915',):
        fn = '/Users/adam/observations/arecibo/2012{date}/W51_{line}_spectra_{date}.fits'.format(line=linename,date=date)
        log.info(" ".join([str(x) for x in (fn, velocityrange, linename, linefreq)]))
        if os.path.exists(fn):
            makecube.add_file_to_cube(fn,
                                      cubename+'.fits',
                                      add_with_kernel=True,
                                      kernel_fwhm=20./3600.,
                                      nhits=cubename+'_nhits.fits',
                                      velocityrange=velocityrange,
                                      excludefitrange=[vmin,vmax],
                                      linefreq=linefreq,
                                      progressbar=True,
                                      chmod=True) # security risk, but too many files!
    #    f = pyfits.open(fn)
    #
    #    fixed=False
    #    for kw in ('TDIM%i' % ii for ii in xrange(5)):
    #        if kw in f[0].header:
    #            fixed = True
    #            del f[0].header[kw]
    #    if fixed:
    #        f[0].writeto(fn.replace(".fits","_fixed.fits"),output_verify='fix')
    #        fn = fn.replace(".fits","_fixed.fits")


    flat_vrange = [20,100]

    if os.path.exists(cubename+".fits"):
        makecube.runscript(cubename)
        makecube.make_flats(cubename,vrange=flat_vrange,noisevrange=[-50,-1])
        #makecube.make_flats(cubename,vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')
    else:
        print "Failure at ",line


    # 6/18/2014: add "superresolution" maps of HELIUM

    helium_freq = {112: 4.62067,
                   111: 4.74612,
                   110: 4.87614,
                   109: 5.01096,
                   108:  5.1508,
                   107: 5.29589,
                   106: 5.44648,}

    transition = int(re.findall('[0-9]+', linename)[0])
    linefreq = helium_freq[transition] * u.GHz
    linename = 'he{0}a'.format(transition)

    makecube.generate_header(49.209553,-0.277137,naxis1=308,naxis2=205,pixsize=15,naxis3=int(naxis3),cd3=cd3,crval3=crval3,clobber=True,
                             restfreq=linefreq, author='Adam Ginsburg')
    cubename = '/Users/adam/work/h2co/maps/w51/W51_%slpha_cube_supersampled' % linename
    makecube.make_blank_images(cubename,clobber=True)

    for date in ('0910','0911','0912','0915',):
        fn = '/Users/adam/observations/arecibo/2012{date}/W51_{line}_spectra_{date}.fits'.format(line=linename_ha,date=date)
        log.info(" ".join([str(x) for x in (fn, velocityrange, linename, linefreq)]))
        if os.path.exists(fn):
            makecube.add_file_to_cube(fn,
                                      cubename+'.fits',
                                      add_with_kernel=True,
                                      kernel_fwhm=20./3600.,
                                      nhits=cubename+'_nhits.fits',
                                      velocityrange=velocityrange,
                                      excludefitrange=[vmin,vmax],
                                      do_runscript=True,
                                      progressbar=True,
                                      linefreq=linefreq, # MAY NEED TO BE MHZ?!  (June 2014)
                                      chmod=True) # security risk, but too many files!
        else:
            log.info("File {0} does not exist".format(fn))
    #    f = pyfits.open(fn)
    #
    #    fixed=False
    #    for kw in ('TDIM%i' % ii for ii in xrange(5)):
    #        if kw in f[0].header:
    #            fixed = True
    #            del f[0].header[kw]
    #    if fixed:
    #        f[0].writeto(fn.replace(".fits","_fixed.fits"),output_verify='fix')
    #        fn = fn.replace(".fits","_fixed.fits")


    flat_vrange = [20,100]

    if os.path.exists(cubename+".fits"):
        makecube.runscript(cubename)
        makecube.make_flats(cubename,vrange=flat_vrange,noisevrange=[-50,-1])
        #makecube.make_flats(cubename,vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')
    else:
        print "Failure at ",line
