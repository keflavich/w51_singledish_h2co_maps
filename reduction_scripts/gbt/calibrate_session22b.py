from pylab import *
import numpy,matplotlib;
from gbtpy import makecube
import time

tau = 0.02

filename = '/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22.raw.acs.fits'
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA
if datapfits.DATA[-1,:].sum() == 0:
    import pdb; pdb.set_trace()
    print "READING USING PFITS"
    t0 = time.time()
    import pfits
    datapfits = pfits.FITS(filename).get_hdus()[1].get_data()
    print "Successfully read in %i seconds" % (time.time() - t0)
    dataarr = numpy.reshape(datapfits['DATA'],datapfits['DATA'].shape[::-1])

figure(1)
for feed in (1,2):
    for sampler in ("A13","A9","C25","C29"):
        OK = (datapfits['FEED'] == feed) * (datapfits['SAMPLER'] == sampler)
        if OK.sum() > 0:
            plot(datapfits['LST'][OK],dataarr[OK,1000:2000].mean(axis=1),',',label='Feed %i Sampler %s' % (feed,sampler))
gca().set_ylim(1,5)
savefig('/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22_continuum.png')
close(1)



scanrange=[184,206]
ref1 = 183
ref2 = 203
refscans = [183,188,193,198,203]
sourcename = "W51WestExtension"
obsmode = "DecLatMap"
mapname = 'W51'
outpath = '/Users/adam/observations/gbt/%smap/' % mapname
makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22.raw.acs.fits',
    outpath+'Session22_%ito%i_A13_F1.fits' %
    (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=1, refscans=refscans,
    sampler='A13', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr, obsmode=obsmode, sourcename=sourcename)
makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22.raw.acs.fits',
    outpath+'Session22_%ito%i_A9_F1.fits' %
    (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=1, refscans=refscans,
    sampler='A9', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr, obsmode=obsmode, sourcename=sourcename)
makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22.raw.acs.fits',
    outpath+'Session22_%ito%i_C25_F2.fits' %
    (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=2, refscans=refscans,
    sampler='C25', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr, obsmode=obsmode, sourcename=sourcename)
makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22.raw.acs.fits',
    outpath+'Session22_%ito%i_C29_F2.fits' %
    (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=2, refscans=refscans,
    sampler='C29', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr, obsmode=obsmode, sourcename=sourcename)

