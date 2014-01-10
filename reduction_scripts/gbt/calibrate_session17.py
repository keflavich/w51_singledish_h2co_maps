
from pylab import *
import numpy,matplotlib;
from gbtpy import makecube

tau = 0.02
refpairs = [[10,41],[42,73],[79,97]]#,[110,186]] # see calibrate_session17b
mapnames = ['W51','W51','W51']#,'W51']

filename = '/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17.raw.acs.fits'
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA
if datapfits.DATA[-1,:].sum() == 0:
    print "READING USING PFITS"
    import pfits
    datapfits = pfits.FITS(filename).get_hdus()[1].get_data()
    dataarr = numpy.reshape(datapfits['DATA'],datapfits['DATA'].shape[::-1])

for feed in (1,2):
    for sampler in ("A13","A9","C25","C29"):
        OK = (datapfits['FEED'] == feed) * (datapfits['SAMPLER'] == sampler)
        if OK.sum() > 0:
            plot(datapfits['LST'][OK],dataarr[OK,1000:2000].mean(axis=1),',',label='Feed %i Sampler %s' % (feed,sampler))
savefig('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17_continuum_nolimits.png')
gca().set_ylim(1,5)
savefig('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17_continuum.png')

for (ref1,ref2),mapname in zip(refpairs,mapnames):
    scanrange = [ref1+1,ref2-1]
    outpath = '/Users/adam/observations/gbt/%smap/' % mapname
    makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17.raw.acs.fits',
        outpath+'Session17_%ito%i_A13_F1.fits' %
        (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=1,
        sampler='A13', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr)
    makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17.raw.acs.fits',
        outpath+'Session17_%ito%i_A9_F1.fits' %
        (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=1,
        sampler='A9', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr)
    makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17.raw.acs.fits',
        outpath+'Session17_%ito%i_C25_F2.fits' %
        (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=2,
        sampler='C25', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr)
    makecube.calibrate_cube_data('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17.raw.acs.fits',
        outpath+'Session17_%ito%i_C29_F2.fits' %
        (ref1,ref2),scanrange=scanrange,refscan1=ref1,refscan2=ref2, feednum=2,
        sampler='C29', filepyfits=filepyfits, datapfits=datapfits, tau=tau, dataarr=dataarr)

