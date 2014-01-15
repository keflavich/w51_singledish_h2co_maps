"""
Note to self: Write more notes to self about creation dates.
"""
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from pylab import *
import numpy,matplotlib;
import sys
import os
from gbtpy import makecube,calibrate_map_scans
import timer

samplers = {
        0: ["A9","A13","C25","C29"],
        1: ["A10","A14","C26","C30"],
        2: ["B17","B21","D33","D37"],
        3: ["B18","B22","D34","D38"],
        }

feeds = {
        0: [1,1,2,2],
        1: [1,1,2,2],
        2: [1,1,2,2],
        3: [1,1,2,2]
        }

def calibrate_session(filename, mapnames, ifnum=0,
        gbtpath='/Users/adam/observations/gbt/', sessionnumber=None,
        refscans_list=None, scanranges=None, tau=0.02, 
        obsmodes=None, sourcenames=None):

    root = filename
    while root != os.path.splitext(root)[0]:
        root = os.path.splitext(root)[0]

    if sessionnumber is None:
        sessionnumber = int(root[-2:])

    filepyfits = pyfits.open(filename,memmap=True)
    datapfits = filepyfits[1].data
    dataarr = datapfits.DATA
    if datapfits.DATA[-1,:].sum() == 0:
        print "READING USING PFITS"
        import pfits
        datapfits = pfits.FITS(filename).get_hdus()[1].get_data()
        dataarr = numpy.reshape(datapfits['DATA'],datapfits['DATA'].shape[::-1])

    for feed in (1,2):
        for sampler in samplers[ifnum]:
            OK = (datapfits['FEED'] == feed) * (datapfits['SAMPLER'] == sampler)
            if OK.sum() > 0:
                plot(datapfits['LST'][OK],dataarr[OK,1000:2000].mean(axis=1),',',label='Feed %i Sampler %s' % (feed,sampler))
    savefig(root+'_secondhalf_continuum_nolimits.png')
    gca().set_ylim(1,5)
    savefig(root+'_secondhalf_continuum.png')

    if scanranges is None:
        scanranges = [(min(r)+1,max(r)-1) for r in refscans_list]
    if obsmodes is None:
        obsmodes = [None for r in refscans_list]
    if sourcenames is None:
        sourcenames = [None for r in refscans_list]

    for refscans,scanrange,mapname,obsmode,sourcename in zip(refscans_list,scanranges,mapnames,obsmodes,sourcenames):
        outpath = '%s%smap/' % (gbtpath,mapname)
        ref1,ref2 = min(refscans),max(refscans)
        for sampler,feednum in zip(samplers[ifnum],feeds[ifnum]):
            calibrate_map_scans.calibrate_cube_data(filename,
                    outpath+'Session%i_%ito%i_%s_F%i.fits' %
                    (sessionnumber,ref1,ref2,sampler,feednum),
                    scanrange=scanrange, feednum=feednum, sampler=sampler,
                    filepyfits=filepyfits, datapfits=datapfits, tau=tau,
                    dataarr=dataarr, refscans=refscans, sourcename=sourcename,
                    obsmode=obsmode)

if __name__ == "__main__":

    refpairs = {
            10:[[32,93]],
            11:[[21,52],[53,78],[89,120],[126,157],[158,189]],
            13:[[78,94]],
            14:[[12,43],[44,75],[81,96],[98,129],[130,161],[167,198]],
            15:[[8,26],[27,41],[48,66],[67,81],[87,102],[103,114]],
            16:[[10,41],[42,73],[74,105]],
            17:[[10,41],[42,73],[79,97]],
            18:[[25,45],[46,66]],
            20:[[38,66],[69,106],[114,142],[145,182]],
            21:[[8,23],[24,34],[35,56],[57,77],[78,96],[97,111]],
            22:[[10,41],[42,73],[74,105],[112,149],[152,180]],
            }
    mapnames = {
            10:['W51'],
            11:['W51','W51','W51','W51','W51'],
            13:['S233IR'],
            14:['W51','W51','W51','W51','W51','W51',],
            15:['GemOB1','GemOB1','GemOB1','GemOB1','S233IR','S233IR'],
            16:['W51','W51','W51'],
            17:['W51','W51','W51'],
            18:['G202','G202'],
            20:['G34','G34','G34'],
            21:['S233IR','S233IR','G202','G202','GemOB1','GemOB1'],
            22:['W51','W51','W51','G34','G34'],
            } 
    tau = {
            10:0.02,
            11:0.02,
            13:0.02,
            14:0.02,
            15:0.02,
            16:0.02,
            17:0.02,
            18:(0.0121710+0.0117295)/2.0,
            20:(0.0211061 + 0.0259094)/2.0,
            21:(0.0121710+0.0117295)/2.0,
            22:0.02,
            }

    for ifnum in [3,2,0,1]:
        for sessionnumber in tau.keys():
            print "Session %i IF %i: " % (sessionnumber,ifnum)
            timer.print_timing(calibrate_session('/Users/adam/observations/gbt/AGBT10B_019_%i/AGBT10B_019_%i.raw.acs.fits' % (sessionnumber,sessionnumber),
                    mapnames[sessionnumber], ifnum=ifnum, tau=tau[sessionnumber],
                    sessionnumber=sessionnumber, refscans_list=refpairs[sessionnumber]))

        timer.print_timing(calibrate_session('/Users/adam/observations/gbt/AGBT10B_019_22/AGBT10B_019_22.raw.acs.fits',
                ["W51"], ifnum=ifnum,scanranges=[[184,206]],
                refscans_list=[[183,188,193,198,203]],
                sourcenames=["W51WestExtension"],
                obsmodes=["DecLatMap"]))

        timer.print_timing(calibrate_session('/Users/adam/observations/gbt/AGBT10B_019_17/AGBT10B_019_17_secondhalf.raw.acs.fits',
                ["W51"], ifnum=ifnum,scanranges=[[112,185]],
                sessionnumber=17,
                refscans_list=[[111,116,121,126,131,136,141,146,151,156,161,166,171,176,181,186]],
                sourcenames=["W51WestExtension"],
                obsmodes=["DecLatMap"]))
