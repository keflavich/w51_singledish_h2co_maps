from astroquery import splatalogue
from astropy.io import fits
from astropy import units as u
import numpy as np
import glob
import re

def rrl(n,dn=1,amu=1.007825):    # compute Radio Recomb Line feqs in GHz
    # from Brown, Lockman & Knapp ARAA 1978 16 445
    nu = 3.289842e6*(1-5.48593e-4/amu)*(1/float(n)**2 - 1/float(n+dn)**2)
    return nu


result = splatalogue.core.Splatalogue.query_lines(4*u.GHz, 6*u.GHz, chemical_name='Hydrogen Recombination')
halphas = result[np.array(['alpha' in x for x in result['Species']])]

# avoid 0910 because of misconfigured backend
files = glob.glob('/Volumes/128gbdisk/arecibo/*0911*02900.fits')

bsgre = re.compile('(b[0-7]s[0-9]g[01])')

lines = []

with open('/Users/adam/work/h2co/arecibo_codes/arecibo_bsg_freqref.txt','w') as f:

    print >>f,"%6s %15s %15s %15s %10s %15s" % ('bsg','minfreq','maxfreq','bandwidth','linename','restfreq')
    print >>f,"%6s %15s %15s %15s %10s %15s" % ('---','-------','-------','---------','--------','--------')

    for fn in files:
        d = fits.getdata(fn)
        freq = d['CRVAL1'][0] + d['CDELT1'][0]*(np.arange(8192)+1-d['CRPIX1'][0])
        # d['data'][0].size  is 16384, but shouldn't be
        #freqTopo=(findgen(nchan) - (hdr[0].crpix1 -1))*hdr[0].cdelt1*1d-6 +$
        #         hdr[0].crval1*1d-6
        minfreq = freq.min() / 1e6 # MHz
        maxfreq = freq.max() / 1e6 # MHz
        bsg = bsgre.findall(fn)[0]
        
        sel = (halphas['Freq-GHz']*1000 < maxfreq) & (halphas['Freq-GHz']*1000 > minfreq)
        if sel.sum() > 0:
            haname = halphas[sel]['Resolved QNs'][0].replace("&alpha;","a").translate(None,'()').lower()
            if haname in lines:
                haname += 'b'
            if haname in lines:
                haname = haname[:-1] + 'c'
            if haname in lines:
                haname = haname[:-1] + 'd'
            lines.append(haname)
            print >>f,"%6s %15f %15f %15f %10s %15f" % (bsg,minfreq,maxfreq,maxfreq-minfreq,haname,halphas[sel]['Freq-GHz'][0]*1000)
        else:
            print >>f,"%6s %15f %15f %15f %10s %15f" % (bsg,minfreq,maxfreq,maxfreq-minfreq,'none',0)
