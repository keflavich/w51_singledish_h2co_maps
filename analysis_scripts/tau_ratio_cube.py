from astropy.io import fits
import numpy as np
import progressbar

h2co11 = fits.getdata('W51_H2CO11_taucube_supersampled.fits')
h2co22 = fits.getdata('W51_H2CO22_pyproc_taucube_lores_supersampled.fits')

noise11 = h2co11[:50,:,:].std(axis=0)
noise22 = h2co22[:50,:,:].std(axis=0)

sn11 = h2co11/noise11
sn22 = h2co22/noise22

mask = (sn11 > 2) & (sn22 > 2)
ratio = h2co11/h2co22
ratio[True-mask] = np.nan

ratioF = fits.open('W51_H2CO11_taucube_supersampled.fits')
ratioF[0].data = ratio
ratioF[0].header['BUNIT'] = 'none'
ratioF.writeto('W51_H2CO11_to_22_tau_ratio_supersampled.fits',clobber=True)

from astropy.convolution import convolve

filt = np.ones([3,3,3],dtype='bool')
filt[1,1,1] = 0
nneighbors = convolve(np.isfinite(ratio), filt)
ratio[(nneighbors<7) + (True-np.isfinite(nneighbors))] = np.nan

ratioF[0].data = ratio
ratioF.writeto('W51_H2CO11_to_22_tau_ratio_supersampled_neighbors.fits',clobber=True)

nfin = np.isfinite(ratio).sum()
nok = (np.isfinite(sn11)*np.isfinite(sn22)).sum()
print "Unfiltered voxels: ",nfin," of ", nok," or ",nfin/float(nok)*100,"%"
nfinpix = np.isfinite(ratio).max(axis=0).sum()
nokpix = (np.isfinite(sn11)*np.isfinite(sn22)).max(axis=0).sum()
print "Unfiltered pixels ",nfinpix," of ", nokpix," or ",nfinpix/float(nokpix)*100,"%"

#import sys
#sys.path.append("/Users/adam/work/h2co/lowdens/code/")
from h2co_modeling import SmoothtauModels
stm = SmoothtauModels('/Users/adam/work/h2co/radex/troscompt_April2013_linearXH2CO/1-1_2-2_XH2CO=1e-9_troscompt.dat')

abund = -9
temperature = 20
opr = 0.1
sigma = 1.0
#tau1,tau2,dens,col = stm.select_data(abund)
trot1,trot2,tex1,tex2,tau1,tau2,dens,col = stm.select_data(abundance=abund,
                                                           opr=opr,
                                                           temperature=temperature)
tau,vtau,vtau_ratio = stm.generate_tau_functions(abundance=abund, opr=opr,
                                                 temperature=temperature)

tauratio = vtau_ratio(dens, line1=tau1, line2=tau2, tex1=tex1, tex2=tex2,
                      sigma=sigma, opr=opr, temperature=temperature)

ok = np.arange(tauratio.size) > np.argmax(tauratio)

def ratio_to_dens(ratio):
    inds = np.argsort(tauratio[ok])
    return np.interp(ratio, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)

dcube = ratio_to_dens(ratio)

ratioF = fits.open('W51_H2CO11_taucube_supersampled.fits')
ratioF[0].data = dcube
ratioF[0].header['BUNIT'] = 'log volume density'
ratioF.writeto('W51_H2CO11_to_22_logdensity_supersampled.fits',clobber=True)

cont11 = fits.getdata('W51_H2CO11_cube_supersampled_continuum.fits')
cont22 = fits.getdata('W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits')

def ratio_to_dens_slow(ratio, c11, c22):
    """
    Shape:
        ratio [z,y,x]
        c11 [y,x]
        c22 [y,x]
    """

    assert c11.size == c22.size == ratio[0,:,:].size

    fshape = [ratio.shape[0], ratio.shape[1]*ratio.shape[2]]
    rrs = ratio.reshape(fshape).T

    outc = (ratio*0).reshape(fshape)

    pb = progressbar.ProgressBar(maxval=np.isfinite(c11*c22).sum())
    pb.start()
    count = 0
    for ii,(r,c1,c2) in enumerate(zip(rrs, c11.flat, c22.flat)):
        #print r.shape,c1,c2
        if np.isfinite(c1) and np.isfinite(c2):
            tauratio = vtau_ratio(dens, line1=tau1, line2=tau2, tex1=tex1,
                                  tex2=tex2,
                                  tbg1=c1+2.73, tbg2=c2+2.73, sigma=sigma,
                                  opr=opr, temperature=temperature)

            ok = np.arange(tauratio.size) > np.argmax(tauratio)

            inds = np.argsort(tauratio[ok])
            outc[:,ii] = np.interp(r, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)
            count += 1
        pb.update(count)
    pb.finish()

    return outc

dcube = ratio_to_dens_slow(ratio,cont11,cont22)

ratioF = fits.open('W51_H2CO11_taucube_supersampled.fits')
ratioF[0].data = dcube
ratioF[0].header['BUNIT'] = 'log volume density'
ratioF.writeto('W51_H2CO11_to_22_logdensity_supersampled_textbg_sigma%0.1f.fits' % sigma,
               clobber=True)
