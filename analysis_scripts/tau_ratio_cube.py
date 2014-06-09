from astropy.io import fits
from FITS_tools.strip_headers import flatten_header
import numpy as np
from astropy.utils.console import ProgressBar
from astropy.convolution import convolve
from h2co_modeling import SmoothtauModels
from common_constants import TCMB
from paths import datapath,dpath,rpath
import pyregion

h2co11 = fits.getdata(dpath('W51_H2CO11_taucube_supersampled.fits'))
h2co22 = fits.getdata(dpath('W51_H2CO22_pyproc_taucube_lores_supersampled.fits'))

noise11 = h2co11[:50,:,:].std(axis=0)
noise22 = h2co22[:50,:,:].std(axis=0)

sn11 = h2co11/noise11
sn22 = h2co22/noise22

mask = (sn11 > 2) & (sn22 > 2)
ratio = h2co11/h2co22
ratio[True-mask] = np.nan

ratioF = fits.open(dpath('W51_H2CO11_taucube_supersampled.fits'))
ratioF[0].data = ratio
ratioF[0].header['BUNIT'] = ''
ratioF.writeto(dpath('W51_H2CO11_to_22_tau_ratio_supersampled.fits'),clobber=True)

filt = np.ones([3,3,3],dtype='bool')
filt[1,1,1] = 0
nneighbors = convolve(np.isfinite(ratio), filt)
ratio[(nneighbors<7) + (True-np.isfinite(nneighbors))] = np.nan

ratioF[0].data = ratio
ratioF.writeto(dpath('W51_H2CO11_to_22_tau_ratio_supersampled_neighbors.fits'),clobber=True)

nfin = np.isfinite(ratio).sum()
nok = (np.isfinite(sn11)*np.isfinite(sn22)).sum()
print "Unfiltered voxels: ",nfin," of ", nok," or ",nfin/float(nok)*100,"%"
nfinpix = np.isfinite(ratio).max(axis=0).sum()
nokpix = (np.isfinite(sn11)*np.isfinite(sn22)).max(axis=0).sum()
print "Unfiltered pixels ",nfinpix," of ", nokpix," or ",nfinpix/float(nokpix)*100,"%"

ratioF = fits.open(datapath+'W51_H2CO11_taucube_supersampled.fits')


cont11 = fits.getdata(datapath+'W51_H2CO11_cube_supersampled_continuum.fits') + TCMB
cont22 = fits.getdata(datapath+'W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits') + TCMB
cont11[cont11<TCMB] = TCMB
cont22[cont22<TCMB] = TCMB

header = fits.getheader(datapath+'W51_H2CO11_cube_supersampled_continuum.fits')

contfrontregions = pyregion.open(rpath('continuum_in_the_front.reg'))
contfrontmask = contfrontregions.get_mask(fits.PrimaryHDU(data=cont11,header=header))
cont11[contfrontmask] = TCMB
cont22[contfrontmask] = TCMB


def get_extrema(x, w=2, threshold=0.01):
    """
    Return the valid range in pixels.

    This function tries to find the first local maximum and minimum and select
    the range between those.  It can run into trouble if the first/last pixel
    are the min/max; most of the code is dedicated to these corner cases
    """
    from astropy.convolution import convolve, Gaussian1DKernel
    sm = convolve(x, Gaussian1DKernel(w), boundary='extend')
    pos = np.diff(sm) > threshold

    # force first point to match second point
    # smoothing can affect edges
    if pos[0] != pos[1]:
        pos[0] = pos[1]
    if pos[0] != pos[2]:
        pos[:2] = pos[2]

    # Find the first point at which the slope is positive
    # (this is very near the 0/0 point)
    pospt = np.argmax(pos)
    if np.count_nonzero(pos) < 5+4:
        # implies there were no extrema found
        return 0,x.size
    elif pospt == 0:
        # find first *negative* extremum
        negpt = np.argmin(pos)
        pos[:negpt] = False
        pospt = np.argmax(pos)
        if pospt < negpt:
            pospt = pos.size - 1
        return negpt,pospt
    else:
        return 0,pospt


modelpath = '/Users/adam/work/h2co/radex/troscompt_April2013_linearXH2CO/'
stm = SmoothtauModels(modelpath+'1-1_2-2_XH2CO=1e-9_troscompt.dat')

abund = -9
temperature = 20
opr = 0.1
for sigma in (0.1,0.5,1.0):
    #tau1,tau2,dens,col = stm.select_data(abund)
    trot1,trot2,tex1,tex2,tau1,tau2,dens,col = stm.select_data(abundance=abund,
                                                               opr=opr,
                                                               temperature=temperature)
    tau,vtau,vtau_ratio = stm.generate_tau_functions(abundance=abund, opr=opr,
                                                     temperature=temperature)

    tauratio = vtau_ratio(dens, line1=tau1, line2=tau2, tex1=tex1, tex2=tex2,
                          sigma=sigma, opr=opr, temperature=temperature)

    tline,vtline,vtline_ratio,tline1,tline2 = stm.generate_tline_functions(abundance=abund,
                                                                           opr=opr,
                                                                           temperature=temperature)

    """
    Build up "grids" of tex/tau for given backgrounds that can then be ratio'd
    This is more efficient that computing a fresh tauratio array each iteration

    5/28/2014:
    We are now computing OBSERVED tau
    """
    # CACHING:
    #if not ('tau1grid' in locals() and 'tau2grid' in locals()):
    kwargs = dict(sigma=sigma, opr=opr, temperature=temperature)
    tbg1grid = np.hstack([np.linspace(TCMB,15,100),
                          np.logspace(np.log10(15),np.log10(450),15)[1:]])
    tau1grid = [vtline(dens, lvg_tau=tau1, tex=tex1, tbg=tbg1, obs_tau=True, **kwargs)
                for tbg1 in ProgressBar(tbg1grid)]
    tbg2grid = np.logspace(np.log10(TCMB),np.log10(45),100)
    tau2grid = [vtline(dens, lvg_tau=tau2, tex=tex2, tbg=tbg2, obs_tau=True, **kwargs)
                for tbg2 in ProgressBar(tbg2grid)]

    okdens = (dens < 7) & (dens > 2)

    def get_tau_ratio(c1,c2,pos=True,
                      okdens=(dens < 7) & (dens > 2)):
        # find nearest match in the grid
        ind1 = np.argmin(np.abs(tbg1grid-c1))
        ind2 = np.argmin(np.abs(tbg2grid-c2))
        # only absorption observed, therefore force models...
        if pos:
            okpos = okdens & (tau1grid[ind1] > 0) & (tau2grid[ind2] > 0) 
        else:
            raise Exception("Why are you accepting non-positive optical depth ratios?")
            okpos = okdens & True
        ratio = tau1grid[ind1]/tau2grid[ind2]
        # extrema are only measured on the "OK" part of the ratio data
        negpt,pospt = get_extrema(ratio[okpos])
        # therefore we modify the OK part by deselecting everything, then selecting
        # the parts that are good
        first_ind = np.argmax(okpos)
        ok = np.zeros(ratio.size, dtype='bool')
        ok[first_ind+negpt:first_ind+pospt] = True

        if np.count_nonzero(ok) < 3:
            import ipdb; ipdb.set_trace()

        return ratio * ok, ok

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

        outc = (ratio*0).reshape(fshape) + np.nan

        # set up a grid...

        pb = ProgressBar((c11.size))
        for ii,(r,c1,c2) in enumerate(zip(rrs, c11.flat, c22.flat)):
            #print r.shape,c1,c2
            if np.isfinite(c1) and np.isfinite(c2) and np.any(np.isfinite(r)):
                tauratio, ok = get_tau_ratio(c1,c2)

                inds = np.argsort(tauratio[ok])
                outc[:,ii] = np.interp(r, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)
            pb.update()
        #pb.finish()

        return outc.reshape(ratio.shape)

    dcube = ratio_to_dens_slow(ratio,cont11,cont22)

    ratioF = fits.open(datapath+'W51_H2CO11_taucube_supersampled.fits')
    ratioF[0].data = dcube.reshape(ratio.shape)
    ratioF[0].header['BUNIT'] = 'log volume density'
    ratioF.writeto(datapath+'W51_H2CO11_to_22_logdensity_supersampled_textbg_sigma%0.1f.fits' % sigma,
                   clobber=True)


    nfin = np.isfinite(ratio).sum(axis=0)
    ratio[:,nfin < 2] = np.nan

    from scipy.stats import nanmedian
    rmax = np.nanmax(ratio, axis=0)[None,:,:]
    rmin = np.nanmin(ratio, axis=0)[None,:,:]
    rmid = nanmedian(ratio, axis=0)[None,:,:]

    densc = fits.getdata(datapath+'W51_H2CO11_to_22_logdensity_supersampled_textbg_sigma%0.1f.fits' % sigma)
    dens_peak = np.nanmax(densc,axis=0)
    dens_peakf = fits.PrimaryHDU(data=dens_peak, header=flatten_header(header))
    dens_peakf.writeto(datapath+'W51_H2CO_logdensity_textbg_peak.fits',clobber=True)

    d_rmaxprojection = ratio_to_dens_slow(rmax, cont11, cont22)
    d_rminprojection = ratio_to_dens_slow(rmin, cont11, cont22)
    d_rmidprojection = ratio_to_dens_slow(rmid, cont11, cont22)

    dens_rmaxf = fits.PrimaryHDU(data=d_rmaxprojection, header=flatten_header(header))
    dens_rmaxf.writeto(datapath+'W51_H2CO_logdensity_textbg_max_ratio_sigma%0.1f.fits' % sigma,clobber=True)
    dens_rminf = fits.PrimaryHDU(data=d_rminprojection, header=flatten_header(header))
    dens_rminf.writeto(datapath+'W51_H2CO_logdensity_textbg_min_ratio_sigma%0.1f.fits' % sigma,clobber=True)
    dens_rmidf = fits.PrimaryHDU(data=d_rmidprojection, header=flatten_header(header))
    dens_rmidf.writeto(datapath+'W51_H2CO_logdensity_textbg_mid_ratio_sigma%0.1f.fits' % sigma,clobber=True)

    dens_rmaxf = fits.PrimaryHDU(data=rmax, header=flatten_header(header))
    dens_rmaxf.writeto(datapath+'W51_H2CO_max_ratio.fits',clobber=True)
    dens_rminf = fits.PrimaryHDU(data=rmin, header=flatten_header(header))
    dens_rminf.writeto(datapath+'W51_H2CO_min_ratio.fits',clobber=True)
    dens_rmidf = fits.PrimaryHDU(data=rmid, header=flatten_header(header))
    dens_rmidf.writeto(datapath+'W51_H2CO_mid_ratio.fits',clobber=True)

    # OBSOLETE ok = np.arange(tauratio.size) > np.argmax(tauratio)
    # OBSOLETE 
    # OBSOLETE def ratio_to_dens(ratio):
    # OBSOLETE     inds = np.argsort(tauratio[ok])
    # OBSOLETE     return np.interp(ratio, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)
    # OBSOLETE 
    # OBSOLETE dcube = ratio_to_dens(ratio)
    # ratioF[0].data = dcube
    # ratioF[0].header['BUNIT'] = 'log volume density'
    # ratioF.writeto(datapath+'W51_H2CO11_to_22_logdensity_supersampled.fits',clobber=True)
    #
    # densc = fits.getdata(datapath+'W51_H2CO11_to_22_logdensity_supersampled.fits')
    # header = fits.getheader(datapath+'W51_H2CO11_cube_supersampled_continuum.fits')
    # dens_peak = np.nanmax(densc,axis=0)
    # dens_peakf = fits.PrimaryHDU(data=dens_peak,header=flatten_header(header))
    # dens_peakf.writeto(datapath+'W51_H2CO_logdensity_peak.fits',clobber=True)
