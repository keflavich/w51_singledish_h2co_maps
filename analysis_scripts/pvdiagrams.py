import pyregion
from astropy.io import fits
from astropy import wcs
import pvextractor
import spectral_cube
from scipy.ndimage import map_coordinates
import numpy as np
from astropy.convolution import convolve
from astropy import convolution
from astropy import coordinates
from astropy import units as u
import scipy.signal
from paths import datapath,datapath_w51
import os

import pylab as pl

def clean(s, g, niter=100, gain=0.1, threshold=1.25):
    """
    A trivial 1D implementation of the CLEAN algorithm
    """
    g = g/g.max()
    result = np.zeros_like(s)
    kern = np.zeros_like(s)
    mid = kern.size/2
    kern[mid] = 1
    sig = np.convolve(kern,g,mode='same')
    for ii in xrange(niter):
        peak = s.argmax()
        mx = s[peak]
        result[peak] += mx*gain
        sig *= mx*gain
        sig = np.roll(sig,peak-mid)
        s = s - sig
        sig /= mx*gain
    return result

cm = pl.cm.gray_r

# number of subplots
nxsp = 3
nysp = 3

titles = ('H$_2$CO 1-1',
          '$^{13}$CO 3-2',
          '$^{12}$CO 3-2',
          'H$_2$CO 2-2',
          '$^{13}$CO 2-1',
          '$^{12}$CO 2-1',
          'C$^{18}$O 3-2',
          '$^{13}$CO 1-0',
          'H$_2$CO n(H$_2$)')
          #'H 112$\\alpha$')
cubefiles = ('W51_H2CO11_taucube.fits',
             '13co_final_cube_c.fits',
             '12co_final_cube_c.fits',
             'W51_H2CO22_pyproc_taucube_lores.fits',
             'w51_bieging_13co32.fits', # these are really CO 2-1
             'w51_bieging_12co32.fits', # these are really CO 2-1
             'c18o_final_cube_c.fits',
             'grs_48and50_cube.fits',
             'W51_h112alpha_modelcube.fits')
cubefiles = ('h2co_singledish/W51_H2CO11_taucube_supersampled.fits',
             '13co_final_cube_c_supersampledh2cogrid.fits',
             '12co_final_cube_c_supersampledh2cogrid.fits',
             'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits',
             'w51_bieging_13co32_supersampledh2cogrid.fits',
             'w51_bieging_12co32_supersampledh2cogrid.fits',
             'c18o_final_cube_c_supersampledh2cogrid.fits',
             'grs_48and50_cube_supersampledh2cogrid.fits',
             'h2co_singledish/W51_H2CO11_to_22_logdensity_supersampled.fits')
             #'W51_h112alpha_modelcube.fits')

def get_pvs(cubefn, endpoints):
    fullcubefn = os.path.join(datapath_w51, cubefn)
    cube = spectral_cube.SpectralCube.read(fullcubefn)
    velo = cube.spectral_axis
    cdelt = pvextractor.utils.wcs_utils.get_spectral_scale(cube.wcs)
    #cube,velo,cdelt = pvextractor.utils.get_cube_info(cubefn)
    pvPath = pvextractor.geometry.path.Path(endpoints, width=60*u.arcsec)
    pv = pvextractor.extract_pv_slice(cube, pvPath)
                                      #respect_nan=False)
    npv = len(endpoints)
    return (pv,npv,velo,cdelt)


if __name__ == "__main__":
    endpoints_wcs = pyregion.open('pvendpoints.reg')

    if not 'colorpvs' in locals():
        colorpvs = {}

    for jj,color in enumerate(('green','red','blue','purple','cyan','yellow','orange')):
        coords = np.array([s.coord_list for s in endpoints_wcs if s.attr[1]['color'] == color])

        pl.figure(jj+1)
        pl.clf()

        if color in colorpvs: #caching^2
            print "Loading from cache: ",color
            pvs = colorpvs[color]
        else:
            print "Processing from scratch: ",color
            pvs = {} # caching
        vmin,vmax = 45,75

        for ii,cubefn in enumerate(cubefiles):
            
            title = titles[ii]
            fullcubefn = os.path.join(datapath_w51, cubefn)
            header = fits.getheader(fullcubefn)
            ww = wcs.WCS(header)

            gcoords = coordinates.Galactic(coords[:,0],coords[:,1], unit=(u.deg,u.deg))
            if 'RA' in ww.wcs.ctype[0]:
                rc = np.array([gcoords.fk5.ra.deg,gcoords.fk5.dec.deg]).T
                #endpoints = ww.wcs_sky2pix(rc,0)
                endpoints = rc
            else:
                #endpoints = ww.wcs_sky2pix(coords,0)
                endpoints = gcoords

            if cubefn in pvs: # caching
                print "Loading ",cubefn," from cache"
                pv, npv, velo, cdelt = pvs[cubefn]
            else:
                print "Extracting ",cubefn
                (pv,npv,velo,cdelt) = pvs[cubefn] = get_pvs(cubefn, endpoints)

            velo = velo.to(u.km/u.s)
            cdelt = cdelt.to(u.km/u.s)
            cdeltv = (velo[1]-velo[0]).value

            vmask = (velo.value < vmax) * (velo.value > vmin)

            #if ii==5:
            #    g = np.exp(-np.linspace(-5,5,101)**2/2.)
            #    g2 = np.exp(-np.linspace(-5,5,101)**2/(2.*0.25**2))
            #    g2 /= g2.sum()
            #    for mm in np.arange(pv.shape[1]):
            #        w = scipy.signal.wiener(pv[:,mm], 9)
            #        cl = clean(w,g,niter=1000,threshold=2)
            #        cn = np.convolve(cl, g2, mode='same')
            #        pv[:,mm] = cn + np.random.randn(cn.size)*1.0

            hdu = pv
            pv = pv.data

            pl.subplot(nxsp,nysp,ii+1)
            aspect = float(pv.shape[1])/float(vmask.sum()) * cdeltv / ((vmax-vmin)/float(vmask.sum()))
            index_order = 1 if cdeltv>0 else -1

            pvtoshow = pv[::index_order,:][vmask[::index_order],:]

            min_val = min([pvtoshow.min(),0])
            max_val = 30 if 'h' not in cubefn.lower() else None

            pl.imshow(pvtoshow,
                      extent=[pv.shape[1]*cdeltv,0,vmin,vmax],
                      aspect=aspect,
                      cmap=cm,
                      vmin=min_val,
                      vmax=max_val,
                      origin='lower')

            pl.vlines(np.cumsum(npv)*cdeltv,vmin,vmax,color='b',linestyle='--',linewidth=2,alpha=0.5)
            pl.hlines(65,0,pv.shape[1]*cdeltv, color='b', linestyle=':',alpha=0.5)
            pl.annotate(title,(0.05,0.9),xycoords='axes fraction',color='k')

        colorpvs[color] = pvs

        # regrid tau22 onto tau11 after smoothing along the velo axis
        #pv1,npv1,velo1,cdelt1 = pvs[datapath_w51+'h2co_singledish/W51_H2CO11_taucube_supersampled.fits']
        #pv2,npv2,velo2,cdelt2 = pvs[datapath_w51+'h2co_singledish/W51_H2CO22_pyproc_taucube_lores_supersampled.fits']
        pv2,npv2,velo2,cdelt2 = pvs['h2co_singledish/W51_H2CO11_to_22_logdensity_supersampled.fits']
        #vmask1 = (velo1.value < vmax) * (velo1.value > vmin)
        #vaxfix = np.interp(velo1[vmask1],velo2,np.arange(velo2.shape[0]))
        #kernelwidth = np.median(np.diff(vaxfix)) / 2.
        #kernel = np.exp(-np.linspace(-5,5,21)**2/(2.*kernelwidth**2))
        #K = np.zeros([21,21]); K[:,10] = kernel
        #pv2s = convolve(pv2,K,normalize_kernel=True,boundary='fill')
        #pv2resample = map_coordinates(pv2s, [np.outer(vaxfix,np.ones(pv1.shape[1])),
        #                                    np.outer(np.ones(vaxfix.size),np.arange(pv1.shape[1]))])

        #pv = ratio = pv1[vmask1,:]/pv2resample


        #pl.subplot(nxsp,nysp,6)
        #aspect = float(pv.shape[1])/float(vmask1.sum()) * cdelt / ((vmax-vmin)/float(vmask1.sum()))
        #pl.imshow(pv,
        #          extent=[0,pv.shape[1]*cdelt,vmin,vmax],
        #          aspect=aspect,
        #          cmap=cm,
        #          vmin=0,vmax=20,
        #          origin='lower')

        velo2 = velo2.to(u.km/u.s)
        vmask2 = (velo2.value < vmax) * (velo2.value > vmin)
        pv2toshow = pv2.data[vmask2[::index_order],:]
        cdelt2 = cdelt2.to(u.km/u.s).value

        # contours are H2CO 2-2
        #kernel = convolution.Gaussian2DKernel(width=0.5,x_size=11,y_size=11)
        #conv_pv2resample = convolve(pv2resample,kernel,normalize_kernel=True,boundary='fill')
        for kk in xrange(1,nxsp*nysp+1):
            pl.subplot(nxsp,nysp,kk)
            # kernel = make_kernel([11,11],kernelwidth=0.5)
            pl.contour(np.linspace(0,pv2toshow.shape[1]*cdelt2,pv2toshow.shape[1])[::-1],
                       np.linspace(vmin,vmax,pv2toshow.shape[0]),
                       pv2toshow,
                       #conv_pv2resample,
                       #levels=[0.02],
                       levels=[4,4.5],
                       colors=[color,'r'])
                       #levels=[1,2,5],
                       #colors=['y','c','r'])

        #pl.subplot(nxsp,nysp,7)
        #pl.imshow(pv2resample,
        #          extent=[0,pv.shape[1]*cdelt,vmin,vmax],
        #          aspect=aspect,
        #          cmap=cm,
        #          origin='lower')
        #pl.colorbar()

        pl.subplots_adjust(hspace=0,wspace=0)
        # for kk in [2,3,5,6]:
        #     pl.subplot(nxsp,nysp,kk).yaxis.set_ticklabels([])
        # for kk in [1,2,3]:
        #     pl.subplot(nxsp,nysp,kk).xaxis.set_ticklabels([])
        for kk in [2,3,5,6,8,9]:
            pl.subplot(nxsp,nysp,kk).yaxis.set_ticklabels([])
        for kk in [4,7]:
            pl.subplot(nxsp,nysp,kk).yaxis.set_ticks(pl.subplot(nxsp,nysp,kk).yaxis.get_ticklocs()[:-1])
        for kk in [1,2,3,4,5,6]:
            pl.subplot(nxsp,nysp,kk).xaxis.set_ticklabels([])
        for kk in xrange(1,nxsp*nysp+1):
            pl.subplot(nxsp,nysp,kk).axis([0,pv.shape[1]*cdeltv,vmin,vmax])
        for kk in [7,8,9]:
            pl.subplot(nxsp,nysp,kk).set_xlabel("Offset (\")",fontsize=18)
        for kk in [1,4,7]:
            pl.subplot(nxsp,nysp,kk).set_ylabel("V$_{LSR}$ (km/s)",fontsize=18)
        pl.savefig(os.path.join(datapath_w51,'pvdiagrams','W51_FilamentPVDiagrams_co32included_regrid_%s.png' % color),
                   bbox_inches='tight')
        pl.savefig(os.path.join(datapath_w51,'pvdiagrams','W51_FilamentPVDiagrams_co32included_regrid_%s.pdf' % color),
                   bbox_inches='tight')

    pl.show()
