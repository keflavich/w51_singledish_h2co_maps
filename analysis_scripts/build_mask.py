from paths import h2co11taufn, h2co22taufn
from spectral_cube import BooleanArrayMask
from astropy import wcs

### copied from tau_ratio_cube
# These are not used for fitting, only for masking!
# The absorption lines themselves (from cube1,cube2 above) are used for the
# fitting
h2co11 = fits.getdata(h2co11taufn)
h2co22 = fits.getdata(h2co22taufn)

noise11 = h2co11[:50,:,:].std(axis=0)
noise22 = h2co22[:50,:,:].std(axis=0)

sn11 = h2co11/noise11
sn22 = h2co22/noise22

finite = (np.isfinite(sn11) & np.isfinite(sn22))
includemask = ((((sn11 > 2) & (sn22 > 2))
                | (sn11 > 4)
                | (sn22 > 4))
               & finite)

log.info("SN22 > 2 & SN22 > 2: {0}".format(np.count_nonzero((sn11 > 2) & (sn22 > 2))))
log.info("SN11 > 4 & SN22 < 2: {0}".format(np.count_nonzero((sn11 > 4) & (sn22 < 2))))
log.info("SN22 > 4 & SN11 < 2: {0}".format(np.count_nonzero((sn22 > 4) & (sn11 < 2))))

ratio = h2co11/h2co22
ratio[True-includemask] = np.nan

filt = np.ones([3,3,3],dtype='bool')
filt[1,1,1] = 0
nneighbors = convolve(np.isfinite(ratio), filt)
includemask[~np.isfinite(nneighbors)] = False
includemask[(nneighbors<7)] = False
includemask[(nneighbors>=10)] = True
nneighbors2 = convolve(includemask, filt)
includemask[(nneighbors2>=5)] = True
includemask[~finite] = False

maskwcs = wcs.WCS(fits.getheader(h2co11taufn))
cubemask = BooleanArrayMask(incluemask, wcs=maskwcs)

# End masking
###
