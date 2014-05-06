
arecibo_datadirs = '/Users/adam/observations/arecibo/{}/'

"""
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0910_offspectra.fits
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0910_offspectra_nomask.fits
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0911_offspectra.fits
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0911_offspectra_nomask.fits
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0912_offspectra.fits
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0912_offspectra_nomask.fits
/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0915_offspectra.fits
"""

"""
Notes 5/2/2014:
    It looks like there are genuine problems with the Arecibo offs having
    primarily to do with the limitations of IDL.  In Python, we take percentiles
    along an axis to get, e.g., the 99th percentile in brightness at each pixel
    for our off position.  In IDL, we're limited to doing some sort of weird
    sort-on-totals.  It may be that getting around that limitation isn't so
    hard, but the data exploration time required to figure that out exceeds its
    usefulness right now.

    The problem is that there appears to be a *very slight* continuum excess
    nearly everywhere that is probably due to bad baseline interpolation.  The
    excess only becomes apparent (i.e. significant) when averaging over large
    areas.
"""
