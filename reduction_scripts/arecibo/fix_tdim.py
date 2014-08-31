from astropy.io import fits

def fix_TDIM_in_header(fn):
    """
    Apparently TDIM that's stored in headers started being checked against the
    data in a recent version of astropy.io.fits, and that leads to files being
    unreadable if they have TDIM parameters.  I apparently discovered this
    during the September 2013 observing run on NGC 1333 and made the
    appropriate changes in the Arecibo reduction files but never re-reduced the
    W51 data
    """
    f = fits.open(fn)

    fixed=False
    for HDU in f:
        for kw in HDU.header.keys():
            if any([x in kw for x in ('TDIM',)]):
                fixed = True
                del HDU.header[kw]

    if fixed:
        outf = fn.replace(".fits","_fixed.fits")
        f.writeto(outf,output_verify='fix',clobber='fixed' in outf)
        print "Fixed file ",fn," -> ",outf
    else:
        outf = fn
        print "Did not fix file ",fn

    return outf


