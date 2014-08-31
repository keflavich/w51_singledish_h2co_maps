from astropy.utils.data import download_file
import paths
import FITS_tools
from astropy.io import fits

url = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/VGPS/MOS_049_contincluded.Tb'

dl = download_file(url, cache=True)

outfn = os.path.split(url)[1]

try:
    os.link(dl, paths.dpath2(outfn))
except OSError as ex:
    if ex.errno == 17 and ex.strerror == 'File exists':
        # All is good.  Continue.
        pass
    else:
        warnings.warn("Data is downloaded onto a different filesystem;"
                      " using symlinks instead of hardlinks")
        os.symlink(dl, paths.dpath2(outfn))

outhead = fits.getheader(paths.h2co11subfn)
hicubefn = paths.dpath('MOS_049_Tb_reprojsupersampled.fits')

FITS_tools.cube_regrid.regrid_fits_cube(outfn, outheader=outheader,
                                        outfilename=hicubefn, order=1)
