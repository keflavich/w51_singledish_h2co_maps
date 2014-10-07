import numpy as np
import astrodendro
import FITS_tools
from astropy import units as u
from paths import datapath
from astropy.io import fits
from allpaths import cont6cm as oneonefn
from allpaths import cont2cm as twotwofn

oneoned = (fits.getdata(oneonefn) *
           u.K).to(u.Jy, equivalencies=u.brightness_temperature(aobeamarea,
                                                                h2co11freq))
twotwod = fits.getdata(twotwofn)
oneoneh = fits.getheader(oneonefn)
twotwoh = fits.getheader(twotwofn)

dendone = astrodendro.dendrogram.Dendrogram.compute(oneoned, 0.1*u.Jy, 0.05*u.Jy, 6)

metadata = {}
#metadata['data_unit'] = u.K
metadata['spatial_scale'] = FITS_tools.header_tools.header_to_platescale(oneoneh) * u.deg
metadata['beam_major'] = aobeam
metadata['beam_minor'] = aobeam

cat = astrodendro.pp_catalog(dendone, metadata)
