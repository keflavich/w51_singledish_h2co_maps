from astropy import units as u
from sdpy import makecube
from paths import dpath

cubename_lores_supersampled = dpath('W51_H2CO22_pyproc_cube_lores_supersampled')
linefreq = 14.488479e9
makecube.make_taucube(cubename_lores_supersampled,
                      cubename_lores_supersampled+"_continuum.fits",
                      etamb=0.886, linefreq=linefreq*u.Hz, tex=0)
makecube.make_flats(cubename_lores_supersampled.replace("cube","taucube"),
                    vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')


linefreq = 4.8296594e9
cubename_supersampled = dpath('W51_H2CO11_cube_supersampled')
makecube.make_taucube(cubename_supersampled,
                      continuum=cubename_supersampled+"_continuum.fits",
                      linefreq=linefreq, tex=0) # etamb accounted for already , etamb=0.51)
makecube.make_flats(cubename_supersampled.replace("cube","taucube"),
                    vrange=flat_vrange,noisevrange=[-50,-1],suffix='.fits')
