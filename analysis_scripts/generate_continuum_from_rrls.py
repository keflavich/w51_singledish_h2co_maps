"""
May 27:
    Attempt to generate a cleaner continuum map by extrapolating from the RRL map
    

While this works, in the end I won't be using them because the residuals are
structured
"""
from paths import dpath
from agpy import fit_a_line
np.seterr(all='ignore')

rrl6cm = fits.getdata(dpath('W51_Halpha_6cm_integrated_supersampled.fits'))
rrl2cm = fits.getdata(dpath('W51_h77a_pyproc_integrated_supersampled.fits'))

cont11filename = dpath('W51_h110alpha_cube_supersampled_continuum.fits')
cont22filename = dpath('W51_h77a_pyproc_cube_supersampled_continuum.fits')
cont11 = fits.getdata(cont11filename)
cont22 = fits.getdata(cont22filename)

rrlerr6cm = 0.0415
rrlerr2cm = 0.0065
rrl6cmmask = rrl6cm > rrlerr6cm * 2
rrl2cmmask = rrl2cm > rrlerr2cm * 2

offset2cm, offset6cm = 0,0

rrldata,contdata,rrlmask = (rrl2cm,cont22,rrl2cmmask)
rrlmask *= (contdata > 0.1)
nok = np.count_nonzero(rrlmask)
factor2cm,offset2cm = fit_a_line.total_least_squares(rrldata[rrlmask],
                                                     contdata[rrlmask],
                                                     data1err=np.ones(nok)*rrlerr2cm,
                                                     data2err=np.ones(nok)*0.05,
                                                     print_results=True,
                                                     intercept=True)
print "2cm: ",offset2cm,factor2cm
artificial_cont_2cm = rrldata * factor2cm + offset2cm


rrldata,contdata,rrlmask = (rrl6cm,cont11,rrl6cmmask)
rrlmask *= (contdata > 0.1)
nok = np.count_nonzero(rrlmask)
factor6cm, offset6cm = fit_a_line.total_least_squares(rrldata[rrlmask],
                                                      contdata[rrlmask],
                                                      data1err=np.ones(nok)*rrlerr6cm,
                                                      data2err=np.ones(nok)*0.05,
                                                      print_results=True,
                                                      intercept=True)
print "6cm: ",offset6cm,factor6cm
artificial_cont_6cm = rrldata * factor6cm + offset6cm

import pylab as pl
pl.figure(1)
pl.clf()
pl.subplot(2,1,1)
pl.imshow(cont11 - artificial_cont_6cm)
pl.subplot(2,1,2)
pl.plot(cont11[rrl6cmmask], artificial_cont_6cm[rrl6cmmask], ',')
pl.plot(cont11[rrl6cmmask], cont11[rrl6cmmask])
#pl.plot(cont11[rrl6cmmask], rrl6cm[rrl6cmmask]*factor6cm + offset6cm, '-')

pl.figure(2)
pl.clf()
pl.subplot(2,1,1)
pl.imshow(cont22 - artificial_cont_2cm)
pl.subplot(2,1,2)
pl.plot(cont22[rrl2cmmask], artificial_cont_2cm[rrl2cmmask], ',')
pl.plot(cont22[rrl2cmmask], cont22[rrl2cmmask])
#pl.plot(rrl2cm[rrl2cmmask], rrl2cm[rrl2cmmask]*factor2cm + offset2cm, '-')

pl.show()
