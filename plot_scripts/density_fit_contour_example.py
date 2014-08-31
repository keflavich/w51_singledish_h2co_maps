from load_pyspeckit_cubes import (both, T, F, cont11, cont22, h2co11filename,
                                  cube1, cube2, texgrid1,  taugrid1,  texgrid2,
                                  taugrid2,  hdr, fit_a_pixel, densityarr,
                                  columnarr, temparr, oprarr, )
from h2co_modeling.grid_fitter import grid_2p_getmatch,grid_getmatch
import pylab as pl
import numpy as np


# FITS -> python convention
xx,yy,zz = (87-1, 76-1, 109-1)
xx,yy,zz = (86-1, 75-1, 107-1)
xx,yy,zz = (82-1, 73-1, 107-1)
xx,yy,zz = (111-1, 110-1, 98-1) # lowdens / highstd?  EXTREME degeneracy!
xx,yy,zz = (90-1, 110-1, 110-1) # low dens, high *column* std but low *dens* std

cont1 = cont11[yy,xx]
cont2 = cont22[yy,xx]
tline1 = cube1.cube[zz,yy,xx]
tline2 = cube2.cube[zz,yy,xx]
etline1 = cube1.errorcube[zz,yy,xx]
etline2 = cube2.errorcube[zz,yy,xx]

pargrid1 = (cont1*np.exp(-taugrid1) + (1-np.exp(-taugrid1))*texgrid1)
pargrid2 = (cont2*np.exp(-taugrid2) + (1-np.exp(-taugrid2))*texgrid2)
#spec = (1.0-np.exp(-np.array(tau_nu_cumul)))*(Tex-Tbackground)

match, indbest, chi2 = grid_2p_getmatch(tline1+cont1, etline1, pargrid1,
                                        tline2+cont2, etline2, pargrid2)
match11, indbest11, chi2_11 = grid_getmatch(tline1+cont1, etline1, pargrid1)
match22, indbest22, chi2_22 = grid_getmatch(tline2+cont2, etline2, pargrid2)

bw,bz,by,bx = np.unravel_index(indbest, pargrid1.shape)

pl.jet()

pl.clf()
extent = [densityarr.flat[0],densityarr.flat[-1],
          columnarr.flat[0],columnarr.flat[-1]]
pl.imshow(pargrid1[bw,bz,:,:], alpha=0.5, extent=extent)
pl.contour(densityarr[bw,bz,:,:],columnarr[bw,bz,:,:],
           pargrid2[bw,bz,:,:], levels=np.linspace(tline2-20*etline2+cont2,
                                                   tline2+20*etline2+cont2,
                                                   10))

H,xe,ye = np.histogram2d(densityarr[match], columnarr[match],
                         bins=[densityarr[0,0,0,::3],
                               columnarr[0,0,::3,0]])
xc = xe[:-1] + 0.5 * (xe[1:] - xe[:-1])
yc = ye[:-1] + 0.5 * (ye[1:] - ye[:-1])
X,Y = np.meshgrid(xc,yc)
pl.contourf(X,Y,H.T, levels=[0.5,1000], colors='g', alpha=0.5)

if False:
    colors = 'rb'

    for ii,match in enumerate((match11,match22)):
        mask = np.where(match, match, np.nan)
        #pl.imshow(pargrid1[bw,bz,:,:]*mask[bw,bz,:,:], alpha=1.0, extent=extent)
        #pl.plot(densityarr[match], columnarr[match], 'x')
        H,xe,ye = np.histogram2d(densityarr[match], columnarr[match],
                                 bins=[densityarr[0,0,0,::3],
                                       columnarr[0,0,::3,0]])
        xc = xe[:-1] + 0.5 * (xe[1:] - xe[:-1])
        yc = ye[:-1] + 0.5 * (ye[1:] - ye[:-1])
        X,Y = np.meshgrid(xc,yc)
        pl.contourf(X,Y,H.T, levels=[0.5], colors=colors[ii])
pl.axis(extent)
pl.show()
