import numpy as np

vmin = 50
vmax = 80
velocityrange = [-50,150] # match GBT exactly
cd3 = 1.0 # match GBT exactly
naxis3 = int((velocityrange[1]-velocityrange[0]) / cd3) + 1
crval3 = 50.0
vels = crval3+cd3*(np.arange(naxis3)+1-naxis3/2-1)

glon, glat = 49.209553, -0.277137
naxis1 = 308
naxis2 = 205
