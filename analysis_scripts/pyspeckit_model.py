import numpy as np
import pyfits
import scipy
from pyspeckit.spectrum.models.formaldehyde import formaldehyde_vtau
#freq_dict = pyspeckit.models.formaldehyde.central_freq_dict
#freq_dict = pyspeckit.models.formaldehyde.hf_freq_dict
ckms = 2.99792458e5
import scipy
scipyOK=True

def formaldehyde_radex_tau_colden(xarr,
        density=4, column=13, xoff_v=0.0, width=1.0, 
        **kwargs):
    "docstring"
    return formaldehyde_radex_tau(xarr,
            (density,column,xoff_v,width),
            **kwargs)

def formaldehyde_radex_tau_coldentemopr(xarr,
        density=4, column=13, orthopara=0.1, temperature=20, xoff_v=0.0, width=1.0, 
        **kwargs):
    "docstring"
    return formaldehyde_radex_tau(xarr,
            (density,column,orthopara,temperature,xoff_v,width),
            **kwargs)

# Create a tau-only fit:
def formaldehyde_radex_tau(xarr, args, grid_vwidth=1.0,
        grid_vwidth_scale=False,
        taugrid=None,
        hdr=None,
        path_to_taugrid='',
        temperature_gridnumber=3,
        opr_gridnumber=10, # or None if there isn't one
        debug=False,
        verbose=False,
        return_hyperfine_components=False,
        **kwargs):
    """
    Use a grid of RADEX-computed models to make a model line spectrum

    The RADEX models have to be available somewhere.
    OR they can be passed as arrays.  If as arrays, the form should be:
    texgrid = ((minfreq1,maxfreq1,texgrid1),(minfreq2,maxfreq2,texgrid2))

    xarr must be a SpectroscopicAxis instance
    xoff_v, width are both in km/s

    grid_vwidth is the velocity assumed when computing the grid in km/s
        this is important because tau = modeltau / width (see, e.g., 
        Draine 2011 textbook pgs 219-230)
    grid_vwidth_scale is True or False: False for LVG, True for Sphere
    """
    if len(args) == 6:
        density,column,orthopara,temperature,xoff_v,width = args
        if np.any(np.array([np.isnan(x) for x in args])):
            raise ValueError("NAN argument")
        # assert orthopara > 0
        # orthopara = np.log10(orthopara)
    elif len(args) == 4:
        density,column,xoff_v,width = args
        orthopara,temperature = None,None
        # old defaults density=4, column=13, xoff_v=0.0, width=1.0, 
    else:
        raise ValueError("Not enough or Too Many arguments")

    if verbose:
        print "Parameters: dens=%f, column=%f, xoff=%f, width=%f" % (density, column, xoff_v, width)

    if taugrid is None:
        if path_to_taugrid=='':
            raise IOError("Must specify model grids to use.")
        else:
            taugrid = [pyfits.getdata(path_to_taugrid)]
            hdr = pyfits.getheader(path_to_taugrid)
            minfreq = (4.8,)
            maxfreq = (5.0,)
    elif hdr is not None:
        minfreq,maxfreq,taugrid = zip(*taugrid)
    else:
        raise Exception
    
    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('Hz', quiet=True)

    if orthopara is not None and temperature is not None:
        tau = taufunc_opr_temp(hdr, density, column, temperature, orthopara, taugrid)
    else:
        tau = taufunc(hdr, density, column, taugrid, temperature_gridnumber, opr_gridnumber=opr_gridnumber)

    if debug:
        print args
        print tau

    #tau_nu = []
    #for ii in xrange(len(tau)):
    #    linename = [k for k,v in freq_dict.items() if (v/1e9>minfreq[ii]) and (v/1e9<maxfreq[ii])][0]
    #    lines = freq_dict[linename]
    #    nuoff = xoff_v/ckms*lines
    #    nuwidth = np.abs(width/ckms*lines)
    #    tauspec = np.array(tau[ii] * np.exp(-(xarr.as_unit('Hz')+nuoff-freq_dict[linename])**2/(2.0*nuwidth**2)))
    #    tau_nu.append(tauspec)

    # let the hyperfine module determine the hyperfine components, and pass all of them here
    spec_components = [(formaldehyde_vtau(xarr.as_unit('Hz', quiet=True),
        tau=float(tau[ii]), xoff_v=xoff_v, width=width,
        return_tau=True, return_hyperfine_components=True, **kwargs) *
            (xarr.as_unit('GHz')>minfreq[ii]) *
            (xarr.as_unit('GHz')<maxfreq[ii])) 
                for ii in  xrange(len(tau))]

    # get an array of [n_lines, n_hyperfine, len(xarr)]
    if return_hyperfine_components:
        return np.array(spec_components).sum(axis=0)
    else:
        return np.sum(spec_components, axis=0).sum(axis=0)

    #if return_hyperfine_components:
    #    return np.array(tau_nu)
    #else:
    #    spec = np.sum(tau_nu, axis=0)

    #    if np.any(np.isnan(spec)) or np.any(np.isinf(spec)):
    #        raise ValueError("nan or infinite detected")
    #  
    #    return spec



def taufunc_opr_temp(hdr, density, column, temperature, orthopara, taugrid):

    densityarr = (np.arange(taugrid[0].shape[3])+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
    columnarr  = (np.arange(taugrid[0].shape[2])+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
    temparr  = (np.arange(taugrid[0].shape[1])+hdr['CRPIX3']-1)*hdr['CDELT3']+hdr['CRVAL3'] # temperature
    oprarr  = (np.arange(taugrid[0].shape[0])+hdr['CRPIX4']-1)*hdr['CDELT4']+hdr['CRVAL4'] # log ortho/para ratio

    gridval1 = np.interp(density,     densityarr,  np.arange(len(densityarr)))
    gridval2 = np.interp(column,      columnarr,   np.arange(len(columnarr)))
    gridval3 = np.interp(temperature, temparr,     np.arange(len(temparr)))
    gridval4 = np.interp(orthopara,   oprarr,      np.arange(len(oprarr)))
    if np.isnan(gridval1) or np.isnan(gridval2):
        raise ValueError("Invalid column/density")

    if scipyOK:
        slices = [slice(int(np.floor(gv)),int(np.floor(gv)+2)) for gv in (gridval4,gridval3,gridval2,gridval1)]
        tau = [scipy.ndimage.map_coordinates(tg[slices],
            np.array([[gridval4%1],[gridval3%1],[gridval2%1],[gridval1%1]],dtype='float'),
            order=1,prefilter=False) 
            for tg in taugrid]
        #tex = [scipy.ndimage.map_coordinates(tg[slices],np.array([[gridval4%1],[gridval3%1],[gridval2%1],[gridval1%1]]),order=1,prefilter=False) for tg in texgrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")

    return tau

def taufunc(hdr, density, column, taugrid, temperature_gridnumber, opr_gridnumber=None):

    densityarr = (np.arange(taugrid[0].shape[3])+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
    columnarr  = (np.arange(taugrid[0].shape[2])+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column

    gridval1 = np.interp(density,     densityarr,  np.arange(len(densityarr)))
    gridval2 = np.interp(column,      columnarr,   np.arange(len(columnarr)))
    if np.isnan(gridval1) or np.isnan(gridval2):
        raise ValueError("Invalid column/density")

    if scipyOK:
        slices = [temperature_gridnumber] + [slice(int(np.floor(gv)),int(np.floor(gv)+2)) for gv in (gridval2,gridval1)]
        if opr_gridnumber is not None:
            slices = [opr_gridnumber] + slices
        tau = [scipy.ndimage.map_coordinates(tg[slices].squeeze(),
            np.array([[gridval2%1],[gridval1%1]],dtype='float'),order=1,prefilter=False) 
            for tg in taugrid]
        #tex = [scipy.ndimage.map_coordinates(tg[slices],np.array([[gridval4%1],[gridval3%1],[gridval2%1],[gridval1%1]]),order=1,prefilter=False) for tg in texgrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")

    return tau
