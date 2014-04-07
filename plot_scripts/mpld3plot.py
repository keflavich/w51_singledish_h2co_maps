import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import mpld3
from mpld3 import plugins, utils
import jinja2
import json
import glob

import pyspeckit
# TODO: replace with library
import sys
sys.path.append('/Users/adam/work/w51_singledish_maps/analysis_scripts/')
import pyspeckit_individual_spectra


class LinkedView(plugins.PluginBase):
    """A simple plugin showing how multiple axes can be linked"""

    JAVASCRIPT = """
    var LinkedViewPlugin = function(fig, prop){
      this.fig = fig;
      this.prop = mpld3.process_props(this, prop, {},
                                      ["idpts", "idline", "data"]);
    }

    LinkedViewPlugin.prototype.draw = function(){
      var pts = mpld3.get_element(this.prop.idpts);
      var line = mpld3.get_element(this.prop.idline);
      var data = this.prop.data;

      function mouseover(d, i){
        line.data = data[i];
        line.elements().transition()
            .attr("d", line.datafunc(line.data))
            .style("stroke", this.style.stroke);
      }
      pts.elements().on("mouseover", mouseover);
    };

    mpld3.register_plugin("linkedview", LinkedViewPlugin);
    """

    def __init__(self, points, line, linedata):
        if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None

        self.dict_ = {"type": "linkedview",
                      "idpts": utils.get_id(points, suffix),
                      "idline": utils.get_id(line),
                      "data": linedata}


#spectra = pyspeckit.Spectra(pyspeckit_individual_spectra.do_all_h2co())
spectra = pyspeckit.Spectra([pyspeckit.Spectrum(fn) for fn in
                             glob.glob('/Users/adam/work/h2co/maps/W51/spectralfits/*.fits')
                             if ('CO' not in fn and 'RRL' not in fn and 'HI'
                                 not in fn)])

for sp in spectra:
    sp.xarr.refX = 14.488
    sp.xarr.refX_units = 'GHz'
    sp.xarr.convert_to_unit('km/s')

P = np.array([sp.header['CRVAL2'] for sp in spectra])
A = np.array([sp.header['CRVAL3'] for sp in spectra])
s = np.array([sp.header['APRADIUS'] for sp in spectra])

# scatter periods and amplitudes
#np.random.seed(0)
#P = 0.2 + np.random.random(size=20)
#A = np.random.random(size=20)
#x = np.linspace(0, 10, 100)
#data = np.array([[x, Ai * np.sin(x / Pi)]
#                 for (Ai, Pi) in zip(A, P)])

fig = plt.figure()
#fig.clf()
#fig,ax = plt.subplots(2)

dowcs=False

from astropy.io import fits
fn = '/Users/adam/work/w51/W51_irac4.fits'
img = fits.getdata(fn)
from astropy.wcs import WCS
hdr = fits.getheader(fn)
wcs = WCS(hdr)
pixll,pixur = [(1,1),(hdr['NAXIS1'],hdr['NAXIS2'])]
pixll,pixur = np.array([(480,720),(1412,1211)])+1

if dowcs:
    from wcsaxes import WCSAxes

    ax = [None,None]
    ax[1] = WCSAxes(fig, [0.1, 0.1, 0.8, 0.4], wcs=wcs)
    fig.add_axes(ax[1])
    ax[1].imshow(img, cmap=plt.cm.gist_heat, origin='lower',
                 norm=matplotlib.colors.LogNorm(), vmin=0.1, vmax=1000,
                 interpolation='nearest')
    ax[0] = matplotlib.axes.Axes(fig, [0.1,0.5,0.8,0.4])
    fig.add_axes(ax[0])

    tr_gal = ax[1].get_transform('galactic')
    points = ax[1].scatter(P, A,
                           s=s*3600*4, alpha=0.5,
                           transform=tr_gal,
                           edgecolor='cyan',
                           facecolor='none',
                           linewidths=1)
else:
    ax = (matplotlib.axes.Axes(fig, [0.1,0.65,0.8,0.3]),
          matplotlib.axes.Axes(fig, [0.1,0.05,0.8,0.55]))
    for a in ax:
        fig.add_axes(a)
    ll = wcs.wcs_pix2world([pixll],1)[0]
    ur = wcs.wcs_pix2world([pixur],1)[0]
    ax[1].imshow(img[pixll[1]-1:pixur[1]-1,
                     pixur[0]-1:pixll[0]-1:-1],
                 extent=[ur[0],ll[0],ll[1],ur[1],],
                 interpolation='nearest', vmin=0.1, cmap=plt.cm.gist_heat,
                 vmax=1000, norm=matplotlib.colors.LogNorm())
    points = ax[1].scatter(P, A,
                           s=s*3600*4, alpha=0.5,
                           edgecolor='black',
                           facecolor='cyan',
                           linewidths=1)
    ax[1].set_xlabel('GLON')
    ax[1].set_ylabel('GLAT')

# create the line object
#lines = ax[0].plot(x, 0 * x, '-w', lw=3, alpha=0.5)
#spectra.ploteach(axis=ax[0],clear=False,xmin=40,xmax=80)
l1 = False

if l1:
    xarr1 = spectra[0].xarr[201:]
    yarr1 = np.array([sp.data[201:]/sp.data[201:].min() for sp in spectra])
    lines = ax[0].plot(xarr1, -10+0*yarr1.T, lw=1, alpha=1.0, color='w')
    data = np.array([(xarr1,y) for y in yarr1])
    yarr = yarr1
else:
    xarr2 = spectra[0].xarr[:201]
    xarr2.convert_to_unit('GHz')
    xarr2.refX = 4.82966
    xarr2.convert_to_unit('km/s')
    yarr2 = np.array([sp.data[:201]/sp.data[:201].min() for sp in spectra])
    lines = ax[0].plot(xarr2, -10+0*yarr2.T, lw=1, alpha=1.0, color='w')
    data = np.array([(xarr2,y) for y in yarr2])
    yarr = yarr2

#lines = ax[0].lines
ax[0].set_xlim(30,80)
ax[0].set_ylim(-0.1,yarr.max())

ax[0].set_title("Hover over points to see lines")

#data = np.array([(np.array(sp.xarr),sp.data) for sp in spectra])

#plt.show()

# transpose line data and add plugin
linedata = data.transpose(0, 2, 1).tolist()
plugins.connect(fig, LinkedView(points, lines[0], linedata))

mpld3.show()
