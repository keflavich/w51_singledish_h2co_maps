import astrodendro
from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy import constants
import pylab as pl
pl.rc('font',size=24)



datapath = '/Users/adam/work/w51/'

def col_conversion_factor(TK=20):
    return 2.19e22 * (np.exp(13.01/TK - 1))

area_per_pixel = (7.2*u.arcsec*5100*u.pc/(u.radian))**2
area_per_pixel_cm = area_per_pixel.to(u.cm**2).value
area_per_pixel_pc = area_per_pixel.to(u.pc**2).value


v2im = fits.getdata(datapath+'v2.0_ds2_l050_13pca_map20.fits')
w51coldens = v2im[80:160,2540:2610] * col_conversion_factor()

dend = astrodendro.Dendrogram.compute(w51coldens,min_value=3e21,min_npix=4,verbose=True)
dend.viewer()

def get(obj, property=None):
    func = getattr(obj, 'get_' + property)
    return func()

def get_mostest_children(item, props=('npix','f_sum'), mostest='f_sum'):
    """
    Return certain properties of all children selecting the mostest
    of something (e.g., brightest)
    """
    item.get_dens = lambda: item.values().sum()/item.get_npix()**(1.5)
    item.get_f_sum = lambda: item.values().sum()
    if item.is_leaf:
        d = dict([(p,[get(item,p)]) for p in props])
        d['branch'] = [item]
        return d
    else:
        brightest = item.children[0]
        brightest.get_dens = lambda: brightest.values().sum()/item.get_npix()**(1.5)
        brightest.get_f_sum = lambda: brightest.values().sum()
        for child in item.children[1:]:
            child.get_dens = lambda: child.values().sum()/item.get_npix()**(1.5)
            child.get_f_sum = lambda: child.values().sum()
            if get(child,mostest) > get(brightest,mostest):
                brightest = child
        brightest_props = get_mostest_children(brightest)
        d = dict([(p,[get(item,p)] + brightest_props[p])
                     for p in props])
        d['branch'] = [item] + brightest_props['branch']
        return d


w51main = dend.structure_at((42,32))
w51irs2 = dend.structure_at((51,32))

cloud = w51main
while cloud.parent:
    cloud = cloud.parent

w51irs2tree = []
c = w51irs2
while c.parent:
    w51irs2tree.append(c)
    c = c.parent

cloud_kids_ALL = get_mostest_children(cloud, mostest='dens')
cloud_kids = cloud_kids_ALL
#cloud_kids = {k:x[:126] for k,x in cloud_kids_ALL.iteritems()}

metadata = {'beam_major':33*u.arcsec,'beam_minor':33*u.arcsec}
cloud_kids['stats'] = [astrodendro.analysis.PPStatistic(i, metadata=metadata) for i in cloud_kids['branch']]
cloud_kids['maj'] = np.array([s.beam_major.value for s in cloud_kids['stats']])
cloud_kids['min'] = np.array([s.beam_minor.value for s in cloud_kids['stats']])

w51irs2dict = {}
w51irs2dict['stats'] = [astrodendro.analysis.PPStatistic(i, metadata=metadata) for i in w51irs2tree]
w51irs2dict['maj'] = np.array([s.beam_major.value for s in w51irs2dict['stats']])
w51irs2dict['min'] = np.array([s.beam_minor.value for s in w51irs2dict['stats']])
w51irs2dict['npix'] = np.array([s.area_exact.value for s in w51irs2dict['stats']])
w51irs2dict['f_sum'] = np.array([c.get_f_sum() for c in w51irs2tree])


cloud_effective_radius_pc = (np.array(cloud_kids['npix'],dtype='float64')*area_per_pixel_pc/np.pi)**0.5 
cloud_mass = np.array(cloud_kids['f_sum'],dtype='float64')*area_per_pixel_cm*1.67e-24*2.72/2e33
cloud_volume_cm3 = cloud_kids['maj'] * cloud_kids['min']**2 * (2*np.pi)**(1.5) * (area_per_pixel_cm)**1.5

w51irs2_effective_radius_pc = (np.array(w51irs2dict['npix'],dtype='float64')*area_per_pixel_pc/np.pi)**0.5 
w51irs2_mass = np.array(w51irs2dict['f_sum'],dtype='float64')*area_per_pixel_cm*1.67e-24*2.72/2e33
w51irs2_volume_cm3 = w51irs2dict['maj'] * w51irs2dict['min']**2 * (2*np.pi)**(1.5) * (area_per_pixel_cm)**1.5


yscale=1e-3

pl.figure(4)
pl.clf()
pl.plot(cloud_effective_radius_pc,cloud_mass*yscale,'s-', color='k', linewidth=3, alpha=0.3)
pl.plot(w51irs2_effective_radius_pc,w51irs2_mass*yscale,'s-', color='k', linewidth=3, alpha=0.3)

mh2 = (constants.m_p * 2.72 / u.solMass).decompose().value 
msun = u.solMass.to('g')
r_pc = np.linspace(0,8,1000)
pl.plot(r_pc, 4/3.*np.pi*(r_pc*3.08e18)**3 * 10   * mh2*yscale, 'g--', linewidth=3, alpha=0.5, label='$n=10  $')
pl.plot(r_pc, 4/3.*np.pi*(r_pc*3.08e18)**3 * 100  * mh2*yscale, 'r--', linewidth=3, alpha=0.5, label='$n=100 $')
pl.plot(r_pc, 4/3.*np.pi*(r_pc*3.08e18)**3 * 1000 * mh2*yscale, 'c--', linewidth=3, alpha=0.5, label='$n=1000$')
pl.plot(r_pc, 4/3.*np.pi*(r_pc*3.08e18)**3 * 1e4  * mh2*yscale, 'b--', linewidth=3, alpha=0.5, label='$n=10^4$')
pl.plot(r_pc, 4/3.*np.pi*(r_pc*3.08e18)**3 * 1e5  * mh2*yscale, 'm--', linewidth=3, alpha=0.5, label='$n=10^5$')
pl.gca().axis([0,6,0,1e5*yscale])
pl.legend(loc='best',prop=dict(size=20))
pl.xlabel("Effective Radius (pc)")
pl.ylabel("Mass (1000 M$_\odot$)")

pl.figure(5)
pl.clf()
pc = u.pc.to('cm')
n_sphere = cloud_mass/mh2/(4/3.*np.pi*(cloud_effective_radius_pc*pc)**3)
#n_pancake = cloud_mass/mh2/(4/3.*np.pi*(cloud_effective_radius_pc*pc)**2*(2*pc))
#cloud_n_ellipsoid = cloud_mass/mh2/cloud_volume_cm3
pl.plot(cloud_effective_radius_pc, n_sphere, 's-', color='r', label='Cloud $n(H_2)$ (sphere)')
#pl.plot(cloud_effective_radius_pc, n_pancake,'o-', color='r', label='Cloud $n(H_2)$ (pancake)')
#pl.plot(cloud_effective_radius_pc,cloud_n_ellipsoid,'^:', color='r', label='Cloud $n(H_2)$ (moments)')
#pl.fill_between(cloud_effective_radius_pc,n_sphere,n_pancake,color='r', zorder=-2, alpha=0.1)

w51irs2_n_sphere = w51irs2_mass/mh2/(4/3.*np.pi*(w51irs2_effective_radius_pc*pc)**3)
#w51irs2_n_pancake = w51irs2_mass/mh2/(4/3.*np.pi*(w51irs2_effective_radius_pc*pc)**2*(2*pc))
#w51irs2_n_ellipsoid = w51irs2_mass/mh2/w51irs2_volume_cm3
pl.plot(w51irs2_effective_radius_pc,w51irs2_n_sphere, 's:', color='g', label='w51irs2 $n(H_2)$ (sphere)' )
#pl.plot(w51irs2_effective_radius_pc,w51irs2_n_pancake,'o:', color='g', label='w51irs2 $n(H_2)$ (pancake)')
#pl.plot(w51irs2_effective_radius_pc,w51irs2_n_ellipsoid,'^:', color='g', label='w51irs2 $n(H_2)$ (moments)')
#pl.fill_between(w51irs2_effective_radius_pc,w51irs2_n_sphere,w51irs2_n_pancake,color='g', zorder=-1, alpha=0.1)

#pl.hlines([8,15,150],*pl.gca().get_xlim(), linestyle='--', color='k')
#pl.hlines([10,200],*pl.gca().get_xlim(), linestyle='--', color='k')
pl.axis([0,6,1e3,5e5])
pl.gca().set_yscale('log')

pl.legend(loc='best',prop=dict(size=20))
pl.xlabel("Effective Radius (pc)")
pl.ylabel("Density ($n(H_2)$)")

