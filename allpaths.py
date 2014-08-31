import os
root = os.path.expanduser('~/work/')
datapath_w51 = os.path.join(root, 'w51/')
datapath = os.path.join(datapath_w51, 'h2co_singledish/')
source_root = os.path.join(root,'w51_singledish_maps/')
analysis_figurepath = os.path.join(source_root, 'figures/')
figurepath = os.path.join(source_root, 'tex/figures/')
datapath_spectra = os.path.join(datapath, 'spectralfits/')
regionpath = os.path.join(source_root, 'regions/')
modelpath = os.path.join(source_root, 'models/')
analysis_path = os.path.split(os.path.realpath(__file__))[0]
rcfilepath = os.path.join(analysis_path, 'plot_scripts/pubfiguresrc')
if not os.path.exists(rcfilepath):
    # probably means we're in a subdirectory
    rcfilepath = os.path.join(analysis_path, '../plot_scripts/pubfiguresrc')

if not os.path.exists(figurepath):
    # The parent directory *should* exist no matter what
    os.mkdir(figurepath)

def dpath(x, datapath=datapath):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)

def dpath2(x, datapath=datapath_w51):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)

def rpath(x, datapath=regionpath):
    """
    Shortcut function
    """
    return os.path.join(datapath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def afpath(x, figurepath=analysis_figurepath):
    return os.path.join(figurepath, x)

def mpath(x, modelpath=modelpath):
    return os.path.join(modelpath, x)

# File paths
cont6cm = dpath('W51_H2CO11_cube_supersampled_continuum.fits')
cont2cm = dpath('W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits')
h2co11subfn = dpath('W51_H2CO11_cube_supersampled_sub.fits')
h2co22subfn = dpath('W51_H2CO22_pyproc_cube_lores_supersampled_sub.fits')
h2co11taufn = dpath('W51_H2CO11_taucube_supersampled.fits')
h2co22taufn = dpath('W51_H2CO22_pyproc_taucube_lores_supersampled.fits')

model11tex = mpath('1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tex1.fits')
model11tau = mpath('1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau1.fits')
model22tex = mpath('1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tex2.fits')
model22tau = mpath('1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau2.fits')

h213co11subfn = dpath('W51_H213CO11_cube_supersampled_sub.fits')

h77asubfn = dpath('W51_h77a_pyproc_cube_supersampled_sub.fits')
h110asubfn = dpath('W51_h110alpha_cube_supersampled_sub.fits')
h110112asubfn = dpath('W51_Halpha_6cm_cube_supersampled_sub.fits')
