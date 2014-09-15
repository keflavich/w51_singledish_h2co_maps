W51 Single Dish maps
====================

Repository for a paper presently in prep.  Let me know if you want access to
the data or want to join in on the paper.

Instructions
------------

There is a fairly large suite of dependencies for this program:

 1. Go here: http://thedata.harvard.edu/dvn/dv/W51_H2CO and download the data.
    It is unfortunately impossible to provide a direct download link, as the
    dataverse obscures the URLs.

    Required Files:

     W51_H2CO11_cube_supersampled_sub.fits
     W51_H2CO11_cube_supersampled_continuum.fits
     W51_H2CO22_pyproc_cube_lores_supersampled_sub.fits
     W51_H2CO22_pyproc_cube_lores_supersampled_continuum.fits
     h77a
     h110a
     h112a

 2. Install python, if not yet available.  `Conda
    <http://continuum.io/downloads>`_ is very efficient for this purpose
 2. Download & install sdpy, FITS_tools, spectral_cube, pyspeckit, astropy...


pip install https://github.com/astropy/astropy/archive/master.zip
pip install https://github.com/astropy/astroquery/archive/master.zip [CURRENTLY REQUIRES AN UNRELEASED VERSION]
pip install https://github.com/aplpy/aplpy/archive/master.zip
pip install https://github.com/keflavich/FITS_tools/archive/master.zip
pip install https://github.com/keflavich/sdpy/archive/master.zip
pip install https://github.com/pyspeckit/pyspeckit/archive/master.zip
pip install https://github.com/keflavich/image_tools/archive/master.zip
pip install https://github.com/keflavich/image_registration/archive/master.zip
pip install https://github.com/radio-astro-tools/spectral-cube/archive/master.zip   [CURRENTLY REQUIRES AN UNRELEASED VERSION]
pip install https://github.com/keflavich/h2co_modeling/archive/master.zip
pip install https://github.com/keflavich/pyregion/archive/master.zip

 3. Set up the paths appropriately in `allpaths.py <allpaths.py>`_

    .. note:: TODO: Add simple example here

 4. Run the `analysis scripts <analysis_scripts/run_all.py>`_
 5. Run the `plot scripts <plot_scripts/run_all.py>`_

.. image:: https://zenodo.org/badge/3821/keflavich/w51_singledish_h2co_maps.png   :target: http://dx.doi.org/10.5281/zenodo.11737


..
    http://thedata.harvard.edu/dvn/dv/W51_H2CO/FileDownload/?fileId=2387750&xff=0&versionNumber=1
    2387749
