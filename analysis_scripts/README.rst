Analysis Scripts and their descriptions
=======================================

In order to run these, you should first modify paths.py

Order matters below!

Setup Scripts
-------------

These are meant to be imported by other scripts, not used directly

 * common_constants.py
 * paths.py
 * load_pyspeckit_cubes.py
 
   Loads appropriate data cubes with `pyspeckit <pyspeckit.bitbucket.org>`_

Standalone Scripts
------------------

These scripts can be run in any order.

 * contcompare_dc_2cm.py

   Compare the 2cm continuum to `Glen Langston's GPA
   <http://www.gb.nrao.edu/~glangsto/GPA/>`_

 * contcompare_dc_6cm.py

   Compare the 6cm continuum to `Urumqi 25m maps <http://zmtt.bao.ac.cn/6cm/>`_
   (requires the file w51.iuq.fits, which is presently private)

 * digging_in_to_G49.27-0.34.py

   A careful examination of a particular HII region

 * examine_one_spectrum.py

   A deep examination of an apparently dense region.  This region is subject to
   worse systematic error than much of the rest of the data.

 * generate_continuum_from_rrls.py

   An abandoned attempt to produce a "cleaner" continuum using the RRLs.
   The fact that there is some structure in the residual maps from this
   approach, hinting that synchrotron emission may indeed be significant, and
   the lower S/N in the RRLs both led me to leave this out fo further analysis.

 * get_simbad_regions.py

   Use `astroquery <astroquery.readthedocs.org>`_ to find all HII regions and
   save the result as a `ds9 <ds9.si.edu>`_ .reg file.

 * integrated_h2co_properties.py

   Compute various properties (optical depth integral, ratio) of the H2CO maps

 * luminosity.py

   Determine W51's total luminosity from HiGal SED fits.  WARNING: the SED fits
   interpolate across NaNs!

 * merge_apertures.py

   Support script.  Merges region files.

 * merge_mips_msx.py

   Support script.  Fill in missing MIPS data with MSX data.  Not used, but
   could be useful to others.

 * pyspeckit_RRL_cube_fit.py

   Single-component fits to each line of sight in the RRL cubes

 * pyspeckit_cube_fit_justlineprops.py

   Multi-component H2CO 1-1 fits to each line of sight.  Used to create
   velocity maps.

 * pyspeckit_cube_fit_textau.py

   Single-component formaldehyde fitting of physical parameters to both lines
   simultaneously.  Requires load_pyspeckit_cubes.py.

Figure-Making Scripts
---------------------

Most figure-making scripts are in the ``plot_scripts`` directory, but some live
here also.

 * pvdiagrams.py
 * pvfigure_h2co13co.py

   Scripts to make PV diagrams with various cuts.


Ordered Scripts
---------------

 * density_fit.py

   Fits densities using a grid and chi-squared matching.  Creates density cubes
   and density maps.

 * extract_spectra.py

   OBSOLETE?
   Use aperture extraction to create individual spectra for a number of regions
   of interest.

 * pyspeckit_individual_spectra.py
 * pyspeckit_model.py
 * radial_profiles.py
 * radial_velocity_dispersion.py
 * regrid_higal.py
 * rrl_analysis.py
 * rrl_average.py
 * set_headers.py
 * tau_ratio_cube.py
 * vpeak.py


 * co_intmaps.py

   Analysis of the CO maps by integrating over regions masked by various
   properties of the H2CO, e.g. density.

Throwaway Scripts
-----------------

TBDeleted

 * continumm_dendro.py

   An attempt to dendrogram the W51 column density map.  Not used in this
   project.

 * pyspeckit_cube_fit_absorption.py

   A failed attempt to directly fit the absorption lines.  The idea behind this
   attempt was moved into another file.

 * pyspeckit_cube_fit.py

   The original fitter; other tools branched from this one.
