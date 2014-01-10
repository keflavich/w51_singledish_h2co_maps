Arecibo Reduction Scripts
=========================

The Arecibo data reduction process requires Phil Perrilat's IDL environment.
Some instructions for its use can be found at
http://code.google.com/p/casaradio/wiki/aoIDL.  The `@masinit` command is
needed to initialize the Mock spectrometer reduction tools.

Getting Started
---------------
It is first necessary to set some environmental variables describing where
the data are stored.

Start by changing the data path in `aoinit.pro` to the appropriate location.
The raw data files must be stored in subdirectories with filenames like
YYYYMMDD (e.g., 20110910 for September 10, 2011).

Running the scripts
-------------------
The various `merge` scripts do all the hard reduction work.  Run them with the idl `@` command::

    @merge_spectra_W51_0910.pro

Creating the Maps
-----------------


