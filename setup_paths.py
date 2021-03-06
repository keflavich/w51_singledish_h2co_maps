import os, errno

def mkdir_p(path):
    """
    http://stackoverflow.com/a/600612/814354
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(os.path.realpath(path)):
            pass
        else:
            raise

from allpaths import (datapath, datapath_w51, analysis_figurepath, figurepath,
                      datapath_spectra, regionpath, modelpath)

for path in (datapath, datapath_w51, analysis_figurepath, figurepath,
             datapath_spectra, regionpath, modelpath):
    mkdir_p(path)
