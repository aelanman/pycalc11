import os
import sys
import numpy as np
from astropy.utils.data import download_file

from .data import DATA_PATH
from . import calc11


# Get JPL ephemeris data
de421_url = f"https://svn.atnf.csiro.au/difx/applications/difxcalc11/trunk/data/DE421_{sys.byteorder}_Endian"
de421_path = download_file(de421_url, cache=True)
calc11.datafiles.jpl_de421 = de421_path.ljust(128) 


def _format(fname):
    pth = os.path.join(DATA_PATH, fname).ljust(128)
    return pth

# Set data file paths
calc11.datafiles.a_tilts = _format("tilt.dat")
calc11.datafiles.oc_file = _format("ocean_load.coef")
calc11.datafiles.optl_file = _format("ocean_pole_tide.coef")
calc11.datafiles.dfleap = _format("ut1ls.dat")

from .interface import Calc
