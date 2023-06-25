import os
import sys
import numpy as np
from .data import DATA_PATH
import ctypes
from . import calc11

def _format(fname):
    pth = os.path.join(DATA_PATH, fname).ljust(128)
    return pth

# Set data file paths
calc11.datafiles.jpl_de421 = _format(f"DE421_{sys.byteorder}_Endian")
calc11.datafiles.a_tilts = _format("tilt.dat")
calc11.datafiles.oc_file = _format("ocean_load.coef")
calc11.datafiles.optl_file = _format("ocean_pole_tide.coef")
calc11.datafiles.dfleap = _format("ut1ls.dat")

from .interface import Calc
