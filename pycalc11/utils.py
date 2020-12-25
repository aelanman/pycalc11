
import numpy as np
import warnings
from datetime import datetime
from astropy.time import Time, TimeDelta
from astropy.utils import data
from astropy.utils import iers
import astropy.coordinates as ac
from astropy import units as un
from astropy.constants import c as speed_of_light
from multiprocessing import Process
import pylab as pl



iers_tab = iers.earth_orientation_table.get()

def get_leap_seconds(tobj):
    # Find current TAI - UTC for a given time.
    lsec_table = iers.LeapSeconds.auto_open().as_array()
    for ti, (yr, mo, tai_utc) in enumerate(lsec_table):
        dtobj = datetime(yr, mo, 1, 0, 0, 0)
        if dtobj > tobj.datetime:
            ti -= 1
            break

    if tobj.datetime.year < 1960:
        return 0.0

    (_, _, tai_utc) = lsec_table[ti]
    return tai_utc
