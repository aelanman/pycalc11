
import numpy as np
import warnings
from datetime import datetime
from astropy.time import Time, TimeDelta
from astropy.utils import data
from astropy.utils import iers
import astropy.coordinates as ac
from astropy import units as un
from astropy.constants import c as speed_of_light
from astropy.coordinates.earth import OMEGA_EARTH
from multiprocessing import Process
import pylab as pl


OMEGA_EARTH_ITRS = np.array([0, 0, 1]) * OMEGA_EARTH

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


def astropy_delay(src, time, ant0, ant1):
    """
    Basic calculation of geometric delay using astropy in ITRS.

    Parameters
    ----------
    src: astropy.coordiantes.SkyCoord
        Pointing center
    time: astropy.time.Time
        Time
    ant0: astropy.coordinates.EarthLocation
        Location of antenna 0
    ant1: astropy.coordinates.EarthLocation
        Location of antenna 1

    Returns
    -------
    float:
        Geometric delay in microseconds
    """
    with ac.solar_system_ephemeris.set('jpl'):
        geo_src = src.transform_to(ac.ITRS(obstime=time))
        svec = geo_src.cartesian.xyz
        itrs0 = ant0.get_itrs(time)
        itrs1 = ant1.get_itrs(time)

        p0vec = itrs0.cartesian.xyz
        p1vec = itrs1.cartesian.xyz

    dx = p1vec - p0vec
    delay0 = np.dot(dx, svec) / speed_of_light

    return delay0.to_value('us')


def astropy_delay_rate(src, time, ant0, ant1):
    """
    Basic calculation of geometric delay rate using astropy in ITRS.

    Parameters
    ----------
    src: astropy.coordiantes.SkyCoord
        Pointing center
    time: astropy.time.Time
        Time
    ant0: astropy.coordinates.EarthLocation
        Location of antenna 0
    ant1: astropy.coordinates.EarthLocation
        Location of antenna 1

    Returns
    -------
    float:
        Geometric delay rate in us Hz
    """
    with ac.solar_system_ephemeris.set('jpl'):
        geo_src = src.transform_to(ac.ITRS(obstime=time))
        svec = geo_src.cartesian.xyz
        itrs0 = ant0.get_itrs(time)
        itrs1 = ant1.get_itrs(time)

        p0vec = itrs0.cartesian.xyz
        p1vec = itrs1.cartesian.xyz

    dx = p1vec - p0vec
    bxw = -np.cross(dx, OMEGA_EARTH_ITRS)
    dr0 = np.dot(bxw, svec) / speed_of_light

    return dr0.to_value("us Hz")
