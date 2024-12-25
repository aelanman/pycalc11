"""Utility functions for testing and line profiling."""

import numpy as np
from datetime import datetime
import atexit
import warnings
from functools import partial
from inspect import isclass, isfunction

from astropy.utils import iers
import astropy.coordinates as ac
from astropy.constants import c as speed_of_light
from astropy.coordinates.earth import OMEGA_EARTH


OMEGA_EARTH_ITRS = np.array([0, 0, 1]) * OMEGA_EARTH

iers_tab = iers.earth_orientation_table.get()


def get_leap_seconds(tobj):
    """Find current TAI - UTC for a given time."""
    lsec_table = iers.LeapSeconds.auto_open().as_array()
    for ti, (yr, mo, _) in enumerate(lsec_table):
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
    Calculate geometric delay using astropy in ITRS.

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
    astropy.Quantity:
        Geometric delay in microseconds
    """
    with ac.solar_system_ephemeris.set("jpl"):
        geo_src = src.transform_to(ac.ITRS(obstime=time))
        svec = geo_src.cartesian.xyz
        itrs0 = ant0.get_itrs(time)
        itrs1 = ant1.get_itrs(time)

        p0vec = itrs0.cartesian.xyz
        p1vec = itrs1.cartesian.xyz

    dx = p1vec - p0vec
    delay0 = np.dot(dx, svec) / speed_of_light

    return delay0.to("us")


def astropy_delay_rate(src, time, ant0, ant1):
    """
    Calculate geometric delay rate using astropy in ITRS.

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
    astropy.Quantity:
        Geometric delay rate in us Hz
    """
    with ac.solar_system_ephemeris.set("jpl"):
        geo_src = src.transform_to(ac.ITRS(obstime=time))
        svec = geo_src.cartesian.xyz
        itrs0 = ant0.get_itrs(time)
        itrs1 = ant1.get_itrs(time)

        p0vec = itrs0.cartesian.xyz
        p1vec = itrs1.cartesian.xyz

    dx = p1vec - p0vec
    bxw = -np.cross(dx, OMEGA_EARTH_ITRS)
    dr0 = np.dot(bxw, svec) / speed_of_light

    return dr0.to("us Hz")


prof = None


def do_profiling(func_list=None, time=True, memory=False):
    """
    Run time or memory profiling on module functions.

    Places a LineProfiler object in the module namespace, and registers its
    dumping/printing functions to run at the end. When the Python environment closes,
    the profiler functions print_stats (and dump_stats, if dump_raw is True) will
    execute, saving profiler data to file.

    Parameters
    ----------
    func_list: list
        List of function names (strings) to profile.
    time: bool
        Enable line by line time profiling (Default True)
    memory: bool
        Enable line by line memory profiling (Default False)
    """
    global prof

    if func_list is None:
        func_list = ["run_driver", "set_stations", "set_sources", "alloc_out_arrays", "check_sites"]

    if memory and time:
        warnings.warn("Cannot run memory and time profilers at the same time. Skipping profiling.")
        return

    if time:
        from line_profiler import LineProfiler

        prof = LineProfiler()
        printfunc = prof.print_stats
        ofname = "time_profile.out"

    if memory:
        from memory_profiler import LineProfiler, show_results

        prof = LineProfiler()
        printfunc = partial(show_results, prof)
        ofname = "memory_profile.out"

    import pycalc11 as _pycalc11

    if prof is None:
        warnings.warn("No profiling requested. Choose memory or time = True to use profiling.")
        return

    # Add module functions to profiler.
    for mod_it in _pycalc11.__dict__.values():
        if isfunction(mod_it) and mod_it.__name__ in func_list:
            prof.add_function(mod_it)
        if isclass(mod_it):
            for item in mod_it.__dict__.values():
                if isfunction(item) and item.__name__ in func_list:
                    prof.add_function(item)

    # Write out profiling report to file.
    ofile = open(ofname, "w")
    atexit.register(ofile.close)
    atexit.register(printfunc, stream=ofile)
    prof.enable_by_count()
