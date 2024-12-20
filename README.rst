*************************
Python bindings to CALC11
*************************

This package provides a Python interface to running the CALC delay modeling tool by providing parameters as
astropy objects.

CALC has a long history, originally written in Fortran 77 with parts updated into Fortran 90 over time. CALC aims to compute
VLBI baseline delays to picosecond precision by incorporating general relativistic deflection from the Sun, Earth, and other planets, and including baseline shifts due to 
solid Earth and ocean tide loading, and deflection due to atmospheric contributions [EUBANKS91]_.

This repository carries a modified version of the CALC source code. The modifications enable easier operation
as a library, rather than as a standalone program, and support dynamic array allocation for some arrays that were previously
fixed-size.


References
----------
    .. [EUBANKS91] Eubanks, Marshall, et al. Proceedings of the U.S. Naval Observatory Workshop on Relativistic Models for Use in Space Geodesy. 1991.


Installation
------------

You will need cmake and a fortran compiler (e.g., gfortran) to build pycalc11.

Install directly from the repository with pip::

    pip install git+https://github.com/aelanman/pycalc11.git

The source code for CALC is in the repository and will be built and added as an
extension to the pycalc11 module.

On a first run, `pycalc11` will download and cache the JPL DE421 ephemeris file from ATNF. This can take around 30s.

Quick Start
-----------

The interface to running CALC11 is provided by the ``Calc`` class. A ``Calc`` instance may be initialized
with lists of stations, sources, a start time, and a duration, with the station locations given as
astropy ``EarthLocation`` instances, source locations given as astropy ``SkyCoord``, start time as an astropy ``Time``
instance, and duration as a float representing the length of the scan in minutes.::

    from pycalc11 import Calc
    import numpy as np
    from astropy import coordinates as ac
    from astropy import units as un
    from astropy.time import Time

    # ------------
    # First, create parameters as astropy objects.
    # ------------

    # Times
    time = Time("2020-10-02T15:30:00.00", format="isot", scale='utc')
    duration_min = 10

    # Names of sites in astropy
    anames = ['gbt', 'vla', 'ALMA', 'mwa', 'chime']
    site_locs = [ac.EarthLocation.of_site(nm) for nm in anames]

    # Names of sites in ocean loading data
    site_names = ['GBT-VLBA', 'VLA_E9', 'ALMA', 'MWA', 'CHIME']

    n_srcs = 500    # Arbitrarily many. Previously limited to 300.
    source_coords = ac.SkyCoord(
        ra=np.random.uniform(0, 360, n_srcs),
        dec=np.random.uniform(-90,90, n_srcs),
        unit='deg',
        frame='icrs',
    )

    # ------------
    # Initialize the Calc object.
    # ------------

    ci = Calc(
        station_names=site_names,
        station_coords=site_locs,
        source_coords=source_coords,
        time=time,
        duration_min=duration_min,
    )

    # ------------
    # Run the driver
    # ------------
    ci.run_driver()


To include ocean loading effects, the provided `station_names` must match entries in the ocean loading (OL) and
ocean pole tide loading (OPTL) data sets. The `OceanFiles` class gives more information on these.

Alternatively, a `.calc` file may be given via the `calc_file` keyword. If both a calc file and keywords are
set, the keywords will override settings from the file.

Once the `Calc` object is initialized, the method `run_driver` will run the main function of CALC11. The results
of this (delays, delay rates, partial derivatives) are then set in the attributes `Calc.delay`, `Calc.delay_rate`, 
and `Calc.partials`, respectively. Those properties' docstrings define the array axes.
