#!/bin/env python
"""
Functions to make a difxcalc input file (.calc) given telescope locations,
times, sources, and other parameters in astropy classes.
"""

import numpy as np
import warnings
from astropy.time import TimeDelta

from .utils import get_leap_seconds, iers_tab

# Up to date EOPs and leap second tables from the IERS


def make_calc(station_coords, station_names, source_coords,
              source_names, time, duration_min, ofile_name=None,
              im_filename=None):
    """
    Make a .calc file as input to difxcalc.

    Parameters
    ----------
    station_coords: list of astropy.coordinates.EarthLocation
        Positions of telescopes on the Earth.
    station_names: list of str
        Corresponding telescope names. These must be unique.
    source_coords: list of astropy.coordinates.SkyCoord
        Source positions to include.
    source_names: list of str
        Names of sources in the list.
    time: astropy.time.Time
        Start time of the observation.
    duration_min: float
        Duration of the observation in minutes.
    ofile_name: str
        Output .calc file name.
        Defaults to "new.calc".
    im_filename: str
        Name of the .im file to be produced by difxcalc.
        Defaults to "new.im"
    """
    # Filename defaults
    if ofile_name is None:
        ofile_name = 'new.calc'
    if im_filename is None:
        ls = ofile_name.split('.')
        im_filename = ofile_name + '.im'
        if ls[-1] == 'calc':
            ls.pop()
            im_filename = ".".join(ls) + '.im'

    # ----------------------------
    # Start time and job params.
    # ----------------------------
    lines = []
    newlines = [
        "JOB ID:             4",
        "JOB START TIME:     {:.8f}".format(time.mjd),
        "JOB STOP TIME:      {:.8f}".format(time.mjd + duration_min / (24 * 60)),
        "DUTY CYCLE:         1.000",
        "OBSCODE:            DUMMY",
        "DIFX VERSION:       DIFX-2.6.2",
        "DIFX LABEL:         VLBADIFX-2.6.2",
        "SUBJOB ID:          0",
        "SUBARRAY ID:        0",
        "VEX FILE:           dummy.vex.obs",
        "START MJD:          {:.8f}".format(time.mjd),
        "START YEAR:         {:.0f}".format(time.datetime.year),
        "START MONTH:        {:.0f}".format(time.datetime.month),
        "START DAY:          {:.0f}".format(time.datetime.day),
        "START HOUR:         {:.0f}".format(time.datetime.hour),
        "START MINUTE:       {:.0f}".format(time.datetime.minute),
        "START SECOND:       {:.0f}".format(time.datetime.second),
        "IM FILENAME:        dummy.im",
        "FLAG FILENAME:      dummy.flag",

    ]
    lines.extend(newlines)

    # ----------------------------
    # Earth Orientation Parameters
    # ----------------------------
    tai_utc = []
    ut1_utc = []
    mjd = []
    xy = []
    times = time + TimeDelta(range(2), format='jd')
    for tt in times:
        mjd.append(np.floor(tt.mjd))
        tai_utc.append(get_leap_seconds(tt))
        ut1_utc.append(tt.delta_ut1_utc.reshape(1)[0])

        # polar motion
        xy.append([z.to_value('arcsec') for z in iers_tab.pm_xy(tt)])

    lines.append("NUM EOPS: {:d}".format(len(times)))
    for ti in range(len(times)):
        newlines = [
            "EOP {:d} TIME (mjd):   {:.0f}".format(ti, mjd[ti]),
            "EOP {:d} TAI_UTC (sec):{:.0f}".format(ti, tai_utc[ti]),
            "EOP {:d} UT1_UTC (sec): {:.10f}".format(ti, ut1_utc[ti]),
            "EOP {:d} XPOLE (arcsec): {:.10f}".format(ti, xy[ti][0]),
            "EOP {:d} YPOLE (arcsec): {:.10f}".format(ti, xy[ti][1]),
        ]
        lines.extend(newlines)

    # ----------------------------
    # Sources
    # ----------------------------
    # CALCODE = calibration code, typically A,B,C for calibrators, G for a gated pulsar, or blank for normal target
    # https://www.atnf.csiro.au/vlbi/dokuwiki/lib/exe/fetch.php/difx/difxuserguide.pdf
    lines.append("NUM SOURCES: {:d}".format(len(source_names)))
    for si, (coord, name) in enumerate(zip(source_coords, source_names)):
        newlines = [
            "SOURCE {:d} NAME:      {}".format(si, name),
            "SOURCE {:d} RA:        {:.8f}".format(si, coord.ra.rad),      # radians
            "SOURCE {:d} DEC:       {:.8f}".format(si, coord.dec.rad),      # radians
            "SOURCE {:d} CALCODE:   B".format(si),
            "SOURCE {:d} QUAL:      0".format(si),
        ]
        lines.extend(newlines)

    # ----------------------------
    # Telescopes
    # ----------------------------
    n_ants = len(station_names)
    # check for lowercase in telescope names
    for char in "".join(station_names):
        if char.islower():
            warnings.warn("Station names should be all uppercase. "
                          "Lowercase detected. Changing to uppercase.")
            break
    lines.append("NUM TELESCOPES:     {}".format(n_ants))
    for ti in range(n_ants):
        newlines=[
            "TELESCOPE {:d} NAME:   {}".format(ti, station_names[ti].upper()),
            "TELESCOPE {:d} MOUNT:  AZEL".format(ti),
            "TELESCOPE {:d} OFFSET (m): 0.0000".format(ti),
            "TELESCOPE {:d} X (m): {:.8f}".format(ti, station_coords[ti].x.to_value('m')),
            "TELESCOPE {:d} Y (m): {:.8f}".format(ti, station_coords[ti].y.to_value('m')),
            "TELESCOPE {:d} Z (m): {:.8f}".format(ti, station_coords[ti].z.to_value('m')),
            "TELESCOPE {:d} SHELF:  None".format(ti),
        ]
        lines.extend(newlines)

    # ----------------------------
    # Scan Info (one scan, multiple pointings centers).
    # ----------------------------

    newlines = [
        "NUM SCANS:          1",
        "NUM SPACECRAFT:     0",
        "SCAN 0 IDENTIFIER:  No0004",
        "SCAN 0 START (S):   0",
        "SCAN 0 DUR (S):     {}".format(duration_min * 60),
        "SCAN 0 OBS MODE NAME:JWST",
        "SCAN 0 UVSHIFT INTERVAL (NS):2000000000",
        "SCAN 0 AC AVG INTERVAL (NS):2000000",
        "SCAN 0 NUM PHS CTRS: {}".format(len(source_names)),
        "SCAN 0 POINTING SRC:0",
    ]
    for si in range(len(source_names)):
        newlines.extend([
            "SCAN 0 PHS CTR {}:   {}".format(si, si),
        ])
    lines.extend(newlines)


    # ----------------------------
    # Other necessary attributes.
    # ----------------------------
    other = [
        "SPECTRAL AVG:       1",
        "TAPER FUNCTION:     UNIFORM",
        "IM FILENAME:        {}".format(im_filename),
        "FLAG FILENAME:      {}".format(im_filename + '.flag'),
    ]
    lines.extend(other)

    lines = [l + '\n' for l in lines]
    with open(ofile_name, 'w') as ofile:
        ofile.writelines(lines)
