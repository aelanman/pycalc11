"""
Functions to make a difxcalc input file (.calc) given telescope locations,
times, sources, and other parameters in astropy classes.
"""

import numpy as np
import warnings
from astropy.time import TimeDelta

from .utils import get_leap_seconds, iers_tab

# Up to date EOPs and leap second tables from the IERS


def make_calc(
    station_coords,
    station_names,
    source_coords,
    start_time,
    duration_min,
    ofile_name=None,
    im_filename=None,
):
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
    start_time: astropy.time.Time
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
        ofile_name = "new.calc"
    if im_filename is None:
        ls = ofile_name.split(".")
        im_filename = ofile_name + ".im"
        if ls[-1] == "calc":
            ls.pop()
            im_filename = ".".join(ls) + ".im"

    # ----------------------------
    # Start time and job params.
    # ----------------------------
    lines = []
    newlines = [
        "JOB ID:             4",
        f"JOB START TIME:     {start_time.mjd:.8f}",
        f"JOB STOP TIME:      {start_time.mjd + duration_min / (24 * 60):.8f}",
        "DUTY CYCLE:         1.000",
        "OBSCODE:            DUMMY",
        "DIFX VERSION:       DIFX-2.6.2",
        "DIFX LABEL:         VLBADIFX-2.6.2",
        "SUBJOB ID:          0",
        "SUBARRAY ID:        0",
        "VEX FILE:           dummy.vex.obs",
        f"START MJD:          {start_time.mjd:.8f}",
        f"START YEAR:         {start_time.datetime.year:.0f}",
        f"START MONTH:        {start_time.datetime.month:.0f}",
        f"START DAY:          {start_time.datetime.day:.0f}",
        f"START HOUR:         {start_time.datetime.hour:.0f}",
        f"START MINUTE:       {start_time.datetime.minute:.0f}",
        f"START SECOND:       {start_time.datetime.second:.0f}",
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
    times = start_time + TimeDelta(range(2), format="jd")
    for tt in times:
        mjd.append(np.floor(tt.mjd))
        tai_utc.append(get_leap_seconds(tt))
        ut1_utc.append(tt.delta_ut1_utc.reshape(1)[0])

        # polar motion
        xy.append([z.to_value("arcsec") for z in iers_tab.pm_xy(tt)])

    lines.append(f"NUM EOPS: {len(times):d}")
    for ti in range(len(times)):
        newlines = [
            f"EOP {ti:d} TIME (mjd):{mjd[ti]:.0f}",
            f"EOP {ti:d} TAI_UTC (sec):{tai_utc[ti]:.0f}",
            f"EOP {ti:d} UT1_UTC (sec):{ut1_utc[ti]:.8f}",
            f"EOP {ti:d} XPOLE (arcsec):{xy[ti][0]:.10f}",
            f"EOP {ti:d} YPOLE (arcsec):{xy[ti][1]:.10f}",
        ]
        lines.extend(newlines)

    # ----------------------------
    # Sources
    # ----------------------------
    # CALCODE = calibration code, typically A,B,C for calibrators,
    #           G for a gated pulsar, or blank for normal target
    # https://www.atnf.csiro.au/vlbi/dokuwiki/lib/exe/fetch.php/difx/difxuserguide.pdf
    lines.append(f"NUM SOURCES: {len(source_coords):d}")
    for si, coord in enumerate(source_coords):
        name = f"src{si:d}"
        newlines = [
            f"SOURCE {si:d} NAME:      {name}",
            f"SOURCE {si:d} RA:        {coord.ra.rad:.10f}",  # radians
            f"SOURCE {si:d} DEC:       {coord.dec.rad:.10f}",  # radians
            f"SOURCE {si:d} CALCODE:   B",
            f"SOURCE {si:d} QUAL:      0",
        ]
        lines.extend(newlines)

    # ----------------------------
    # Telescopes
    # ----------------------------
    n_ants = len(station_names)
    # check for lowercase in telescope names
    for char in "".join(station_names):
        if char.islower():
            warnings.warn(
                "Station names should be all uppercase. "
                "Lowercase detected. Changing to uppercase."
            )
            break
    lines.append(f"NUM TELESCOPES:     {n_ants}")
    for ti in range(n_ants):
        newlines = [
            f"TELESCOPE {ti:d} NAME:   {station_names[ti].upper()}",
            f"TELESCOPE {ti:d} MOUNT:  AZEL",
            f"TELESCOPE {ti:d} OFFSET (m): 0.0000",
            "TELESCOPE {:d} X (m): {:.8f}".format(ti, station_coords[ti].x.to_value("m")),
            "TELESCOPE {:d} Y (m): {:.8f}".format(ti, station_coords[ti].y.to_value("m")),
            "TELESCOPE {:d} Z (m): {:.8f}".format(ti, station_coords[ti].z.to_value("m")),
            f"TELESCOPE {ti:d} SHELF:  None",
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
        f"SCAN 0 DUR (S):     {duration_min * 60}",
        "SCAN 0 OBS MODE NAME:JWST",
        "SCAN 0 UVSHIFT INTERVAL (NS):2000000000",
        "SCAN 0 AC AVG INTERVAL (NS):2000000",
        f"SCAN 0 NUM PHS CTRS: {len(source_coords)}",
        "SCAN 0 POINTING SRC:0",
    ]
    for si in range(len(source_coords)):
        newlines.extend([f"SCAN 0 PHS CTR {si}:   {si}"])
    lines.extend(newlines)

    # ----------------------------
    # Other necessary attributes.
    # ----------------------------
    other = [
        "SPECTRAL AVG:       1",
        "TAPER FUNCTION:     UNIFORM",
        f"IM FILENAME:        {im_filename}",
        "FLAG FILENAME:      {}".format(im_filename + ".flag"),
    ]
    lines.extend(other)

    lines = [ll + "\n" for ll in lines]
    with open(ofile_name, "w") as ofile:
        ofile.writelines(lines)
