#!/bin/env python
"""
Module to parse difxcalc output files (.im).
"""

import numpy as np
from datetime import datetime
from astropy.time import Time, TimeDelta


class TimeRange:
    """
    Helper class to quickly represent a range of times.

    Parameters
    ----------
    start: astropy.time.Time
        Start time
    end: astropy.time.Time
        End time
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __contains__(self, other):
        """
        Check if start <= other <= end

        Parameters
        ----------
        other: astropy.time.Time
            Time to check.
        """
        return (
            (other.mjd > self.start.mjd or np.isclose(other.mjd, self.start.mjd, atol=1e-6, rtol=0.0)) and
            (other.mjd <= self.end.mjd or np.isclose(other.mjd, self.end.mjd, atol=1e-6, rtol=0.0))
        )

    def __repr__(self):
        """String representation."""
        return "{} -- {}".format(self.start.isot, self.end.isot)


class CalcReader:
    """
    Reads difxcalc output .im files and evaluates polynomials for delay.

    Parameters
    ----------
    filename: str
        Path to file to load. Optional.
    """

    poly_order = None
    interval = None
    antnames = None
    antnums = None
    num_scans = None
    params = None
    filename = None

    def __init__(self, filename=None):
        if filename is not None:
            self.read_im(filename)


    def read_im(self, filename):
        """
        Read in .im file.

        Parameters
        ----------
        filename: str
            Path to file to load. Optional.
        """
        self.filename = filename
        lines = open(filename, 'r').readlines()

        # Form a dictionary of data in the file.
        # This needs to be done line by line, because some line prefixes are the same,
        # and their meaning changes based on the order in the file.
        dat_dict = {}

        cur_scan = -1
        cur_poly = -1
        for line in lines:

            pref, dat = line.split(':')
            dat = dat.rstrip()
            pref = pref.replace(' ', '_')

            try:
                dat = int(dat)
            except ValueError:
                try:
                    dat = float(dat)
                except ValueError:
                    dat = dat.strip()
                    pass

            if pref.startswith("SCAN"):
                spl = pref.split('_')
                cur_scan = int(spl[1])
                cur_poly = -1   # reset

                if spl[2] == 'POLY':
                    cur_poly = int(spl[3])

            if cur_poly > -1 and not pref.startswith("SCAN"):
                pref = f"POLY_{cur_poly}_" + pref

            if cur_scan > -1 and not pref.startswith("SCAN"):
                pref = f"SCAN_{cur_scan}_" + pref

            if any(
                key in pref
                for key in ['DELAY', 'DRY', 'WET', 'AZ', 'EL_GEOM', 'U_(m)', 'V_(m)', 'W_(m)']
            ):
                lst = np.array(dat.split()).astype(float)
                dat = np.poly1d(lst[::-1])

            dat_dict[pref] = dat

            self.params = dat_dict

        # Create a dictionary identifying valid time ranges for each set of polynomials.
        self.poly_ranges = {}
        for scan_i in range(self.params['NUM_SCANS']):
            npoly = self.params[f"SCAN_{scan_i}_NUM_POLY"]
            for pol_i in range(npoly):
                key = (scan_i, pol_i)
                ti = Time(
                        self.params[f"SCAN_{scan_i}_POLY_{pol_i}_MJD"]
                        + self.params[f"SCAN_{scan_i}_POLY_{pol_i}_SEC"] / (24 * 3600),
                format='mjd')
                tf = ti + TimeDelta(self.params['INTERVAL_(SECS)'], format='sec')

                self.poly_ranges[key] = TimeRange(ti, tf)

        # Collect the start time from other keywords.
        start_date = datetime(
            self.params['START_YEAR'],
            self.params['START_MONTH'],
            self.params['START_DAY'],
            self.params['START_HOUR'],
            self.params['START_MINUTE'],
            self.params['START_SECOND'],
        )
        self.start_date = Time(start_date, format='datetime')

        # Get the telescope names and numbers:
        n_ants = self.params['NUM_TELESCOPES']
        self.antnames = [None] * n_ants
        self.antnums = list(range(n_ants))
        for ai in range(n_ants):
            self.antnames[ai] = self.params[f"TELESCOPE_{ai}_NAME"]

    def _get_polykey(self, time, ant_num, src_num, scan):
        # Find polynomial for time/ant/scan.
        polyind = None
        for pi, ran in self.poly_ranges.items():
            if (time in ran) and (pi[0] == scan):
                polyind = pi
                break
        if polyind is None:
            raise ValueError(f"Time {time} is not covered by current polynomials for scan {scan}.")

        pscan, pnum = polyind
        key = f"SCAN_{pscan}_POLY_{pnum}_SRC_{src_num}_ANT_{ant_num}"

        return polyind, key

    def delay(self, ant_num, time, src_num, scan=0):
        """
        Geocentric delay in microseconds.

        Parameters
        ----------
        ant_num: int
            Index of antenna.
        time: astropy.time.Time
            Time of observation.
        src_num: int
            Source number.
        scan: int
            Scan index number (default 0)

        Returns
        -------
        float
            Delay in microseconds.
        """
        if self.poly_ranges is None:
            raise ValueError("Need to read a file first.")

        polyind, key = self._get_polykey(time, ant_num, src_num, scan)
        poly = self.params[key + "_DELAY_(us)"]
        polystart = self.poly_ranges[polyind].start
        # Evaluate polynomial:
        dt = (time - polystart).sec
        delay = np.polyval(poly, dt)

        return delay

    def baseline_delay(self, ant1, ant2, time, src_num, scan=0):
        """
        Relative delay between two antennas, in microseconds.

        Returns light travel time from ant1 to ant2.

        Parameters
        ----------
        ant1: int
            Index of antenna 1.
        ant2: int
            Index of antenna 2
        time: astropy.time.Time
            Time of observation.
        src_num: int
            Source number.
        scan: int
            Scan index number (default 0)

        Returns
        -------
        float
            Delay in microseconds.
        """
        return (
            self.delay(ant2, time, src_num, scan=scan)
            - self.delay(ant1, time, src_num, scan=scan)
        )

