
import numpy as np
import warnings
from astropy.time import Time, TimeDelta
from astropy import coordinates as ac
from astropy.constants import R_earth
from astropy.units import Quantity
from datetime import datetime
from scipy.interpolate import Akima1DInterpolator as Akima1DInterpolator
from .utils import get_leap_seconds, iers_tab
from . import calc11 as calc

class Calc:
    """
    A class to handle the current state of the calc11 fortran module and run its functions.

    A unique state is described by:
        - A set of telescope positions and names
        - A set of times [or a single time to start]
        - A set of source positions and names
        - A set of command line options:
            * Add atmosphere contributions or not
            * Baseline or geocenter mode
            * uvw mode
            * Near or far field
    When the state changes, the driver subroutine `adrvr` needs to be rerun before
    accessing results.

    Parameters
    ----------
    station_names: array_like of str
        Names of VLBI stations.
    station_coords: array_like of astropy.coordinates.EarthLocation
        Positions of VLBI stations
    source_coords: array_like of astropy.coordinates.SkyCoord
        Positions of sources in ICRS
    start_time: astropy.time.Time
        Start time of scan.
    duration_min: float
        Duration of scan in minutes
    base_mode: str
        CALC baseline mode.
        Either 'baseline' or 'geocenter'.
        Default 'geocenter'
    uvw_mode: str
        CALC uvw calculation mode
        exact [include partials], uncorr, approx, noatmo [exact with no atmo]
        Default 'exact.
    dry_atm: bool
        Include dry troposphere contribution to delay.
        Default True
    wet_atm: bool
        Include wet troposphere contribution to delay.
        Default True
    calc_file: str
        Path to a .calc configuration file.
        If provided, all other parameters are optional.

    Notes
    -----
    Keyword arguments supersede options in a loaded .calc file. For instance,
    if a set of sources is provided via the `source_names` and `source_coords` options,
    then the stations used in CALC will match those and not the ones in the provided .calc file.
    """

    _rerun = True       # Rerun driver before accessing results

    _delay = None
    _delay_rate = None
    _partials = None

    def _state_dep(func):
        """Decorator ensuring that state is consistent with products."""
        def require_rerun(self):
            if self._rerun:
                raise ValueError("Need to rerun adrivr before accessing data. Use Calc.run_driver().")
            else:
                return func(self)
        require_rerun.__doc__ = func.__doc__
        return require_rerun

    def __init__(self, station_names=None, station_coords=None, source_coords=None,
            start_time=None, duration_min=None,
            base_mode='geocenter', uvw_mode='exact', dry_atm=True, wet_atm=True,
            calc_file=None, d_interval=24,
        ):
        # Setting defaults
        self._reset()        # Clear if there's another instance.
        calc.mode.c_mode = 'difx  '
        calc.contrl.near_far = 'Far-field '
        calc.contrl.l_time = 'dont-solve'   # Optionally solve for LTT if using near-field
        calc.nfo.spoffset = 'NoOffset'      # Use near-field ephemeris epoch offset (set nfo.t_offset accordingly)
        calc.contrl.verbose = 0             # Print out pointing source name

        self.base_mode = base_mode
        self.dry_atm = dry_atm
        self.wet_atm = wet_atm
        self.uvw_mode = uvw_mode            # Options are exact [include partials], uncorr, approx, noatmo [exact with no atmo]

        # Steps within 2 min chunks
        calc.contrl.d_interval = d_interval        # Step size in s throughout 2 min epoch
        calc.contrl.epoch2m = (120.0001/calc.contrl.d_interval) + 1 # Number of steps in 2 min epoch

        # Check for required parameters
        kwargs = {
                "station_names" : station_names, "station_coords" : station_coords, "source_coords": source_coords,
                "start_time" : start_time, "duration_min" : duration_min, "base_mode" : base_mode,
        }

        # Behavior:
        #   - If calc_file is provided, set up from it.
        #   - If other parameters are provided, they override calc_file settings.
        if calc_file is not None:
            self.calc_file = calc_file
        elif not all(v is not None for k,v in kwargs.items()):
            raise ValueError(
                "If calc_file is not set, then all of the following must be set:"
                "\n \t " + ", ".join(kwargs.keys()) +
                "\nMissing: " + ", ".join([k for k,v in kwargs.items() if v is None])
            )
        else:
            del self.calc_file

        # Initialize stations, sources, scan, and eops
        self.station_names = station_names
        self.station_coords = station_coords
        self.source_coords = source_coords
        self._start_time = start_time
        self._duration_min = duration_min

        if start_time is not None:
            self.set_eops(start_time)
            self.set_scan(start_time, duration_min)

        # Add geocenter station
        self._add_geocenter()

        calc.dinitl(1)

    def parse_calcfile(self, calcfile):
        """
        Get axis sizes from calcfile.

        (Number of sources and number of stations)
        """
        nsrcs = None
        nstat = None
        with open(calcfile, 'r') as cfile:
            ln0 = cfile.readline()
            while ln0:
                if ln0.startswith("NUM SOURCES"):
                    nsrcs = int(ln0.split(":")[-1])
                elif ln0.startswith("NUM TELESCOPES"):
                    nstat = int(ln0.split(":")[-1])
                ln0 = cfile.readline()

        return nsrcs, nstat

    def run_driver(self):
        """Run adrivr to get results."""
        calc.dsiti(1)      # In case stations changed, need to reload ocean loading etc.
        calc.dwobi()       # Recalc wobble coefs in case EOPs changed.
        self.alloc_out_arrays()
        e2m = calc.contrl.epoch2m - 1
        for ii in range(calc.calc_input.intrvls2min):
            calc.adrivr(1,ii+1)
            slc = np.s_[ii*e2m: ii*e2m + e2m, :, :, :]
            self._delay[slc] = calc.outputs.delay_f[:-1, :, :, 1:]   # Skip pointing source
            self._delay_rate[slc] = calc.outputs.rate_f[:-1, :, :, 1:]
            self._partials["dD_dRA"][slc] = calc.outputs.partials_f[0, 0, :-1, :, :, 1:]
            self._partials["dDR_dRA"][slc] = calc.outputs.partials_f[0, 1, :-1, :, :, 1:]
            self._partials["dD_dDEC"][slc] = calc.outputs.partials_f[1, 0, :-1, :, :, 1:]
            self._partials["dDR_dDEC"][slc] = calc.outputs.partials_f[1, 1, :-1, :, :, 1:]
            self._times[ii*e2m: ii*e2m + e2m] = [
                    np.datetime64(datetime(*t)) for t in calc.outputs.iymdhms_f[:-1]
            ]

        self._rerun = False

    def interpolate_delays(self, t_0):
        """ Evaluates the delays at an absolute time t_0 by calculating an akima spline and interpolating delays along the time axis.
        The inputs to the akima spline is the time in seconds between self.times[0] and the time the delays should be evaluated.

        Parameters
        -----------------
        t_0 : astropy.time.Time
            Absolute time at which delays are to be evaluated

        Returns
        -----------
        np.ndarray
            The delays in seconds evaluated at t_0
        """
        x=(t_0-self.times[0]).sec
        return self.delays_dt(x)

    @property
    def delays_dt(self):
        return Akima1DInterpolator((self.times - self.times[0]).sec, self.delay.to_value('s'), axis=0)

    def _reset(self):
        """Reset all common block items that were not initialized at startup."""
        self._delay = None
        self._delay_rate = None
        self._partials = None
        calc.srcmod.numstr = 0
        calc.units.ipoint = 40              # Reset units counter
        calc.units.iutot = 3

        for key, part in calc.__dict__.items():
            # Select only common blocks
            if (key.startswith('_')
                    or all(dk.startswith('_') for
                           dk in part.__dict__)):
                continue
            for item, value in part.__dict__.items():
                if not is_initialized(key, item):
                    value[...] = b'' if value.dtype.kind == 'S' else 0
                    part.__dict__[item] = value
            calc.__dict__[key] = part

    def _add_geocenter(self):
        """Insert a station at index 0 for geocenter."""
        if calc.calc_input.sites[0] == 'GEOCENTR':
            raise ValueError("Site 0 is already geocenter")

        calc.calc_input.sites[0] = "GEOCENTR"
        calc.calc_input.axis[0] = "AZEL"
        calc.sitcm.sitaxo[0] = 0.0    # Axis offset
        calc.sitcm.sitxyz[:, 0] = [0, 0, 0]

    def set_scan(self, time, duration_min):
        """
        Set start time and duration of scan.

        Only one scan is supported.

        Parameters
        ----------
        time : astropy.time.Time
            Start time of scan
        duration_min : float
            Duration of scan in minutes.
        """
        jstart = time
        jstop = time + TimeDelta(duration_min * 60, format='sec')
        calc.con.kctic = 1      # Ensures CT = AT. Shouldn't make a difference, since CT isn't used.
        calc.ut1cm.xintv[0] = jstart.jd
        calc.ut1cm.xintv[1] = jstop.jd
        calc.ut1cm.intrvl[:, 0] = list(jstart.ymdhms)[:5]    # intrvl[:, 1] unused
        calc.calc_input.numscans = 1
        calc.calc_input.scanid = "Scan001".rjust(10)
        calc.calc_input.scandur = duration_min * 60  # Seconds
        calc.calc_input.pointingsrc = 1

        # As defined in dscan.f
        calc.calc_input.intrvls2min = np.ceil(duration_min / 2 + 0.001) + 1     # Number of 2 min intervals
        self._rerun = True

    def set_eops(self, time):
        """Set Earth orientation parameters.

        Parameters
        ----------
        time : astropy.time.Time
            Start time.
            EOPs will be computed for (time, time+1d, time+2d).
        """
        _times = time + TimeDelta(range(2), format='jd')
        mjd = [np.floor(tt.mjd) for tt in _times]
        tai_utc = [get_leap_seconds(tt) for tt in _times]
        ut1_utc = [tt.delta_ut1_utc.reshape(1)[0] for tt in _times]
        xy = [[z.to_value('milliarcsecond') for z in iers_tab.pm_xy(tt)] for tt in _times]

        eoptag = mjd
        numeop = len(eoptag)
        calc.ut1cm.ut1if = [eoptag[0] + 2400000.5, eoptag[1] - eoptag[0], numeop, 1]
        calc.wobcm.wobif = [eoptag[0] + 2400000.5, eoptag[1] - eoptag[0], numeop]
        calc.wobcm.xywob[:, :2] = np.transpose(xy)
        # Rounded to match difxcalc's internal value
        calc.ut1cm.ut1pt[:2] = np.round([tai_utc[ti] - ut1_utc[ti] for ti in range(2)], decimals=8)
        calc.calc_input.xleap_sec = tai_utc[0]
        self._rerun = True

    def alloc_out_arrays(self):
        """Allocate output arrays."""
        Nstat1 = self.nants
        Nstat2 = self.nants
        Nsrcs = self.nsrcs

        Ntimes = (calc.contrl.epoch2m - 1) * (calc.calc_input.intrvls2min)

        if self.base_mode == 'geocenter':
            Nstat1 = 1
        # These arrays are output for each epoch
        calc.alloc_out_arrays(calc.contrl.epoch2m, Nstat1, Nstat2, Nsrcs + 1)

        # These arrays collect results from all 2 min epochs.
        self._delay = np.zeros((Ntimes, Nstat1, Nstat2, Nsrcs))
        self._delay_rate = np.zeros((Ntimes, Nstat1, Nstat2, Nsrcs))
        self._partials = np.zeros(
            (Ntimes, Nstat1, Nstat2, Nsrcs),
            dtype=np.dtype([("dD_dRA", "f8"), ("dDR_dRA", "f8"), ("dD_dDEC", "f8"), ("dDR_dDEC", "f8")])
        )
        self._times = np.full(Ntimes, fill_value='NaT', dtype="<M8[us]")
        self._rerun = True

    @property
    def base_mode(self):
        """Baseline mode."""
        return (calc.contrl.base_mode[()]).decode('utf-8').strip()

    @base_mode.setter
    def base_mode(self, value):
        valid = ['baseline', 'geocenter']
        if not isinstance(value, str) or value.strip() not in valid:
            raise ValueError("base_mode must be one of the following: "  + ", ".join(valid))
        # In calc, base_mode must be a 10 character lowercase string.
        calc.contrl.base_mode[()] = value.ljust(10)
        self._rerun = True

    @property
    def dry_atm(self):
        """Include dry troposphere."""
        return calc.contrl.atmdr[()] == b'Add-dry   '

    @dry_atm.setter
    def dry_atm(self, value):
        if not isinstance(value, bool):
            raise ValueError("dry_atm must be a boolean type.")
        calc.contrl.atmdr[()] = b'Add-dry   ' if value else b''
        self._rerun = True

    @property
    def wet_atm(self):
        """Include wet troposphere."""
        return calc.contrl.atmwt[()] == b'Add-wet   '

    @wet_atm.setter
    def wet_atm(self, value):
        if not isinstance(value, bool):
            raise ValueError("wet_atm must be a boolean type.")
        calc.contrl.atmwt[()] = b'Add-wet   ' if value else b''
        self._rerun = True

    @property
    def uvw_mode(self):
        """UVW calculation mode."""
        return calc.contrl.uvw

    @uvw_mode.setter
    def uvw_mode(self, value):
        valid = ['exact', 'uncorr', 'approx', 'noatmo']
        if not value.strip() in valid:
            raise ValueError("uvw_mode must be one of the following: " + ", ".join(valid))
        calc.contrl.uvw = value.ljust(6)
        self._rerun = True

    @property
    def calc_file(self):
        """Path to .calc file."""
        return calc.contrl.calc_file_name[()]

    @calc_file.setter
    def calc_file(self, fname):
        calc.contrl.calc_file_name[()] = fname.ljust(128)[:128]
        nsrcs, nstat = self.parse_calcfile(fname)
        calc.alloc_source_arrays(nsrcs)
        calc.dget_input(1)
        calc.sitcm.numsit += 1      # Numsit in calc must include the geocenter.
        calc.dscan(1,1)
        self._rerun = True

    @calc_file.deleter
    def calc_file(self):
        calc.contrl.calc_file_name[()] = b""

    @property
    def nsrcs(self):
        """Number of sources."""
        return calc.calc_input.numphcntr

    @property
    def nants(self):
        """Number of stations."""
        return calc.sitcm.numsit -  1  # Drop geocenter

    @property
    def start_time(self):
        return self._start_time

    @start_time.setter
    def start_time(self, time):
        if time is None:
            return
        self._start_time = time
        self.set_eops(time)
        self.set_scan(time, self.duration_min)
        self._rerun = True

    @property
    def duration_min(self):
        return self._duration_min

    @duration_min.setter
    def duration_min(self, dur):
        if dur is None:
            return
        self._duration_min = dur
        self.set_scan(self._start_time, dur)

    @property
    def ant1_ind(self):
        """Index in stations array corresponding with ant1 axis of output arrays."""
        if self.base_mode == 'geocenter':
            return 0
        else:
            return list(range(self.nants - 1))

    @property
    def ant2_ind(self):
        """Index in stations array corresponding with ant2 axis of output arrays."""
        return list(range(self.nants))

    @property
    def station_names(self):
        """Station names."""
        return calc.calc_input.sites[1:self.nants+1].astype("U8")

    @station_names.setter
    def station_names(self, station_names):
        if station_names is None:
            return
        calc.sitcm.numsit = len(station_names) + 1
        tnames = [s.upper() for s in station_names]
        for ti in range(self.nants):
            # Offset of 1 to account for zeroth site being the geocenter
            calc.calc_input.sites[ti+1] = np.bytes_(tnames[ti].ljust(10))

        # Check ocean loading params are available
        OceanFiles.check_sites(self.station_names)
        self._rerun = True

    @property
    def station_coords(self):
        """Station coordinates in ITRF."""
        return Quantity(
            calc.sitcm.sitxyz[:, 1:self.nants+1],
            unit='m',
            copy=False,
        )

    @station_coords.setter
    def station_coords(self, tpos):
        if tpos is None:
            return
        calc.sitcm.numsit = len(tpos) + 1
        for ti in range(self.nants):
            # Offset of 1 to account for zeroth site being the geocenter
            calc.calc_input.axis[ti+1] = "AZEL"
            calc.sitcm.sitaxo[ti+1] = 0.0    # Axis offset
            calc.sitcm.sitxyz[0, ti+1] = tpos[ti].x.to_value('m')
            calc.sitcm.sitxyz[1, ti+1] = tpos[ti].y.to_value('m')
            calc.sitcm.sitxyz[2, ti+1] = tpos[ti].z.to_value('m')
        self._rerun = True

    @property
    def source_coords(self):
        """Source coordinates in ICRS ra/dec"""
        return ac.SkyCoord(
            ra=calc.srcmod.radec[0],
            dec=calc.srcmod.radec[1],
            unit='rad',
            frame='icrs',
            copy=False,     # This might not be respected...
        )

    @source_coords.setter
    def source_coords(self, src_coords):
        """Set source positions."""
        if src_coords is None:
            return
        Nsrc = 0 if src_coords.isscalar else len(src_coords)
        calc.calc_input.numphcntr = Nsrc
        calc.alloc_source_arrays(self.nsrcs)
        calc.srcmod.numstr = Nsrc
        calc.srcmod.radec[0, :] = np.round(src_coords.ra.rad, decimals=10)
        calc.srcmod.radec[1, :] = np.round(src_coords.dec.rad, decimals=10)

        calc.srcmod.phcntr[:Nsrc] = range(1, Nsrc+1)
        calc.srcmod.phcntr[Nsrc:] = -1
        self._rerun = True

    @property
    def src_ra(self):
        """Source right ascension in ICRS"""
        return ac.Longitude(
            calc.srcmod.radec[0],
            unit='rad',
            copy=False,
        )

    @property
    def src_dec(self):
        """Source declination in ICRS"""
        return ac.Latitude(
            calc.srcmod.radec[1],
            unit='rad',
            copy=False,
        )

    # ---------------
    # Results of running the driver
    # ---------------

    @property
    @_state_dep
    def times(self):
        """List of times for time axis of delay/delay_rate/partials arrays."""
        return Time(self._times, format='datetime64', scale='utc', copy=False)

    @property
    @_state_dep
    def delay(self):
        """
        Delay in seconds.

        Axes: (time, ant1, ant2, source)
        """
        # In baseline mode, the geocenter station is still in the arrays.
        if self.base_mode == 'baseline':
            val = self._delay[:, 1:, :, :]
        else:
            val = self._delay
        return Quantity(val, unit="s", copy=False)

    @property
    @_state_dep
    def delay_rate(self):
        """
        Delay rate in s Hz.

        Axes: (time, ant1, ant2, source)
        """
        if self.base_mode == 'baseline':
            val = self._delay_rate[:, 1:, :, :]
        else:
            val = self._delay_rate
        return Quantity(val, unit="s Hz", copy=False)

    @property
    @_state_dep
    def partials(self):
        """
        All partial derivatives of delay.

        Structured array, with units:
            dD_dRA = partial of delay wrt right ascension
            dDR_dRA = partial of delay rate wrt right ascension
            dD_dDEC = partial of delay wrt declination
            dDR_dDEC = partial of delay rate wrt declination
        """
        if self.base_mode == 'baseline':
            vals = self._partials[..., 1:, :, :]
        else:
            vals = self._partials
        return Quantity(
            vals,
            unit=("s/rad", "s Hz / rad", "s/rad", "s Hz/rad"),
            copy=False,
        )


class OceanFiles:
    """Interface to ocean loading and pole tide coef data."""

    # Ocean pole tide loading coefficients
    OPTL_file = calc.datafiles.optl_file.item().decode('utf-8').strip()
    # Ocean loading catalog
    OC_file = calc.datafiles.oc_file.item().decode('utf-8').strip()

    optl_data = None
    oc_data = None
    stat_pos = None     # Lons/lats/heights of stations in degrees/meters, from file.

    def __new__(cls):
        cls.read_optl()
        cls.read_oc()

    @classmethod
    def read_optl(cls):
        """Read ocean pole tide loading coefficients from file."""
        nms = ['name', 'code', 'lat', 'lon', 'urr', 'uri', 'unr', 'uni', 'uer', 'uei']
        dt = np.dtype(list(zip(nms, ['U8', 'U4'] + ['f8'] * 8)))
        with open(cls.OPTL_file, 'r') as optl_file:
            lns = optl_file.readlines()[5:]
            dat = []
            for li, ln in enumerate(lns):
                name = ln[1:9]
                code = ln[10:14]
                vals = list(map(float, ln[15:].split()))
                if len(vals) != 8:
                    continue
                dat.append(tuple([name, code] + vals))
        cls.optl_data = np.array(dat, dtype=dt)

    @classmethod
    def read_oc(cls):
        """Read ocean loading coefficients from file."""
        dt = np.dtype([('name', 'U8'), ('code', 'U4'),
                       ('lonlat_deg', 'f8', (3,)),('coefs', 'f8', (6,11))])
        with open(cls.OC_file, 'r') as oc_file:
            ocean_load = []
            stat = ""
            coef = []
            code = ""
            lonlat = None
            for ln in oc_file.readlines():
                spl = ln.strip().split()
                if "lon/lat" in ln:
                    # Get site position
                    stat = ln[3:11].replace(',', '')
                    if stat.lower() == 'geocentr':
                        lonlat = (0,0,0)
                    else:
                        lon, lat, height= tuple(map(float, spl[-3:]))
                        lonlat = (lon, lat, height)
                    continue
                if ln.startswith("$$"):
                    continue
                if len(spl) < 6:
                    if len(coef) != 0:
                        ocean_load.append(
                            tuple([stat.ljust(10), code.ljust(4)] + [lonlat] + [coef])
                        )
                    stat = ln[2:11]
                    code = ln[12:15] if len(ln) >= 11 else ""
                    code = code.replace('\n', '')
                    coef = []
                else:
                    coef.append(list(map(float, spl)))
        cls.oc_data = np.asarray(ocean_load, dtype=dt)

    @classmethod
    def list_optl_sites(cls):
        return cls.optl_data["name"]

    @classmethod
    def list_oc_sites(cls):
        return cls.oc_data["name"]

    @classmethod
    def _get_ind(cls, arr, name, code=None):
        name = name.ljust(8)
        if code is not None:
            code = code.ljust(4)
            return np.where(arr["code"] == code)[0]
        return np.where(arr["name"] == name)[0]

    @classmethod
    def get_oc_ind(cls, name, code=None):
        return cls._get_ind(cls.oc_data, name, code=code)

    @classmethod
    def get_optl_ind(cls, name, code=None):
        return cls._get_ind(cls.optl_data, name, code=code)

    @classmethod
    def check_sites(cls, site_names, site_pos=None, sep_tol=10):
        """
        Verify that the sites given are in the ocean data files.

        Parameters
        ----------
        site_names: list of str
            List of strings giving station names.
            Each name, stripped of trailing and preceding whitespace, should be in
            the `name` column of the oc_data and optl_data arrays.
        site_pos: list of astropy.coordinates.EarthLocation
            Optional. Positions corresponding with stations in site_names.
            Verify that the given lon/lat positions are close to positions listed in
            ocean pole tide loading file.
        sep_tol: float
            Tolerance for comparing site positions, in meters.
            Optional. Default = 10 meters
        """
        names = [s.ljust(8) for s in site_names]
        oc_names = [s for s in cls.oc_data['name']]

        oc_not_found = list(set(names) - set(cls.oc_data['name']))
        optl_not_found = list(set(names) - set(cls.optl_data['name']))

        if len(optl_not_found) > 0:
            warnings.warn(
                "No ocean pole tide loading coefficients found for " +
                ", ".join([s.strip() for s in optl_not_found])
            )

        if len(oc_not_found) > 0:
            warnings.warn(
                "No ocean loading coefficients found for " +
                ", ".join([s.strip() for s in oc_not_found])
            )

        if site_pos is not None:
            sndists = []
            for si, sn in enumerate(site_names):
                sn = sn.ljust(8)
                if sn in oc_not_found or sn in optl_not_found:
                    continue
                ind = oc_names.index(sn)
                lon, lat, height = cls.oc_data['lonlat_deg'][ind]
                ref_pos = ac.SkyCoord(ac.EarthLocation.from_geodetic(
                    lon=lon, lat=lat).itrs
                )
                dist = ref_pos.separation(ac.SkyCoord(site_pos[si].itrs))
                if dist.rad * R_earth.to_value('m') > sep_tol:
                    d = dist.rad * R_earth
                    sndists.append(f"\t {sn} : {d}")

            if len(sndists) > 0:
                warnings.warn(
                  f"Station positions given are farther than {sep_tol} m from positions in "
                  "ocean loading coefficients file: \n"
                  + "\n".join(sndists)
                )

# Init an OceanFiles
OceanFiles()


def is_initialized(cb, item):
    """True if the item from common block cb is initialized on start-up.

    Also true if cb is module.

    Handcoded using output of get_initialized, checking where the data
    assignment actually happens.
    """
    return (
        cb == 'cmath'
        or cb in ('outputs', 'srcmod','datafiles')  # Modules
        or cb in 'units'    # File unit vars
        or cb == 'cticm' and item != 'a1tai'  # cctiu.f
        or cb == 'ut1cm' and item in ('centj', 'dj1900', 'dj2000')  # cut1m.f
        or cb == 'axocm' and item == 'richm'  # caxom.f
        or cb == 'wobcm' and item == 'len_wob_table'  # dinit.f
        or cb == 'stacm' and item == 'nflag'  # dstrt.f
        or cb == 'atmcm' and item in ('rf', 'rtrace', 'wetcon')  # catmm.f
        or cb == 'hmf2_coef'  # catmm.f
        or cb == 'stcomx' and item == 'km'  # cpepu.f
        or cb == 'xwahr'  # cnutu.f
        or cb == 'nutcmw'  # cnutu.f
        or cb == 'nutcm'  # cnutm.f
        or cb == 'tide_speed'  # cocem.f
    )
