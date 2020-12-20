
import numpy as np
import warnings
from datetime import datetime
from astropy.time import Time, TimeDelta
from astropy.utils import data
from astropy.utils import iers
import astropy.coordinates as ac
from astropy import units as un
from astropy.constants import c as speed_of_light

import calc11 as _calc11


iers_tab = iers.earth_orientation_table.get()

def _get_leap_seconds(tobj):
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


def _reset_calc_params():
    ## TODO 
    pass


def set_calc_params(telescope_positions, telescope_names, source_coords,
                    source_names, time, duration_min):
    """
    Set the common parameters in _calc11 that would be read from a .calc file.
    Also set command line parameters
    """

    _calc11.mode.c_mode = 'difx  '
    _calc11.contrl.i_out = 0
    _calc11.contrl.im_out = 1
    _calc11.contrl.near_far = 'Far-field '
    _calc11.nfo.spoffset = 'NoOffset'
    _calc11.contrl.base_mode = 'geocenter '
    _calc11.contrl.l_time = 'dont-solve'
    _calc11.contrl.atmdr = 'Add-dry   '
    _calc11.contrl.atmwt = 'Add-wet   '
    _calc11.contrl.verbose = 1 
    _calc11.contrl.overwrite = 'no  '
    _calc11.contrl.uvw = 'exact '
    _calc11.contrl.d_interval = 24.
    _calc11.contrl.epoch2m = (120.0001/_calc11.contrl.d_interval) + 1
    _calc11.contrl.numjobs = 0 


    min_per_day = 60 * 24
    jstart, jstop = time.jd, time.jd + duration_min / min_per_day

    _calc11.ut1cm.xintv[0] = time.jd
    _calc11.ut1cm.xintv[1] = time.jd + duration_min / min_per_day
    _calc11.ut1cm.intrvl[:, 0] = list(time.ymdhms)[:5]    # intrvl[:, 1] doesn't get used.

    # EOP  -- will need to review this.
    # Can handle at most 20 EOPs
    _times = time + TimeDelta(range(2), format='jd')
    jd = [np.floor(tt.jd) for tt in _times]
    tai_utc = [_get_leap_seconds(tt) for tt in _times]
    ut1_utc = [tt.delta_ut1_utc.reshape(1)[0] for tt in _times]
    xy = [[z.to_value('milliarcsecond') for z in iers_tab.pm_xy(tt)] for tt in _times]

    _calc11.wobcm.wobif = [jd[0], jd[1] - jd[0], 2]     # Wobble info
    _calc11.wobcm.xywob[:, :2] = np.transpose(xy)

    _calc11.ut1cm.ut1if = [jd[0], jd[1] - jd[0], 2, 1.0]  # Leap seconds etc.
    _calc11.ut1cm.ut1pt[:2] = [tai_utc[ti] - ut1_utc[ti] for ti in range(2)]


    # Sources
    if len(source_names) > 300:
        raise ValueError("No more than 300 sources may be used at once.")
    _calc11.strcm.numstr = len(source_names)
    for si, (coord, name) in enumerate(zip(source_coords, source_names)):
        _calc11.strcm.radec[:, si] = [coord.ra.rad, coord.dec.rad]


    # Telescopes
    n_ants = len(telescope_names)
    _calc11.sitcm.numsit = n_ants + 1
    # Zeroth site is geocenter
    telescope_names.insert(0, "GEOCENTR")
    telescope_positions.insert(0, ac.EarthLocation.from_geocentric(0, 0, 0, unit='m'))
    for ti in range(n_ants + 1):
        _calc11.calc_input.sites[ti] = telescope_names[ti].upper()
        _calc11.calc_input.axis[ti] = "AZEL"
        _calc11.sitcm.sitaxo[ti] = 0.0    # Axis offset
        _calc11.sitcm.sitxyz[0, ti] = telescope_positions[ti].x.to_value('m')
        _calc11.sitcm.sitxyz[1, ti] = telescope_positions[ti].y.to_value('m')
        _calc11.sitcm.sitxyz[2, ti] = telescope_positions[ti].z.to_value('m')

    # Scans
    _calc11.calc_input.numscans = 1
    _calc11.calc_input.scanid = "Scan001".rjust(10)
    _calc11.calc_input.scanstr = 0.0
    _calc11.calc_input.scandur = duration_min * 60  # Seconds
    _calc11.calc_input.pointingsrc = 0
    _calc11.calc_input.numphcntr = 1
    _calc11.calc_input.phcntr[0] = 1
    _calc11.calc_input.phcntr[-1] = 1



_calc11.start(1,1)

time = Time.now()
gbo_loc = ac.EarthLocation.of_site("GBT")
chime_loc = ac.EarthLocation.from_geodetic(lat=ac.Latitude('49d19m15.6s'), lon=ac.Longitude('119d37m26.4s'))
crab = ac.SkyCoord('05:34:31.9383', '22:00:52.175', unit=(un.hourangle, un.deg), frame='icrs')

duration_min = 60 * 8


telescope_positions = [chime_loc, gbo_loc]
telescope_names = ['chime', 'gbo']
source_coords = [crab]
source_names = ['src']
duration_min = 60 * 8

set_calc_params(telescope_positions, telescope_names, source_coords,
          source_names, time, duration_min)

#import IPython; IPython.embed()
_calc11.dinitl(1,)
_calc11.drivr(0,0)
