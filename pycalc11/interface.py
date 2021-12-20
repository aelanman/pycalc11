
import numpy as np
import copy
from astropy.time import TimeDelta
import astropy.coordinates as ac

from .utils import get_leap_seconds, iers_tab
from . import calc11 as calc


def is_initialized(cb, item):
    """True if the item from common block cb is initialized on start-up.

    Handcoded using output of get_initialized, checking where the data
    assignment actually happens.
    """
    return (
        cb == 'cmath'
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


def run(**kwargs):
    """Run calc11 subroutines."""
    for routine, args in kwargs.items():
        getattr(calc, routine)(*args)

def set_geocenter_telescope():
    calc.sitcm.numsit += 1
    calc.calc_input.sites[0] = "GEOCENTR"
    calc.calc_input.axis[0] = "AZEL"
    calc.sitcm.sitaxo[0] = 0.0    # Axis offset
    calc.sitcm.sitxyz[0, 0] = 0
    calc.sitcm.sitxyz[1, 0] = 0
    calc.sitcm.sitxyz[2, 0] = 0

def set_telescopes(telescope_positions, telescope_names):
    n_ants = len(telescope_names)
    calc.sitcm.numsit += n_ants

    # Zeroth site is geocenter
    # TODO -- Only set it to geocenter if in "geocenter" mode.
    #         (i.e., Allow for baseline mode.)
    set_geocenter_telescope()
    tnames = telescope_names
    tpos = telescope_positions
    for ti in range(n_ants):
        calc.calc_input.sites[ti+1] = tnames[ti].upper()
        calc.calc_input.axis[ti+1] = "AZEL"
        calc.sitcm.sitaxo[ti+1] = 0.0    # Axis offset
        calc.sitcm.sitxyz[0, ti+1] = tpos[ti].x.to_value('m')
        calc.sitcm.sitxyz[1, ti+1] = tpos[ti].y.to_value('m')
        calc.sitcm.sitxyz[2, ti+1] = tpos[ti].z.to_value('m')


def set_sources(source_coords, source_names):
    Nsrc = len(source_names)
    calc.alloc_source_arrays(Nsrc+1)
    calc.srcmod.numstr = len(source_names)
    for si, (coord, name) in enumerate(zip(source_coords, source_names)):
        calc.srcmod.radec[:, si] = [coord.ra.rad, coord.dec.rad]
        calc.srcmod.lnstar[:, si] = np.frombuffer(bytes(name.ljust(20), encoding='utf-8'), dtype=np.int16)
    calc.calc_input.numphcntr = Nsrc
    calc.calc_input.phcntr[:Nsrc] = range(1, Nsrc+1)
    calc.calc_input.phcntr[Nsrc:] = -1

    # The call below causes a variety of different errors on shutdown:
    #   > corrupted_size vs. prev_size
    #   > segmentation fault
    #   > double free or corruption (out)
    #   > (partway through source list) realloc(): invalid old size
    # calc doesn't actually use source names in difx mode, except in dscan when handling spacecraft.

    # Fixing this is a TODO for the future.
    #calc.set_srcname(Nsrc)


def set_scan(time, duration_min):
    """Only one scan is supported."""
    jstart = time
    jstop = time + TimeDelta(duration_min * 60, format='sec')
    calc.ut1cm.xintv[0] = jstart.jd
    calc.ut1cm.xintv[1] = jstop.jd
    calc.ut1cm.intrvl[:, 0] = list(jstart.ymdhms)[:5]    # intrvl[:, 1] unused
    calc.calc_input.numscans = 1
    calc.calc_input.scanid = "Scan001".rjust(10)
    calc.calc_input.scandur = duration_min * 60  # Seconds
    calc.calc_input.pointingsrc = 1


def set_eops(time):
    _times = time + TimeDelta(range(2), format='jd')
    jd = [np.floor(tt.jd) for tt in _times]
    tai_utc = [get_leap_seconds(tt) for tt in _times]
    ut1_utc = [tt.delta_ut1_utc.reshape(1)[0] for tt in _times]
    xy = [[z.to_value('milliarcsecond') for z in iers_tab.pm_xy(tt)] for tt in _times]

    Nerp = 2
    Uintv = 1
    xjds = jd[0] - (Nerp * Uintv) / 2
    calc.wobcm.wobif = [xjds, Uintv, Nerp]     # Wobble info
    calc.wobcm.xywob[:, :2] = np.transpose(xy)
    calc.ut1cm.ut1if = [xjds, Uintv, Nerp, 1.0]
    calc.ut1cm.ut1pt[:2] = [tai_utc[ti] - ut1_utc[ti] for ti in range(2)]
    calc.calc_input.xleap_sec = tai_utc[0]


def set_calcfile(calc_file_name):
    calc.contrl.calc_file_name = calc_file_name.ljust(128)


def set_cl():
    """
    Set command line arguments

    No options for now.

    TODO -- Enable baseline mode.
    """
    calc.mode.c_mode = 'difx  '
    calc.contrl.i_out = 0
    calc.contrl.im_out = 1
    calc.contrl.near_far = 'Far-field '
    calc.nfo.spoffset = 'NoOffset'
    calc.contrl.base_mode = 'geocenter '
    calc.contrl.l_time = 'dont-solve'
    calc.contrl.atmdr = 'Add-dry   '
    calc.contrl.atmwt = 'Add-wet   '
    calc.contrl.verbose = 1
    calc.contrl.overwrite = 'no  '
    calc.contrl.uvw = 'exact '
    calc.contrl.d_interval = 24.
    calc.contrl.epoch2m = (120.0001/calc.contrl.d_interval) + 1
    calc.contrl.numjobs = 0


def reset():
    """Reset all common block items that were not initialized at startup."""
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


def delay():
    return calc.outputs.delay_f

def delay_rate():
    return calc.outputs.rate_f

def partials():
    return calc.outputs.partials_f

def nsrcs():
    return calc.calc_input.numphcntr

def nants():
    return calc.sitcm.numsit - 1
