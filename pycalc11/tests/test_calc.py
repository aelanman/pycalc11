import os
from collections import defaultdict

import numpy as np
from numpy.testing import assert_allclose
from astropy import constants as const
from astropy import coordinates as ac
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename

from pycalc11.io import parse_im
from pycalc11.funcs import get_delay
from pycalc11.calcfile import make_calc
from pycalc11 import calc11


CRAB_CALC = get_pkg_data_filename(os.path.join(
    'data', 'B0531+21_CHIME_ARO10m_59153.45269019096.calc'))
CRAB_IM = get_pkg_data_filename(os.path.join(
    'data', 'B0531+21_CHIME_ARO10m_59153.45269019096.im'))


def get_initialized(calc):
    """Get values from common blocks which are non zero on startup.

    Obviously, for a useful result, this should be run right after
    loading calc11.
    """
    initialized = defaultdict(dict)
    for key, part in calc.__dict__.items():
        if (key.startswith('_')
                or all(dk.startswith('_') for
                       dk in part.__dict__)):
            continue

        for pk, value in part.__dict__.items():
            if not np.all(value == (b'' if value.dtype.kind == 'S' else 0)):
                initialized[key][pk] = value

    return initialized


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


def reset_calc(calc):
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


def check_initialized(calc):
    """Check only common-block items initialized at startup are non-zero."""
    for key, part in calc.__dict__.items():
        # Select only common blocks
        if (key.startswith('_')
                or all(dk.startswith('_') for
                       dk in part.__dict__)):
            continue
        for item, value in part.__dict__.items():
            zero = b'' if value.dtype.kind == 'S' else 0
            all_zero = np.all(value == zero)
            initialized = is_initialized(key, item)
            if (initialized and all_zero
                    or not initialized and not all_zero):
                return False

    return True


def test_initialization():
    """This should be the first test run on calc!!

    Check that get_initialized returns items which are all caught
    by is_initialized.
    """
    initial = get_initialized(calc11)
    for key, cm in initial.items():
        for item, value in cm.items():
            assert is_initialized(key, item)
            zero = b'' if value.dtype.kind == 'S' else 0
            assert not np.all(value == zero)

    # Now check is_initialized itself, both is and is not.
    assert check_initialized(calc11)
    # And check that reset_calc does not change that initialized state.
    reset_calc(calc11)
    assert check_initialized(calc11)
    # Check that no initialized values were changed by reset.
    for key, cm in initial.items():
        for item, value in cm.items():
            assert np.all(getattr(getattr(calc11, key), item) == value)


class CALCTestBase:
    @classmethod
    def setup_class(cls):
        reset_calc(calc11)

    @classmethod
    def teardown_class(cls):
        reset_calc(calc11)


class TestFromCalcFile(CALCTestBase):
    @classmethod
    def setup_class(cls):
        super().setup_class()
        # Load comparison data from direct run of executable.
        cls.im = parse_im(CRAB_IM)

    @classmethod
    def teardown_class(cls):
        # Test wiping in the process.
        super().teardown_class()
        # Explicit tests.
        assert calc11.mode.c_mode == b''
        assert calc11.contrl.im_out == 0
        assert np.all(calc11.sitcm.sitxyz == 0)
        assert np.all(calc11.out_c.delay_f == 0.)
        # Should cover above, but just in case.
        assert check_initialized(calc11)

    def do_calc(self, **kwargs):
        reset_calc(calc11)
        # Initialize following dmain.
        calc11.mode.c_mode = 'difx  '
        # Follow get_cl in dstrt.
        calc11.contrl.i_out = 0
        calc11.contrl.im_out = 1
        calc11.nfo.spoffset = 'NoOffset'
        calc11.contrl.base_mode = 'geocenter '
        calc11.contrl.l_time = 'dont-solve'
        calc11.contrl.atmdr = 'Add-dry   '
        calc11.contrl.atmwt = 'Add-wet   '
        calc11.contrl.verbose = 0
        calc11.contrl.overwrite = 'no  '
        calc11.contrl.uvw = 'exact '
        calc11.contrl.d_interval = 24.
        calc11.contrl.epoch2m = (120.0001/calc11.contrl.d_interval) + 1
        calc11.contrl.numjobs = 0
        # Read inputs.
        calc11.contrl.calc_file_name = f"{CRAB_CALC:<128s}"
        for routine, args in kwargs.items():
            getattr(calc11, routine)(*args)

    def test_setup(self):
        self.do_calc()
        assert calc11.mode.c_mode == b'difx  '
        assert calc11.contrl.d_interval == 24.

    def test_dstart(self):
        self.do_calc(dstart=(1, 1), dget_input=(1,))
        assert np.all(calc11.sitcm.sitxyz[:, 0] == 0.)  # Center of Earth
        assert_allclose(calc11.sitcm.sitxyz[:, 1],
                        [-2059159.52756903,
                         -3621260.06650439,
                         4814325.37572648])
        assert_allclose(calc11.strcm.radec[:, 0], [1.45967254, 0.38422539])

    def test_dinitl(self):
        # Initialize constants.
        self.do_calc(dstart=(1, 1), dinitl=(1,))
        # CALC11 does not have IAU value.
        assert_allclose(calc11.cphys.gmplanet[3], const.GM_jup.value,
                        rtol=0.0003)

    def test_dscan(self):
        self.do_calc(dstart=(1, 1), dinitl=(1,), dscan=(1, 1))
        # initializes near_far = 'Far-field'
        assert calc11.calc_input.phcntr[0] == 1
        assert calc11.calc_input.phcntr[1] == -1

    def test_calc(self):
        self.do_calc(dstart=(1, 1), dinitl=(1,), dscan=(1, 1), ddrivr=(1, 1))
        assert np.all(calc11.out_c.iymdhms_f[0] == [2020, 10, 31, 10, 50, 0])
        assert np.all(calc11.out_c.iymdhms_f[:, -1] == [0, 24, 48, 12, 36, 0])
        # Note: can only compare first and last point for both telescopes
        # with zero point of first and second polynomial.
        assert_allclose(calc11.out_c.delay_f[0, 0, :2, 0],
                        -1e-6*self.im['scan'][0]['delay'][0, 0, :, 0],
                        atol=0, rtol=4e-16)
        assert_allclose(calc11.out_c.delay_f[-1, 0, :2, 0],
                        -1e-6*self.im['scan'][0]['delay'][1, 0, :, 0],
                        atol=0, rtol=4e-16)


def test_delay_calc(tmpdir):
    # Run with and without using a .calc file as intermediary.
    # Compare delays

    time = Time("2020-10-02T15:30:00.00", format="isot", scale='utc')
    gbo_loc = ac.EarthLocation.of_site("GBT")
    chime_loc = ac.EarthLocation.from_geodetic(lat=ac.Latitude('49d19m15.6s'), lon=ac.Longitude('119d37m26.4s'))

    nsrcs = 10
    srcs = ac.SkyCoord(
        az=np.random.uniform(0, 2 * np.pi, nsrcs),
        alt=np.random.uniform(0, np.pi / 2, nsrcs),
        unit='rad',
        frame='altaz',
        obstime=time,
        location=gbo_loc
    )
    srcs = srcs.transform_to(ac.ICRS)
    
    duration_min = 60 * 8
    
    telescope_positions = [chime_loc, gbo_loc]
    telescope_names = ['chime', 'gbo']
    source_coords = [s for s in srcs]
    source_names = [f"src{si}" for si in range(nsrcs)]
    duration_min = 60 * 8

    calcfile = str(tmpdir.join('temp.calc'))
    make_calc(telescope_positions, telescope_names, source_coords,
              source_names, time, duration_min, ofile_name=calcfile)
    
    
    delays0 = get_delay(telescope_positions, telescope_names, source_coords, source_names, time, duration_min).copy()
    delays1 = get_delay(calc_file_name=calcfile)
    
    assert np.allclose(delays0, delays1)
