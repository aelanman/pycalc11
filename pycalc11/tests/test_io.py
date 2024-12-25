import os

from numpy.testing import assert_allclose
from astropy import units as u
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename

from pycalc11.file_io import parse_calc, parse_im


CRAB_CALC = get_pkg_data_filename(
    os.path.join("data", "B0531+21_CHIME_ARO10m_59153.45269019096.calc")
)
CRAB_IM = get_pkg_data_filename(
    os.path.join("data", "B0531+21_CHIME_ARO10m_59153.45269019096.im")
)


def test_read_calc():
    # For now, just read the meta data.
    calc = parse_calc(CRAB_CALC)
    assert calc["vex file"] == "dummy.vex.obs"
    assert calc["eop"][0]["time"] == Time(59153, format="mjd")
    assert calc["eop"][1]["tai_utc"] == 37 * u.s


def test_read_im():
    """Simple test that delays get read in correctly."""
    im = parse_im(CRAB_IM)
    assert im["telescope"][0]["name"].item().strip() == "CHIME"
    assert im["calc program"] == "DIFXCALC"
    assert len(im["scan"]) == 1
    scan0 = im["scan"][0]
    assert len(scan0) == 3
    assert set(scan0.meta.keys()) == {"pointing src", "phs ctr", "main"}
    assert_allclose(
        scan0["delay"][0, 0, 0],
        [
            1.890259772582982e04,
            1.428524766492706e-02,
            -3.424386838841481e-05,
            -1.265439323865132e-11,
            1.509492047497244e-14,
            2.741316895347732e-19,
        ],
    )
    assert_allclose(
        scan0["delay"][2, 1, 1],
        [
            1.596502423550063e04,
            -6.661073449020709e-01,
            -2.727294852756184e-05,
            5.903191697810408e-10,
            1.213093754392361e-14,
            -3.629969553109230e-20,
        ],
    )
    assert scan0["delay"].info.meta["unit"] == u.us
