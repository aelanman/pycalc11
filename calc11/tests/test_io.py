import os

import numpy as np
from numpy.testing import assert_allclose
from astropy import units as u
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename

from calc11.io import read_im, calc2dict


CRAB_CALC = get_pkg_data_filename(os.path.join(
    'data', 'B0531+21_CHIME_ARO10m_59153.45269019096.calc'))
CRAB_IM = get_pkg_data_filename(os.path.join(
    'data', 'B0531+21_CHIME_ARO10m_59153.45269019096.im'))


def test_read_calc():
    # For now, just read the meta data.
    calc = calc2dict(CRAB_CALC)
    assert calc['vex file'] == 'dummy.vex.obs'
    assert calc['eop'][0]['time'] == Time(59153, format='mjd')
    assert calc['eop'][1]['tai_utc'] == 37 * u.s


def test_read_im():
    """Simple test that delays get read in correctly."""
    im = read_im(CRAB_IM)
    assert len(im) == 6
    assert im.meta['telescope 0 name'] == 'CHIME'
    assert im.meta['calc program'] == 'DIFXCALC'
    assert_allclose(im['delay'][0, 0],
                    [1.890259772582982e+04,
                     1.428524766492706e-02,
                     -3.424386838841481e-05,
                     -1.265439323865132e-11,
                     1.509492047497244e-14,
                     2.741316895347732e-19])
    assert_allclose(im['delay'][-2, 1],
                    [1.596502423550063e+04,
                     -6.661073449020709e-01,
                     -2.727294852756184e-05,
                     5.903191697810408e-10,
                     1.213093754392361e-14,
                     -3.629969553109230e-20])
    assert np.all(im['scan'] == 0)
