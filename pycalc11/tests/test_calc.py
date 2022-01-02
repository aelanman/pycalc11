
import numpy as np
from astropy import constants as const
from astropy import coordinates as ac
from astropy.time import Time
from numpy.testing import assert_allclose, assert_array_equal
from copy import deepcopy
import warnings

import pytest

from pycalc11.calcfile import make_calc
from pycalc11.interface import Calc, calc, OceanFiles
from pycalc11.utils import astropy_delay, astropy_delay_rate


def get_mod_state(fmod):
    odict = {}
    for k,v in fmod.__dict__.items():
        if (k.startswith('_')
                or all([dk.startswith('_') for
                       dk in v.__dict__])):
            continue
        for item, value in v.__dict__.items():
            odict[k + "_" + item] = deepcopy(value)
    return odict 

def compare_dicts(c0,c1, quiet=False):
    diffs = []
    for k in c1.keys():
        v0 = c0[k]
        v1 = c1[k]
        if np.atleast_1d(np.isreal(v0)).all():
            if not np.allclose(v0, v1):
                diffs.append(k)
        else:
            if not np.atleast_1d(v0 == v1).all():
                diffs.append(k)
    if len(diffs) == 0:
        return True
    else:
        if not quiet:
            print(diffs)
        return False

#-----------------------------
# Set up station locations and source/scan data
#-----------------------------
def make_params(nsrcs=305, duration_min=10):
    time = Time("2020-10-02T15:30:00.00", format="isot", scale='utc')
    gbo_loc = ac.EarthLocation.of_site("GBT")
    chime_loc = ac.EarthLocation.from_geodetic(lat=ac.Latitude('49d19m15.6s'), lon=ac.Longitude('119d37m26.4s'))
    wf_loc = ac.EarthLocation.from_geocentric(      # Haystack westford antenna
        x= 1492206.5970,
        y=-4458130.5170,
        z= 4296015.5320,
        unit='m',
    )

    ggao_loc = ac.EarthLocation.from_geocentric(    # Goddard, Greenbelt, Maryland
         x=  1130794.76936,
         y= -4831233.80170,
         z=  3994217.03883,
        unit='m',
    )

    srcs = ac.SkyCoord(
        az=np.random.uniform(0, 2 * np.pi, nsrcs),
        alt=np.random.uniform(np.radians(70), np.pi / 2, nsrcs),
        unit='rad',
        frame='altaz',
        obstime=time,
        location=gbo_loc
    )
    srcs = srcs.transform_to(ac.ICRS())

    telescope_positions = [chime_loc, gbo_loc, wf_loc, ggao_loc]
    telescope_names = ['chime', 'gbo', 'WESTFORD', 'GGAO7108']
    source_coords = [s for s in srcs]
    source_names = [f"src{si}" for si in range(nsrcs)]

    kwargs = {
        "station_coords": telescope_positions,
        "station_names": telescope_names,
        "source_coords": source_coords,
        "source_names": source_names,
        "time": time,
        "duration_min": duration_min,
    }

    return kwargs

@pytest.fixture(scope='session')
def params():
    # One realization of parameters
    return make_params()


def run_calc_2x(calcf_kwargs={}, calcfile_name="temp.calc", base_mode="geocenter",
               dry_atm=False, wet_atm=False):

    # 0 = calc file, 1 = keywords
    quantities = ['delay', 'delay_rate', 'partials', 'times']

    # Run with keywords
    ci = Calc(**calcf_kwargs, base_mode=base_mode, dry_atm=dry_atm, wet_atm=wet_atm)
    ci.run_driver()
    c1 = get_mod_state(calc)
    quants1 = {q: getattr(ci, q).copy() for q in quantities}

    # Run with calc file
    make_calc(**calcf_kwargs, ofile_name=calcfile_name)
    ci = Calc(calc_file=calcfile_name, base_mode=base_mode, dry_atm=dry_atm, wet_atm=wet_atm)
    ci.run_driver()
    del ci.calc_file
    c0 = get_mod_state(calc)
    quants0 = {q: getattr(ci, q).copy() for q in quantities}

    return quants0, quants1, c0, c1


def test_reset(params):
#   > Reset works properly
    stat0 = get_mod_state(calc)
    ci = Calc(**params)
    assert not compare_dicts(stat0, get_mod_state(calc), quiet=True)
    ci.reset()
    stat1 = get_mod_state(calc)

    res = compare_dicts(stat0, stat1)
    assert res


@pytest.mark.parametrize("atmo", [True, False])
def test_file_vs_kwds(atmo, tmpdir):
#   > Module attributes match when run from calcfile or from keywords
#   > Delays, delay rates, partials match between calcfile and keywords
    calcname = str(tmpdir.join("temp.calc"))
    params = make_params(nsrcs=5)
    quants0, quants1, c0, c1 = run_calc_2x(calcf_kwargs=params,
                                           calcfile_name=calcname,
                                           dry_atm=atmo, wet_atm=atmo)
    # TODO -- Fix so these compare properly
    #assert compare_dicts(c0, c1)
    for k, qi in quants1.items():
        if k.startswith('times'):
            assert np.all(np.abs(quants0[k] - qi) < np.timedelta64('1', 'us'))
        else:
            assert_allclose(quants0[k], qi, rtol=1e-1)


def test_calc_props(params):
#   > Calc attributes match passed params
    ci = Calc(**params)
    ra, dec = ci.src_coords
    srcs = ac.SkyCoord(ra, dec, unit='rad', frame='icrs')
    psrcs = ac.SkyCoord(params['source_coords'])
    # Source positions
    assert_allclose(srcs.ra.rad, psrcs.ra.rad)
    assert_allclose(srcs.dec.rad, psrcs.dec.rad)

    # Station positions
    stats = [ac.EarthLocation.from_geocentric(*c, unit='m') for c in ci.stat_coords.T]
    assert all([stats[ci] == loc for ci, loc in enumerate(params['station_coords'])])


def test_kwd_override(params):
#   > Providing calc file and keywords --> Keywords override .calc file
    calcfile_name = 'temp.calc'
    make_calc(**params, ofile_name=calcfile_name)

    # Remove a station
    stat_locs = params['station_coords'][:-1]
    stat_nmes = params['station_names'][:-1]

    ci = Calc(calc_file=calcfile_name, station_coords=stat_locs, station_names=stat_nmes)

    assert ci.stat_coords.shape[1] == len(stat_locs) == 3


@pytest.mark.skip("Can't set tolerance due to some systematic with astropy."
                  " To be explored further.")
@pytest.mark.parametrize("base_mode", ['geocenter', 'baseline'])
def test_ddr_vals(base_mode, params):
#   > D, DR are close to expectation from a rough astropy calculation.
#   (Eventually, set tolerances based on expected astropy accuracy.)
#  NOTE -- Currently failing due to some conditions with some source positions.
    pars = deepcopy(params)
    # remove CHIME, because the one ~3000 km NW -- SE baseline is creating some
    # weird conditions with sources too near the horizon.
    pars['station_names'] = pars['station_names'][1:]
    pars['station_coords'] = pars['station_coords'][1:]
    ci = Calc(**pars, base_mode=base_mode)
    ci.run_driver()
    st_nmes = pars['station_names']
    st_locs = pars['station_coords']
    srcs = ac.SkyCoord(params['source_coords'])
    nstat = len(st_nmes)
    t0 = params['time']
    gc = ac.EarthLocation.from_geocentric(0,0,0, unit='m')
    for ti, tt in enumerate(ci.times):
        if base_mode == 'geocenter':
            for si in range(nstat):
                ap_del = astropy_delay(srcs, tt, st_locs[si], gc) / 1e6
                ap_dlr = astropy_delay_rate(srcs, tt, st_locs[si], gc) / 1e6
                ra,dec = ci.src_coords
                assert_allclose(ci.delay[ti, 0, si, :], ap_del, rtol=0.05)
                assert_allclose(ap_dlr, ci.delay_rate[ti, 0, si, :], rtol=0.05)
        if base_mode == 'baseline':
            for s1i in range(nstat):
                for s2i in range(s1i, nstat-1):
                    ap_del = astropy_delay(srcs, tt, st_locs[s2i], st_locs[s1i]) / 1e6
                    ap_dlr = astropy_delay_rate(srcs, tt, st_locs[s2i], st_locs[s1i]) / 1e6
                    assert_allclose(ci.delay[ti, s1i, s2i, :], ap_del, rtol=0.05)
                    assert_allclose(ap_dlr, ci.delay_rate[ti, s1i, s2i, :], rtol=0.05)

def test_oc_warnings(params):
#   > Check that OceanFiles.check_sites warns about the correct sites being missing
#   > from ocean loading and ocean pole tide loading files.
#   > TODO I shouldn't need to use the full setup to run this test. Rewrite to just use dinit
#          to run the file parser in CALC (filling oc_coefs_found) before comparing with
#          OceanFiles.check_sites.

    wf_loc = ac.EarthLocation.from_geocentric(      # Haystack westford antenna
        x= 1492206.5970,
        y=-4458130.5170,
        z= 4296015.5320,
        unit='m',
    )

    ggao_loc = ac.EarthLocation.from_geocentric(    # Goddard, Greenbelt, Maryland
         x=  1130794.76936,
         y= -4831233.80170,
         z=  3994217.03883,
        unit='m',
    )
    mwa_loc = ac.EarthLocation.of_site("MWA")
    locs = [wf_loc, ggao_loc, mwa_loc]
    nmes = ["WESTFORD", "GGAO7108", "mwa"]
    params['station_names'] = nmes
    params['station_coords'] = locs
    with warnings.catch_warnings(record=True) as warning_list:
        ci = Calc(**params)

    msgs = [str(w.message) for w in warning_list]

    # lists of antenna names found by OceanFiles.check_sites
    # Start with full list and remove if in warning messages.
    oc_bad0 = []
    optl_bad0 = []
    for msg in msgs:
        ant = msg.split()[-1]
        if "ocean loading" in msg:
            oc_bad0.append(ant.ljust(8))
        elif "ocean pole tide" in msg:
            optl_bad0.append(ant.ljust(8))

    oc_bad0 = np.asarray(oc_bad0, dtype='<U8')
    optl_bad0 = np.asarray(optl_bad0, dtype='<U8')
    # calc.sitcm.oc_coefs_found[0, ii] == 1 if calc.calc_input.sites[ii] found in
    #   ocean pole tide loading coefficients file
    # calc.sitcm.oc_coefs_found[1, ii] == 1 if calc.calc_input.sites[ii] found in
    #   ocean loading data file.
    oc_bad = ci.stat_names[~calc.sitcm.oc_coefs_found[1][1:ci.nants+1].astype(bool)]
    optl_bad = ci.stat_names[~calc.sitcm.oc_coefs_found[0][1:ci.nants+1].astype(bool)]

    for l in [oc_bad0, oc_bad, optl_bad0, optl_bad]:
        l.sort()

    assert_array_equal(oc_bad0, oc_bad)
    assert_array_equal(optl_bad0, optl_bad)


def test_oc_dist_warnings():
    anmes = ['gbt', 'mwa', 'vla', 'ALMA', 'CHIME']
    site_locs = [ac.EarthLocation.of_site(nm) for nm in anmes]
    # Last one is assumed to be at DRAO, but given to calc as algo. Should give big offset.
    site_names = ['GBT-VLBA', 'mwa', 'VLA', 'ALMA', 'ALGOPARK']
    with warnings.catch_warnings(record=True) as warning_list:
        OceanFiles.check_sites(site_names=site_names, site_pos=site_locs)

    for w in (str(wn.message) for wn in warning_list):
        if "Station positions" in w:
            assert all(st in w for st in ['ALMA', 'ALGOPARK'])
            assert all(st not in w for st in ['GBT-VLBA', 'VLA'])

#   > Correct error is raised when trying to access without running driver
#   > Correct error is raised when initialized without params
#   > Time steps make sense 
