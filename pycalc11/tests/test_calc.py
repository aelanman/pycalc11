
import numpy as np
from astropy import coordinates as ac
from astropy.time import Time, TimeDelta
from astropy import units as un
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose, assert_array_equal
from copy import deepcopy
import warnings

import pytest

from pycalc11.calcfile import make_calc
from pycalc11.interface import Calc, calc, OceanFiles
from pycalc11.utils import astropy_delay, astropy_delay_rate
from pycalc11.runner import run_difxcalc, difxcalc
from pycalc11.imfile import CalcReader


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


    anames = ['gbt', 'vla', 'ALMA', 'mwa', 'chime']
    site_locs = [ac.EarthLocation.of_site(nm) for nm in anames]
    site_names = ['GBT-VLBA', 'VLA_E9', 'ALMA', 'MWA', 'CHIME 0']


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

    # Add the site NRAO85, since this has multiple unique entries in ocean loading
    # that are distinguished by unique codes, while all having the same name.
    nrao85_3_loc = ac.EarthLocation.from_geodetic(
        lon= 280.1566,
        lat= 38.4296,
        height=785.847,
    )
    nrao85_1_loc = ac.EarthLocation.from_geodetic(
        lon= 280.1718,
        lat= 38.4359,
        height=807.687,
    )

    site_names = ["WESTFORD", "GGAO7108", "NRAO85 3", "NRAO85 1"] + site_names
    site_locs = [wf_loc, ggao_loc, nrao85_3_loc, nrao85_1_loc] + site_locs

    srcs = ac.SkyCoord(
        az=np.random.uniform(0, 2 * np.pi, nsrcs),
        alt=np.random.uniform(np.radians(70), np.pi / 2, nsrcs),
        unit='rad',
        frame='altaz',
        obstime=time,
        location=gbo_loc
    )
    srcs = srcs.transform_to(ac.ICRS())

    kwargs = {
        "station_coords": site_locs,
        "station_names": site_names,
        "source_coords": srcs,
        "start_time": time,
        "duration_min": duration_min,
    }

    return kwargs

@pytest.fixture(scope='session')
def params_all():
    # One realization of parameters
    return make_params()

@pytest.fixture(scope='session')
def params_vlbi():
    # Keep stations in CALC database
    # (i.e., remove MWA and CHIME from lists).
    par = make_params()
    par['station_coords'] = par['station_coords'][:-2]
    par['station_names'] = par['station_names'][:-2]
    return par


def run_calc_2x(calcf_kwargs={}, calcfile_name="temp.calc", base_mode="geocenter",
               dry_atm=False, wet_atm=False):

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


def test_reset(params_vlbi):
#   Reset works properly
    stat0 = get_mod_state(calc)
    ci = Calc(**params_vlbi)
    assert not compare_dicts(stat0, get_mod_state(calc), quiet=True)
    ci._reset()
    stat1 = get_mod_state(calc)

    res = compare_dicts(stat0, stat1)
    assert res


@pytest.mark.parametrize("atmo", [True, False])
def test_file_vs_kwds(atmo, params_vlbi, tmpdir):
#   Module attributes match when run from calcfile or from keywords
#   Delays, delay rates, partials match between calcfile and keywords
    calcname = str(tmpdir.join("temp.calc"))
    quants0, quants1, c0, c1 = run_calc_2x(calcf_kwargs=params_vlbi,
                                           calcfile_name=calcname,
                                           dry_atm=atmo, wet_atm=atmo)

    # TODO -- Fix so these compare properly
    #assert compare_dicts(c0, c1)

    quantities = ['delay', 'delay_rate', 'partials', 'times']
    atols = {
        'delay': 0.01 * un.ns,
        'delay_rate': 1e-11 * un.s * un.Hz,
    }

    assert_quantity_allclose(quants0['delay'], quants1['delay'], atol=0.01*un.ns)
    assert_quantity_allclose(quants0['delay_rate'], quants1['delay_rate'], atol=1e-11*un.s*un.Hz)

    for nm in quants0['partials'].dtype.names:
        if 'dDR' in nm:
            atol = atols['delay_rate'] / un.rad
        else:
            atol = atols['delay'] / un.rad
        assert_quantity_allclose(quants0['partials'][nm], quants1['partials'][nm], atol=atol)

    assert np.all(quants0['times'] == quants1['times'])


def test_calc_props(params_vlbi):
#   Calc attributes match passed params
    ci = Calc(**params_vlbi)
    srcs = ci.source_coords
    psrcs = ac.SkyCoord(params_vlbi['source_coords'])

    # Source positions
    assert_allclose(srcs.ra.rad, psrcs.ra.rad)
    assert_allclose(srcs.dec.rad, psrcs.dec.rad)

    # Station positions
    stats = [ac.EarthLocation.from_geocentric(*c, unit='m') for c in ci.station_coords.T]
    assert all([stats[ci] == loc for ci, loc in enumerate(params_vlbi['station_coords'])])


@pytest.mark.filterwarnings("ignore:No ocean")
def test_kwd_override(params_all, tmpdir):
#   Providing calc file and keywords --> Keywords override .calc file
    calcfile_name = str(tmpdir.join("temp.calc"))
    make_calc(**params_all, ofile_name=calcfile_name)

    # Remove a station
    stat_locs = params_all['station_coords'][:-1]
    stat_nmes = params_all['station_names'][:-1]

    ci = Calc(calc_file=calcfile_name, station_coords=stat_locs, station_names=stat_nmes)
    assert ci.station_coords.shape[1] == len(stat_locs) == 8


@pytest.mark.skip("Can't set tolerance due to some systematic with astropy."
                  " To be explored further.")
@pytest.mark.parametrize("base_mode", ['geocenter', 'baseline'])
def test_ddr_vals(base_mode, params_vlbi):
#   > D, DR are close to expectation from a rough astropy calculation.
#   (Eventually, set tolerances based on expected astropy accuracy.)
    pars = deepcopy(params_vlbi)
    # remove CHIME, because the one ~3000 km NW -- SE baseline is creating some
    # weird conditions with sources too near the horizon.
    pars['station_names'] = pars['station_names'][1:]
    pars['station_coords'] = pars['station_coords'][1:]
    ci = Calc(**pars, base_mode=base_mode, dry_atm=False, wet_atm=False)
    ci.run_driver()
    st_nmes = pars['station_names']
    st_locs = pars['station_coords']
    srcs = ac.SkyCoord(params_vlbi['source_coords'])
    nstat = len(st_nmes)
    gc = ac.EarthLocation.from_geocentric(0,0,0, unit='m')
    for ti, tt in enumerate(ci.times):
        if base_mode == 'geocenter':
            for si in range(nstat):
                ap_del = astropy_delay(srcs, tt, st_locs[si], gc)
                ap_dlr = astropy_delay_rate(srcs, tt, st_locs[si], gc)
                ra,dec = ci.source_coords
                assert_quantity_allclose(ci.delay[ti, 0, si, :], ap_del, rtol=0.05)
                assert_quantity_allclose(ap_dlr, ci.delay_rate[ti, 0, si, :], rtol=0.05)
        if base_mode == 'baseline':
            for s1i in range(nstat):
                for s2i in range(s1i, nstat-1):
                    ap_del = astropy_delay(srcs, tt, st_locs[s2i], st_locs[s1i])
                    ap_dlr = astropy_delay_rate(srcs, tt, st_locs[s2i], st_locs[s1i])
                    assert_quantity_allclose(ci.delay[ti, s1i, s2i, :], ap_del, rtol=0.05)
                    assert_quantity_allclose(ap_dlr, ci.delay_rate[ti, s1i, s2i, :], rtol=0.05)


def test_oc_warnings(params_all):
#   Check that OceanFiles.check_sites warns about the correct sites being missing
#   from ocean loading and ocean pole tide loading files.

#   TODO I shouldn't need to use the full setup to run this test. Rewrite to just use dinit
#          to run the file parser in CALC (filling oc_coefs_found) before comparing with
#          OceanFiles.check_sites.

    with warnings.catch_warnings(record=True) as warning_list:
        ci = Calc(**params_all)

    msgs = [str(w.message) for w in warning_list]

    # lists of antenna names not found by OceanFiles.check_sites
    oc_bad0 = []
    optl_bad0 = []
    for msg in msgs:
        ants = [a.strip().ljust(8)[:8] for a in msg.split('for ')[1].split(',')]
        if "ocean loading" in msg:
            oc_bad0 += ants
        if "ocean pole tide" in msg:
            optl_bad0 += ants

    oc_bad0 = np.asarray(oc_bad0, dtype='<U8')
    optl_bad0 = np.asarray(optl_bad0, dtype='<U8')
    # calc.sitcm.oc_coefs_found[0, ii] == 1 if calc.calc_input.sites[ii] found in
    #   ocean pole tide loading coefficients file
    # calc.sitcm.oc_coefs_found[1, ii] == 1 if calc.calc_input.sites[ii] found in
    #   ocean loading data file.
    oc_bad = ci.station_names[~calc.sitcm.oc_coefs_found[1][1:ci.nants+1].astype(bool)]
    optl_bad = ci.station_names[~calc.sitcm.oc_coefs_found[0][1:ci.nants+1].astype(bool)]

    for l in [oc_bad0, oc_bad, optl_bad0, optl_bad]:
        l.sort()

    assert_array_equal(oc_bad0, oc_bad)
    assert_array_equal(optl_bad0, optl_bad)


def test_oc_dist_warnings(params_all):
    # Last one is assumed to be at DRAO, but given to calc as algo. Should give big offset.
    with warnings.catch_warnings(record=True) as warning_list:
        OceanFiles.check_sites(
            site_names=params_all['station_names'], site_pos=params_all['station_coords']
        )

    print([str(wn.message) for wn in warning_list])
    for w in (str(wn.message) for wn in warning_list):
        if "Station positions" in w:
            assert all([st in w for st in ['VLA_E9', 'ALMA']])
            assert all([st not in w for st in ['GBT-VLBA']])


def test_partials_calc(params_vlbi):
    # Get partial derivatives. Compare to approximate calculation from delays.

    dth = 1 * un.arcsec

    c_ra = np.pi / 6 * un.rad
    c_dec = np.pi / 4 * un.rad

    srcs = ac.SkyCoord(
        ra=[c_ra, c_ra, c_ra, c_ra - dth, c_ra + dth],
        dec=[c_dec, c_dec - dth, c_dec + dth, c_dec, c_dec],
        unit='rad', frame='icrs'
    )

    pars = deepcopy(params_vlbi)
    pars['source_coords'] = srcs

    ci = Calc(**pars, dry_atm=False, wet_atm=False)
    ci.run_driver()

    delays = ci.delay
    partials = ci.partials

    # Partial derivatives at the zeroth source position (our test point)
    partials = partials[..., 0]

    dtau_dra = (delays[..., 4] - delays[..., 3]) / (2 * dth)
    dtau_ddec = (delays[..., 2] - delays[..., 1]) / (2 * dth)

    # dtau_dra and dtau_ddec axes = (time, ant1, ant2) where ant1 is just geocenter
    atol = 1e-6 * un.s / un.rad
    assert_quantity_allclose(dtau_dra[:, 0, :], partials['dD_dRA'][:, 0, :], atol=atol)
    assert_quantity_allclose(dtau_ddec[:, 0, :], partials['dD_dDEC'][:, 0, :], atol=atol)


def test_errors(params_all):
    # Correct error is raised when trying to access without running driver
    # Correct error is raised when initialized without params
    with pytest.raises(ValueError, match="If calc_file is not set"):
        Calc()

    with pytest.raises(ValueError, match="Need to rerun adrivr"):
        ci = Calc(**params_all)
        ci.delay


@pytest.mark.skipif(difxcalc is None, reason='difxcalc must be installed for this test')
def test_compare_to_difxcalc(params_vlbi, tmpdir):
    # Compare running Calc here with results from a separately-installed difxcalc.

    # Difxcalc can't handle more than 300 sources
    # Just keep ten sources for this test.
    params_vlbi['source_coords'] = params_vlbi['source_coords'][:10]

    params_vlbi['duration_min'] = 20

    quantities = ['delay', 'delay_rate', 'partials', 'times']
    ci = Calc(**params_vlbi, base_mode='geocenter', dry_atm=False, wet_atm=False)
    ci.run_driver()
    quants1 = {q: getattr(ci, q).copy() for q in quantities}

    calcname = str(tmpdir.join("temp.calc"))
    make_calc(**params_vlbi, ofile_name=calcname)

    impath = run_difxcalc(calcname, verbose=False)

    im_obj = CalcReader(impath)
    times = Time(quants1['times'])

    delays = -quants1['delay'].to_value('us')# * (-1e6)
    im_dels = np.zeros_like(delays)
    for ti, tt in enumerate(times):
        for si in range(10):
            for ai in range(7):
                # si + 1 because the first source (in difxcalc) is the same as 0th
                # This is because the zeroth source is treated as the pointing center.
                im_dels[ti, 0, ai, si] = im_obj.delay(ai, tt, si + 1)

    # Tolerance of picoseconds
    assert_allclose(delays, im_dels, atol=1e-6)


def test_coef_vals(params_vlbi):
    # Check that coefficient values found by OceanFiles match those loaded by calc.
    with warnings.catch_warnings():
        ci = Calc(**params_vlbi)

    nants = len(params_vlbi['station_coords'])
    names = params_vlbi['station_names']
    inds = [OceanFiles.get_oc_ind(n)[0] for n in names]

    # Phases are given as deg in the file, but converted to rad in CALC.
    coefs = np.rollaxis(np.stack((
        calc.sitcm.sitoam[..., 1:nants + 1],    # Vertical amplitude [m]
        calc.sitcm.sithoa[:, 0, 1:nants + 1],   # horizontal amplitudes [m]
        calc.sitcm.sithoa[:, 1, 1:nants + 1],
        np.degrees(calc.sitcm.sitoph[..., 1:nants + 1]),  # Vertical phase [deg]
        np.degrees(calc.sitcm.sithop[:, 0, 1:nants + 1]), # horizontal phases [deg]
        np.degrees(calc.sitcm.sithop[:, 1, 1:nants + 1]),
    )), 2, 0)

    # Compare the above to the corresponding values in OceanFiles.oc_data['coefs']
    assert_allclose(coefs, OceanFiles.oc_data['coefs'][inds])

    # Now compare ocean pole tide loading coefs
    keys = ['urr', 'uri', 'unr', 'uni', 'uer', 'uei']
    inds = [OceanFiles.get_optl_ind(n)[0] for n in names]
    coefs = calc.sitcm.optl6[:, 1:nants+1]
    comp = np.zeros_like(coefs)
    for ii in range(nants):
        comp[:, ii] = OceanFiles.optl_data[keys][inds][ii].tolist()

    assert_allclose(comp, coefs)


@pytest.mark.parametrize("kv", make_params(nsrcs=30).items())
def test_change_quantities(params_vlbi, kv):
    ci = Calc(**params_vlbi)
    ci.run_driver()
    key, val = kv
    setattr(ci, key, val)
    assert ci._rerun

    with pytest.raises(ValueError, match="Need to rerun adrivr before accessing data. "):
        print(ci.delay)

    if key == 'source_coords':
        ci.run_driver()
        assert ci.nsrcs == 30
        assert ci.delay.shape[-1] == 30


def test_compare_baseline_geocenter(params_vlbi):
    ci = Calc(**params_vlbi)
    ci.run_driver()
    delays_geoc = ci.delay.copy()
    ci.base_mode = "baseline"
    ci.run_driver()
    delays_bl = ci.delay.copy()
    delays_bl_fromgeo = delays_geoc[:, 0, [ci.ant2_ind], :] - delays_geoc[:, 0, ci.ant1_ind, None, :]
    # delays_bl has zeros for ant2 < ant1
    for a1, a2 in [(a1, a2) for a1 in ci.ant1_ind for a2 in ci.ant2_ind]:
        if a2 < a1:
            delays_bl_fromgeo[:, a1, a2, :] = 0.0

    # !! NOTE Seeing differences on scales of 1 to 12 ns between two methods. Investigate!!


## TODO
#   Test that epochs cover the requested scan times
#   Fix tests with astropy comparison
