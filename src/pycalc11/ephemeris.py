""Ephemeris computation using jplephem SPK kernels."""

import numpy as np


# JPL planet indices used by PEP's SPLANET/GPLANET arrays:
#   1=Mercury, 2=Venus, 3=Mars, 4=Jupiter, 5=Saturn, 6=Uranus, 7=Neptune
# Corresponding SPK segment center-target pairs (SSB -> body):
_PLANET_SPK_TARGETS = [1, 2, 4, 5, 6, 7, 8]

# Conversion factor: km/day -> m/s
_KM_PER_DAY_TO_M_PER_S = 1.0e3 / 86400.0


def _compute_earth_ssb(kernel, tdb_jd):
    """Compute Earth SSB position (km) and velocity (km/day)."""
    emb_pos, emb_vel = kernel[0, 3].compute_and_differentiate(tdb_jd)
    e_from_emb_pos, e_from_emb_vel = kernel[3, 399].compute_and_differentiate(tdb_jd)
    return emb_pos + e_from_emb_pos, emb_vel + e_from_emb_vel


def fill_ephem_arrays(kernel, tdb_jd_array, calc_module):
    """Pre-compute ephemeris and fill Fortran ephcom arrays for a set of time steps.

    Computes all the quantities that the Fortran PEP subroutine normally produces
    (Earth, Sun, Moon, and planet positions/velocities) using a jplephem SPK kernel,
    and writes them into the calc11.ephcom module arrays.

    Parameters
    ----------
    kernel : jplephem.spk.SPK
        An open JPL SPK ephemeris kernel.
    tdb_jd_array : array_like
        TDB Julian dates for each time step. Length must not exceed
        calc_module.ephcom.max_eph_steps.
    calc_module : module
        The compiled calc11 Fortran module (typically ``pycalc11.calc11``).
    """
    tdb_jd = np.atleast_1d(np.asarray(tdb_jd_array, dtype=np.float64))
    n_steps = len(tdb_jd)

    eph = calc_module.ephcom
    if n_steps > eph.ext_earth.shape[2]:
        raise ValueError(
            f"Number of time steps ({n_steps}) exceeds Fortran array limit "
            f"({eph.ext_earth.shape[2]})"
        )

    # --- Earth (barycentric) ---
    earth_pos, earth_vel = _compute_earth_ssb(kernel, tdb_jd)
    # Earth acceleration: numerical diff of velocity at +/-1 second (matching PEP)
    dt = 1.0 / 86400.0  # 1 second in days
    _, vel_m1 = _compute_earth_ssb(kernel, tdb_jd - dt)
    _, vel_p1 = _compute_earth_ssb(kernel, tdb_jd + dt)
    # PEP does: accel = (v_plus_km_s - v_minus_km_s) * 1e3 / 2
    # jplephem gives km/day, so convert to km/s first
    earth_accel_m_s2 = (vel_p1 - vel_m1) / 86400.0 * 1.0e3 / 2.0

    earth_pos_m = earth_pos * 1.0e3  # km -> m
    earth_vel_m_s = earth_vel * _KM_PER_DAY_TO_M_PER_S

    # --- Sun (barycentric) ---
    sun_pos, sun_vel = kernel[0, 10].compute_and_differentiate(tdb_jd)
    sun_pos_m = sun_pos * 1.0e3
    sun_vel_m_s = sun_vel * _KM_PER_DAY_TO_M_PER_S

    # --- Moon (barycentric) ---
    emb_pos, emb_vel = kernel[0, 3].compute_and_differentiate(tdb_jd)
    moon_from_emb_pos, moon_from_emb_vel = kernel[3, 301].compute_and_differentiate(tdb_jd)
    moon_pos_m = (emb_pos + moon_from_emb_pos) * 1.0e3
    moon_vel_m_s = (emb_vel + moon_from_emb_vel) * _KM_PER_DAY_TO_M_PER_S

    # --- Geocentric Sun and Moon ---
    sun_geo_pos = sun_pos_m - earth_pos_m
    sun_geo_vel = sun_vel_m_s - earth_vel_m_s
    moon_geo_pos = moon_pos_m - earth_pos_m
    moon_geo_vel = moon_vel_m_s - earth_vel_m_s

    # --- Fill Fortran arrays ---
    # f2py exposes Fortran arrays in column-major order matching the Fortran shape.
    # ext_earth(3,3,MAX_EPH_STEPS) -> Python shape (3,3,150)
    for t in range(n_steps):
        eph.ext_earth[:, 0, t] = earth_pos_m[:, t]
        eph.ext_earth[:, 1, t] = earth_vel_m_s[:, t]
        eph.ext_earth[:, 2, t] = earth_accel_m_s2[:, t]

        eph.ext_sun[:, 0, t] = sun_geo_pos[:, t]
        eph.ext_sun[:, 1, t] = sun_geo_vel[:, t]

        eph.ext_xmoon[:, 0, t] = moon_geo_pos[:, t]
        eph.ext_xmoon[:, 1, t] = moon_geo_vel[:, t]

        eph.ext_sunb[:, 0, t] = sun_pos_m[:, t]
        eph.ext_sunb[:, 1, t] = sun_vel_m_s[:, t]

        eph.ext_moonb[:, 0, t] = moon_pos_m[:, t]
        eph.ext_moonb[:, 1, t] = moon_vel_m_s[:, t]

    # --- Planets (7 planets: Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune) ---
    for pidx, spk_targ in enumerate(_PLANET_SPK_TARGETS):
        p_pos, p_vel = kernel[0, spk_targ].compute_and_differentiate(tdb_jd)
        p_pos_m = p_pos * 1.0e3
        p_vel_m_s = p_vel * _KM_PER_DAY_TO_M_PER_S
        for t in range(n_steps):
            eph.ext_splanet[:, 0, pidx, t] = p_pos_m[:, t]
            eph.ext_splanet[:, 1, pidx, t] = p_vel_m_s[:, t]
            eph.ext_gplanet[:, 0, pidx, t] = p_pos_m[:, t] - earth_pos_m[:, t]
            eph.ext_gplanet[:, 1, pidx, t] = p_vel_m_s[:, t] - earth_vel_m_s[:, t]

    # Reset step counter for Fortran PEP to read from index 1
    eph.eph_step_idx = 1
