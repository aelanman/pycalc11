import os
import sys
from astropy.utils.data import download_file

from .data import DATA_PATH
from . import calc11
from . import runner


# Default ephemeris used when none is specified. This is a JPL SPK kernel
# read in Python via jplephem (see get_spk / Calc._setup_ephemeris).
DEFAULT_EPHEMERIS = "de440s"

# Legacy DE421 binary URL for the Fortran ephemeris reader (ephemeris="legacy").
# Served from the difx GitHub mirror (the original ATNF SVN host is defunct).
# Override this value if the file location changes.
_DE421_URL = f"https://github.com/difx/difx/raw/refs/heads/main/applications/difxcalc11/data/DE421_{sys.byteorder}_Endian"
_de421_loaded = False


def _ensure_de421():
    """Download and set the legacy DE421 binary path if not already done.

    Only used by the legacy Fortran ephemeris reader (``ephemeris="legacy"``).
    Called lazily so that ``import pycalc11`` and the default JPL SPK path do
    not require this (potentially unavailable) download.
    """
    global _de421_loaded
    if _de421_loaded:
        return
    try:
        de421_path = download_file(_DE421_URL, cache=True)
    except Exception as err:
        raise RuntimeError(
            "Could not download the legacy DE421 binary from\n"
            f"    {_DE421_URL}\n"
            "Use the default JPL ephemeris instead (omit the `ephemeris` "
            "argument, or pass a name like 'de421'/'de440s'), or override "
            "pycalc11._DE421_URL to point at an available mirror."
        ) from err
    calc11.datafiles.jpl_de421 = de421_path.ljust(128)
    _de421_loaded = True


def _format(fname):
    pth = os.path.join(DATA_PATH, fname).ljust(128)
    return pth


# Set data file paths
calc11.datafiles.a_tilts = _format("tilt.dat")
calc11.datafiles.oc_file = _format("ocean_load.coef")
calc11.datafiles.optl_file = _format("ocean_pole_tide.coef")
calc11.datafiles.dfleap = _format("ut1ls.dat")


# Known SPK ephemeris URLs from NAIF/JPL
_SPK_URLS = {
    "de421": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp",
    "de430": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp",
    "de440": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
    "de440s": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
    "de441": "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp",
}


def get_spk(name="de440s"):
    """Download a JPL SPK ephemeris file.

    Parameters
    ----------
    name : str
        Ephemeris name
        Options: SPK_AVAIL
        The 's' variants (e.g., 'de440s') are smaller files covering a shorter
        time span and are recommended for typical use.

    Returns
    -------
    str
        Path to the cached SPK file.
    """
    name = name.lower()
    if name not in _SPK_URLS:
        raise ValueError(f"Unknown ephemeris '{name}'. Available: {list(_SPK_URLS.keys())}")
    return download_file(_SPK_URLS[name], cache=True)


get_spk.__doc__ = get_spk.__doc__.replace("SPK_AVAIL", "'" + "', '".join(list(_SPK_URLS)) + "'")

from .interface import Calc

__all__ = ["Calc", "DATA_PATH", "calc11", "runner", "get_spk", "DEFAULT_EPHEMERIS"]
