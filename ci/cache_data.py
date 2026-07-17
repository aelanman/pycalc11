"""
Cache ephemeris and IERS data for CI.
"""

from astropy.utils.data import (
    download_file,
    export_download_cache,
    import_download_cache,
    get_cached_urls,
)
from astropy.time import Time
import os
import sys

home = os.path.expanduser("~")
cache_file = os.path.join(home, "astropy_cache.zip")

# JPL SPK kernels read in Python via jplephem. 'de440s' is the default
# ephemeris; 'de421' is used by the ephemeris comparison tests.
_SPK_URLS = [
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp",
]

# Legacy Fortran DE421 binary reader (ephemeris="legacy").
_DE421_URL = (
    "https://github.com/difx/difx/raw/refs/heads/main/"
    f"applications/difxcalc11/data/DE421_{sys.byteorder}_Endian"
)

# Download iers data
Time.now().ut1

if sys.argv[1] == "save":
    if not os.path.exists(cache_file):
        for url in _SPK_URLS:
            download_file(url, cache=True)
        download_file(_DE421_URL, cache=True)
        urls = get_cached_urls()
        export_download_cache(cache_file, urls=urls, overwrite=True)

if sys.argv[1] == "load":
    import_download_cache(cache_file)
