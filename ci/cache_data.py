"""
Cache DE421 kernel data and
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
path = os.path.join(home, "ap_cache")
cache_file = os.path.join(path, "astropy_cache.zip")
os.makedirs(path, exist_ok=True)

# Download iers data
Time.now().ut1

if sys.argv[1] == "save":
    if not os.path.exists(cache_file):
        de421_url = f"https://svn.atnf.csiro.au/difx/applications/difxcalc11/trunk/data/DE421_{sys.byteorder}_Endian"
        de421_path = download_file(de421_url, cache=True)
        urls = get_cached_urls()
        export_download_cache(cache_file, urls=urls, overwrite=True)

if sys.argv[1] == "load":
    import_download_cache(cache_file)
