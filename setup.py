from skbuild import setup
import os
import requests
from tqdm.auto import tqdm
import sys

de421_url = f"https://svn.atnf.csiro.au/difx/applications/difxcalc11/trunk/data/DE421_{sys.byteorder}_Endian"
de421_file = os.path.join(os.path.dirname(__file__), "pycalc11", "data", os.path.basename(de421_url))

if not os.path.exists(de421_file):
    with requests.get(de421_url, stream=True) as r:
        r.raise_for_status()
        file_size = int(r.headers.get('Content-Length', 0))
        pbar = tqdm(total=file_size, unit='iB', unit_scale=True, desc=f"{os.path.basename(de421_file)}")
        with open(de421_file, 'wb') as ofile:
            for data in r.iter_content(4096):
                pbar.update(len(data))
                ofile.write(data)

setup(
    packages=['pycalc11'],
    python_requires=">=3.7",
    include_package_data=True,
)
