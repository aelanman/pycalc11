from skbuild import setup
import os
import requests
from tqdm.auto import tqdm
import sys

setup(
    packages=['pycalc11'],
    python_requires=">=3.7",
    include_package_data=True,
)
