# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import sys

import setuptools
from setuptools.config import read_configuration


#INSTALL_HELP = """
#To install pycalc11, you should either have difxcalc11 installed in the
#calc11 directory or have defined an environment variable DIFXCALC11
#that contains the path.
#"""

#BASE = os.getenv('DIFXCALC11', default=os.path.join(
#    os.path.abspath(os.path.dirname(__file__)), 'difxcalc11'))
BASE = os.path.join(os.path.abspath(os.path.dirname(__file__)), "calc11")
SRCDIR = os.path.join(BASE, 'src')
print(BASE)
DATADIR = os.path.join(BASE, 'data')
if not all(os.path.isdir(p) for p in (SRCDIR, DATADIR)):
    print(f"Failed to find source or data directories {SRCDIR} and {DATADIR}")
    sys.exit(1)

MODNAME = 'calc11'
F90_COMBINED = os.path.join(SRCDIR, MODNAME+'.f90')

SRC_FILES = os.listdir(SRCDIR)
# Move the file data_module.f to the start of SRC_FILES
fn = "data_module.f"
SRC_FILES.remove(fn)
SRC_FILES.insert(0, fn)

F_SOURCES = [os.path.join(SRCDIR, f) for f in SRC_FILES if f.endswith('.f')
             and not f.startswith('dmain')]
C_SOURCES = [os.path.join(SRCDIR, f) for f in SRC_FILES if f.endswith('.c')]
INCLUDES = [os.path.join(SRCDIR, f) for f in SRC_FILES if f.endswith('.i')
            and not f.startswith('param11')]
PARAM11_IN = os.path.join(SRCDIR, 'param11.i.in')
PARAM11 = os.path.join(SRCDIR, 'param11.i')

DATA_FILES = [os.path.join(BASE, 'data', f) for f in os.listdir(DATADIR)
              if f.endswith(('.dat', '.coef')) or sys.byteorder in f]
DE421_FILE_NAME = [f for f in os.listdir(DATADIR) if sys.byteorder in f][0]

def get_extensions():
    from distutils.dep_util import newer
    from numpy.distutils.core import Extension

    # TODO: build inside a build directory rather than inside difxcalc11!!
#    if not (os.path.isfile(PARAM11) and newer(PARAM11, PARAM11_IN)):
    with open(PARAM11_IN, 'rt') as fr, open(PARAM11, 'wt') as fw:
        # TODO: point to installed data files, make link safe for windows.
        fw.writelines(
            [line.replace('@prefix@/share/difxcalc', str(DATADIR))
             .replace('@de421_file@', DE421_FILE_NAME)
             for line in fr.readlines()])

    if any(newer(fn, F90_COMBINED) for fn in F_SOURCES + INCLUDES):
        print(f'Creating combined f90 file {F90_COMBINED}')
        # Files do not get automatically recognized as f90, so be blunt and
        # get all of them.
        with open(F90_COMBINED, 'wt') as f90file:
            f90file.write('! -*- f90 -*-' + os.linesep)
            for f in F_SOURCES:
                with open(os.path.join(SRCDIR, f)) as fh:
                    f90file.writelines(fh.readlines())
    calc_ext = Extension(
        name="pycalc11.calc11",
        sources=[F90_COMBINED] + C_SOURCES,
        include_dirs=[SRCDIR],
        extra_f90_compile_args=['-Wargument-mismatch', '-fcheck=bounds', '-DF2PY_REPORT_ATEXIT'],
    )

    return [calc_ext]


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    # TODO: install data files as well and use them in library.
    #      data_files=[(os.path.join('difxcalc', 'data'), DATA_FILES)])
    config = Configuration('', parent_package, top_path,
                           ext_modules=get_extensions())
    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup_config = read_configuration('setup.cfg')
    metadata, options = setup_config['metadata'], setup_config['options']
    if 'python_requires' in options:
        options['python_requires'] = str(options['python_requires'])
    setup(configuration=configuration, **metadata, **options)
