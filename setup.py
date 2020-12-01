# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import setuptools
from distutils.dep_util import newer

from numpy.distutils.core import Extension, setup


BASE = os.path.abspath('difxcalc11')
SRCDIR = os.path.join(BASE, 'src')
DATADIR = os.path.join(BASE, 'data')
MODNAME = 'calc11'
F90_COMBINED = os.path.join(SRCDIR, MODNAME+'.f90')

SRC_FILES = os.listdir(SRCDIR)

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


# https://mail.python.org/pipermail/distutils-sig/2007-September/008253.html
class NumpyExtension(setuptools.Extension):
    """Extension type that adds the NumPy include directory to include_dirs."""

    @property
    def include_dirs(self):
        from numpy import get_include
        return self._include_dirs + [get_include()]

    @include_dirs.setter
    def include_dirs(self, include_dirs):
        self._include_dirs = include_dirs


def get_extensions():
    if not (os.path.isfile(PARAM11) and newer(PARAM11, PARAM11_IN)):
        with open(PARAM11_IN, 'rt') as fr, open(PARAM11, 'wt') as fw:
            # TODO: point to installed data files, make link safe for windows.
            fw.writelines(
                [line.replace('@prefix@/share/difxcalc', str(DATADIR))
                 .replace('@de421_file@', DE421_FILE_NAME)
                 for line in fr.readlines()])

    if not (os.path.isfile(F90_COMBINED)
            and all(newer(F90_COMBINED, fn) for fn in F_SOURCES + INCLUDES)):
        print(f'Creating combined f90 file {F90_COMBINED}')
        # Files do not get automatically recognized as f90, so be blunt and
        # get all of them.
        with open(F90_COMBINED, 'wt') as f90file:
            f90file.write('! -*- f90 -*-' + os.linesep)
            for f in F_SOURCES:
                with open(os.path.join(SRCDIR, f)) as fh:
                    f90file.writelines(fh.readlines())

    calc_ext = Extension(
        name="difxcalc.calc11",
        sources=[F90_COMBINED] + C_SOURCES,
        include_dirs=[SRCDIR],
        libraries=['gsl'],
        extra_f90_compile_args=['-fallow-argument-mismatch'],
    )

    return [calc_ext]


setup(name='calc11',
      zip_safe=False,
      version='0.1',
      ext_modules=get_extensions())
# TODO: install data files as well and use them in library.
#      data_files=[(os.path.join('difxcalc', 'data'), DATA_FILES)])
