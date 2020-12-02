*************************
Python bindings to CALC11
*************************

To install, you also need the ``difxcalc11`` source code, and have ``numpy``
installed (as well as ``pytest`` for running tests).

Assuming you are at the top of the directory of ``pycalc11``, you can do the
following::

  git svn clone https://svn.atnf.csiro.au/difx/applications/difxcalc11/trunk difxcalc11
  pip install -e .

If you already have the ``difxcalc11`` source and data directories somewhere,
either create a symbolic link to its location, or set the environment variable
``DIFXCALC11`` with the path to it.
