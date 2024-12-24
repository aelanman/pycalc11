"""Runner for difxcalc, if installed."""

import shutil
import os
from subprocess import check_output

difxcalc = shutil.which("difxcalc")


def run_difxcalc(
    calcfile,
    dry_tropo=False,
    wet_tropo=False,
    force=False,
    verbose=True,
):
    if difxcalc is None:
        raise ValueError("difxcalc is not installed.")

    imfile = calcfile.replace(".calc", ".im")
    comm = [difxcalc]

    if not dry_tropo:
        comm.append("-dry")
    if not wet_tropo:
        comm.append("-wet")
    if force:
        comm.append("-f")
    if verbose:
        comm.append("-v")

    path, cfile = os.path.split(calcfile)
    cwd = os.getcwd()

    comm.append(cfile)

    try:
        os.chdir(path)
        out = check_output(comm).decode("utf-8")
    finally:
        os.chdir(cwd)

    if verbose:
        print(out)

    return imfile
