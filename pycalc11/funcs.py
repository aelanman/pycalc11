from .interface import CalcInterface, OceanFiles

ci = CalcInterface()

def _get_quantity(quantity,
    telescope_positions=None, telescope_names=None, source_coords=None,
    source_names=None, time=None, duration_min=None, calc_file_name=None
):
    """
    If calc_file_name is given, the other keywords are not required.
    """
    # TODO Option to only rerun if something changed
    ci.reset()
    if calc_file_name is not None:
        # Do setup with calcfile.
        ci.set_cl()
        ci.set_calcfile(calc_file_name)
        # TODO -- Parse calcfile for number of sources first.
        ci.run(alloc_out_arrays=(6, 1, 41, 301), alloc_source_arrays=(300,))    # NOTE -- Treating these like old fixed-size arrays for now
        ci.run(dget_input=(1,))
        ci.set_geocenter_telescope()
        ci.calc.dscan(1,1)
        ci.calc.dinitl(1)
        ci.run(adrivr=(1, 1))    # TODO After dscan, Intrvls2min is set. Loop over this to run adrivr for each 2 min chunk.
    else:
        ci.set_cl()
        ci.run(alloc_out_arrays=(6, 1, len(telescope_names)+1, len(source_names)+1))  # TODO Remove distinction between pointing and phase src
#        ci.run(alloc_source_arrays=(len(source_names)+1))
        ci.set_sources(source_coords, source_names)
        ci.set_eops(time)
        ci.set_scan(time, duration_min)
        ci.set_telescopes(telescope_positions, telescope_names)
        ci.run(dinitl=(1,), adrivr=(1, 1))

    return getattr(ci, quantity)()


def get_delay(**kwargs):
    # Returned axes: (time, ant1, ant2, src)
    # Same keyword args as _get_quantity
    delay = _get_quantity('delay', **kwargs)
    nsrcs = ci.nsrcs()
    nants = ci.nants()

    # Cut back on long axes to just those used.
    # 0th on source axis is pointing center
    # 1st and on are phase centers
    # TODO -- With allocatable arrays, is this still necessary?
    delay = delay[:, :, :nants, 1:nsrcs+1]
    return delay


def get_delay_rate(**kwargs):
    # Same keyword args as _get_quantity
    delay_rate = _get_quantity('delay_rate', **kwargs)
    nsrcs = ci.nsrcs()
    nants = ci.nants()
    # TODO -- With allocatable arrays, is this still necessary?
    delay_rate = delay_rate[:, :, :nants, 1:nsrcs+1]
    return delay_rate


def get_partials(**kwargs):
    # Returned axes: (ra/dec, d/dr, time, ant1, ant2, src))
    # Zeroth axis: (0) d/dra, (1) d/ddec
    # First axis: (0) delay, (1) delay rate
    partials = _get_quantity('partials', **kwargs)
    nsrcs = ci.nsrcs()
    nants = ci.nants()
    return partials[..., :nants, 1:nsrcs+1]
