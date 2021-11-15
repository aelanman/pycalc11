from . import interface as ci


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
        ci.run(dstart=(1, 1), dget_input=(1,), dinitl=(1,), dscan=(1, 1), ddrivr=(1, 1))
    else:
        ci.set_cl()
        ci.run(dstart=(1, 1))
        ci.set_sources(source_coords, source_names)
        ci.set_eops(time)
        ci.set_scan(time, duration_min)
        ci.set_telescopes(telescope_positions, telescope_names)
        ci.run(dinitl=(1,), ddrivr=(1, 1))

    return getattr(ci, quantity)


def get_delay(**kwargs):
    # Returned axes: (time, ant1, ant2, src)
    # Same keyword args as _get_quantity
    delay = _get_quantity('delay', **kwargs)
    nsrcs = ci.get_nsrcs()
    nants = ci.get_nants()

    # Cut back on long axes to just those used.
    # 0th on source axis is pointing center
    # 1st and on are phase centers
    delay = delay[:, :, :nants, 1:nsrcs+1]
    return delay


def get_delay_rate(**kwargs):
    # Same keyword args as _get_quantity
    delay_rate = _get_quantity('delay_rate', **kwargs)
    nsrcs = ci.get_nsrcs()
    nants = ci.get_nants()

    delay_rate = delay_rate[:, :, :nants, 1:nsrcs+1]
    return delay_rate


def get_partials(**kwargs):
    # Returned axes: (ra/dec, d/dr, time, ant1, ant2, src))
    # Zeroth axis: (0) d/dra, (1) d/ddec
    # First axis: (0) delay, (1) delay rate
    partials = _get_quantity('partials', **kwargs)
    nsrcs = ci.get_nsrcs()
    nants = ci.get_nants()
    return partials[..., :nants, 1:nsrcs+1]
