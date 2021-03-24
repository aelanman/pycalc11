from . import interface as ci


def _get_quantity(quantity,
    telescope_positions=None, telescope_names=None, source_coords=None,
    source_names=None, time=None, duration_min=None, calc_file_name=None
):
    """
    If calc_file_name is given, the other keywords are not required.
    """
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
    # Same keyword args as _get_quantity
    delay = _get_quantity('delay', **kwargs)
    nsrcs = ci.calc.calc_input.numphcntr
    nants = len([s for s in ci.calc.calc_input.sites if s != b''])

    # Cut back on long axes to just those used.
    delay = delay[:, :nants, :nants, :nsrcs+1]
    return delay


def get_delay_rate(**kwargs):
    # Same keyword args as _get_quantity
    delay_rate = _get_quantity('delay_rate', **kwargs)

    nsrcs = ci.calc.calc_input.numphcntr
    nants = len([s for s in ci.calc.calc_input.sites if s != b''])

    delay_rate = delay_rate[:, :nants, :nants, :nsrcs+1]
    return delay_rate
