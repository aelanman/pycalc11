
from . import interface as ci



def get_delay(telescope_positions=None, telescope_names=None, source_coords=None,
             source_names=None, time=None, duration_min=None, calc_file_name=None):
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

    return ci.delay
