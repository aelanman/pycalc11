import numpy as np
from astropy import units as u
from astropy.table import QTable, vstack
from astropy.time import Time


__all__ = ['parse_calc', 'parse_im']


def key_value(line):
    """Get key parts and the value from a line.

    Interpret the value as a number if possible, including attaching
    a unit if given by the last part of the key.
    """
    k, v = line.strip().split(':')
    k_list = k.split()
    # Interpret value as a number, Quantity, or Time if possible.
    try:
        v = int(v)
    except Exception:
        try:
            v = float(v)
        except Exception:
            v = v.strip()
            # Not simply numeric, so don't try units.
            return k_list, v

    # Check whether a unit is given -- as the last item of the key.
    if k_list[-1][0] == '(' and k_list[-1][-1] == ')':
        unit = k_list[-1][1:-1].lower()
        if unit.startswith('sec'):
            unit = 's'
        if 'TIME' in k:  # unit "(mjd)"
            v = Time(v, format=unit)
        else:
            v = u.Quantity(v, unit)

        k_list = k_list[:-1]

    return k_list, v


def add_to_table(t, field, item, v, col_shape):
    # If we already created the column in the table, just set
    # the appropriate item, otherwise create the column, making
    # sure that for the first column, we use the correct length.
    print(field, item, v)
    if isinstance(v, str) and len(v) < 8:
        v = f"{v:<8s}"
    if field in t.colnames:
        t[field][item] = v
    else:
        # First access to this column; get right shape
        t[field] = np.broadcast_to(v, col_shape, subok=True)


def parse_calc(name):
    """Represent the content of a .calc file as a dictionary.

    The dictionary will be keyed following entries in the file, except that
    they are in lower case, and that any unit is stripped from the key
    and added to the value.  Furthermore, sequences are converted to small
    tables. For instance::

        NUM TELESCOPES:     2
        TELESCOPE 0 NAME:   CHIME
        TELESCOPE 0 MOUNT:  AZEL
        TELESCOPE 0 OFFSET (m): 0.0000
        TELESCOPE 0 X (m): -2059159.52756903
        TELESCOPE 0 Y (m): -3621260.06650439
        TELESCOPE 0 Z (m): 4814325.37572648
        TELESCOPE 0 SHELF:  None
        TELESCOPE 1 NAME:   ARO10m
        TELESCOPE 1 MOUNT:  AZEL
        TELESCOPE 1 OFFSET (m): 0.0000
        TELESCOPE 1 X (m): 918239.85303214
        TELESCOPE 1 Y (m): -4346109.57976893
        TELESCOPE 1 Z (m): 4562002.27426509
        TELESCOPE 1 SHELF:  None

    Will be converted to a table with columns 'name', 'mount', etc.,
    and the whole will be stored under entry 'telescope'.
    """
    # Initialize output directory
    meta = {}
    # Many items such as TELESCOPE come in sequences; these
    # get gathered together in small tables.
    kind = None  # kind of items (e.g., TELESCOPE).
    n_kind = None  # number of elements to expect.
    t = None  # table being created.
    with open(name) as im:
        while line := im.readline():
            # Separate out keyword and its value.
            k_list, v = key_value(line)

            # If this is one of a list of items we are gathering, add it.
            # E.g., "EOP 0 TIME (mjd)"
            if k_list[0] == kind:
                assert t is not None
                # The second part is always the sequence number (e.g., 0).
                # While the rest describes what item it is (e.g., TIME).
                add_to_table(t, ' '.join(k_list[2:]).lower(), int(k_list[1]),
                             v, col_shape=(n_kind,))
                continue

            # If we were creating a table but got here, the table is
            # complete, so just store it.
            if kind:
                meta[kind.lower()] = t
                kind = None

            # If the key starts with "NUM", a list will follow
            # (e.g., NUM TELESCOPES).  Items themselves do not have
            # the trailing "S", logically, but not all lists have
            # plurals (e.g., NUM SPACECRAFT).  This is obviously fragile.
            if k_list[0] == 'NUM':
                kind = k_list[1]
                if kind.endswith('S'):
                    kind = kind[:-1]
                t = QTable()
                n_kind = v
                continue

            # And if we are not making a list, just create an entry.
            meta[' '.join(k_list).lower()] = v

    return meta


def parse_im(name):
    """Convert the content of a .im file to a dictionary."""
    # Initialize outputs
    scans = {}
    main = meta = {}
    # Many items such as TELESCOPE come in sequences; these
    # get gathered together in small tables.
    kind = None  # kind of items (e.g., TELESCOPE).
    poly = None
    n_kind = None  # number of elements to expect.
    t = None  # table being created.
    n_src = n_ant = None
    with open(name) as im:
        while line := im.readline():
            k_list, v = key_value(line)

            if k_list[0] == 'SRC' and k_list[2] == 'ANT':
                assert kind[0] == 'SCAN' and kind[2] == 'POLY'
                field = k_list[4].lower()
                src, ant = int(k_list[1]), int(k_list[3])
                v = np.array([float(v_) for v_ in v.split()])
                add_to_table(t, field, (poly, src, ant), v,
                             col_shape=(n_kind, n_src, n_ant, 6))
                t[field].info.meta.setdefault(
                    'unit', ('deg' if field == 'az' or field == 'el'
                             else k_list[5][1:-1]))
                continue

            if kind and k_list[:len(kind)] == kind:
                item = int(k_list[len(kind)])
                add_to_table(t, k_list[len(kind)+1].lower(), item, v,
                             col_shape=(n_kind,))
                if kind[0] == 'SCAN' and kind[2] == 'POLY':
                    poly = item
                continue

            if kind:
                if kind[0] == 'SCAN':
                    if kind[2] == 'POLY':
                        t.meta.update(meta)
                        scans[int(kind[1])] = t
                    else:
                        meta[' '.join(kind[2:]).lower()] = t
                else:
                    meta[' '.join(kind).lower()] = t
                kind = poly = None

            if 'NUM' in k_list:
                if k_list == ['NUM', 'SCANS']:
                    # For scans, we start a new meta, which we will store
                    # in the SCAN POLY table once that comes along.
                    meta = {'main': main}
                    continue

                i_num = k_list.index('NUM')
                kind = k_list[:i_num] + k_list[i_num+1:]
                if kind[-1].endswith('S'):
                    kind[-1] = kind[-1][:-1]
                t = QTable()
                n_kind = v
                if kind == ['TELESCOPE']:
                    n_ant = v
                elif kind[0] == 'SCAN' and kind[2:4] == ['PHS', 'CTR']:
                    n_src = v+1
                continue

            if k_list[0] == 'SCAN':
                k_list = k_list[2:]
            meta[' '.join(k_list).lower()] = v

    if kind:
        if kind[0] == 'SCAN' and kind[-1] == 'POLY':
            t.meta.update(meta)
            scans[int(kind[1])] = t
        else:
            meta[' '.join(kind).lower()] = t

    result = main.copy()
    result['scan'] = scans
    return result
