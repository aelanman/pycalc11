import numpy as np
from astropy import units as u
from astropy.table import QTable, vstack
from astropy.time import Time


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
    if isinstance(v, str) and len(v) < 8:
        v = f"{v:<8s}"
    if not t.colnames:
        # First access to this table; get right length
        t[field] = np.broadcast_to(v, col_shape, subok=True)
    elif field in t.colnames:
        t[field][item] = v
    else:
        t[field] = v


def calc2dict(name):
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


def im2dict(name):
    """Convert the content of a .calc file to a dictionary."""
    meta = {}
    results = {}
    scan_id = None
    poly = None
    with open(name) as im:
        while line := im.readline():
            k_list, v = key_value(line)

            if k_list[0] == 'SRC' and k_list[2] == 'ANT':
                src, ant = int(k_list[1]), int(k_list[3])
                item = k_list[4].lower()
                if item == 'az' or item == 'el':
                    unit = 'deg'
                else:
                    unit = k_list[5][1:-1]
                assert scan_id is not None and poly is not None
                results[scan_id][poly].setdefault(item, {})[src, ant] = (
                    v.split(), unit)
            elif k_list[0] == 'SCAN':
                if k_list[2] == 'POLY':
                    assert int(k_list[1]) == scan_id
                    poly = int(k_list[3])
                    item = k_list[4].lower()
                    results[scan_id].setdefault(poly, {})[item] = v
                else:
                    scan_id = int(k_list[1])
                    item = ' '.join(k_list[2:]).lower()
                    results.setdefault(scan_id, {})[item] = v
            else:
                meta[' '.join(k_list).lower()] = v

        return results, meta


def read_im(name):
    results, meta = im2dict(name)
    tables = []
    n_tel = meta['num telescopes']
    order = meta['polynomial order']
    for scan in range(meta['num scans']):
        s = results[scan]
        n_poly = s['num poly']
        n_src = s['num phs ctrs'] + 1
        n = n_poly * n_src
        t = QTable()
        t['scan'] = np.full((n,), scan)
        t['poly'] = np.full((n,), 0)
        t['pointing'] = list(range(n_src)) * n_poly
        t['source'] = ([s['pointing src']]
                       + [s[f'phs ctr {i} src']
                          for i in range(s['num phs ctrs'])]) * n_poly
        # Set up actual result rows.
        for k, v in s[0].items():
            if not isinstance(v, dict):
                t[k] = np.full((n,), v)
            else:
                t[k] = np.zeros((n, n_tel, order+1))
                t[k].info.meta['unit'] = v[0, 0][1]
        for poly in range(n_poly):
            slc = slice(poly*n_src, (poly+1)*n_src)
            p = s[poly]
            t['poly'][slc] = poly
            for k, v in p.items():
                if not isinstance(v, dict):
                    t[k][slc] = v
                else:
                    t[k][slc] = np.stack([
                        np.stack([v[(src, ant)][0] for ant in range(n_tel)])
                        for src in range(n_src)])
        tables.append(t)

    result = vstack(tables)
    result.meta = meta
    return result
