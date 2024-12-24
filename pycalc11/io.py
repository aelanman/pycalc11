import numpy as np
from astropy import units as u
from astropy.table import QTable
from astropy.time import Time


__all__ = ["parse_calc", "parse_im"]


def key_value(line):
    """Get key parts and the value from a line.

    Interpret the value as a number if possible, including attaching
    a unit if given by the last part of the key.
    """
    k, v = line.strip().split(":")
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
    if k_list[-1][0] == "(" and k_list[-1][-1] == ")":
        unit = k_list[-1][1:-1].lower()
        if unit.startswith("sec"):
            unit = "s"
        if "TIME" in k:  # unit "(mjd)"
            v = Time(v, format=unit)
        else:
            v = u.Quantity(v, unit)

        k_list = k_list[:-1]

    return k_list, v


def add_to_table(t, field, item, v, col_shape, **meta):
    # If we already created the column in the table, just set
    # the appropriate item, otherwise create the column, making
    # sure that for the first column, we use the correct length.
    if isinstance(v, str) and len(v) < 8:
        v = f"{v:<8s}"
    if field in t.colnames:
        t[field][item] = v
    else:
        # First access to this column; get right shape and
        # set possible metadata.
        t[field] = np.broadcast_to(v, col_shape, subok=True)
        if meta:
            t[field].info.meta = meta.copy()


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
    # Initialize output.
    meta = {}
    # Many items such as TELESCOPE come in sequences; these
    # get gathered together in small tables.
    kind = None  # kind of items (e.g., TELESCOPE).
    t = None  # table being created.
    n_kind = None  # number of elements to expect.
    with open(name) as im:
        for line in im.readlines():
            # Separate out keyword and its value.
            k_list, v = key_value(line)

            if kind:
                # We're in the process of making a table.
                if k_list[0] == kind:
                    # Got a match.  Add item to table.
                    assert t is not None
                    # The second part is always the sequence number (e.g., 0).
                    # While the rest describes what item it is (e.g., TIME).
                    item = int(k_list[1])
                    field = " ".join(k_list[2:]).lower()
                    add_to_table(t, field, item, v, col_shape=(n_kind,))
                    continue

                # No match, so signal that we are no longer matching,
                # and pass through to do something else with this line.
                kind = None

            # If the key starts with "NUM", a list will follow
            # (e.g., NUM TELESCOPES).  Items themselves do not have
            # the trailing "S", logically, but not all lists have
            # plurals (e.g., NUM SPACECRAFT).  This is obviously fragile.
            if k_list[0] == "NUM":
                kind = k_list[1]
                if kind.endswith("S"):
                    kind = kind[:-1]
                t = meta[kind.lower()] = QTable()
                n_kind = v
                continue

            # And if we are not making a list, just create an entry.
            meta[" ".join(k_list).lower()] = v

    return meta


def parse_im(name):
    """Represent the contents of a .im file as a dictionary.

    The dictionary will be keyed following entries in the file, except that
    they are in lower case, and that any unit is stripped from the key
    and added to the value.  Furthermore, sequences are converted to small
    tables. For instance::

        NUM TELESCOPES:     2
        TELESCOPE 0 NAME:   CHIME
        TELESCOPE 1 NAME:   ARO10m

    Will be converted to a table with column 'name', and the whole will be
    stored under entry 'telescope'.

    Scans are similarly stored as tables under entry 'scan', but are
    diffferent in that each scan has its own meta-data.  For instance::

        NUM SCANS:          1
        SCAN 0 POINTING SRC:B0531+21
        SCAN 0 NUM PHS CTRS:1
        SCAN 0 PHS CTR 0 SRC:B0531+21
        SCAN 0 NUM POLY:    3
        SCAN 0 POLY 0 MJD:  59153
        SCAN 0 POLY 0 SEC:  39000
        SRC 0 ANT 0 DELAY (us): 1.890259772582982e+04    1.428524766492706e-02  -3.424386838841481e-05  -1.265439323865132e-11   1.509492047497244e-14   2.741316895347732e-19
        SRC 0 ANT 0 DRY (us): 8.118429973475821e-03     -6.103368456164730e-09   1.463522056046739e-11  -1.655725632336728e-17   1.985277088714070e-20  -3.632645144125028e-27

    Will have a table as entry ['scans'][0], with meta-data 'pointing src' and
    'phs ctr' (as well as, for convenience, an 'main' entry that links back to
    the main metadata), and columns 'mjd', 'sec', 'delay', 'dry', etc.

    One gets a particular polynomial representation by indexing the result
    columns by a tuple of three integers representing poly, source, and
    antenna.  For the polynomial result columns, the unit is stored in
    ``col.info.meta['unit']``.

    """  # noqa
    # Initialize outputs
    main = meta = {}
    scans = []
    # Many items such as TELESCOPE come in sequences; these
    # get gathered together in small tables.
    kind = None  # kind of items (e.g., TELESCOPE).
    t = None  # table being created.
    n_kind = n_src = n_ant = order = poly = None
    scan_poly = False
    with open(name) as im:
        for line in im.readlines():
            # Separate out keyword and its value.
            k_list, v = key_value(line)

            if kind:
                # We're in the process of making a table.
                if k_list[: len(kind)] == kind:
                    # Got a match.  Add item to table.
                    item = int(k_list[len(kind)])
                    field = k_list[len(kind) + 1].lower()
                    add_to_table(t, field, item, v, col_shape=(n_kind,))
                    if scan_poly:
                        # For use below, remember polynomial we are doing.
                        poly = item
                    continue

                if scan_poly and k_list[0] == "SRC" and k_list[2] == "ANT":
                    # For a SCAN POLY, SRC and ANT items belong to the table.
                    field = k_list[4].lower()
                    src, ant = int(k_list[1]), int(k_list[3])
                    v = np.array([float(v_) for v_ in v.split()])
                    unit = u.Unit(
                        "deg" if field == "az" or field == "el" else k_list[5][1:-1]
                    )
                    add_to_table(
                        t,
                        field,
                        (poly, src, ant),
                        v,
                        col_shape=(n_kind, n_src, n_ant, order + 1),
                        unit=unit,
                    )
                    continue

                # Really not a match.  Signal that we are no longer matching,
                # and pass through to do something else with this line.
                kind = None
                assert not (
                    scan_poly and meta
                ), "Meta data found after start of SCAN POLY."

            if "NUM" in k_list:
                if k_list == ["NUM", "SCANS"]:
                    # For scans, we start a list instead of a table, since we
                    # will want to store the individual scans in tables.
                    # We special-case this, so do not set "kind".
                    scans = [None] * v
                    # We start a new meta for the first table.
                    meta = {}
                    continue

                # Line can be like "NUM TELESCOPES" or "SCAN 0 NUM POLY".
                i_num = k_list.index("NUM")
                # Key used by subsequent lines to indicate items in the table.
                # E.g., TELESCOPE, SCAN 0 PHS CTR, or SCAN 0 POLY.
                kind = k_list[:i_num] + k_list[i_num + 1 :]
                if kind[-1].endswith("S"):
                    kind[-1] = kind[-1][:-1]
                # Keep number for when creating the first column later.
                n_kind = v

                scan_poly = kind[0] == "SCAN" and kind[2] == "POLY"
                if scan_poly:
                    # Start a SCAN POLY table, with table-specific meta-data
                    # gathered so far.
                    t = scans[int(kind[1])] = QTable(meta=meta, copy=False)
                    # Reset meta for next scan poly table.
                    meta = {}
                else:
                    # Start a regular table in meta.
                    key = " ".join(kind[i_num:]).lower()
                    t = meta[key] = QTable()
                    # Store some relevant numbers for the scan poly tables.
                    if key == "telescope":
                        n_ant = n_kind
                    elif key == "phs ctr":
                        n_src = n_kind + 1

                continue

            # Just a regular meta item. Store it.
            key = " ".join(k_list[(0 if k_list[0] != "SCAN" else 2) :]).lower()
            meta[key] = v
            if key == "polynomial order":
                order = v

    # Keep a link to the main meta in all the scans.
    for scan in scans:
        scan.meta["main"] = main
    result = main.copy()
    result["scan"] = scans
    return result
