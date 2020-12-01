import numpy as np
from astropy import units as u
from astropy.table import QTable, vstack
from astropy.time import Time


def calc2dict(name):
    meta = {}
    kind = None
    n_kind = None
    t = None
    with open(name) as im:
        while line := im.readline():
            k, v = line.strip().split(':')
            k_list = k.split()
            v_list = v.split()
            if len(v_list) == 1:
                try:
                    v = int(v)
                except Exception:
                    try:
                        v = float(v)
                    except Exception:
                        v = v.strip()
            if k_list[-1][0] == '(' and k_list[-1][-1] == ')':
                unit = k_list[-1][1:-1].lower()
                if unit.startswith('sec'):
                    unit = 's'
                if 'TIME' in k:
                    v = Time(v, format=unit)
                else:
                    v = u.Quantity(v, unit)

                k_list = k_list[:-1]

            if k_list[0] == kind:
                j = int(k_list[1])
                field = ' '.join(k_list[2:]).lower()
                if field in t.colnames:
                    t[field][j] = v
                elif len(t) == n_kind:
                    t[field] = v
                else:
                    if isinstance(v, str) and len(v) < 8:
                        v = f"{v:<8s}"
                    t[field] = np.full((n_kind,), v)
                continue

            if kind:
                meta[kind.lower()] = t
                kind = None

            if k_list[0] == 'NUM':
                kind = k_list[1]
                if kind.endswith('S'):
                    kind = kind[:-1]
                n_kind = v
                t = QTable()
                continue

            else:
                meta[' '.join(k_list).lower()] = v

    return meta


def im2dict(name):
    meta = {}
    results = {}
    scan_id = None
    poly = None
    with open(name) as im:
        while line := im.readline():
            k, v = line.strip().split(':')
            k_list = k.split()
            v_list = v.split()
            if len(v_list) == 1:
                try:
                    v = int(v)
                except Exception:
                    v = v.strip()
            if k_list[0] == 'SRC' and k_list[2] == 'ANT':
                src, ant = int(k_list[1]), int(k_list[3])
                item = k_list[4].lower()
                if item == 'az' or item == 'el':
                    unit = 'deg'
                else:
                    unit = k_list[5][1:-1]
                assert scan_id is not None and poly is not None
                results[scan_id][poly].setdefault(item, {})[src, ant] = (
                    v_list, unit)
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
                if k_list[-1][0] == '(' and k_list[-1][-1] == ')':
                    unit = k_list[-1][1:-1].lower()
                    if unit.startswith('sec'):
                        unit = 's'
                    v = u.Quantity(v, unit)

                    k_list = k_list[:-1]
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
