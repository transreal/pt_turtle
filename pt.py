"""
pt.py -- Penrose rhomb tilings (pentagrid dual and exact substitution) and
turtle walks, mirroring the conventions of CellularAutomata.wl / TurtleTiling.wl.

Turtle state: (tile c, edge e in 0..3).  Edges of a tile are listed in
counter-clockwise vertex order; edge e joins vertex e and vertex (e+1)%4.
  S : cross edge e into the neighbour c', new edge = opposite of the shared
      edge on the c' side  (opp(e) = (e+2)%4)
  R : e -> (e+3)%4   (clockwise, CCW winding)
  L : e -> (e+1)%4
"""
import numpy as np
from collections import defaultdict

PHI = (1 + 5 ** 0.5) / 2
EVEC = np.array([[np.cos(2 * np.pi * i / 5), np.sin(2 * np.pi * i / 5)] for i in range(5)])


# ----------------------------------------------------------------------------
# Pentagrid (de Bruijn) generator -- same conventions as GeneratePenroseRhombs
# ----------------------------------------------------------------------------
def pentagrid_tiles(xmin, ymin, width, height, offsets):
    """Return dict with 'keys' (T,4,5) int, 'verts' (T,4,2) float, 'fat' (T,) bool,
    'label' (T,7) int.  Vertex order is CCW, identical to computeVertices[]."""
    off = np.asarray(offsets, float)
    xmax, ymax = xmin + width, ymin + height
    corners = np.array([[xmin, ymin], [xmax, ymin], [xmin, ymax], [xmax, ymax]], float)
    proj = corners @ EVEC.T - off  # (4,5)
    lo = np.floor(proj.min(0)).astype(int) - 1
    hi = np.ceil(proj.max(0)).astype(int) + 1
    tiles = {}
    for i in range(4):
        for j in range(i + 1, 5):
            ms = np.arange(lo[i], hi[i] + 1)
            ns = np.arange(lo[j], hi[j] + 1)
            M, N = np.meshgrid(ms, ns, indexing="ij")
            M = M.ravel(); N = N.ravel()
            # solve  x.e_i = m+g_i ,  x.e_j = n+g_j
            A = np.array([EVEC[i], EVEC[j]])
            Ainv = np.linalg.inv(A)
            rhs = np.stack([M + off[i], N + off[j]], axis=0)  # (2,K)
            P = (Ainv @ rhs).T  # (K,2)
            inside = (P[:, 0] >= xmin) & (P[:, 0] <= xmax) & (P[:, 1] >= ymin) & (P[:, 1] <= ymax)
            P = P[inside]; M = M[inside]; N = N[inside]
            K = np.ceil(P @ EVEC.T - off).astype(int)  # (K,5)
            K[:, i] = M; K[:, j] = N
            r, s = i, j
            if s - r in (3, 4):
                r, s = j, i
            for kk in K:
                label = (*kk, r, s)
                if label not in tiles:
                    tiles[label] = None
    labels = np.array(list(tiles.keys()), dtype=int)  # (T,7)
    T = len(labels)
    K = labels[:, :5]
    r = labels[:, 5]; s = labels[:, 6]
    er = np.eye(5, dtype=int)[r]; es = np.eye(5, dtype=int)[s]
    keys = np.stack([K, K + er, K + er + es, K + es], axis=1)  # (T,4,5)
    rot = ((K.sum(1) + 1) % 5) == 2
    keys[rot] = np.roll(keys[rot], -2, axis=1)  # start from v2: {v2,v3,v0,v1}
    verts = keys @ EVEC
    fat = ((s - r) % 5) == 1
    return {"keys": keys, "verts": verts, "fat": fat, "label": labels}


# ----------------------------------------------------------------------------
# Navigation tables from vertex keys
# ----------------------------------------------------------------------------
def build_nav(keys):
    """keys: (T,4,5) int (any consistent integer vertex labels).
    Returns edgeNbr (T,4) int32 (-1 = boundary), edgeEntry (T,4) int32, winding check."""
    T = keys.shape[0]
    edges = defaultdict(list)
    kt = [tuple(map(tuple, keys[t])) for t in range(T)]
    for t in range(T):
        vk = kt[t]
        for e in range(4):
            a, b = vk[e], vk[(e + 1) % 4]
            key = (a, b) if a <= b else (b, a)
            edges[key].append((t, e))
    edgeNbr = -np.ones((T, 4), dtype=np.int32)
    edgeEntry = -np.ones((T, 4), dtype=np.int32)
    bad = 0
    for key, lst in edges.items():
        if len(lst) == 2:
            (t1, e1), (t2, e2) = lst
            edgeNbr[t1, e1] = t2; edgeEntry[t1, e1] = e2
            edgeNbr[t2, e2] = t1; edgeEntry[t2, e2] = e1
        elif len(lst) > 2:
            bad += 1
    if bad:
        raise RuntimeError(f"{bad} edges shared by >2 tiles")
    return edgeNbr, edgeEntry


def windings(verts):
    v = verts
    x = v[:, :, 0]; y = v[:, :, 1]
    xn = np.roll(x, -1, axis=1); yn = np.roll(y, -1, axis=1)
    return np.sign((x * yn - xn * y).sum(1))


def parse_prog(prog):
    return np.array([{"S": 0, "R": 1, "L": 2}[ch] for ch in prog], dtype=np.int32)


# ----------------------------------------------------------------------------
# Vectorised turtle scan
# ----------------------------------------------------------------------------
def turtle_scan(edgeNbr, edgeEntry, prog, start_cells, max_cycles=100000, edges=(0, 1, 2, 3)):
    """Run the cyclic program from every (cell, edge) start.
    Returns dict: period (-1 if not closed), flag (0 closed, 1 escaped, 2 timeout),
    orbit_id (min state index c*4+e over cycle boundaries; equal for starts on the same orbit)."""
    cmds = parse_prog(prog)
    sc = np.repeat(np.asarray(start_cells, np.int64), len(edges))
    se = np.tile(np.asarray(edges, np.int64), len(start_cells))
    n = sc.size
    c = sc.copy(); e = se.copy()
    period = -np.ones(n, np.int64)
    flag = np.full(n, 2, np.int8)
    orbit = sc * 4 + se
    active = np.ones(n, bool)
    eN = edgeNbr.reshape(-1); eE = edgeEntry.reshape(-1)
    idx_all = np.arange(n)
    for cyc in range(1, max_cycles + 1):
        ia = idx_all[active]
        if ia.size == 0:
            break
        ca = c[ia]; ea = e[ia]
        alive = np.ones(ia.size, bool)
        for cmd in cmds:
            if cmd == 0:
                idx = ca * 4 + ea
                nb = eN[idx]
                esc = (nb < 0) & alive
                if esc.any():
                    flag[ia[esc]] = 1
                    alive &= ~esc
                ea = np.where(alive, (eE[idx] + 2) & 3, ea)
                ca = np.where(alive, nb, ca)
            elif cmd == 1:
                ea = (ea + 3) & 3
            else:
                ea = (ea + 1) & 3
        c[ia] = ca; e[ia] = ea
        st = ca * 4 + ea
        orbit[ia] = np.minimum(orbit[ia], np.where(alive, st, orbit[ia]))
        closed = alive & (ca == sc[ia]) & (ea == se[ia])
        if closed.any():
            period[ia[closed]] = cyc
            flag[ia[closed]] = 0
        active[ia[~alive | closed]] = False
    return {"start_cell": sc, "start_edge": se, "period": period, "flag": flag, "orbit_id": orbit}


def turtle_run(edgeNbr, edgeEntry, prog, cell, edge, max_cycles=100000):
    """Scalar run returning (period, flag, visited_cells list in order of first visit, states)."""
    cmds = parse_prog(prog)
    c, e = int(cell), int(edge)
    seen = {c: 0}
    order = [c]
    states = []
    for cyc in range(1, max_cycles + 1):
        for cmd in cmds:
            if cmd == 0:
                nb = edgeNbr[c, e]
                if nb < 0:
                    return -1, 1, order, states
                e = (edgeEntry[c, e] + 2) & 3
                c = int(nb)
                if c not in seen:
                    seen[c] = len(order); order.append(c)
            elif cmd == 1:
                e = (e + 3) & 3
            else:
                e = (e + 1) & 3
        states.append((c, e))
        if c == cell and e == edge:
            return cyc, 0, order, states
    return -1, 2, order, states


# ----------------------------------------------------------------------------
# Geometric D10 canonical form of a set of tiles (centroid based, as turtledeep.wl)
# ----------------------------------------------------------------------------
def d10_canon(cents, cells, ndig=2):
    """isometry-invariant signature: sorted distances (rounded) from the centroid of the swept
    tiles' centroids (same signature as the CUDA survey)."""
    cs = cents[cells]
    r = np.sqrt(((cs - cs.mean(0)) ** 2).sum(1))
    return np.round(r, ndig).tobytes() if False else tuple(np.sort(np.round(r * 10 ** ndig).astype(np.int64)).tolist())


def classify_orbits(edgeNbr, edgeEntry, cents, prog, scan, max_cycles=100000):
    """From a turtle_scan result, list distinct closed orbits and D10 classes.
    Returns list of dict(period, ncells, count, rep_cell, rep_edge, radius, center)."""
    ok = scan["flag"] == 0
    oid = scan["orbit_id"][ok]
    uniq, first, counts = np.unique(oid, return_index=True, return_counts=True)
    sc = scan["start_cell"][ok][first]; se = scan["start_edge"][ok][first]
    per = scan["period"][ok][first]
    classes = {}
    orbits = []
    for u, c0, e0, p, cnt in zip(uniq, sc, se, per, counts):
        pp, fl, order, _ = turtle_run(edgeNbr, edgeEntry, prog, c0, e0, max_cycles)
        assert fl == 0 and pp == p
        cells = np.array(order)
        ctr = cents[cells].mean(0)
        rad = np.sqrt(((cents[cells] - ctr) ** 2).sum(1)).max()
        key = (int(p), len(cells), d10_canon(cents, cells))
        orbits.append({"orbit_id": int(u), "period": int(p), "ncells": len(cells), "nstarts": int(cnt),
                       "rep": (int(c0), int(e0)), "radius": float(rad), "center": ctr})
        cl = classes.setdefault(key, {"period": int(p), "ncells": len(cells), "norbits": 0, "nstarts": 0,
                                       "rep": (int(c0), int(e0)), "radius": float(rad), "cells": cells})
        cl["norbits"] += 1; cl["nstarts"] += int(cnt)
        cl["radius"] = max(cl["radius"], float(rad))
    cls = sorted(classes.values(), key=lambda d: (d["period"], d["ncells"]))
    return orbits, cls


def centroids(verts):
    return verts.mean(1)


def depth_from_bbox(cents):
    x, y = cents[:, 0], cents[:, 1]
    return np.minimum.reduce([x - x.min(), x.max() - x, y - y.min(), y.max() - y])
