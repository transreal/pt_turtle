"""
subst.py -- exact Penrose rhomb substitution via Robinson triangles in Z[zeta]
(zeta = exp(2 pi i/5)).  Points are integer 5-vectors K with embedding
sum_j K_j zeta^j ; canonical representative has sum(K) in 0..4
(the all-ones vector maps to 0).

Triangle = (color, A, B, C): color 0 = red = acute golden triangle
(apex A, angle 36 deg, legs AB = AC = 1, base BC = 1/phi) -- half of a THIN rhomb;
color 1 = blue = obtuse golden triangle (apex A, 108 deg, legs 1, base BC = phi)
-- half of a FAT rhomb.  The order (B, C) encodes the Robinson marking (chirality).

One substitution step = inflate by phi (exact in Z[zeta]) then subdivide:
  red  (A,B,C): P = phiA + (B-A);        -> red(C',P,B'), blue(P,C',A')
  blue (A,B,C): Q = phiB + (A-B), R = phiB + (C-B);
                                          -> blue(R,C',A'), blue(Q,R,B'), red(R,Q,A')
(this is the standard subdivision, cf. preshing.com 'Penrose tiling explained').
Every triangle stays congruent to the level-0 ones; the patch grows by phi.
"""
import numpy as np
from collections import defaultdict
from pt import EVEC

ONES = np.ones(5, dtype=np.int64)


def canon(K):
    s = sum(K)
    q = s // 5
    if q == 0:
        return tuple(K)
    return (K[0] - q, K[1] - q, K[2] - q, K[3] - q, K[4] - q)


def mulphi(K):
    return canon((K[0] + K[4] + K[1], K[1] + K[0] + K[2], K[2] + K[1] + K[3],
                  K[3] + K[2] + K[4], K[4] + K[3] + K[0]))


def add(K, L):
    return canon((K[0] + L[0], K[1] + L[1], K[2] + L[2], K[3] + L[3], K[4] + L[4]))


def sub(K, L):
    return canon((K[0] - L[0], K[1] - L[1], K[2] - L[2], K[3] - L[3], K[4] - L[4]))


def rot36(K):
    """multiply by -zeta^3 = exp(i pi/5): rotation by 36 degrees.  (roll by 3: out[j] = -K[j-3])"""
    return canon((-K[2], -K[3], -K[4], -K[0], -K[1]))


def embed(K):
    return np.asarray(K, float) @ EVEC


def unit(k):
    """direction exp(i pi k/5) as a Z^5 vector (python ints)."""
    v = [0, 0, 0, 0, 0]
    v[(3 * k) % 5] = 1 if k % 2 == 0 else -1
    return canon(tuple(v))


def seed_sun():
    """10 red triangles around the origin (preshing's wheel)."""
    O = (0, 0, 0, 0, 0)
    tris = []
    for i in range(10):
        B, C = unit(i), unit(i + 1)
        if i % 2 == 0:
            B, C = C, B
        tris.append((0, O, B, C))
    return tris


def subdivide(tris):
    out = []
    for col, A, B, C in tris:
        A2, B2, C2 = mulphi(A), mulphi(B), mulphi(C)
        if col == 0:
            P = add(A2, sub(B, A))
            out.append((0, C2, P, B2))
            out.append((1, P, C2, A2))
        else:
            Q = add(B2, sub(A, B))
            R = add(B2, sub(C, B))
            out.append((1, R, C2, A2))
            out.append((1, Q, R, B2))
            out.append((0, R, Q, A2))
    return out


def subdivide_tracked(tris):
    """like subdivide but returns (child, parent_index) pairs."""
    out = []
    for pi, (col, A, B, C) in enumerate(tris):
        for ch in subdivide([(col, A, B, C)]):
            out.append((ch, pi))
    return out


def merge_rhombs(tris, want_parents=None):
    """Pair triangles sharing their base BC (same colour) into rhombs.
    Returns list of rhombs: dict(tris=(i,j), verts=[4 keys CCW], fat=bool, tri_idx=(i,j)).
    Triangles whose partner is missing (patch boundary) are dropped."""
    bases = defaultdict(list)
    for i, (col, A, B, C) in enumerate(tris):
        key = (min(B, C), max(B, C), col)
        bases[key].append(i)
    rhombs = []
    for key, lst in bases.items():
        if len(lst) != 2:
            if len(lst) > 2:
                raise RuntimeError("base edge shared by >2 triangles")
            continue
        i, j = lst
        col, A1, B1, C1 = tris[i]
        _, A2, B2, C2 = tris[j]
        assert A1 != A2
        # rhomb vertices A1, B, A2, C (a quadrilateral); orient CCW
        quad = [A1, B1, A2, C1]
        P = np.array([embed(q) for q in quad])
        area = 0.0
        for k in range(4):
            x1, y1 = P[k]; x2, y2 = P[(k + 1) % 4]
            area += x1 * y2 - x2 * y1
        if area < 0:
            quad = [A1, C1, A2, B1]
        rhombs.append({"verts": quad, "fat": col == 1, "tri_idx": (i, j)})
    return rhombs


def rhomb_arrays(rhombs):
    keys = np.array([[list(v) for v in r["verts"]] for r in rhombs], dtype=np.int64)
    verts = keys @ EVEC
    fat = np.array([r["fat"] for r in rhombs])
    return keys, verts, fat


def patch_from_seed(n):
    tris = seed_sun()
    for _ in range(n):
        tris = subdivide(tris)
    return tris


# ----------------------------------------------------------------------------
# coronas
# ----------------------------------------------------------------------------
def vertex_tiles(rhombs):
    vt = defaultdict(list)
    for t, r in enumerate(rhombs):
        for v in r["verts"]:
            vt[v].append(t)
    return vt


def tile_angle_units(rhomb, v):
    """interior angle of rhomb at vertex v in units of 36 degrees."""
    q = rhomb["verts"]
    i = q.index(v)
    P = np.array([embed(x) for x in q])
    a = P[(i - 1) % 4] - P[i]; b = P[(i + 1) % 4] - P[i]
    ang = np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))
    u = int(round(ang / 36))
    assert abs(u * 36 - ang) < 1e-6
    return u


def full_vertex(rhombs, vt, v):
    return sum(tile_angle_units(rhombs[t], v) for t in vt[v]) == 10


def corona(rhombs, vt, t):
    """indices of tiles sharing a vertex with t (including t); None if incomplete."""
    r = rhombs[t]
    for v in r["verts"]:
        if not full_vertex(rhombs, vt, v):
            return None
    s = set()
    for v in r["verts"]:
        s.update(vt[v])
    return sorted(s)


def marked_corona_canon(tris, rhombs, tiles, center):
    """Canonical form of the marked corona up to translation and rotation by 36 deg.
    Representative-independent: minimum over all (origin vertex, rotation) of the sorted
    tuple of triangles (color, A-o, B-o, C-o, is_center) with differences canonicalised."""
    items = []
    for t in tiles:
        for ti in rhombs[t]["tri_idx"]:
            col, A, B, C = tris[ti]
            items.append((col, A, B, C, int(t == center)))
    return _canon_items(items, 10)


def _canon_items(items, nrot):
    """items: list of (col, v1, v2, v3, flag) or (v1..v4 as tuple, flag) generic: handled by
    treating every element that is a 5-tuple of ints as a vertex."""
    best = None
    cur = items
    step = 10 // nrot
    for k in range(nrot):
        if k:
            for _ in range(step):
                cur = [tuple(rot36(x) if _isvert(x) else x for x in it) for it in cur]
        verts = sorted({x for it in cur for x in it if _isvert(x)})
        for o in verts:
            O = o
            tr = tuple(sorted(tuple(sub(x, O) if _isvert(x) else x for x in it) for it in cur))
            if best is None or tr < best:
                best = tr
    return best


def _isvert(x):
    return isinstance(x, tuple) and len(x) == 5 and all(isinstance(a, (int, np.integer)) for a in x)


def corona_atlas(tris, rhombs=None):
    """all complete marked coronas in the patch, as dict canon -> (tile index, tiles)."""
    if rhombs is None:
        rhombs = merge_rhombs(tris)
    vt = vertex_tiles(rhombs)
    atlas = {}
    for t in range(len(rhombs)):
        cor = corona(rhombs, vt, t)
        if cor is None:
            continue
        key = marked_corona_canon(tris, rhombs, cor, t)
        if key not in atlas:
            atlas[key] = (t, cor)
    return atlas, rhombs, vt


def corona_triangles(tris, rhombs, tiles):
    """list of (triangle, is_center_tile_triangle) for the given tiles."""
    out = []
    for t in tiles:
        for ti in rhombs[t]["tri_idx"]:
            out.append(tris[ti])
    return out


def substitute_patch(tris_list, n):
    """apply n substitution steps to a list of triangles, tracking ancestry to the
    original triangles.  Returns (tris, ancestor index list)."""
    cur = list(tris_list)
    anc = list(range(len(cur)))
    for _ in range(n):
        nxt = []; nanc = []
        for (t, pi) in subdivide_tracked(cur):
            nxt.append(t); nanc.append(anc[pi])
        cur, anc = nxt, nanc
    return cur, anc


def unmarked_corona_canon(rhombs, tiles, center, nrot=10):
    """canonical form of the unmarked corona up to translation+rotation (36 deg, or 72 if nrot=5)."""
    items = []
    for t in tiles:
        q = rhombs[t]["verts"]
        items.append((q[0], q[1], q[2], q[3], int(t == center)))
    # order-free rhomb: sort the 4 translated vertices inside _canon_items? do it here by
    # canonicalising each rhomb as a sorted 4-tuple after translation -> use helper
    best = None
    cur = items
    step = 10 // nrot
    for k in range(nrot):
        if k:
            for _ in range(step):
                cur = [(rot36(a), rot36(b), rot36(c), rot36(d), f) for (a, b, c, d, f) in cur]
        verts = sorted({x for it in cur for x in it[:4]})
        for o in verts:
            O = o
            tr = tuple(sorted((tuple(sorted(sub(x, O) for x in it[:4])), it[4]) for it in cur))
            if best is None or tr < best:
                best = tr
    return best
