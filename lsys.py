"""
lsys.py -- scattering of the turtle through standard level-m supertriangles.

Standard marked triangles (type = (col, sign)):  col 0 = acute (red), 1 = obtuse (blue);
sign = orientation of (A,B,C) (+1 counter-clockwise).  Standard position: A = 0, B = 1 (= unit(0)).
Level-m supertile: substitute the standard triangle m times (coordinates scale by phi^m, no rotation).

Rhomb ownership: every rhomb consists of two mirror-image fine triangles; it is OWNED by the
supertile containing its positively oriented half.  R+(X) = rhombs owned by X.  Boundary triangles
whose base lies on the boundary of X get a ghost partner (mirror image across the base), so the
geometry of every rhomb meeting X is determined by X alone.

A crossing is (tail, head, i): the turtle executes its i-th command (an S) and crosses the rhomb edge
tail->head, written in the counter-clockwise order of the rhomb being left.
A pass of X is a maximal run of turtle states in rhombs owned by X: entry crossing -> exit crossing.
"""
import numpy as np
from subst import canon, mulphi, add, sub, rot36, unit, embed
from subst_fast import tris_to_arrays, subdivide_v, vcanon, vmulphi, _lexsort_rows, _rowmin_rowmax
from pt import EVEC, parse_prog

PHI = (1 + 5 ** 0.5) / 2
TYPES = [(0, 1), (0, -1), (1, 1), (1, -1)]


def std_triangle(col, sign):
    O = (0, 0, 0, 0, 0)
    b = unit(0)
    c = unit(1) if col == 0 else unit(3)
    return (col, O, b, c) if sign > 0 else (col, O, c, b)


def tri_sign(A, B, C):
    a, b, c = embed(A), embed(B), embed(C)
    return 1 if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0 else -1


def children(tri):
    """level-1 children of a marked triangle, in the coordinates of the inflated parent."""
    col, A, B, C = tris_to_arrays([tri])
    c, a, b, cc, anc = subdivide_v(col, A, B, C)
    return [(int(c[i]), tuple(int(x) for x in a[i]), tuple(int(x) for x in b[i]), tuple(int(x) for x in cc[i])) for i in range(len(c))]


def frame_of(tri):
    """orientation-preserving isometry g with g(standard triangle of the same type) = tri.
    returns (k, t): x -> rot36^k(x) + t, plus the type."""
    col, A, B, C = tri
    s = tri_sign(A, B, C)
    std = std_triangle(col, s)
    d = sub(B, A)
    # find k with rot36^k(std.B) == d   (|B-A| = 1 for level-0 triangles)
    v = std[2]
    for k in range(10):
        if v == d:
            # check C as well
            w = std[3]
            for _ in range(k):
                w = rot36(w)
            assert add(w, A) == C, "frame mismatch"
            return (col, s), k, A
        v = rot36(v)
    raise RuntimeError("no frame")


def apply_frame(k, t, scale_m, K):
    """g_m(K) = rot36^k(K) + phi^scale_m * t   (t given at unit scale)."""
    v = K
    for _ in range(k % 10):
        v = rot36(v)
    tt = t
    for _ in range(scale_m):
        tt = mulphi(tt)
    return add(v, tt)


def inv_frame(k, t, scale_m, K):
    tt = t
    for _ in range(scale_m):
        tt = mulphi(tt)
    v = sub(K, tt)
    for _ in range((10 - k) % 10):
        v = rot36(v)
    return v


class Supertile:
    """fine structure of the standard level-m supertile of a given type."""

    def __init__(self, typ, m):
        self.typ = typ; self.m = m
        tri = std_triangle(*typ)
        if m == 0:
            col, A, B, C = tris_to_arrays([tri]); slot = np.zeros(1, int)
        else:
            ch = children(tri)
            col, A, B, C = tris_to_arrays(ch)
            slot = np.arange(len(ch))
            for _ in range(m - 1):
                col, A, B, C, anc = subdivide_v(col, A, B, C)
                slot = slot[anc]
        self.col, self.A, self.B, self.C, self.slot = col, A, B, C, slot
        n = len(col)
        PA, PB, PC = A @ EVEC, B @ EVEC, C @ EVEC
        self.sign = np.where((PB[:, 0] - PA[:, 0]) * (PC[:, 1] - PA[:, 1]) - (PB[:, 1] - PA[:, 1]) * (PC[:, 0] - PA[:, 0]) > 0, 1, -1)
        # pair triangles by base
        lo, hi = _rowmin_rowmax(B, C)
        Mx = np.concatenate([lo, hi], axis=1)
        order = _lexsort_rows(Mx); Ms = Mx[order]
        same = np.all(Ms[1:] == Ms[:-1], axis=1)
        partner = -np.ones(n, int)
        f = np.where(same)[0]
        partner[order[f]] = order[f + 1]; partner[order[f + 1]] = order[f]
        # rhombs: one per positive triangle (owned).  partner inside or ghost.
        pos = np.where(self.sign > 0)[0]
        self.owner_tri = pos                        # rhomb r <-> positive triangle pos[r]
        Ap = A[pos]; Bp = B[pos]; Cp = C[pos]
        Aq = vcanon(Bp + Cp - Ap)                   # apex of the partner (mirror image across the base)
        quad = np.stack([Ap, Bp, Aq, Cp], axis=1)   # positive triangle (A,B,C) is CCW -> A,B,A*,C is CCW
        self.quad = quad
        self.nr = len(pos)
        self.rslot = slot[pos]
        # non-owned rhombs meeting D (negative triangle in D with ghost positive partner) are not needed:
        # a move into them is an exit.  Navigation among owned rhombs only:
        from subst_fast import nav_from_keys
        self.eN, self.eE = nav_from_keys(quad)
        self.cent = (quad @ EVEC).mean(1)

    def key(self, r, v):
        return tuple(int(x) for x in self.quad[r, v % 4])

    def entries(self, nS):
        """all entry states: (rhomb, heading edge, next command index, crossing)"""
        out = []
        br, be = np.where(self.eN < 0)
        for r, e in zip(br.tolist(), be.tolist()):
            for i in range(nS):
                cross = (self.key(r, e + 1), self.key(r, e), i)
                out.append((r, (e + 2) & 3, i + 1, cross))
        return out

    def run_pass(self, prog, r, e, ci, max_steps=10 ** 8, record=False):
        """run from state (rhomb r, heading e, next command index ci) until exit.
        returns (exit crossing, number of S moves, list of (rhomb) if record) or (None, ...) if trapped."""
        cmds = prog; n = len(cmds)
        eN, eE = self.eN, self.eE
        steps = 0; path = [r] if record else None
        s0 = (r, e, ci % n); first = True
        while steps < max_steps:
            c = cmds[ci % n]
            if c == 0:
                nb = eN[r, e]
                if nb < 0:
                    return (self.key(r, e), self.key(r, e + 1), ci % n), steps, path
                e = (eE[r, e] + 2) & 3; r = int(nb); steps += 1
                if record: path.append(r)
            elif c == 1:
                e = (e + 3) & 3
            else:
                e = (e + 1) & 3
            ci += 1
            if (r, e, ci % n) == s0:
                return None, steps, path
        return None, steps, path


def scattering(typ, m, prog_str):
    st = Supertile(typ, m)
    prog = parse_prog(prog_str).tolist()
    nS = [i for i, c in enumerate(prog) if c == 0]
    res = []
    br, be = np.where(st.eN < 0)
    for r, e in zip(br.tolist(), be.tolist()):
        for i in nS:
            cross_in = (st.key(r, e + 1), st.key(r, e), i)
            ex, steps, _ = st.run_pass(prog, r, (e + 2) & 3, i + 1)
            res.append((cross_in, ex, steps, (r, (e + 2) & 3, i + 1)))
    return st, res


if __name__ == "__main__":
    import sys
    prog = sys.argv[1]; m0 = int(sys.argv[2]); m1 = int(sys.argv[3])
    for typ in TYPES:
        for m in range(m0, m1 + 1):
            st, res = scattering(typ, m, prog)
            res.sort(key=lambda t: -t[2])
            top = res[:4]
            def fmt(cr):
                if cr is None: return "trapped"
                p = (embed(cr[0]) + embed(cr[1])) / 2 / PHI ** m
                return f"({p[0]:+.3f},{p[1]:+.3f};S{cr[2]})"
            print(f"type {typ} level {m}: tris {len(st.col):6d} owned rhombs {st.nr:6d} entries {len(res):5d} | longest: " +
                  "  ".join(f"{t[2]}: {fmt(t[0])}->{fmt(t[1])}" for t in top))
