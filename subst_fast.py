"""
subst_fast.py -- vectorised (numpy int64) version of the exact Robinson-triangle substitution
of subst.py.  Triangles: col (n,), A,B,C (n,5) canonical Z^5 keys (sum in 0..4).
"""
import numpy as np
from subst import ONES

ONESV = np.ones(5, dtype=np.int64)


def vcanon(K):
    K = np.asarray(K, np.int64)
    s = K.sum(-1, keepdims=True)
    return K - ONESV * (s // 5)


def vmulphi(K):
    K = np.asarray(K, np.int64)
    return vcanon(K + np.roll(K, 1, axis=-1) + np.roll(K, -1, axis=-1))


def tris_to_arrays(tris):
    col = np.array([t[0] for t in tris], np.int64)
    A = np.array([t[1] for t in tris], np.int64)
    B = np.array([t[2] for t in tris], np.int64)
    C = np.array([t[3] for t in tris], np.int64)
    return col, vcanon(A), vcanon(B), vcanon(C)


def subdivide_v(col, A, B, C, anc=None):
    """one inflate+subdivide step on arrays; returns new arrays (+ ancestry)."""
    A2, B2, C2 = vmulphi(A), vmulphi(B), vmulphi(C)
    red = col == 0
    blue = ~red
    if anc is None:
        anc = np.arange(len(col))
    # red: P = A2 + (B - A) ; children red(C2,P,B2), blue(P,C2,A2)
    P = vcanon(A2[red] + (B[red] - A[red]))
    r_col = np.concatenate([np.zeros(red.sum(), np.int64), np.ones(red.sum(), np.int64)])
    r_A = np.concatenate([C2[red], P]); r_B = np.concatenate([P, C2[red]]); r_C = np.concatenate([B2[red], A2[red]])
    r_anc = np.concatenate([anc[red], anc[red]])
    # blue: Q = B2 + (A-B), R = B2 + (C-B); children blue(R,C2,A2), blue(Q,R,B2), red(R,Q,A2)
    Q = vcanon(B2[blue] + (A[blue] - B[blue])); R = vcanon(B2[blue] + (C[blue] - B[blue]))
    b_col = np.concatenate([np.ones(blue.sum(), np.int64), np.ones(blue.sum(), np.int64), np.zeros(blue.sum(), np.int64)])
    b_A = np.concatenate([R, Q, R]); b_B = np.concatenate([C2[blue], R, Q]); b_C = np.concatenate([A2[blue], B2[blue], A2[blue]])
    b_anc = np.concatenate([anc[blue], anc[blue], anc[blue]])
    return (np.concatenate([r_col, b_col]), np.concatenate([r_A, b_A]), np.concatenate([r_B, b_B]),
            np.concatenate([r_C, b_C]), np.concatenate([r_anc, b_anc]))


def substitute_v(col, A, B, C, n):
    anc = np.arange(len(col))
    for _ in range(n):
        col, A, B, C, anc = subdivide_v(col, A, B, C, anc)
    return col, A, B, C, anc


def _rowkey(X):
    """pack a (n,5) int64 array into one int64 key per row (12 bits per coordinate).
    Only used where the range is asserted; the matching functions below use exact
    column-wise lexsort and do not depend on packing."""
    X = np.asarray(X, np.int64)
    if X.size and (X.min() < -2048 or X.max() > 2047):
        raise OverflowError("coordinate outside the 12-bit packing range")
    X = X + (1 << 11)
    k = np.zeros(len(X), np.int64)
    for j in range(5):
        k = (k << 12) | X[:, j]
    return k


def _lexsort_rows(M):
    """argsort of the rows of an integer (n,m) array in lexicographic order (exact)."""
    return np.lexsort(tuple(M[:, j] for j in range(M.shape[1] - 1, -1, -1)))


def _rowmin_rowmax(P, Q):
    """elementwise lexicographic min/max of two (n,5) integer arrays (row-wise)."""
    # compare rows lexicographically
    less = np.zeros(len(P), bool); decided = np.zeros(len(P), bool)
    for j in range(P.shape[1]):
        lt = (P[:, j] < Q[:, j]) & ~decided
        gt = (P[:, j] > Q[:, j]) & ~decided
        less |= lt; decided |= lt | gt
    lo = np.where(less[:, None], P, Q); hi = np.where(less[:, None], Q, P)
    return lo, hi


def merge_v(col, A, B, C, anc=None):
    """pair triangles sharing base (B,C) and colour into rhombs (exact lexsort matching on
    the full integer coordinates; no packing).  Returns keys (m,4,5) CCW vertex keys,
    fat (m,), tri_idx (m,2)."""
    lo, hi = _rowmin_rowmax(np.asarray(B, np.int64), np.asarray(C, np.int64))
    M = np.concatenate([lo, hi, col[:, None]], axis=1)  # (n,11)
    order = _lexsort_rows(M)
    Ms = M[order]
    same = np.all(Ms[1:] == Ms[:-1], axis=1)
    first = np.where(same)[0]
    if len(first) > 1 and np.any(first[1:] == first[:-1] + 1):
        raise RuntimeError("base shared by >2 triangles")
    i1 = order[first]; i2 = order[first + 1]
    quad = np.stack([A[i1], B[i1], A[i2], C[i1]], axis=1)  # (m,4,5)
    from pt import EVEC
    P = quad @ EVEC
    x = P[:, :, 0]; y = P[:, :, 1]
    area = (x * np.roll(y, -1, 1) - np.roll(x, -1, 1) * y).sum(1)
    neg = area < 0
    quad[neg] = quad[neg][:, [0, 3, 2, 1]]
    return quad, col[i1] == 1, np.stack([i1, i2], 1)


def nav_from_keys(quad):
    """edgeNbr/edgeEntry from (m,4,5) integer keys (exact column-wise lexsort; no packing)."""
    m = quad.shape[0]
    a = quad.reshape(-1, 5)                       # vertex e of tile t at row 4t+e
    b = np.roll(quad, -1, axis=1).reshape(-1, 5)  # vertex e+1
    lo, hi = _rowmin_rowmax(a, b)
    M = np.concatenate([lo, hi], axis=1)          # (4m,10)
    order = _lexsort_rows(M)
    Ms = M[order]
    same = np.all(Ms[1:] == Ms[:-1], axis=1)
    first = np.where(same)[0]
    if len(first) > 1 and np.any((first[1:] == first[:-1] + 1)):
        raise RuntimeError("edge shared by >2 tiles")
    p1 = order[first]; p2 = order[first + 1]
    t1, e1 = p1 // 4, p1 % 4; t2, e2 = p2 // 4, p2 % 4
    eN = -np.ones((m, 4), np.int32); eE = -np.ones((m, 4), np.int32)
    eN[t1, e1] = t2; eE[t1, e1] = e2; eN[t2, e2] = t1; eE[t2, e2] = e1
    return eN, eE
