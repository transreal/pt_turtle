"""
verify_cert.py -- independent checker of a pass-substitution certificate (cert2_PROG.json).

It re-implements everything it needs (exact Z[zeta] arithmetic, Robinson-triangle substitution,
rhomb ownership, turtle simulation) with plain Python integers and tuples; it shares no code with
the extractor.

Statement H(m) proved by induction on m >= B for every pass family pi whose parity equals m mod 2:
  Let X be the standard level-m supertriangle of type tau(pi).  The crossing port(pi.in, m) leads
  into a rhomb owned by X; starting there the turtle performs exactly N(pi, m) S-moves, the first
  N-1 of them between rhombs owned by X, and the N-th is the crossing port(pi.out, m), which leads
  to a rhomb NOT owned by X.
Obligations checked here:
  O1  H(B) by direct simulation (families of parity B mod 2).
  O2  port identities (a),(b),(c) of every decomposition at the two lowest levels of each parity
      class above the base (they then hold at all levels because both sides obey u_{m+2}=phi^2 u_m+const).
  O3  closure and type/parity consistency of the pass table.
  O5  exit property at all higher levels: the scaled limit point of the exit port lies on an edge of
      the last child which is contained in a side of the parent, with sufficient margins.
  O4  growth: a root family occurs in its own two-level expansion together with other passes.
Extra (redundant) simulations at higher levels are run as a consistency check.
"""
import sys, json, math
from collections import Counter

# ---------------------------------------------------------------- exact arithmetic in Z[zeta]
def canon(K):
    q = sum(K) // 5
    return tuple(k - q for k in K)
def add(K, L): return canon(tuple(a + b for a, b in zip(K, L)))
def sub(K, L): return canon(tuple(a - b for a, b in zip(K, L)))
def neg(K): return canon(tuple(-a for a in K))
def mulphi(K): return canon(tuple(K[j] + K[(j - 1) % 5] + K[(j + 1) % 5] for j in range(5)))   # phi = 1 + zeta + zeta^4
def invphi(K): return sub(mulphi(K), K)                                                      # 1/phi = phi - 1
def rot36(K): return canon(tuple(-K[(j - 3) % 5] for j in range(5)))                          # multiply by -zeta^3
def zmul(K, L):
    out = [0] * 5
    for i in range(5):
        for j in range(5):
            out[(i + j) % 5] += K[i] * L[j]
    return canon(tuple(out))
def conj(K): return canon(tuple(K[(-j) % 5] for j in range(5)))
def is_real_multiple(a, b):
    """True iff a and b (as complex numbers) are parallel: Im(conj(a) b) = 0."""
    return zmul(conj(a), b) == zmul(a, conj(b))
C5 = [math.cos(2 * math.pi * j / 5) for j in range(5)]; S5 = [math.sin(2 * math.pi * j / 5) for j in range(5)]
def emb(K): return (sum(k * c for k, c in zip(K, C5)), sum(k * s for k, s in zip(K, S5)))
def norm(K): x, y = emb(K); return math.hypot(x, y)
def unit(k):
    v = [0] * 5; v[(3 * k) % 5] = 1 if k % 2 == 0 else -1
    return canon(tuple(v))
ZERO = (0, 0, 0, 0, 0)
PHI = (1 + 5 ** 0.5) / 2

# ---------------------------------------------------------------- triangles and substitution
def std_triangle(col, sign):
    b = unit(0); c = unit(1) if col == 0 else unit(3)
    return (col, ZERO, b, c) if sign > 0 else (col, ZERO, c, b)
def orientation(A, B, C):
    a, b, c = emb(A), emb(B), emb(C)
    return 1 if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0 else -1
def subdivide_one(tri):
    col, A, B, C = tri
    A2, B2, C2 = mulphi(A), mulphi(B), mulphi(C)
    if col == 0:
        P = add(A2, sub(B, A))
        return [(0, C2, P, B2), (1, P, C2, A2)]
    Q = add(B2, sub(A, B)); R = add(B2, sub(C, B))
    return [(1, R, C2, A2), (1, Q, R, B2), (0, R, Q, A2)]
def child_frames(typ):
    """for each child slot: (child type, k, t) with child = rot36^k(standard child type) + t (unit scale)."""
    out = []
    for (col, A, B, C) in subdivide_one(std_triangle(*typ)):
        s = orientation(A, B, C); std = std_triangle(col, s)
        d = sub(B, A); v = std[2]; k = None
        for kk in range(10):
            if v == d: k = kk; break
            v = rot36(v)
        assert k is not None
        w = std[3]
        for _ in range(k): w = rot36(w)
        assert add(w, A) == C
        out.append(((col, s), k, A, (A, B, C)))
    return out
TYPES = [(0, 1), (0, -1), (1, 1), (1, -1)]
FR = {t: child_frames(t) for t in TYPES}
def G(typ, slot, m, K):
    """position in the level-m parent frame of the point K given in the frame of child `slot` (level m-1)."""
    _, k, t, _ = FR[typ][slot]
    v = K
    for _ in range(k): v = rot36(v)
    tt = t
    for _ in range(m - 1): tt = mulphi(tt)
    return add(v, tt)

class Supertile:
    def __init__(self, typ, m):
        tris = [std_triangle(*typ)]
        slots = [0]
        for lev in range(m):
            nt, ns = [], []
            for t, s in zip(tris, slots):
                ch = subdivide_one(t)
                for j, c in enumerate(ch):
                    nt.append(c); ns.append(j if lev == 0 else s)
            tris, slots = nt, ns
        self.rhombs = []; self.slot = []
        for (col, A, B, C), s in zip(tris, slots):
            if orientation(A, B, C) > 0:
                Astar = sub(add(B, C), A)
                self.rhombs.append((A, B, Astar, C)); self.slot.append(s)
        self.edge = {}
        for r, q in enumerate(self.rhombs):
            for e in range(4):
                key = (q[e], q[(e + 1) % 4])
                assert key not in self.edge
                self.edge[key] = (r, e)
    def entry_state(self, cross):
        tail, head, i = cross
        if (head, tail) not in self.edge: return None          # no owned rhomb on the far side
        if (tail, head) in self.edge: return None              # near side is owned too: not an entry
        r, eb = self.edge[(head, tail)]
        return (r, (eb + 2) % 4, i + 1)
    def run(self, prog, state, limit=10 ** 8):
        r, e, ci = state; n = len(prog); N = 0; seen0 = (r, e, ci % n)
        while N < limit:
            c = prog[ci % n]
            if c == "S":
                q = self.rhombs[r]; a, b = q[e], q[(e + 1) % 4]
                N += 1
                if (b, a) not in self.edge:
                    return (a, b, ci % n), N
                r, e2 = self.edge[(b, a)]; e = (e2 + 2) % 4
            elif c == "R": e = (e + 3) % 4
            else: e = (e + 1) % 4
            ci += 1
            if (r, e, ci % n) == seen0: return None, N
        return None, N

# ---------------------------------------------------------------- certificate
cert = json.load(open(sys.argv[1]))
prog = cert["prog"]; B = cert["B"]
assert isinstance(prog, str) and len(prog) > 0 and set(prog) <= set("SRL") and "S" in prog, "program must be a word over S, R, L"
print(f"program ({prog})*   base level B = {B}")
def Rq(q): return B if B % 2 == q else B + 1
PORT = {int(k): v for k, v in cert["ports"].items()}
PASS = {int(k): v for k, v in cert["passes"].items()}
def port_at(pid, m):
    f = PORT[pid]; q = f["parity"]; mm = Rq(q)
    assert m >= mm and (m - mm) % 2 == 0, (pid, m)
    tail, head, i = tuple(f["ref"][0]), tuple(f["ref"][1]), f["ref"][2]
    dt, dh = tuple(f["delta"][0]), tuple(f["delta"][1])
    while mm < m:
        tail = add(mulphi(mulphi(tail)), dt); head = add(mulphi(mulphi(head)), dh); mm += 2
    return (tail, head, i)
def ptype(pid): return tuple(PORT[pid]["typ"])
errors = []
def check(cond, msg):
    if not cond: errors.append(msg)

# O3: closure / consistency
for pin, p in PASS.items():
    typ = ptype(pin); q = PORT[pin]["parity"]
    check(ptype(p["out"]) == typ and PORT[p["out"]]["parity"] == q, f"pass {pin}: exit port type/parity")
    check(len(p["children"]) >= 1, f"pass {pin}: no children")
    prev = None
    for slot, a, b in p["children"]:
        ctyp = FR[typ][slot][0]
        check(a in PASS, f"pass {pin}: child {a} not in table")
        if a in PASS: check(PASS[a]["out"] == b, f"pass {pin}: child exit mismatch")
        check(ptype(a) == ctyp and ptype(b) == ctyp, f"pass {pin}: child type mismatch")
        check(PORT[a]["parity"] == 1 - q and PORT[b]["parity"] == 1 - q, f"pass {pin}: child parity")
        check(slot != prev, f"pass {pin}: consecutive child passes in the same slot")
        prev = slot
    for f in (PORT[pin], PORT[p["out"]]):
        t, h = tuple(f["ref"][0]), tuple(f["ref"][1])
        check(norm(sub(h, t)) > 0.999 and norm(sub(h, t)) < 1.001, "port is not a unit leg")
        check(sub(tuple(f["delta"][1]), tuple(f["delta"][0])) == sub(sub(h, t), mulphi(mulphi(sub(h, t)))), "delta of head and tail inconsistent")
print(f"O3 closure/consistency: {len(PASS)} pass families, {len(PORT)} port families, errors so far {len(errors)}")

# O1: base simulation at level B
cache = {}
def ST(typ, m):
    if (typ, m) not in cache: cache[(typ, m)] = Supertile(typ, m)
    return cache[(typ, m)]
N = {}
nbase = 0
for pin, p in PASS.items():
    if PORT[pin]["parity"] != B % 2: continue
    st = ST(ptype(pin), B)
    s = st.entry_state(port_at(pin, B))
    check(s is not None, f"O1 pass {pin}: entry port is not an entry crossing at level {B}")
    if s is None: continue
    ex, n = st.run(prog, s)
    check(ex == port_at(p["out"], B), f"O1 pass {pin}: exit mismatch at level {B}")
    N[(pin, B)] = n; nbase += 1
print(f"O1 base simulations at level {B}: {nbase} passes, errors so far {len(errors)}")

# O2: port identities at the two lowest levels of each parity class above the base
nid = 0
for pin, p in PASS.items():
    typ = ptype(pin); q = PORT[pin]["parity"]
    m1 = B + 1 if (B + 1) % 2 == q else B + 2
    for m in (m1, m1 + 2):
        ch = p["children"]
        def GP(slot, pid):
            t, h, i = port_at(pid, m - 1)
            return (G(typ, slot, m, t), G(typ, slot, m, h), i)
        check(port_at(pin, m) == GP(ch[0][0], ch[0][1]), f"O2a pass {pin} level {m}")
        for (s1, a1, b1), (s2, a2, b2) in zip(ch, ch[1:]):
            check(GP(s1, b1) == GP(s2, a2), f"O2b pass {pin} level {m}")
            nid += 1
        check(port_at(p["out"], m) == GP(ch[-1][0], ch[-1][2]), f"O2c pass {pin} level {m}")
        nid += 2
print(f"O2 port identities checked: {nid}, errors so far {len(errors)}")

# O5: exit property at all levels above the base
SIN36 = math.sin(math.pi / 5)
worst = 1e9
for pin, p in PASS.items():
    typ = ptype(pin); q = PORT[pin]["parity"]
    slot, a, b = p["children"][-1]
    ctyp, k, t, (cA, cB, cC) = FR[typ][slot]
    n0 = (B + 1 if (B + 1) % 2 == q else B + 2) - 1          # lowest child level used in an inductive step
    f = PORT[b]
    tail0 = port_at(b, n0)[0]
    c = neg(invphi(tuple(f["delta"][0])))                   # fixed point part: p_n = c + phi^(n-n0) (p_n0 - c)
    y = sub(tail0, c)
    # sanity: the affine law really gives p_{n0+2} = c + phi^2 y
    check(port_at(b, n0 + 2)[0] == add(c, mulphi(mulphi(y))), f"O5 pass {pin}: affine form")
    std = std_triangle(*ctyp); V = [std[1], std[2], std[3]]
    def scale(K, n):
        for _ in range(n): K = mulphi(K)
        return K
    ok = False
    parent = std_triangle(*typ); PV = [mulphi(parent[1]), mulphi(parent[2]), mulphi(parent[3])]   # inflated parent (level-1 coordinates)
    for iu in range(3):
        U, W = V[iu], V[(iu + 1) % 3]
        Un, Wn = scale(U, n0), scale(W, n0)
        if not is_real_multiple(sub(y, Un), sub(Wn, Un)): continue
        yu, wu = emb(sub(y, Un)), emb(sub(Wn, Un))
        s = (yu[0] * wu[0] + yu[1] * wu[1]) / (wu[0] ** 2 + wu[1] ** 2)
        if not (0 < s < 1): continue
        dU, dW = norm(sub(y, Un)), norm(sub(y, Wn)); cc = norm(c)
        need = (1.903 + cc) / SIN36 + 0.05
        # child edge (U,W) must lie on a side of the parent (exact collinearity with a parent side, within it)
        def place(K):
            v = K
            for _ in range(k): v = rot36(v)
            return add(v, t)
        gU, gW = place(U), place(W)
        onside = False
        for ip in range(3):
            P1, P2 = PV[ip], PV[(ip + 1) % 3]
            if is_real_multiple(sub(gU, P1), sub(P2, P1)) and is_real_multiple(sub(gW, P1), sub(P2, P1)):
                e1, e2, ee = emb(sub(gU, P1)), emb(sub(gW, P1)), emb(sub(P2, P1))
                l2 = ee[0] ** 2 + ee[1] ** 2
                s1 = (e1[0] * ee[0] + e1[1] * ee[1]) / l2; s2 = (e2[0] * ee[0] + e2[1] * ee[1]) / l2
                if -1e-9 <= s1 <= 1 + 1e-9 and -1e-9 <= s2 <= 1 + 1e-9: onside = True
        if onside and min(dU, dW) >= need:
            ok = True; worst = min(worst, min(dU, dW) - need)
    check(ok, f"O5 pass {pin}: exit port of last child is not safely on a boundary edge (n0={n0})")
print(f"O5 exit property certificates: errors so far {len(errors)}; smallest margin surplus {worst:.2f}")

# lengths by recursion, redundant simulations
def Nrec(pin, m):
    if (pin, m) in N: return N[(pin, m)]
    assert m > B
    v = sum(Nrec(a, m - 1) for _, a, _ in PASS[pin]["children"])
    N[(pin, m)] = v; return v
nsim = 0
for extra in ((B + 1, B + 2) if '--fast' in sys.argv else (B + 1, B + 2, B + 3, B + 4)):
    for pin, p in PASS.items():
        if PORT[pin]["parity"] != extra % 2: continue
        st = ST(ptype(pin), extra)
        s = st.entry_state(port_at(pin, extra))
        check(s is not None, f"sim pass {pin} level {extra}: not an entry")
        if s is None: continue
        ex, n = st.run(prog, s)
        check(ex == port_at(p["out"], extra), f"sim pass {pin} level {extra}: exit mismatch")
        check(n == Nrec(pin, extra), f"sim pass {pin} level {extra}: length {n} != recursion {Nrec(pin, extra)}")
        nsim += 1
print(f"redundant simulations at levels {B+1}..{B+4}: {nsim}, errors so far {len(errors)}")

# O4: growth
def expand2(pin):
    out = Counter()
    for _, a, _ in PASS[pin]["children"]:
        for _, a2, _ in PASS[a]["children"]: out[a2] += 1
    return out
grow = []
selfrep = []
for pin in PASS:
    ex = expand2(pin); mult = ex[pin]; total = sum(ex.values())
    if mult >= 1 and total >= mult + 1: selfrep.append((pin, mult, total))
check(len(selfrep) > 0, "O4: no self-reproducing pass family (no growth certificate)")
def first_level(pin):
    q = PORT[pin]["parity"]; return B if q == B % 2 else B + 1
best = sorted(selfrep, key=lambda t: -Nrec(t[0], first_level(t[0]) + 8))[:4]
for pin, mult, total in best:
    m0_ = first_level(pin)
    print(f"O4 self-reproducing pass {pin} type {ptype(pin)} parity {PORT[pin]['parity']}: occurs {mult}x in its own two-level expansion of {total} passes; "
          f"N(m) for m={m0_},{m0_+2},... = {[Nrec(pin, m) for m in range(m0_, m0_ + 14, 2)]}")
print(f"O4: {len(selfrep)} self-reproducing pass families (N(m+2) >= mult*N(m) + (total-mult) -> unbounded)")
# O6: exact recurrence N(m+2) = a N(m) + b.  The vector of lengths of the families of one parity
# satisfies N(m+2) = M N(m) with a fixed nonnegative integer matrix M (two-level expansion), so
# y_s = N(pi, m0+2s+2) - a N(pi, m0+2s) - b obeys a linear recurrence of order <= dim(M)+1; if it
# vanishes for dim(M)+2 consecutive s it vanishes for all s (Cayley-Hamilton).
def recurrence(pin):
    q = PORT[pin]["parity"]; m0_ = first_level(pin)
    dim = sum(1 for p in PASS if PORT[p]["parity"] == q)
    K = dim + 3
    seq = [Nrec(pin, m0_ + 2 * s) for s in range(K + 1)]
    d0, d1 = seq[1] - seq[0], seq[2] - seq[1]
    if d0 == 0 or d1 % d0 != 0: return None
    a_ = d1 // d0; b_ = seq[1] - a_ * seq[0]
    if all(seq[s + 1] == a_ * seq[s] + b_ for s in range(K)): return (a_, b_, m0_, seq[0], K)
    return None
for r in cert["roots"]:
    check(r in PASS, f"root {r} is not a pass family")
    if r not in PASS: continue
    ex = expand2(r); mult = ex[r]; total = sum(ex.values())
    rec = recurrence(r)
    m0_ = first_level(r)
    print(f"   root {r} type {ptype(r)} parity {PORT[r]['parity']}: self-multiplicity {mult} in a two-level expansion of {total} passes; "
          f"N = {[Nrec(r, m) for m in range(m0_, m0_ + 14, 2)]}" + (f"; O6 exact recurrence N(m+2) = {rec[0]} N(m) + {rec[1]} (N({rec[2]}) = {rec[3]}), verified for {rec[4]} steps" if rec else "; no affine recurrence"))
for pin, mult, total in best:
    rec = recurrence(pin)
    if rec: print(f"O6 pass {pin} type {ptype(pin)} parity {PORT[pin]['parity']}: N(m+2) = {rec[0]} N(m) + {rec[1]}, N({rec[2]}) = {rec[3]}  [verified for {rec[4]} steps > number of families of that parity + 1]")
print("RESULT:", "CERTIFICATE VALID" if not errors else f"INVALID ({len(errors)} errors)")
for e in errors[:20]: print("  ", e)
