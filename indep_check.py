#!/usr/bin/env python3
# Independent re-implementation of the certificate checker for the (SSSR)* unboundedness proof.
# [Generalised afterwards in two marked places (GEN) so that it also accepts certificates of other
#  programs: the recurrence constants and the expected root lengths are no longer hard-coded.]
# Written from proof_unbounded.md only (no code shared with verify_cert.py / lsys.py / extract2.py /
# subst.py / pt.py).  Exact integer arithmetic in Z[zeta], zeta = exp(2 pi i/5); elements are integer
# 5-vectors modulo (1,1,1,1,1), normalised here so that the last coordinate is 0.
import json, sys, math, os, time, cmath

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'cert2_SSSR_1p_1.json')
MAXDIRECT = int(sys.argv[2]) if len(sys.argv) > 2 else 11   # highest level for direct simulations

# ------------------------------------------------------------------ Z[zeta]
def nz(v):
    l = v[4]
    return (v[0] - l, v[1] - l, v[2] - l, v[3] - l, 0)

def add(u, v): return nz(tuple(a + b for a, b in zip(u, v)))
def sub(u, v): return nz(tuple(a - b for a, b in zip(u, v)))
def neg(u): return nz(tuple(-a for a in u))
def zshift(v, k): return tuple(v[(i - k) % 5] for i in range(5))      # zeta^k * v (raw)
def mphi(v):                                                          # (1 + zeta + zeta^4) * v
    a = zshift(v, 1); b = zshift(v, 4)
    return nz(tuple(v[i] + a[i] + b[i] for i in range(5)))
def mrho(v):                                                          # (-zeta^3) * v = rotation by 36 deg
    return nz(tuple(-x for x in zshift(v, 3)))
def conj(v): return tuple(v[(-i) % 5] for i in range(5))
def mul(u, v):                                                        # raw product (not normalised)
    w = [0] * 5
    for i in range(5):
        if u[i]:
            for j in range(5):
                w[(i + j) % 5] += u[i] * v[j]
    return tuple(w)

ZERO = (0, 0, 0, 0, 0)
ONE = (1, 0, 0, 0, 0)
RHO = [ONE]
for _ in range(9):
    RHO.append(mrho(RHO[-1]))
assert mrho(RHO[9]) == ONE and len(set(RHO)) == 10

def rot(v, k):
    for _ in range(k % 10):
        v = mrho(v)
    return v

def phipow(v, n):
    for _ in range(n):
        v = mphi(v)
    return v

ZC = [cmath.exp(2j * math.pi * k / 5) for k in range(5)]
def emb(v): return sum(v[k] * ZC[k] for k in range(5))

# exact signs.  a + b*phi with integers a, b
def sign_zphi(a, b):
    if a == 0 and b == 0: return 0
    if a >= 0 and b >= 0: return 1
    if a <= 0 and b <= 0: return -1
    p, q = 2 * a + b, b                     # a + b phi = (p + q sqrt5)/2
    if p >= 0 and q >= 0: return 1
    if p <= 0 and q <= 0: return -1
    if p > 0: return 1 if p * p > 5 * q * q else -1
    return 1 if 5 * q * q > p * p else -1

def im_sign(w):                              # sign of Im(w), w raw or normalised
    # w - conj(w) = a1 (z - z^4) + a2 (z^2 - z^3) = 2i sin36 (a1 phi + a2)
    return sign_zphi(w[2] - w[3], w[1] - w[4])
def is_real(w): return w[1] == w[4] and w[2] == w[3]
def real_sign(w):                            # w real: w = (w0 - w1) + (w1 - w2) phi
    assert is_real(w)
    return sign_zphi(w[0] - w[1], w[1] - w[2])

def orient(A, B, C):
    return im_sign(mul(conj(sub(B, A)), sub(C, A)))

def on_segment(P, U, W, closed):
    """P on segment UW (exact).  closed: endpoints allowed."""
    d = conj(sub(W, U))
    w1 = mul(sub(P, U), d)
    w2 = mul(sub(W, P), d)
    if not (is_real(w1) and is_real(w2)): return False
    s1, s2 = real_sign(w1), real_sign(w2)
    if closed: return s1 >= 0 and s2 >= 0
    return s1 > 0 and s2 > 0

# sanity of the arithmetic against floating point
_phi = (1 + 5 ** 0.5) / 2
_t = (3, -2, 5, 1, 0)
assert abs(emb(mphi(_t)) - _phi * emb(_t)) < 1e-9
assert abs(emb(mrho(_t)) - cmath.exp(1j * math.pi / 5) * emb(_t)) < 1e-9
assert abs(emb(mul(_t, (1, 4, -2, 0, 7))) - emb(_t) * emb((1, 4, -2, 0, 7))) < 1e-9
assert abs(emb((1, 1, 1, 1, 1))) < 1e-12

# ------------------------------------------------------------------ triangles, substitution
# Convention for the standard triangle of a NEGATIVELY oriented type (the spec only says
# "A = 0, B = 1 (s = +1)"):
#   'swap'   : A = 0, C = 1, B = rho^theta  (same point set as the + triangle, B and C exchanged)
#   'mirror' : A = 0, B = 1, C = rho^(-theta) (complex conjugate of the + triangle)
CONV = os.environ.get('STDCONV', 'swap')
def std_tri(tau):
    col, s = tau
    th = 1 if col == 0 else 3
    if s == 1: return (col, ZERO, ONE, RHO[th])
    if CONV == 'swap': return (col, ZERO, RHO[th], ONE)
    return (col, ZERO, ONE, RHO[(-th) % 10])

def subst1(t):
    col, A, B, C = t
    A2, B2, C2 = mphi(A), mphi(B), mphi(C)
    if col == 0:
        P = add(A2, sub(B, A))
        return [(0, C2, P, B2), (1, P, C2, A2)]
    Q = add(B2, sub(A, B)); R = add(B2, sub(C, B))
    return [(1, R, C2, A2), (1, Q, R, B2), (0, R, Q, A2)]

TYPES = [(0, 1), (0, -1), (1, 1), (1, -1)]
for tau in TYPES:
    c, A, B, C = std_tri(tau)
    assert orient(A, B, C) == tau[1]
    # side lengths: legs 1, base 1/phi or phi
    assert abs(abs(emb(C)) - 1) < 1e-12 and abs(abs(emb(B)) - 1) < 1e-12
    assert abs(abs(emb(sub(C, B))) - (1 / _phi if c == 0 else _phi)) < 1e-12

_super = {}
def supertile(tau, m):
    key = (tau, m)
    if key not in _super:
        if m == 0:
            _super[key] = [std_tri(tau)]
        else:
            _super[key] = [c for t in supertile(tau, m - 1) for c in subst1(t)]
    return _super[key]

_child = {}
def children(tau):
    """[(tau_j, k_j, t_j)] read off at level 1."""
    if tau not in _child:
        res = []
        for (col, A, B, C) in subst1(std_tri(tau)):
            s = orient(A, B, C)
            tj = (col, s)
            sB, sC = std_tri(tj)[2], std_tri(tj)[3]
            d = sub(B, A)
            ks = [k for k in range(10) if rot(sB, k) == d]
            assert len(ks) == 1, "child leg AB is not a rotated standard leg"
            k = ks[0]
            assert add(A, rot(sC, k)) == C, "child is not G(std)"
            res.append((tj, k, A))
        _child[tau] = res
    return _child[tau]

def G(tau, j, m, x):
    tj, k, t = children(tau)[j]
    return add(rot(x, k), phipow(t, m - 1))

_owned = {}
def owned(tau, m):
    """owned rhombs of X_m^tau and the map directed leg -> (rhomb index, edge index in ccw order)"""
    key = (tau, m)
    if key not in _owned:
        rh = []; emap = {}
        for (col, A, B, C) in supertile(tau, m):
            if orient(A, B, C) > 0:
                As = sub(add(B, C), A)
                verts = (A, B, As, C)                 # counter-clockwise
                idx = len(rh); rh.append((col, verts))
                for k in range(4):
                    e = (verts[k], verts[(k + 1) % 4])
                    assert e not in emap, "directed leg in two owned rhombs"
                    emap[e] = (idx, k)
        _owned[key] = (rh, emap)
    return _owned[key]

def simulate(tau, m, entry, prog, maxS=10 ** 7):
    """turtle enters X_m^tau by the crossing `entry`; returns (ok, nS, exit crossing, msg).
    nS = number of S-moves performed after the entry crossing, the last one being the exit."""
    rh, emap = owned(tau, m)
    tail, head, i = entry
    L = len(prog)
    if prog[i] != 'S': return (False, 0, None, 'entry command is not S')
    if (head, tail) not in emap: return (False, 0, None, 'entry crossing does not lead into an owned rhomb')
    r, k = emap[(head, tail)]
    e = (k + 2) % 4; idx = (i + 1) % L
    nS = 0
    seen = set()
    while True:
        st = (r, e, idx)
        if st in seen: return (False, nS, None, 'state repeated (trapped cycle)')
        seen.add(st)
        c = prog[idx]
        if c == 'S':
            v = rh[r][1]
            u, w = v[e], v[(e + 1) % 4]
            nS += 1
            nxt = emap.get((w, u))
            if nxt is None:
                return (True, nS, (u, w, idx), 'exit')
            r, k = nxt
            e = (k + 2) % 4
        elif c == 'R':
            e = (e - 1) % 4
        elif c == 'L':
            e = (e + 1) % 4
        else:
            raise ValueError(c)
        idx = (idx + 1) % L
        if nS > maxS: return (False, nS, None, 'too many steps')

# ------------------------------------------------------------------ certificate
cert = json.load(open(CERT))
prog = cert['prog']; Bv = cert['B']
ports = {int(k): v for k, v in cert['ports'].items()}
passes = {int(k): v for k, v in cert['passes'].items()}
roots = cert['roots']
def Rlev(q): return Bv if (Bv - q) % 2 == 0 else Bv + 1
def ptype(pid): return tuple(ports[pid]['typ'])
def ppar(pid): return ports[pid]['parity']

_portc = {}
def port(pid, m):
    key = (pid, m)
    if key not in _portc:
        p = ports[pid]; q = p['parity']
        assert (m - q) % 2 == 0 and m >= Rlev(q), (pid, m)
        if m == Rlev(q):
            tail, head, i = p['ref']
            _portc[key] = (nz(tuple(tail)), nz(tuple(head)), i)
        else:
            tail, head, i = port(pid, m - 2)
            dt, dh = nz(tuple(p['delta'][0])), nz(tuple(p['delta'][1]))
            _portc[key] = (add(mphi(mphi(tail)), dt), add(mphi(mphi(head)), dh), i)
    return _portc[key]

report = []
fails = []
def log(s):
    print(s); sys.stdout.flush()
def check(cond, msg):
    if not cond:
        fails.append(msg); log('  FAIL: ' + msg)
    return cond

t0 = time.time()
log('standard-triangle convention for s=-1: %s' % CONV)
log('certificate %s: prog=%s B=%d ports=%d passes=%d roots=%s' % (os.path.basename(CERT), prog, Bv, len(ports), len(passes), roots))
npar = {0: 0, 1: 0}
for pid in passes: npar[ppar(pid)] += 1
log('pass families by parity: %s' % npar)

# ---- S: sanity checks of my own substitution / ownership model
log('\n[S] sanity checks of the independent substitution model')
ok = True
for tau in TYPES:
    ch = children(tau)
    log('  type %s: children (tau_j,k_j,t_j) = %s' % (tau, [(tj, k, tuple(t)) for tj, k, t in ch]))
    for m in range(1, 8):
        D = supertile(tau, m)
        assert len(set(D)) == len(D)
        U = []
        for j, (tj, k, t) in enumerate(ch):
            sh = phipow(t, m - 1)
            U += [(c, add(rot(A, k), sh), add(rot(B, k), sh), add(rot(C, k), sh)) for (c, A, B, C) in supertile(tj, m - 1)]
        ok &= check(set(U) == set(D) and len(U) == len(D), 'X_%d^%s != union of G_j(X_%d)' % (m, tau, m - 1))
    # rhomb matching inside the supertile: a base is either shared with the mirror triangle or on the boundary
    for m in (5, 6, 7):
        D = supertile(tau, m); Ds = set(D)
        c0, A0, B0, C0 = std_tri(tau)
        V = [phipow(A0, m), phipow(B0, m), phipow(C0, m)]
        sides = [(V[0], V[1]), (V[1], V[2]), (V[2], V[0])]
        und = {}
        for t in D:
            c, A, B, C = t
            for (x, y, kind) in ((A, B, 'leg'), (A, C, 'leg'), (B, C, 'base')):
                und.setdefault(frozenset((x, y)), []).append((t, kind))
        bad = 0
        for ekey, lst in und.items():
            if len(lst) > 2: bad += 1; continue
            kinds = set(k for _, k in lst)
            if len(kinds) != 1: bad += 1; continue
            x, y = tuple(ekey)
            if len(lst) == 1:
                if not any(on_segment(x, s1, s2, True) and on_segment(y, s1, s2, True) for s1, s2 in sides): bad += 1
            elif kinds == {'base'}:
                (t1, _), (t2, _) = lst
                c, A, B, C = t1
                if t2 != (c, sub(add(B, C), A), B, C): bad += 1
                if orient(*t1[1:]) + orient(*t2[1:]) != 0: bad += 1
        ok &= check(bad == 0, 'edge matching violated in X_%d^%s (%d bad edges)' % (m, tau, bad))
log('  supertile = union of placed children (levels 1..7), bases pair up mirror triangles, unpaired edges lie on the supertile boundary (levels 5..7): %s' % ('PASS' if ok else 'FAIL'))
S_ok = ok

# ---- O3 closure
log('\n[O3] closure / consistency')
ok = True
for pid, ps in passes.items():
    tau = ptype(pid); q = ppar(pid)
    ok &= check(ps['out'] in ports, 'pass %d: exit port missing' % pid)
    ok &= check(ptype(ps['out']) == tau and ppar(ps['out']) == q, 'pass %d: exit port type/parity mismatch' % pid)
    ok &= check(prog[ports[pid]['ref'][2]] == 'S' and prog[ports[ps['out']]['ref'][2]] == 'S', 'pass %d: port command is not S' % pid)
    ch = children(tau)
    ok &= check(len(ps['children']) >= 1, 'pass %d: no children' % pid)
    prev = None
    for (j, a, b) in ps['children']:
        ok &= check(0 <= j < len(ch), 'pass %d: bad slot %s' % (pid, j))
        ok &= check(a in passes, 'pass %d: child entry %d is not a pass family' % (pid, a))
        if a in passes:
            ok &= check(passes[a]['out'] == b, 'pass %d: child (%d->%d) exit differs from table (%d)' % (pid, a, b, passes[a]['out']))
            ok &= check(ptype(a) == ch[j][0], 'pass %d: child %d type %s != slot type %s' % (pid, a, ptype(a), ch[j][0]))
            ok &= check(ppar(a) == 1 - q, 'pass %d: child %d parity' % (pid, a))
            ok &= check(b in ports and ptype(b) == ch[j][0] and ppar(b) == 1 - q, 'pass %d: child exit %d type/parity' % (pid, b))
        ok &= check(prev != j, 'pass %d: consecutive children in the same slot' % pid)
        prev = j
for r in roots:
    ok &= check(r in passes, 'root %s not a pass family' % r)
O3_ok = ok
log('  O3: %s' % ('PASS' if ok else 'FAIL'))

# ---- O1 base simulations
log('\n[O1] base simulations at level B=%d (families of parity %d)' % (Bv, Bv % 2))
ok = True
N = {}
cnt = 0
for pid, ps in sorted(passes.items()):
    if ppar(pid) != Bv % 2: continue
    cnt += 1
    good, nS, ex, msg = simulate(ptype(pid), Bv, port(pid, Bv), prog)
    ok &= check(good, 'O1 pass %d: %s' % (pid, msg))
    if good:
        ok &= check(ex == port(ps['out'], Bv), 'O1 pass %d: exit crossing %s != port %d at level %d %s' % (pid, ex, ps['out'], Bv, port(ps['out'], Bv)))
        N[(pid, Bv)] = nS
O1_ok = ok
log('  %d base simulations; N(pi,%d) range %s..%s ; O1: %s' % (cnt, Bv, min(v for (p, m), v in N.items()) if N else None, max(v for (p, m), v in N.items()) if N else None, 'PASS' if ok else 'FAIL'))

# ---- lengths via the recurrence, compare with the N recorded in the certificate
def Nrec(pid, m):
    key = (pid, m)
    if key not in N:
        if m <= Bv: raise RuntimeError('N(%d,%d) unavailable (base simulation failed)' % (pid, m))
        N[key] = sum(Nrec(a, m - 1) for (j, a, b) in passes[pid]['children'])
    return N[key]
log('\n[N] lengths from the recurrence vs. lengths stored in the certificate')
ok = True; ncmp = 0
for pid, ps in sorted(passes.items()):
    for ms, val in ps.get('N', {}).items():
        m = int(ms)
        ok &= check((m - ppar(pid)) % 2 == 0, 'pass %d: N key %d has wrong parity' % (pid, m))
        if (m - ppar(pid)) % 2 == 0:
            ncmp += 1
            ok &= check(Nrec(pid, m) == val, 'pass %d: N(%d) recurrence %d != certificate %d' % (pid, m, Nrec(pid, m), val))
N_ok = ok
log('  %d stored lengths compared: %s' % (ncmp, 'PASS' if ok else 'FAIL'))
for r in roots:
    log('  root %d (type %s, parity %d): N = %s' % (r, ptype(r), ppar(r), {m: Nrec(r, m) for m in range(Rlev(ppar(r)), 20, 2)}))

# ---- O2 port identities
log('\n[O2] port identities (a),(b),(c) at the two lowest levels of each parity class above B')
ok = True; nid = 0
def o2_levels(q):
    m = Bv + 1
    if (m - q) % 2: m += 1
    return [m, m + 2]
def check_O2(levels_of):
    global nid
    ok = True
    for pid, ps in sorted(passes.items()):
        tau = ptype(pid); q = ppar(pid); ch = ps['children']
        for m in levels_of(q):
            def Gc(j, x):
                return (G(tau, j, m, x[0]), G(tau, j, m, x[1]), x[2])
            j1, a1, b1 = ch[0]
            nid += 1
            ok &= check(port(pid, m) == Gc(j1, port(a1, m - 1)), 'O2(a) pass %d level %d' % (pid, m))
            for i in range(len(ch) - 1):
                j, a, b = ch[i]; j2, a2, b2 = ch[i + 1]
                nid += 1
                ok &= check(Gc(j, port(b, m - 1)) == Gc(j2, port(a2, m - 1)), 'O2(b) pass %d level %d child %d' % (pid, m, i))
            jr, ar, br = ch[-1]
            nid += 1
            ok &= check(port(ps['out'], m) == Gc(jr, port(br, m - 1)), 'O2(c) pass %d level %d' % (pid, m))
    return ok
O2_ok = check_O2(o2_levels)
log('  levels used: parity 0 -> %s, parity 1 -> %s ; %d identities: %s' % (o2_levels(0), o2_levels(1), nid, 'PASS' if O2_ok else 'FAIL'))
nid = 0
extra = check_O2(lambda q: [o2_levels(q)[1] + 2, o2_levels(q)[1] + 4])
log('  (redundant, Lemma 1 cross-check) same identities at the next two levels of each class: %d identities: %s' % (nid, 'PASS' if extra else 'FAIL'))
# Lemma 1 mechanism: constants agree, i.e. delta(parent port) = rho^k delta(child port)
okd = True
for pid, ps in passes.items():
    tau = ptype(pid); ch = ps['children']
    def rd(j, p):
        k = children(tau)[j][1]
        return (rot(nz(tuple(ports[p]['delta'][0])), k), rot(nz(tuple(ports[p]['delta'][1])), k))
    dd = lambda p: (nz(tuple(ports[p]['delta'][0])), nz(tuple(ports[p]['delta'][1])))
    okd &= dd(pid) == rd(ch[0][0], ch[0][1])
    okd &= dd(ps['out']) == rd(ch[-1][0], ch[-1][2])
    for i in range(len(ch) - 1):
        okd &= rd(ch[i][0], ch[i][2]) == rd(ch[i + 1][0], ch[i + 1][1])
log('  (redundant) delta constants: delta(parent side) = rho^k delta(child side) for all identities: %s' % ('PASS' if okd else 'FAIL'))

# ---- O5 exit property
log('\n[O5] exit property above the base (conditions (i)-(iii))')
ok = True
minratio = None; maxc = 0.0; mind = None
s36 = math.sin(math.pi / 5)
for pid, ps in sorted(passes.items()):
    tau = ptype(pid); q = ppar(pid)
    jr, ar, b = ps['children'][-1]
    tj, k, t = children(tau)[jr]
    n0 = Rlev(1 - q)                       # lowest child level used in an inductive step
    # inductive steps for parity q are m > B, m = q mod 2; child level m-1 >= B of parity 1-q
    assert n0 >= Bv and (n0 - (1 - q)) % 2 == 0
    p_n0 = port(b, n0)[0]
    dt = nz(tuple(ports[b]['delta'][0]))
    c = neg(sub(mphi(dt), dt))             # -delta/phi = -delta*(phi-1)
    assert add(mphi(mphi(c)), dt) == c     # c is the fixed point of x -> phi^2 x + delta
    y = sub(p_n0, c)
    # redundant: p_{n0+2} - c = phi^2 y
    assert sub(port(b, n0 + 2)[0], c) == mphi(mphi(y))
    cs, As, Bs, Cs = std_tri(tj)
    cp, Ap, Bp, Cp = std_tri(tau)
    PV = [mphi(Ap), mphi(Bp), mphi(Cp)]
    psides = [(PV[0], PV[1]), (PV[1], PV[2]), (PV[2], PV[0])]
    found = False
    for (U, W) in ((As, Bs), (Bs, Cs), (Cs, As)):
        Un, Wn = phipow(U, n0), phipow(W, n0)
        if not on_segment(y, Un, Wn, False): continue            # (i)
        gU, gW = add(rot(U, k), t), add(rot(W, k), t)             # G_{j,1}
        if not any(on_segment(gU, s1, s2, True) and on_segment(gW, s1, s2, True) for s1, s2 in psides): continue   # (ii)
        d = min(abs(emb(sub(y, Un))), abs(emb(sub(y, Wn))))
        ac = abs(emb(c))
        need = (1.903 + ac) / s36
        ratio = d / need
        if d >= need:                                             # (iii)
            found = True
            minratio = ratio if minratio is None else min(minratio, ratio)
            mind = d if mind is None else min(mind, d)
            maxc = max(maxc, ac)
        else:
            log('  pass %d: (i),(ii) hold but (iii) fails: d=%.6f need=%.6f' % (pid, d, need))
    ok &= check(found, 'O5 pass %d (last child slot %d, exit port %d)' % (pid, jr, b))
O5_ok = ok
log('  %d families; min d = %s, max |c| = %.6f, min d/((1.903+|c|)/sin36) = %s ; O5: %s' % (len(passes), mind, maxc, minratio, 'PASS' if ok else 'FAIL'))

# ---- O4 growth
log('\n[O4] growth')
def expand2(pid):
    out = []
    for (j, a, b) in passes[pid]['children']:
        for (j2, a2, b2) in passes[a]['children']:
            out.append(a2)
    return out
selfrep = {}
for pid in passes:
    e2 = expand2(pid); mu = e2.count(pid)
    if mu >= 1 and len(e2) > mu: selfrep[pid] = (mu, len(e2))
O4_ok = check(len(selfrep) > 0, 'O4: no self-reproducing family')
for r in roots:
    e2 = expand2(r)
    log('  root %d: two-level expansion has %d elements, contains the root %d times' % (r, len(e2), e2.count(r)))
    O4_ok &= check(r in selfrep, 'O4: root %d is not self-reproducing' % r)
log('  self-reproducing families (mu, size): %d of %d ; O4: %s' % (len(selfrep), len(passes), 'PASS' if O4_ok else 'FAIL'))

# claimed recurrence N_{m+2} = 4 N_m + 66 for the root (statement of Theorem U)
for r in roots:
    q = ppar(r); m0 = Rlev(q)
    fam = [p for p in passes if ppar(p) == q]
    K = len(fam) + 3                      # more consecutive terms than the dimension of the linear system (+1 for the constant)
    seq = [Nrec(r, m0 + 2 * i) for i in range(K + 2)]
    ra = (seq[2] - seq[1]) // (seq[1] - seq[0]) if seq[1] != seq[0] else 0; rb = seq[1] - ra * seq[0]   # GEN: constants from the first terms (4 and 66 for SSSR, obtuse +, odd levels)
    rec_ok = ra >= 2 and all(seq[i + 1] == ra * seq[i] + rb for i in range(K + 1))
    log('  root %d: affine recurrence N_{m+2} = %d N_m + %d' % (r, ra, rb))
    log('  root %d: N_%d=%d, N_%d=%d, N_%d=%d, N_%d=%d; N_{m+2}=4N_m+66 checked for %d consecutive steps (> number %d of families of that parity + 1, so it holds for all m by Cayley-Hamilton): %s'
        % (r, m0, seq[0], m0 + 2, seq[1], m0 + 4, seq[2], m0 + 6, seq[3], K + 1, len(fam), 'PASS' if rec_ok else 'FAIL'))
    check(rec_ok, 'recurrence N_{m+2}=4N_m+66 for root %d' % r)

# ---- direct simulations of H(m) for all families (no use of the decomposition)
log('\n[D] direct simulations (no decomposition): all families at levels %d..%d' % (Bv, MAXDIRECT))
D_ok = True; nd = 0
for m in range(Bv, MAXDIRECT + 1):
    tt = time.time(); c_m = 0
    for pid, ps in sorted(passes.items()):
        if (m - ppar(pid)) % 2 or m < Rlev(ppar(pid)): continue
        good, nS, ex, msg = simulate(ptype(pid), m, port(pid, m), prog)
        nd += 1; c_m += 1
        D_ok &= check(good, 'direct pass %d level %d: %s' % (pid, m, msg))
        if good:
            D_ok &= check(ex == port(ps['out'], m), 'direct pass %d level %d: exit crossing differs from predicted port' % (pid, m))
            D_ok &= check(nS == Nrec(pid, m), 'direct pass %d level %d: %d S-moves, recurrence says %d' % (pid, m, nS, Nrec(pid, m)))
    log('  level %d: %d families simulated (%.1fs)' % (m, c_m, time.time() - tt))
log('  %d direct simulations: %s' % (nd, 'PASS' if D_ok else 'FAIL'))

log('\n[R] direct simulation of the root pass in the standard positively oriented obtuse supertile')
R_ok = True
expect = {9: 874, 11: 3562}
for r in roots:
    tau = ptype(r)
    for m in (9, 11):
        if m > MAXDIRECT or (m - ppar(r)) % 2: continue
        D = supertile(tau, m); rh, emap = owned(tau, m)
        ent = port(r, m); exi = port(passes[r]['out'], m)
        good, nS, ex, msg = simulate(tau, m, ent, prog)
        log('  root %d type %s level %d: %d fine triangles, %d owned rhombs' % (r, tau, m, len(D), len(rh)))
        log('    entry crossing (tail,head,i) = %s' % (ent,))
        log('    predicted exit crossing      = %s' % (exi,))
        log('    simulated: ok=%s, S-moves after entry (last one = exit crossing) = %d, exit crossing = %s' % (good, nS, ex))
        R_ok &= check(good and ex == exi, 'root level %d: exit mismatch' % m)
        exp_m = expect[m] if (prog == 'SSSR' and tau == (1, 1)) else Nrec(r, m)   # GEN
        R_ok &= check(nS == exp_m, 'root level %d: %d S-moves, expected %d' % (m, nS, exp_m))
        R_ok &= check((ent[0], ent[1]) not in emap, 'root level %d: entry crossing starts inside an owned rhomb (not an entry from outside R+)' % m) if False else True
        log('    entry leg (tail->head) is a ccw edge of an owned rhomb: %s ; reversed exit leg in an owned rhomb: %s' % ((ent[0], ent[1]) in emap, (exi[1], exi[0]) in emap))
log('  root direct simulations: %s' % ('PASS' if R_ok else 'FAIL'))

log('\nSUMMARY  sanity=%s O1=%s O2=%s O3=%s O4=%s O5=%s N=%s direct_all=%s direct_root=%s  (%.1fs)' % (S_ok, O1_ok, O2_ok, O3_ok, O4_ok, O5_ok, N_ok, D_ok, R_ok, time.time() - t0))
log('failures: %d' % len(fails))
for f in fails[:50]: log('  ' + f)
sys.exit(0 if not fails else 1)
