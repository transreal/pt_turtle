"""verify_cert_exact.py -- exact re-verification of the exit property (obligation O5 / condition (V3)) of a
pass-substitution certificate, with no floating-point arithmetic.

usage: python verify_cert_exact.py results/unbounded_certs/cert2_SSSR_1p_1.json [more certificates ...]

verify_cert.py checks O5 with the collinearity tests exact (in Z[zeta]) but the segment parameters and the distance
margins in double precision.  This script repeats O5 with
  * exact arithmetic in Q(sqrt 5) for all squared distances, dot products and segment parameters (a point of
    Z[zeta] has x-coordinate in Q(sqrt 5) and y-coordinate p*r1 + q*r2 with p, q in Q(sqrt 5),
    r1 = sin 72 deg, r2 = sin 144 deg; r1^2, r2^2 and r1*r2 are in Q(sqrt 5)), and
  * rational interval bounds for the square roots in the margin condition
        min(dU, dW) >= (1.903 + |c|) / sin 36 deg + 0.05
    (lower bounds for the distances, upper bounds for the right-hand side), which is therefore proved whenever the
    script reports it.
It prints, for every certificate, whether O5 holds exactly for every pass family and a rigorous lower bound of the
smallest surplus of the margin.
"""
import sys, json, os, re
from fractions import Fraction as F
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "verify_cert.py"), encoding="utf-8").read()
helpers = src.split("# ---------------------------------------------------------------- certificate")[0]
helpers = helpers.split('"""', 2)[2] if helpers.startswith('"""') else helpers      # drop the module docstring
ns = {}
exec(helpers, ns)                                                                    # Z[zeta] arithmetic, std_triangle, FR, ...
canon, add, sub, neg, mulphi, invphi, rot36, is_real_multiple = (ns[k] for k in
    ("canon", "add", "sub", "neg", "mulphi", "invphi", "rot36", "is_real_multiple"))
std_triangle, FR = ns["std_triangle"], ns["FR"]

# ---------------------------------------------------------------- Q(sqrt 5)
class Q5:
    __slots__ = ("a", "b")
    def __init__(self, a, b=0): self.a = F(a); self.b = F(b)          # a + b sqrt5
    def __add__(s, o): o = Q5.of(o); return Q5(s.a + o.a, s.b + o.b)
    def __sub__(s, o): o = Q5.of(o); return Q5(s.a - o.a, s.b - o.b)
    def __mul__(s, o): o = Q5.of(o); return Q5(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)
    def __neg__(s): return Q5(-s.a, -s.b)
    def inv(s):
        n = s.a * s.a - 5 * s.b * s.b
        return Q5(s.a / n, -s.b / n)
    def __truediv__(s, o): return s * Q5.of(o).inv()
    @staticmethod
    def of(x): return x if isinstance(x, Q5) else Q5(x)
    def sign(s):
        if s.b == 0: return (s.a > 0) - (s.a < 0)
        if s.a == 0: return (s.b > 0) - (s.b < 0)
        sa, sb = (s.a > 0) - (s.a < 0), (s.b > 0) - (s.b < 0)
        if sa == sb: return sa
        return sa if s.a * s.a > 5 * s.b * s.b else sb
    def bounds(s):
        """rational interval containing a + b sqrt5"""
        lo5, hi5 = SQRT5_LO, SQRT5_HI
        if s.b >= 0: return (s.a + s.b * lo5, s.a + s.b * hi5)
        return (s.a + s.b * hi5, s.a + s.b * lo5)
    def __repr__(s): return f"({s.a} + {s.b} sqrt5)"

def sqrt_bounds(q, digits=40):
    """rational [lo, hi] with lo <= sqrt(q) <= hi for a rational q >= 0"""
    assert q >= 0
    m = 10 ** digits
    r = isqrt(q.numerator * m * m // q.denominator)
    return (F(r, m), F(r + 2, m))
SQRT5_LO, SQRT5_HI = sqrt_bounds(F(5), 50)

# exact embedding: cos(2 pi j/5) in Q(sqrt5); sin(2 pi j/5) = p r1 + q r2 with r1 = sin 72, r2 = sin 144
COS = [Q5(1), Q5(F(-1, 4), F(1, 4)), Q5(F(-1, 4), F(-1, 4)), Q5(F(-1, 4), F(-1, 4)), Q5(F(-1, 4), F(1, 4))]
SIN = [(Q5(0), Q5(0)), (Q5(1), Q5(0)), (Q5(0), Q5(1)), (Q5(0), Q5(-1)), (Q5(-1), Q5(0))]
R1SQ, R2SQ, R1R2 = Q5(F(5, 8), F(1, 8)), Q5(F(5, 8), F(-1, 8)), Q5(0, F(1, 4))
def embx(K):
    x, p, q = Q5(0), Q5(0), Q5(0)
    for k, c, (sp, sq) in zip(K, COS, SIN):
        x = x + c * k; p = p + sp * k; q = q + sq * k
    return (x, p, q)
def dot(u, v):
    xu, pu, qu = embx(u); xv, pv, qv = embx(v)
    return xu * xv + pu * pv * R1SQ + qu * qv * R2SQ + (pu * qv + qu * pv) * R1R2
def norm2(u): return dot(u, u)
SIN36SQ = Q5(F(5, 8), F(-1, 8))                 # sin^2 36 = (5 - sqrt5)/8

def check_cert(path):
    cert = json.load(open(path))
    prog, B = cert["prog"], cert["B"]
    PORT = {int(k): v for k, v in cert["ports"].items()}
    PASS = {int(k): v for k, v in cert["passes"].items()}
    def Rq(q): return B if B % 2 == q else B + 1
    def port_at(pid, m):
        f = PORT[pid]; q = f["parity"]; mm = Rq(q)
        tail, head, i = tuple(f["ref"][0]), tuple(f["ref"][1]), f["ref"][2]
        dt, dh = tuple(f["delta"][0]), tuple(f["delta"][1])
        while mm < m:
            tail = add(mulphi(mulphi(tail)), dt); head = add(mulphi(mulphi(head)), dh); mm += 2
        return (tail, head, i)
    def ptype(pid): return tuple(PORT[pid]["typ"])
    sin36_lo = sqrt_bounds(SIN36SQ.bounds()[0])[0]
    worst = None; bad = []
    for pin, p in PASS.items():
        typ = ptype(pin); q = PORT[pin]["parity"]
        slot, a, b = p["children"][-1]
        ctyp, k, t, _ = FR[typ][slot]
        n0 = (B + 1 if (B + 1) % 2 == q else B + 2) - 1
        f = PORT[b]
        tail0 = port_at(b, n0)[0]
        c = neg(invphi(tuple(f["delta"][0])))
        y = sub(tail0, c)
        if port_at(b, n0 + 2)[0] != add(c, mulphi(mulphi(y))): bad.append((pin, "affine form")); continue
        std = std_triangle(*ctyp); V = [std[1], std[2], std[3]]
        def scale(K, n):
            for _ in range(n): K = mulphi(K)
            return K
        parent = std_triangle(*typ); PV = [mulphi(parent[1]), mulphi(parent[2]), mulphi(parent[3])]
        def place(K):
            v = K
            for _ in range(k): v = rot36(v)
            return add(v, t)
        cc2 = norm2(c)
        cc_hi = sqrt_bounds(cc2.bounds()[1])[1]
        need_hi = (F(1903, 1000) + cc_hi) / sin36_lo + F(1, 20)        # upper bound of (1.903 + |c|)/sin36 + 0.05
        ok = False
        for iu in range(3):
            U, W = V[iu], V[(iu + 1) % 3]
            Un, Wn = scale(U, n0), scale(W, n0)
            if not is_real_multiple(sub(y, Un), sub(Wn, Un)): continue
            s = dot(sub(y, Un), sub(Wn, Un)) / norm2(sub(Wn, Un))
            if not (s.sign() > 0 and (Q5(1) - s).sign() > 0): continue
            # child edge (U, W) placed in the parent lies on a side of the parent, exactly
            gU, gW = place(U), place(W)
            onside = False
            for ip in range(3):
                P1, P2 = PV[ip], PV[(ip + 1) % 3]
                if is_real_multiple(sub(gU, P1), sub(P2, P1)) and is_real_multiple(sub(gW, P1), sub(P2, P1)):
                    l2 = norm2(sub(P2, P1))
                    s1 = dot(sub(gU, P1), sub(P2, P1)) / l2; s2 = dot(sub(gW, P1), sub(P2, P1)) / l2
                    if s1.sign() >= 0 and (Q5(1) - s1).sign() >= 0 and s2.sign() >= 0 and (Q5(1) - s2).sign() >= 0:
                        onside = True
            if not onside: continue
            dU_lo = sqrt_bounds(max(F(0), norm2(sub(y, Un)).bounds()[0]))[0]
            dW_lo = sqrt_bounds(max(F(0), norm2(sub(y, Wn)).bounds()[0]))[0]
            surplus = min(dU_lo, dW_lo) - need_hi
            if surplus >= 0:
                ok = True
                worst = surplus if worst is None else min(worst, surplus)
        if not ok: bad.append((pin, "exit port of last child not safely on a boundary edge (exact)"))
    return prog, B, len(PASS), bad, worst

if __name__ == "__main__":
    allok = True
    for path in sys.argv[1:]:
        prog, B, npass, bad, worst = check_cert(path)
        ok = not bad
        allok &= ok
        print(f"{os.path.basename(path)}: ({prog})* B={B} passes={npass} O5 exact: {'OK' if ok else 'FAIL'}"
              + (f"; rigorous lower bound of the smallest margin surplus = {float(worst):.6f}" if worst is not None else "")
              + (f"; failures: {bad}" if bad else ""))
    print("ALL CERTIFICATES: O5 EXACT " + ("OK" if allok else "FAIL"))
    sys.exit(0 if allok else 1)
