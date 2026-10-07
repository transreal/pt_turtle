"""
extract2.py -- constructive extraction of a closed pass-substitution system.

A port family of type tau and parity q is a sequence of crossings chi_m (m = R(q), R(q)+2, ...) with
    chi_{m+2} = phi^2 * chi_m + delta        (tail and head separately; exact in Z[zeta]).
A pass family = (tau, q, entry port family); its exit port family and its decomposition into child
pass families (parity 1-q, one level down) are read off from simulations at the two levels
R(q)+2 and R(q)+4 and must agree.

usage: python extract2.py PROG B [root_type_col root_sign]
output: cert2_PROG.json
"""
import sys, json, time
import numpy as np
from lsys import *
from subst import canon, mulphi, add, sub, rot36, embed
from pt import parse_prog

prog_str = sys.argv[1]; B = int(sys.argv[2])
import os
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else "."
TOPK = int(sys.argv[4]) if len(sys.argv) > 4 else 1
NCAND = int(sys.argv[5]) if len(sys.argv) > 5 else 4
PROG = parse_prog(prog_str).tolist()
S_IDX = [i for i, c in enumerate(PROG) if c == 0]
def R(q): return B if B % 2 == q else B + 1

CH = {typ: [frame_of(c) for c in children(std_triangle(*typ))] for typ in TYPES}
_st = {}
def ST(typ, m):
    if (typ, m) not in _st:
        s = Supertile(typ, m); s.bmap = {}
        br, be = np.where(s.eN < 0)
        for r, e in zip(br.tolist(), be.tolist()):
            s.bmap[(s.key(r, e + 1), s.key(r, e))] = (r, e)
        _st[(typ, m)] = s
    return _st[(typ, m)]

def phi2(K): return mulphi(mulphi(K))
def invphi(K): return sub(mulphi(K), K)
def invphi2(K): return invphi(invphi(K))
def up(port, d):   # port at m -> port at m+2
    return (add(phi2(port[0]), d[0]), add(phi2(port[1]), d[1]), port[2])
def down(port, d):
    return (invphi2(sub(port[0], d[0])), invphi2(sub(port[1], d[1])), port[2])
def delta(hi, lo):
    assert hi[2] == lo[2]
    return (sub(hi[0], phi2(lo[0])), sub(hi[1], phi2(lo[1])))

def run_detail(st, entry):
    tail, head, i = entry
    if (tail, head) not in st.bmap: return None
    r, eb = st.bmap[(tail, head)]
    e = (eb + 2) & 3; ci = i + 1; n = len(PROG)
    eN, eE, rslot = st.eN, st.eE, st.rslot
    cur_slot = int(rslot[r]); cur_entry = entry; cnt = 0; total = 0; cps = []
    s0 = (r, e, ci % n)
    while True:
        c = PROG[ci % n]
        if c == 0:
            cross = (st.key(r, e), st.key(r, e + 1), ci % n)
            nb = eN[r, e]
            if nb < 0:
                cps.append((cur_slot, cur_entry, cross, cnt + 1)); return cross, total + 1, cps
            if rslot[nb] != cur_slot:
                cps.append((cur_slot, cur_entry, cross, cnt + 1)); cur_slot = int(rslot[nb]); cur_entry = cross; cnt = 0
            else:
                cnt += 1
            e = (eE[r, e] + 2) & 3; r = int(nb); total += 1
        elif c == 1: e = (e + 3) & 3
        else: e = (e + 1) & 3
        ci += 1
        if (r, e, ci % n) == s0: return None

def to_child(typ, m, slot, cross):
    tj, k, t = CH[typ][slot]
    return tj, (inv_frame(k, t, m - 1, cross[0]), inv_frame(k, t, m - 1, cross[1]), cross[2])

PORTS = {}     # port family id -> dict(typ, parity, ref (port at R(parity)), delta)
PID = {}       # (typ, parity, ref port) -> id
def port_family(typ, q, port_lo, m_lo, port_hi):
    """port_lo at level m_lo, port_hi at level m_lo+2 (same parity q)."""
    d = delta(port_hi, port_lo)
    ref = port_lo; m = m_lo
    while m > R(q):
        ref = down(ref, d); m -= 2
    while m < R(q):
        ref = up(ref, d); m += 2
    key = (typ, q, ref, d)      # two families may share a port at one level but differ in delta
    if key in PID:
        return PID[key]
    pid = len(PORTS); PID[key] = pid
    PORTS[pid] = {"typ": typ, "parity": q, "ref": ref, "delta": d}
    return pid
def port_at(pid, m):
    f = PORTS[pid]; p = f["ref"]; mm = R(f["parity"])
    assert (m - mm) % 2 == 0 and m >= mm
    while mm < m:
        p = up(p, f["delta"]); mm += 2
    return p

PASS = {}      # entry port family id -> dict(out, children [(slot, child entry pid)], N {m: n})
queue = []
def process(pid_in):
    f = PORTS[pid_in]; typ, q = f["typ"], f["parity"]
    L1, L2 = R(q) + 2, R(q) + 4
    e1, e2 = port_at(pid_in, L1), port_at(pid_in, L2)
    r1 = run_detail(ST(typ, L1), e1); r2 = run_detail(ST(typ, L2), e2)
    if r1 is None or r2 is None:
        return f"entry not on boundary or trapped (levels {L1},{L2})"
    x1, N1, c1 = r1; x2, N2, c2 = r2
    if [c[0] for c in c1] != [c[0] for c in c2]:
        return f"decompositions differ between levels {L1} and {L2}: {len(c1)} vs {len(c2)} child passes"
    if x1[2] != x2[2]: return "exit command index differs"
    pid_out = port_family(typ, q, x1, L1, x2)
    kids = []
    for (s1, a1, b1, n1), (s2, a2, b2, n2) in zip(c1, c2):
        tj, ca1 = to_child(typ, L1, s1, a1); _, ca2 = to_child(typ, L2, s2, a2)
        _, cb1 = to_child(typ, L1, s1, b1); _, cb2 = to_child(typ, L2, s2, b2)
        if ca1[2] != ca2[2] or cb1[2] != cb2[2]: return "child command index differs"
        cin = port_family(tj, 1 - q, ca1, L1 - 1, ca2)
        cout = port_family(tj, 1 - q, cb1, L1 - 1, cb2)
        kids.append((s1, cin, cout))
        if cin not in PASS and cin not in queue: queue.append(cin)
    PASS[pid_in] = {"out": pid_out, "children": kids, "N": {L1: N1, L2: N2}}
    return None

def find_roots(typ, q, topk=4, ncand=4):
    """top-k longest passes at level R(q)+2; partner at R(q)+4 chosen among the nearest scaled entries
    such that the decomposition at both levels agrees."""
    L1, L2 = R(q) + 2, R(q) + 4
    best = {}
    for L in (L1, L2):
        st = ST(typ, L); cand = []
        for (tail, head), (r, eb) in st.bmap.items():
            for i in S_IDX:
                ex, steps, _ = st.run_pass(PROG, r, (eb + 2) & 3, i + 1)
                if ex is not None: cand.append((steps, (tail, head, i)))
        best[L] = cand
    out = []
    for n1, p1 in sorted(best[L1], reverse=True)[:topk]:
        tgt = embed(p1[0]) * PHI ** 2
        c2 = [(n, p) for n, p in best[L2] if p[2] == p1[2] and sub(p[1], p[0]) == sub(p1[1], p1[0])]
        c2.sort(key=lambda t: np.linalg.norm(embed(t[1][0]) - tgt))
        for n2, p2 in c2[:ncand]:
            key = (typ, q, p1)
            d = delta(p2, p1)
            # tentative family; test that decompositions agree before registering
            r1 = run_detail(ST(typ, L1), p1); r2 = run_detail(ST(typ, L2), p2)
            if r1 is None or r2 is None: continue
            if [c[0] for c in r1[2]] != [c[0] for c in r2[2]] or r1[0][2] != r2[0][2]: continue
            out.append((port_family(typ, q, p1, L1, p2), (n1, n2))); break
    return out


def reachable(root):
    seen = set(); stack = [root]
    while stack:
        p = stack.pop()
        if p in seen: continue
        if p not in PASS: return None
        seen.add(p)
        for s, a_, b_ in PASS[p]["children"]: stack.append(a_)
    return seen


def write_cert(root, good, fname):
    used = set()
    for p in good:
        used.add(p); used.add(PASS[p]["out"])
        for s, a_, b_ in PASS[p]["children"]: used.add(a_); used.add(b_)
    cert = {"prog": prog_str, "B": B,
            "ports": {str(pid): {"typ": list(PORTS[pid]["typ"]), "parity": PORTS[pid]["parity"],
                                 "ref": [list(PORTS[pid]["ref"][0]), list(PORTS[pid]["ref"][1]), PORTS[pid]["ref"][2]],
                                 "delta": [list(PORTS[pid]["delta"][0]), list(PORTS[pid]["delta"][1])]} for pid in sorted(used)},
            "passes": {str(p): {"out": PASS[p]["out"], "children": [list(c) for c in PASS[p]["children"]],
                                "N": {str(k): v for k, v in PASS[p]["N"].items()}} for p in sorted(good)},
            "roots": [root]}
    json.dump(cert, open(fname, "w"))


if __name__ == "__main__":
    t0 = time.time()
    roots = []
    for typ in TYPES:
        for q in (0, 1):
            for pid, ns in find_roots(typ, q, topk=TOPK, ncand=NCAND):
                roots.append((pid, typ, q, ns))
                if pid not in PASS and pid not in queue: queue.append(pid)
    failures = {}
    while queue:
        pid = queue.pop(0)
        if pid in PASS: continue
        if len(PASS) > 20000: print("too many pass families -- giving up"); break
        err = process(pid)
        if err: failures[pid] = err
    print(f"{prog_str}: pass families {len(PASS)}  port families {len(PORTS)}  failures {len(failures)}  [{time.time()-t0:.0f}s]")
    written = []
    for pid, typ, q, ns in roots:
        good = reachable(pid)
        if good is None or any(p in failures for p in good):
            print(f"  root type {typ} parity {q} lengths {ns}: NOT closed"); continue
        fname = os.path.join(OUTDIR, f"cert2_{prog_str}_{typ[0]}{'p' if typ[1] > 0 else 'm'}_{q}_{pid}.json")
        write_cert(pid, good, fname); written.append(fname)
        print(f"  root type {typ} parity {q} lengths {ns}: closed system of {len(good)} pass families -> {fname}")
