"""Second-implementation check (the repository's Python code, numpy path only) of the extra bounded programs.
Uses the functions of prove2.py / subst_fast.py / pt.py without modifying them."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)          # this file lives in <repo>/notebook_src
sys.path.insert(0, REPO)
import numpy as np
from collections import Counter
from prove2 import build_tri_atlas, closure_check
from subst_fast import tris_to_arrays, substitute_v, merge_v, nav_from_keys
from pt import turtle_scan

jobs = [("SR" * (k - 3) + "SSR", 7) for k in (35, 40, 45, 50)] + [(p, 10) for p in (
    "SRSSRSRSRSSSR", "SRSRSSRSRSSSR", "SRSRSSRSRSRSSRSRSRSRSSRSRSRSRSSR", "SRSRSRSSRSSRSRSRSRSRSSSR", "SRSRSRSRSSRSSRSRSRSRSSSR",
    "SSRSRSSRSRSRSRSRSSRSRSRSRSRSSR", "SSRSRSRSRSRSSRSRSRSRSRSSRSRSSR", "SRSRSRSSRSRSRSRSRSSSSR", "SRSRSRSRSSRSRSRSRSSSSR",
    "SRSSRSRSRSSRSRSRSRSRSRSSRSRSRSRSSR", "SRSSRSRSRSRSSRSRSRSRSRSRSSRSRSRSSR")]
out_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "py_extra.json")
out = json.load(open(out_file)) if os.path.exists(out_file) else {}
entries, keys = build_tri_atlas(8)
new = closure_check(entries, keys)
assert len(new) == 0 and len(entries) == 30
cache = {}
def patches(n):
    if n not in cache:
        lst = []
        for e in entries:
            col, A, B, C = tris_to_arrays(e["tris"])
            col, A, B, C, anc = substitute_v(col, A, B, C, n)
            quad, fat, ti = merge_v(col, A, B, C)
            eN, eE = nav_from_keys(quad)
            inside = (anc[ti[:, 0]] == e["center"]) | (anc[ti[:, 1]] == e["center"])
            lst.append((eN, eE, np.where(inside)[0]))
        cache.clear(); cache[n] = lst
    return cache[n]
for prog, n in jobs:
    if prog in out: continue
    t0 = time.time(); tot = esc = to = 0; spec = Counter()
    for eN, eE, starts in patches(n):
        sc = turtle_scan(eN, eE, prog, starts, max_cycles=20000)
        fl = sc["flag"]
        tot += 4 * len(starts); esc += int((fl == 1).sum()); to += int((fl == 2).sum())
        spec.update(sc["period"][fl == 0].tolist())
    out[prog] = {"n": n, "starts": tot, "escape": esc, "timeout": to, "spectrum": {str(k): v for k, v in sorted(spec.items())}}
    print(f"({prog})* n={n} starts {tot} esc {esc} to {to} periods {sorted(spec)}  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(out, open(out_file, "w"))
print("finished")
