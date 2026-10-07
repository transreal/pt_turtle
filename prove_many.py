"""prove_many.py -- run the proof pipeline (atlas built once) for many programs.
usage: python prove_many.py PROG:n PROG:n ...   or   python prove_many.py --auto  (all survey_300 bounded candidates)
"""
import sys, os, json, glob, time
from prove2 import build_tri_atlas, closure_check, margin, verify_program, HERE
from pt import PHI

def auto_list():
    out = []
    for f in glob.glob(os.path.join(HERE, "survey_300_*.json")):
        r = json.load(open(f)); cls = r["classes"]
        if not cls: continue
        mr = max(x["radius"] for x in cls); c = 100 * r["closed"] / (4 * r["nStarts"])
        if mr < 25 and c > 95 and len(r["prog"]) <= 60:
            out.append((r["prog"], mr))
    return sorted(set(out), key=lambda t: (len(t[0]), t[0]))

def level_for(mr, d0):
    need = 2 * mr + 3.0
    n = 7
    while PHI ** n * d0 < need: n += 1
    return n

if __name__ == "__main__":
    entries, keys = build_tri_atlas(8)
    new = closure_check(entries, keys); assert not new
    d0 = margin(entries); print("delta0", d0)
    if sys.argv[1] == "--auto":
        jobs = [(p, level_for(mr, d0)) for p, mr in auto_list()]
    else:
        jobs = [(a.split(":")[0], int(a.split(":")[1])) for a in sys.argv[1:]]
    done = {}
    if os.path.exists(os.path.join(HERE, "proof_many.json")):
        done = json.load(open(os.path.join(HERE, "proof_many.json")))
    for prog, n in jobs:
        if prog in done and done[prog]["ok"]: continue
        t0 = time.time()
        ok, per_entry, spectra = verify_program(entries, prog, n, use_cuda=(n >= 9), verbose=False)
        starts = sum(p["starts"] for p in per_entry); esc = sum(p["escape"] for p in per_entry); to = sum(p["timeout"] for p in per_entry)
        done[prog] = {"n": n, "ok": ok, "starts": starts, "escape": esc, "timeout": to, "spectrum": spectra, "time": time.time() - t0}
        print(f"{prog:44s} n={n} starts={starts:8d} esc={esc} to={to} -> {'VERIFIED' if ok else 'FAILED'} periods={sorted(int(k) for k in spectra)}  [{time.time()-t0:.0f}s]", flush=True)
        json.dump(done, open(os.path.join(HERE, "proof_many.json"), "w"), indent=1)
