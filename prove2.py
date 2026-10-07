"""
prove2.py -- computer-assisted proof at the Robinson-triangle level that every orbit of a
turtle program P is closed on every Penrose rhomb tiling (see the writeup for the argument).

usage: python prove2.py PROG n [N_atlas] [--cuda]
"""
import sys, time, json, os, subprocess
import numpy as np
from collections import defaultdict, Counter
from pt import *
from subst import *
from subst import _canon_items
from subst_fast import *
import survey

HERE = os.path.dirname(os.path.abspath(__file__))
ANG = {0: (1, 2, 2), 1: (3, 1, 1)}  # angle units (36 deg) at A,B,C for red/blue


def tri_vertex_map(tris):
    vt = defaultdict(list)
    for i, (col, A, B, C) in enumerate(tris):
        vt[A].append((i, 0)); vt[B].append((i, 1)); vt[C].append((i, 2))
    return vt


def tri_corona(tris, vt, i):
    col, A, B, C = tris[i]
    idx = set()
    for v in (A, B, C):
        lst = vt[v]
        if sum(ANG[tris[j][0]][role] for j, role in lst) != 10:
            return None
        idx.update(j for j, _ in lst)
    return sorted(idx)


def tri_corona_canon(tris, idxs, center):
    items = [(tris[j][0], tris[j][1], tris[j][2], tris[j][3], int(j == center)) for j in idxs]
    return _canon_items(items, 10)


def build_tri_atlas(N=8, verbose=True):
    t0 = time.time()
    tris = patch_from_seed(N)
    vt = tri_vertex_map(tris)
    atlas = {}
    counts = Counter()
    for i in range(len(tris)):
        cor = tri_corona(tris, vt, i)
        if cor is None:
            continue
        key = tri_corona_canon(tris, cor, i)
        counts[key] += 1
        if key not in atlas:
            atlas[key] = {"tris": [tris[j] for j in cor], "center": cor.index(i)}
    entries = list(atlas.values())
    keys = list(atlas.keys())
    for e, k in zip(entries, keys):
        e["ntris"] = len(e["tris"]); e["count"] = counts[k]
    if verbose:
        print(f"triangle atlas from sigma^{N}(sun): {len(tris)} triangles, {len(entries)} marked coronas (rot36 classes), {time.time()-t0:.0f}s")
    return entries, keys


def closure_check(entries, keys, verbose=True):
    keyset = set(keys)
    new = {}
    for ei, e in enumerate(entries):
        tris = e["tris"]; c = e["center"]
        sub_tris, anc = substitute_patch(tris, 1)
        vt = tri_vertex_map(sub_tris)
        for i in range(len(sub_tris)):
            if anc[i] != c:
                continue
            cor = tri_corona(sub_tris, vt, i)
            if cor is None:
                raise RuntimeError(f"incomplete corona inside sigma(center) (entry {ei})")
            k = tri_corona_canon(sub_tris, cor, i)
            if k not in keyset and k not in new:
                new[k] = (ei, i)
    if verbose:
        print(f"closure check: {len(new)} coronas outside atlas -> " + ("CLOSED" if not new else "NOT CLOSED"))
    return new


def margin(entries):
    best = np.inf
    for e in entries:
        tris = e["tris"]; c = e["center"]
        cnt = Counter()
        for (col, A, B, C) in tris:
            for a, b in ((A, B), (B, C), (C, A)):
                cnt[(min(a, b), max(a, b))] += 1
        bnd = [k for k, v in cnt.items() if v == 1]
        col, A, B, C = tris[c]
        P = np.array([embed(A), embed(B), embed(C)])
        pts = []
        for s in np.linspace(0, 1, 21):
            for t in np.linspace(0, 1 - s, 21):
                pts.append(P[0] + s * (P[1] - P[0]) + t * (P[2] - P[0]))
        pts = np.array(pts)
        for a, b in bnd:
            Aa, Bb = embed(a), embed(b); d = Bb - Aa
            tt = np.clip(((pts - Aa) @ d) / (d @ d), 0, 1)
            proj = Aa + tt[:, None] * d
            best = min(best, np.sqrt(((pts - proj) ** 2).sum(1)).min())
    return best


def verify_program(entries, prog, n, max_cycles=200000, use_cuda=False, verbose=True):
    t0 = time.time()
    total = 0; esc_total = 0; to_total = 0; spectra = Counter(); per_entry = []
    for ei, e in enumerate(entries):
        col, A, B, C = tris_to_arrays(e["tris"])
        col, A, B, C, anc = substitute_v(col, A, B, C, n)
        quad, fat, ti = merge_v(col, A, B, C)
        eN, eE = nav_from_keys(quad)
        inside = (anc[ti[:, 0]] == e["center"]) | (anc[ti[:, 1]] == e["center"])
        starts = np.where(inside)[0]
        cents = (quad @ EVEC).mean(1).astype(np.float32)
        if use_cuda or len(starts) * 4 > 40000:
            inp = os.path.join(HERE, "prove_in.bin"); outp = os.path.join(HERE, "prove_out.bin")
            survey.write_input(inp, eN, eE, prog, starts.astype(np.int32), cents, max_cycles)
            r = subprocess.run([survey.EXE, inp, outp], capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(r.stderr)
            res = json.loads(r.stdout)
            per, fl, orb, mr = survey.read_raw(outp)
            ncl, nesc, nto = res["closed"], res["escape"], res["timeout"]
            spec = Counter(per[fl == 0].tolist())
            os.remove(inp); os.remove(outp)
        else:
            sc = turtle_scan(eN, eE, prog, starts, max_cycles=max_cycles)
            fl = sc["flag"]
            ncl, nesc, nto = int((fl == 0).sum()), int((fl == 1).sum()), int((fl == 2).sum())
            spec = Counter(sc["period"][fl == 0].tolist())
        # certificate sanity: non-empty start set, and every start accounted for exactly once
        if len(starts) == 0:
            raise RuntimeError(f"entry {ei}: empty start set (ancestry/centre mismatch)")
        if ncl + nesc + nto != 4 * len(starts):
            raise RuntimeError(f"entry {ei}: closed+escape+timeout != 4*starts")
        total += 4 * len(starts); esc_total += nesc; to_total += nto; spectra.update(spec)
        per_entry.append({"entry": ei, "ntris_corona": e["ntris"], "rhombs_patch": int(len(quad)), "starts": int(4 * len(starts)),
                          "closed": int(ncl), "escape": int(nesc), "timeout": int(nto)})
        if verbose:
            print(f"  corona {ei:3d} ({e['ntris']:2d} tris): patch {len(quad):7d} rhombs, starts {4*len(starts):7d}, closed {ncl:7d}, escape {nesc:5d}, timeout {nto:4d}  [{time.time()-t0:.0f}s]", flush=True)
    ok = esc_total == 0 and to_total == 0
    if verbose:
        print(f"program ({prog})*, n={n}: total starts {total}, escapes {esc_total}, timeouts {to_total} -> "
              + ("ALL CLOSED: theorem verified" if ok else "FAILED"))
        print("period spectrum:", dict(sorted(spectra.items())))
    return ok, per_entry, dict(sorted(spectra.items()))


if __name__ == "__main__":
    prog = sys.argv[1]; n = int(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else 8
    use_cuda = "--cuda" in sys.argv
    entries, keys = build_tri_atlas(N)
    new = closure_check(entries, keys)
    d0 = margin(entries)
    print("corona sizes (triangles):", sorted(e["ntris"] for e in entries))
    print(f"margin delta0 = {d0:.4f}; at level n={n}: phi^n*delta0 = {PHI**n*d0:.2f}")
    ok, per_entry, spectra = verify_program(entries, prog, n, use_cuda=use_cuda)
    # the certificate requires BOTH the atlas closure (Step 1) and the exhaustive verification (Step 2)
    ok = ok and len(new) == 0 and len(entries) > 0
    print("CERTIFICATE:", "VALID (atlas closed, all starts closed)" if ok else "INVALID")
    json.dump({"prog": prog, "n": n, "N_atlas": N, "n_coronas": len(entries), "closure_new": len(new), "margin0": d0,
               "ok": ok, "per_entry": per_entry, "spectrum": spectra},
              open(os.path.join(HERE, f"proof2_{prog}_n{n}.json"), "w"), indent=1)
