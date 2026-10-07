"""
survey.py -- GPU survey driver.
  python survey.py gen N            -> builds pentagrid window N x N, caches nav arrays in win_N.npz
  python survey.py run N PROG [maxCycles] [depth]  -> runs turtle_survey.exe, prints JSON summary,
                                       saves survey_N_PROG.json
"""
import sys, os, json, struct, subprocess, time
import numpy as np
from pt import *

HERE = os.path.dirname(os.path.abspath(__file__))
OFFS = (0, 0.1, -0.1, 0.2, -0.2)
EXE = os.path.join(HERE, "turtle_survey.exe")


def gen(N, offs=OFFS, tag=None):
    t0 = time.time()
    tl = pentagrid_tiles(-3, -3, N, N, offs)
    keys, verts = tl["keys"], tl["verts"]
    eN, eE = build_nav(keys)
    cents = centroids(verts).astype(np.float32)
    dep = depth_from_bbox(cents)
    path = os.path.join(HERE, f"win_{tag or N}.npz")
    np.savez_compressed(path, eN=eN, eE=eE, cents=cents, dep=dep, fat=tl["fat"], N=N, offs=np.array(offs))
    print(f"window {N}: {len(keys)} tiles, interior {(eN>=0).all(1).sum()}, {time.time()-t0:.1f}s -> {path}")
    return path


def load(N):
    return np.load(os.path.join(HERE, f"win_{N}.npz"))


def write_input(path, eN, eE, prog, starts, cents, maxCycles):
    with open(path, "wb") as f:
        f.write(struct.pack("<iiii", eN.shape[0], len(starts), len(prog), maxCycles))
        f.write(np.ascontiguousarray(eN, dtype=np.int32).tobytes())
        f.write(np.ascontiguousarray(eE, dtype=np.int32).tobytes())
        f.write(parse_prog(prog).astype(np.int32).tobytes())
        f.write(np.ascontiguousarray(starts, dtype=np.int32).tobytes())
        f.write(np.ascontiguousarray(cents, dtype=np.float32).tobytes())


def read_raw(path):
    with open(path, "rb") as f:
        n = struct.unpack("<i", f.read(4))[0]
        per = np.frombuffer(f.read(4 * n), dtype=np.int32)
        fl = np.frombuffer(f.read(4 * n), dtype=np.int32)
        orb = np.frombuffer(f.read(4 * n), dtype=np.int32)
        mr = np.frombuffer(f.read(4 * n), dtype=np.float32)
    return per, fl, orb, mr


def run(N, prog, maxCycles=200000, depth=0.0, quiet=False, win=None, keep_raw=False):
    w = win if win is not None else load(N)
    eN, eE, cents, dep = w["eN"], w["eE"], w["cents"], w["dep"]
    interior = (eN >= 0).all(1)
    starts = np.where(interior & (dep >= depth))[0].astype(np.int32)
    tag = f"{N}_{prog}"
    inp = os.path.join(HERE, f"in_{tag}.bin"); outp = os.path.join(HERE, f"raw_{tag}.bin")
    write_input(inp, eN, eE, prog, starts, cents, maxCycles)
    t0 = time.time()
    r = subprocess.run([EXE, inp, outp], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr); raise RuntimeError("cuda failed")
    res = json.loads(r.stdout)
    res["prog"] = prog; res["N"] = int(N); res["depth"] = depth; res["time"] = time.time() - t0
    res["starts"] = starts  # keep for raw analysis
    per, fl, orb, mr = read_raw(outp)
    res["raw"] = {"period": per, "flag": fl, "orbit": orb, "maxr": mr}
    os.remove(inp)
    if not keep_raw:
        os.remove(outp)
    if not quiet:
        summarize(res)
    return res


def summarize(res, maxlines=60):
    print(f"[{res['prog']}] N={res['N']} tiles={res['nTiles']} starts={res['nStarts']} closed={res['closed']} "
          f"escape={res['escape']} timeout={res['timeout']} (maxTOradius {res['maxTimeoutRadius']}) "
          f"orbits={res['nOrbits']} classes={res['nClasses']} t={res['time']:.1f}s")
    cls = res["classes"]
    for c in cls[:maxlines]:
        print(f"   period={c['period']:7d} ncells={c['ncells']:7d} norbits={c['norbits']:6d} nstarts={c['nstarts']:8d} radius={c['radius']:8.2f} rep={c['rep']}")
    if len(cls) > maxlines:
        print(f"   ... {len(cls)-maxlines} more classes")


def save(res, name=None):
    d = {k: v for k, v in res.items() if k not in ("raw", "starts")}
    name = name or os.path.join(HERE, f"survey_{res['N']}_{res['prog']}.json")
    with open(name, "w") as f:
        json.dump(d, f)
    return name


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "gen":
        gen(int(sys.argv[2]))
    elif cmd == "run":
        N = int(sys.argv[2]); prog = sys.argv[3]
        mc = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
        dp = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
        res = run(N, prog, mc, dp)
        save(res)
