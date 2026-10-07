"""figures for the report (matplotlib)."""
import json, os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from pt import *
from survey import load

HERE = os.path.dirname(os.path.abspath(__file__))
OFFS = (0, 0.1, -0.1, 0.2, -0.2)
OUT = os.path.join(HERE, "figs"); os.makedirs(OUT, exist_ok=True)


def window_geometry(N):
    tl = pentagrid_tiles(-3, -3, N, N, OFFS)
    keys, verts = tl["keys"], tl["verts"]
    eN, eE = build_nav(keys)
    return verts, eN, eE, centroids(verts)


def gallery(prog, N, ncols=5, margin=1.5, fname=None, title=None, maxpanels=40):
    res = json.load(open(os.path.join(HERE, f"survey_{N}_{prog}.json")))
    verts, eN, eE, cents = window_geometry(N)
    cls = sorted(res["classes"], key=lambda c: (c["period"], c["ncells"]))[:maxpanels]
    n = len(cls); nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows))
    axes = np.array(axes).reshape(-1)
    for ax in axes: ax.axis("off")
    for ax, c in zip(axes, cls):
        c0, e0 = c["rep"]
        per, fl, order, states = turtle_run(eN, eE, prog, c0, e0)
        cells = np.array(order)
        pts = verts[cells].reshape(-1, 2)
        x0, x1 = pts[:, 0].min() - margin, pts[:, 0].max() + margin
        y0, y1 = pts[:, 1].min() - margin, pts[:, 1].max() + margin
        bg = np.where((cents[:, 0] > x0 - 2) & (cents[:, 0] < x1 + 2) & (cents[:, 1] > y0 - 2) & (cents[:, 1] < y1 + 2))[0]
        ax.add_collection(PolyCollection(verts[bg], facecolors="white", edgecolors="0.8", linewidths=0.4))
        ax.add_collection(PolyCollection(verts[cells], facecolors="#4c8be0", edgecolors="0.15", linewidths=0.5))
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal")
        ax.set_title(f"period {c['period']}, {c['ncells']} tiles", fontsize=9)
    fig.suptitle(title or f"({prog})*  orbit classes (window {N})", fontsize=12)
    fig.tight_layout()
    fname = fname or os.path.join(OUT, f"gallery_{prog}.png")
    fig.savefig(fname, dpi=130); plt.close(fig)
    return fname


def trajectory(prog, N, seed=1, depth=140, max_cycles=20000, fname=None):
    verts, eN, eE, cents = window_geometry(N)
    dep = depth_from_bbox(cents)
    rng = np.random.default_rng(seed)
    c0 = int(rng.choice(np.where(dep >= depth)[0])); e0 = int(rng.integers(4))
    per, fl, order, states = turtle_run(eN, eE, prog, c0, e0, max_cycles=max_cycles)
    cells = np.array(order)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.add_collection(PolyCollection(verts[cells], facecolors="#4c8be0", edgecolors="none"))
    pos = cents[[s[0] for s in states]]
    ax.plot(pos[:, 0], pos[:, 1], "-", color="#d0342c", lw=0.4, alpha=0.7)
    ax.plot(cents[c0, 0], cents[c0, 1], "o", color="k", ms=5)
    ax.set_aspect("equal"); ax.autoscale(); ax.set_title(f"({prog})*: {len(states)} cycles, {'closed' if fl == 0 else 'escaped' if fl == 1 else 'not closed'}, {len(cells)} tiles")
    fig.tight_layout(); fname = fname or os.path.join(OUT, f"traj_{prog}.png"); fig.savefig(fname, dpi=130); plt.close(fig)
    return fname


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "gallery":
        print(gallery(sys.argv[2], int(sys.argv[3]), maxpanels=int(sys.argv[4]) if len(sys.argv) > 4 else 40))
    elif what == "traj":
        print(trajectory(sys.argv[2], int(sys.argv[3]), max_cycles=int(sys.argv[4]) if len(sys.argv) > 4 else 20000))
