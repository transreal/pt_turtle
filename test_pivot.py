"""Verify the pivot-walk (vertex-ray walk) reformulation against the turtle.
State (c,e) <-> directed edge d = (v_e -> v_{e+1}) of tile c (c on the left). Pivot p = tail(d).
Rays at a vertex ordered clockwise. Claims:
  S: exit along the ray 1 CW from rho (rho = ray of d), new rho = 1 CW from the back ray
  L: exit along rho itself, new rho = 1 CW from the back ray
  R: exit along the ray 1 CCW from rho, new rho = the back ray
"""
import numpy as np
from pt import *
from collections import defaultdict
OFFS=(0,0.1,-0.1,0.2,-0.2)
tl=pentagrid_tiles(-3,-3,40,40,OFFS); keys,verts=tl["keys"],tl["verts"]; eN,eE=build_nav(keys)
T=len(keys)
kt=[tuple(map(tuple,keys[t])) for t in range(T)]
# rays at each vertex: neighbours via tile edges, sorted clockwise by angle
nbrs=defaultdict(set); vpos={}
for t in range(T):
    for e in range(4):
        a,b=kt[t][e],kt[t][(e+1)%4]; nbrs[a].add(b); nbrs[b].add(a); vpos[a]=verts[t][e]; vpos[b]=verts[t][(e+1)%4]
rays={}
for v,ns in nbrs.items():
    lst=sorted(ns,key=lambda w:-np.arctan2(*(vpos[w]-vpos[v])[::-1]))  # clockwise = decreasing angle
    rays[v]=lst
def cw(v,w,k):
    """ray at v that is k sectors clockwise from ray v->w (k may be negative)"""
    lst=rays[v]; i=lst.index(w); return lst[(i+k)%len(lst)]
def state_to_pivot(c,e): return kt[c][e], kt[c][(e+1)%4]
def pivot_to_state(p,q):
    # tile on the left of p->q : find tile with consecutive vertices p,q
    for t in range(T):
        pass
# build lookup for directed edges
dmap={}
for t in range(T):
    for e in range(4): dmap[(kt[t][e],kt[t][(e+1)%4])]=(t,e)
rng=np.random.default_rng(0)
interior=(eN>=0).all(1); cand=np.where(interior)[0]
bad=0; n=0
for trial in range(4000):
    c=int(rng.choice(cand)); e=int(rng.integers(4)); cmd=rng.choice(["S","L","R"])
    p,q=state_to_pivot(c,e)
    if len(rays[p])==0: continue
    # need full vertex stars at p and at the new pivot: skip near boundary via checking all rays' tiles exist
    if cmd=="S":
        nb=eN[c,e]
        if nb<0: continue
        c2,e2=int(nb),(eE[c,e]+2)&3
        p2,q2=state_to_pivot(c2,e2)
        exp_p2=cw(p,q,1)                       # exit 1 CW from rho
        exp_q2=cw(p2,p,1) if p2 in rays and p in rays[p2] else None   # new rho = 1 CW from back ray
    elif cmd=="L":
        c2,e2=c,(e+1)&3; p2,q2=state_to_pivot(c2,e2)
        exp_p2=q; exp_q2=cw(p2,p,1)
    else:
        c2,e2=c,(e+3)&3; p2,q2=state_to_pivot(c2,e2)
        exp_p2=cw(p,q,-1); exp_q2=p
    n+=1
    if (p2,q2)!=(exp_p2,exp_q2): bad+=1
print("checked",n,"moves; mismatches:",bad)
# exit-sequence claim: consecutive letters (x,y) -> exit offset e(x,y) CW from back ray
tab={("S","S"):2,("L","S"):2,("R","S"):1,("S","L"):1,("L","L"):1,("R","L"):0,("S","R"):0,("L","R"):0,("R","R"):-1}
bad=0;n=0
for trial in range(4000):
    c=int(rng.choice(cand)); e=int(rng.integers(4)); x=rng.choice(["S","L","R"]); y=rng.choice(["S","L","R"])
    # apply x then y, check pivot moves
    def step(c,e,cmd):
        if cmd=="S":
            nb=eN[c,e]; return (int(nb),(eE[c,e]+2)&3) if nb>=0 else None
        return (c,(e+1)&3) if cmd=="L" else (c,(e+3)&3)
    s1=step(c,e,x); 
    if s1 is None: continue
    s2=step(*s1,y)
    if s2 is None: continue
    p0,_=state_to_pivot(c,e); p1,_=state_to_pivot(*s1); p2,_=state_to_pivot(*s2)
    if p1==p0: continue   # need an actual arrival at p1 to define the back ray (never happens: every move moves the pivot)
    n+=1
    if p2!=cw(p1,p0,tab[(x,y)]): bad+=1
print("exit-offset table checked",n,"pairs; mismatches:",bad)
