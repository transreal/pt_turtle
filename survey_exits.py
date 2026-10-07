"""survey R-only programs parametrised by cyclic net-exit sequences over {2,3,...}:
   block (k | 2^r)  ->  (SR)^(k-3) S^(r+2) R  ;  k>=3."""
import sys, json, itertools, numpy as np
from survey import run, save, load
def prog_from_exits(seq):
    m=len(seq); i0=next(i for i,k in enumerate(seq) if k>=3); s=seq[i0:]+seq[:i0]
    blocks=[]; 
    for k in s:
        if k>=3: blocks.append([k,0])
        else: blocks[-1][1]+=1
    return "".join("SR"*(k-3)+"S"*(r+2)+"R" for k,r in blocks)
def canon_cyc(seq):
    m=len(seq); return min(seq[i:]+seq[:i] for i in range(m))
if __name__=="__main__":
    N=300; mc=100000; win=load(N); rows=[]; seen=set()
    fams=[]
    for m in range(1,5): fams+= [tuple(s) for s in itertools.product(range(2,9),repeat=m)]
    for m in (5,6): fams+= [tuple(s) for s in itertools.product(range(2,5),repeat=m)]
    for seq in fams:
        if max(seq)<3: continue
        c=canon_cyc(seq); m=len(c)
        if c in seen or any(c==c[:d]*(m//d) for d in range(1,m) if m%d==0): continue
        seen.add(c); prog=prog_from_exits(list(c))
        r=run(N,prog,mc,0.0,quiet=True,win=win); cls=r["classes"]
        maxper=max((x["period"] for x in cls),default=0); maxrad=max((x["radius"] for x in cls),default=0)
        row=dict(exits=c,prog=prog,N=N,starts=r["nStarts"]*4,closed=r["closed"],escape=r["escape"],timeout=r["timeout"],
                 nclasses=r["nClasses"],maxperiod=maxper,maxradius=maxrad,classes=[(x["period"],x["ncells"],x["norbits"]) for x in cls][:40])
        rows.append(row); save(r)
        print(f"{str(c):18s} {prog:28s} closed={100*row['closed']/row['starts']:5.1f}% to={row['timeout']:6d} classes={row['nclasses']:4d} maxP={maxper:7d} maxR={maxrad:7.1f}",flush=True)
        json.dump(rows,open("surveyexits_300.json","w"))
