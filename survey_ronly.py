"""survey of R-only programs (S^a R)(S^b R)(S^c R)... : exit-sequence family"""
import sys, json, itertools, numpy as np
from survey import run, save, load
N=int(sys.argv[1]); maxm=int(sys.argv[2]); maxk=int(sys.argv[3]); mc=int(sys.argv[4]) if len(sys.argv)>4 else 100000
win=load(N)
seen=set(); rows=[]
for m in range(1,maxm+1):
    for ks in itertools.product(range(1,maxk+1),repeat=m):
        # canonical under cyclic rotation of the block sequence
        c=min(ks[i:]+ks[:i] for i in range(m))
        if c in seen: continue
        seen.add(c)
        # skip if periodic repetition of a shorter block sequence
        if any(c==c[:d]*(m//d) for d in range(1,m) if m%d==0): continue
        prog="".join("S"*k+"R" for k in c)
        r=run(N,prog,mc,0.0,quiet=True,win=win); cls=r["classes"]
        maxper=max((x["period"] for x in cls),default=0); maxrad=max((x["radius"] for x in cls),default=0)
        row=dict(prog=prog,blocks=c,N=N,starts=r["nStarts"]*4,closed=r["closed"],escape=r["escape"],timeout=r["timeout"],
                 nclasses=r["nClasses"],maxperiod=maxper,maxradius=maxrad,classes=[(x["period"],x["ncells"],x["norbits"]) for x in cls])
        rows.append(row); save(r)
        print(f"{str(c):14s} {prog:22s} closed={100*row['closed']/row['starts']:5.1f}% to={row['timeout']:4d} classes={row['nclasses']:4d} maxP={maxper:7d} maxR={maxrad:7.1f}",flush=True)
json.dump(rows,open(f"surveyronly_{N}_{maxm}_{maxk}.json","w"))
