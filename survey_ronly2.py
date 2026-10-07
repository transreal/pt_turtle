"""2-block R-only family (a,b) with a,b<=10 and 4-block families with small k, on the 300 window."""
import sys, json, itertools, numpy as np
from survey import run, save, load
N=300; mc=100000; win=load(N); rows=[]; seen=set()
fam=[(a,b) for a in range(1,11) for b in range(a,11)]
fam+= [c for c in itertools.product(range(1,4),repeat=4)]
for ks in fam:
    m=len(ks); c=min(ks[i:]+ks[:i] for i in range(m))
    if c in seen or any(c==c[:d]*(m//d) for d in range(1,m) if m%d==0): continue
    seen.add(c); prog="".join("S"*k+"R" for k in c)
    r=run(N,prog,mc,0.0,quiet=True,win=win); cls=r["classes"]
    maxper=max((x["period"] for x in cls),default=0); maxrad=max((x["radius"] for x in cls),default=0)
    row=dict(prog=prog,blocks=c,N=N,starts=r["nStarts"]*4,closed=r["closed"],escape=r["escape"],timeout=r["timeout"],
             nclasses=r["nClasses"],maxperiod=maxper,maxradius=maxrad,classes=[(x["period"],x["ncells"],x["norbits"]) for x in cls])
    rows.append(row); save(r)
    print(f"{str(c):16s} {prog:30s} closed={100*row['closed']/row['starts']:5.1f}% to={row['timeout']:5d} classes={row['nclasses']:4d} maxP={maxper:7d} maxR={maxrad:7.1f}",flush=True)
json.dump(rows,open("surveyronly2_300.json","w"))
