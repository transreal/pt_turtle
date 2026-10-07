import sys, json, os, numpy as np
from survey import run, save, load
from progs import programs
N=int(sys.argv[1]); maxlen=int(sys.argv[2]); mc=int(sys.argv[3]) if len(sys.argv)>3 else 50000
win=load(N)
rows=[]
for p in programs(maxlen):
    r=run(N,p,mc,0.0,quiet=True,win=win)
    cls=r["classes"]
    maxper=max((c["period"] for c in cls),default=0); maxrad=max((c["radius"] for c in cls),default=0)
    row=dict(prog=p,N=N,starts=r["nStarts"]*4,closed=r["closed"],escape=r["escape"],timeout=r["timeout"],maxTOr=r["maxTimeoutRadius"],
             nclasses=r["nClasses"],norbits=r["nOrbits"],maxperiod=maxper,maxradius=maxrad,
             classes=[(c["period"],c["ncells"],c["norbits"]) for c in cls])
    rows.append(row); save(r)
    print(f"{p:10s} esc={row['escape']:7d}/{row['starts']:8d} to={row['timeout']:5d} classes={row['nclasses']:4d} maxP={maxper:7d} maxR={maxrad:7.1f}  {row['classes'][:8]}{'...' if len(cls)>8 else ''}",flush=True)
json.dump(rows,open(f"surveyall_{N}_{maxlen}.json","w"))
