import numpy as np, sys
from pt import *
from survey import load
w=load(300); eN,eE,cents,dep=w["eN"],w["eE"],w["cents"],w["dep"]
rng=np.random.default_rng(7)
center=np.where(dep>=145)[0]
checkpoints=[50,100,200,400,800,1600,3200,6400]
for prog in sys.argv[1:]:
    D={t:[] for t in checkpoints}; nesc=0; ncl=0
    for trial in range(60):
        c0=int(rng.choice(center)); e0=int(rng.integers(4))
        per,fl,order,states=turtle_run(eN,eE,prog,c0,e0,max_cycles=6400)
        if fl==0: ncl+=1; continue
        if fl==1: nesc+=1
        pos=cents[[s[0] for s in states]]; d=np.linalg.norm(pos-cents[c0],axis=1)
        for t in checkpoints:
            if len(d)>=t: D[t].append(d[t-1])
    print("==",prog,"closed",ncl,"escaped",nesc)
    rms={t:np.sqrt(np.mean(np.square(v))) for t,v in D.items() if len(v)>=5}
    ts=np.array(sorted(rms)); r=np.array([rms[t] for t in ts])
    print("   t:",ts.tolist()); print("   rms d:",np.round(r,1).tolist(), " n:",[len(D[t]) for t in ts])
    if len(ts)>=3:
        a,b=np.polyfit(np.log(ts),np.log(r),1); print(f"   fit rms ~ t^{a:.2f}")
