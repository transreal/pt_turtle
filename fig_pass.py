import json, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from lsys import *
from subst import mulphi, add, embed
from pt import parse_prog
cert=json.load(open(sys.argv[1])); prog=parse_prog(cert["prog"]).tolist(); B=cert["B"]
root=cert["roots"][0]; f=cert["ports"][str(root)]; typ=tuple(f["typ"]); q=f["parity"]
def port_at(pid,m):
    g=cert["ports"][str(pid)]; mm=B if B%2==g["parity"] else B+1
    t,h,i=tuple(g["ref"][0]),tuple(g["ref"][1]),g["ref"][2]; dt,dh=tuple(g["delta"][0]),tuple(g["delta"][1])
    while mm<m: t=add(mulphi(mulphi(t)),dt); h=add(mulphi(mulphi(h)),dh); mm+=2
    return (t,h,i)
levels=[int(x) for x in sys.argv[2].split(",")]
fig,axes=plt.subplots(1,len(levels),figsize=(9*len(levels),8))
axes=np.atleast_1d(axes)
for ax,m in zip(axes,levels):
    st=Supertile(typ,m)
    bmap={}
    br,be=np.where(st.eN<0)
    for r,e in zip(br.tolist(),be.tolist()): bmap[(st.key(r,e+1),st.key(r,e))]=(r,e)
    t,h,i=port_at(root,m); r,eb=bmap[(t,h)]
    ex,steps,path=st.run_pass(prog,r,(eb+2)&3,i+1,record=True)
    V=st.quad@EVEC
    ax.add_collection(PolyCollection(V,facecolors="white",edgecolors="0.85",linewidths=0.2))
    ax.add_collection(PolyCollection(V[sorted(set(path))],facecolors="#4c8be0",edgecolors="none",alpha=0.8))
    tri=std_triangle(*typ); P=np.array([embed(tri[1]),embed(tri[2]),embed(tri[3]),embed(tri[1])])*PHI**m
    ax.plot(P[:,0],P[:,1],"k-",lw=1.2)
    c=st.cent[path]; ax.plot(c[:,0],c[:,1],"-",color="#d0342c",lw=0.5)
    ax.plot(*c[0],"o",color="green",ms=7); ax.plot(*c[-1],"s",color="k",ms=7)
    ax.set_aspect("equal"); ax.autoscale(); ax.axis("off")
    ax.set_title(f"level {m}: root pass of ({cert['prog']})* through the obtuse supertriangle, {steps+1} S-moves, {len(set(path))} rhombs")
fig.tight_layout(); fig.savefig(sys.argv[3],dpi=110); print("saved",sys.argv[3])
