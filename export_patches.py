"""export the sigma^n(C) rhomb patches used in the proof to JSON for an independent WL check."""
import json, sys, numpy as np
from pt import *
from subst import *
from subst_fast import *
from prove2 import build_tri_atlas, closure_check
prog=sys.argv[1]; n=int(sys.argv[2]); N=int(sys.argv[3]) if len(sys.argv)>3 else 8
entries,keys=build_tri_atlas(N, verbose=True)
out=[]
for ei,e in enumerate(entries):
    col,A,B,C=tris_to_arrays(e["tris"]); col,A,B,C,anc=substitute_v(col,A,B,C,n)
    quad,fat,ti=merge_v(col,A,B,C)
    inside=(anc[ti[:,0]]==e["center"])|(anc[ti[:,1]]==e["center"])
    out.append({"corona":ei,"keys":quad.tolist(),"starts":(np.where(inside)[0]+1).tolist()})  # 1-based starts for WL
json.dump({"prog":prog,"n":n,"patches":out},open(f"patches_{prog}_n{n}.json","w"))
print("exported",len(out),"patches; total rhombs",sum(len(p["keys"]) for p in out))
