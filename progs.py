"""enumerate reduced cyclic programs: words over S,R,L with >=1 S, >=1 turn, no two adjacent turns
(cyclically), canonical under cyclic rotation and R<->L mirror."""
import itertools
def canon(w):
    n=len(w); cands=[]
    for m in (w, w.translate(str.maketrans("RL","LR"))):
        for i in range(n): cands.append(m[i:]+m[:i])
    return min(cands)
def programs(maxlen, minlen=2):
    out=set()
    for n in range(minlen,maxlen+1):
        for t in itertools.product("SRL",repeat=n):
            w="".join(t)
            if "S" not in w or w.count("S")==n: continue
            ok=True
            for i in range(n):
                if w[i]!="S" and w[(i+1)%n]!="S": ok=False;break
            if not ok: continue
            out.add(canon(w))
    return sorted(out,key=lambda s:(len(s),s))
if __name__=="__main__":
    import sys
    ps=programs(int(sys.argv[1]) if len(sys.argv)>1 else 7)
    print(len(ps)); print(" ".join(ps))
