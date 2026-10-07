"""collect one valid certificate per program into certs_final/ and make a summary table"""
import json, glob, os, re, shutil, subprocess, sys
PY=sys.executable
srcs=[("batch_lsys.json",None),("batch3_B6.json",None),("batch3_B8.json",None),("batch4_B8.json",None)]
best={}
for f,_ in srcs:
    if not os.path.exists(f): continue
    d=json.load(open(f))
    for p,r in d.items():
        for v in r.get("valid",[]):
            c=v["cert"].replace("\\","/")
            if not os.path.exists(c): continue
            B=json.load(open(c))["B"]
            key=(v["npass"],B)
            if p not in best or key<best[p][0]: best[p]=(key,c)
os.makedirs("certs_final",exist_ok=True)
GROW=re.compile(r"O4 self-reproducing pass (\d+) type (\(\d, -?\d\)) parity (\d): occurs (\d+)x in its own two-level expansion of (\d+) passes; N\(m\) for m=(\d+),\d+,... = \[([^\]]*)\]")
rows=[]
for p,(key,c) in sorted(best.items(),key=lambda kv:(len(kv[0]),kv[0])):
    dst=f"certs_final/cert_{p}.json"; shutil.copyfile(c,dst)
    v=subprocess.run([PY,"verify_cert.py",dst,"--fast"],capture_output=True,text=True)
    ok="CERTIFICATE VALID" in v.stdout
    g=GROW.findall(v.stdout)
    cert=json.load(open(dst))
    o1=re.search(r"O1 base simulations at level (\d+): (\d+) passes",v.stdout); o2=re.search(r"O2 port identities checked: (\d+)",v.stdout); o5=re.search(r"smallest margin surplus ([\d.]+)",v.stdout)
    row={"prog":p,"valid":ok,"B":cert["B"],"npass":len(cert["passes"]),"nports":len(cert["ports"]),"O1":int(o1.group(2)) if o1 else None,"O2":int(o2.group(1)) if o2 else None,"margin":float(o5.group(1)) if o5 else None,
         "growth":[{"type":x[1],"parity":int(x[2]),"mult":int(x[3]),"total":int(x[4]),"m0":int(x[5]),"N":[int(t) for t in x[6].split(",")]} for x in g[:2]]}
    o6=re.findall(r"O6 pass (\d+) type (\(\d, -?\d\)) parity (\d): N\(m\+2\) = (\d+) N\(m\) \+ (-?\d+), N\((\d+)\) = (\d+)  \[verified for (\d+) steps",v.stdout)
    row["O6"]=[{"type":x[1],"parity":int(x[2]),"a":int(x[3]),"b":int(x[4]),"m0":int(x[5]),"N0":int(x[6]),"steps":int(x[7])} for x in o6[:2]]
    rows.append(row)
    print(p,ok,row["B"],row["npass"],row["growth"][0]["N"][:5] if row["growth"] else None,flush=True)
json.dump(rows,open("certs_final/summary.json","w"),indent=1)
