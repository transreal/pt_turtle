import subprocess, sys, json, glob, os, re
PY=sys.executable
B=sys.argv[1]; progs=sys.argv[2:]
outdir=f"batch_certs_B{B}"; os.makedirs(outdir,exist_ok=True)
GROW=re.compile(r"O4 self-reproducing pass \d+ type (\(\d, -?\d\)) parity (\d): occurs (\d+)x in its own two-level expansion of (\d+) passes; N\(m\) for m=(\d+),\d+,... = \[([^\]]*)\]")
out={}
for p in progs:
    res={"prog":p,"valid":[],"invalid":[]}
    try:
        r=subprocess.run([PY,"extract2.py",p,B,outdir,"3","8"],capture_output=True,text=True,timeout=3600)
        m=re.search(r"pass families (\d+)  port families (\d+)  failures (\d+)",r.stdout)
        res["extract"]=m.groups() if m else (r.stdout[-150:]+r.stderr[-200:]).replace(chr(10)," ")
        for f in sorted(glob.glob(f"{outdir}/cert2_{p}_*.json")):
            v=subprocess.run([PY,"verify_cert.py",f,"--fast"],capture_output=True,text=True,timeout=7200)
            if "CERTIFICATE VALID" in v.stdout: res["valid"].append({"cert":f,"npass":len(json.load(open(f))["passes"]),"growth":GROW.findall(v.stdout)[:2]})
            else: res["invalid"].append((f,v.stdout[-400:].replace(chr(10)," | ")))
    except subprocess.TimeoutExpired: res["extract"]="timeout"
    out[p]=res
    best=res["valid"][0] if res["valid"] else None
    print(f"{p:26s} B={B} extract={res['extract']} valid={len(res['valid'])} invalid={len(res['invalid'])} "+(f"passes={best['npass']} growth={best['growth']}" if best else (res['invalid'][0][1][-200:] if res['invalid'] else "")),flush=True)
    json.dump(out,open(f"batch4_B{B}.json","w"),indent=1)
