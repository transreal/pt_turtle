"""collect_cert_roots.py -- run verify_cert.py (and verify_cert_exact.py) on every certificate and record, per file,
the program, the base level, the sizes, the designated root pass families with their type, lengths and exact
recurrence, and the checker results.  Output: results/unbounded_certs/certs_roots.json (used by
notebook_src/programs_md.wl to generate PROGRAMS.md so that every row refers to one certificate file).

usage: python collect_cert_roots.py
"""
import json, os, re, subprocess, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
CD = os.path.join(HERE, "results", "unbounded_certs")
summary = {e["prog"]: e for e in json.load(open(os.path.join(CD, "summary.json")))}
out = []
for path in sorted(glob.glob(os.path.join(CD, "cert*.json"))):
    cert = json.load(open(path)); fn = os.path.basename(path)
    r = subprocess.run([sys.executable, os.path.join(HERE, "verify_cert.py"), path], capture_output=True, text=True)
    txt = r.stdout
    valid = "RESULT: CERTIFICATE VALID" in txt
    ident = re.search(r"O2 port identities checked: (\d+)", txt)
    surplus = re.search(r"smallest margin surplus ([-\d.]+)", txt)
    roots = []
    for m in re.finditer(r"^\s+root (\d+) type \((-?\d+), (-?\d+)\) parity (\d): .*?N = \[([^\]]*)\]; O6 exact recurrence N\(m\+2\) = (-?\d+) N\(m\) ([+-]) (\d+) \(N\((\d+)\) = (\d+)\)", txt, flags=re.M):
        roots.append({"id": int(m.group(1)), "type": [int(m.group(2)), int(m.group(3))], "parity": int(m.group(4)),
                      "N": [int(x) for x in m.group(5).split(",")], "a": int(m.group(6)),
                      "b": int(m.group(8)) * (1 if m.group(7) == "+" else -1), "m0": int(m.group(9))})
    ri = subprocess.run([sys.executable, os.path.join(HERE, "indep_check.py"), path], capture_output=True, text=True)
    indep_ok = ("failures: 0" in ri.stdout) and ("SUMMARY" in ri.stdout) and ("False" not in ri.stdout.split("SUMMARY")[-1].splitlines()[0])
    rx = subprocess.run([sys.executable, os.path.join(HERE, "verify_cert_exact.py"), path], capture_output=True, text=True)
    exact_ok = "O5 EXACT OK" in rx.stdout
    ex_m = re.search(r"smallest margin surplus = ([-\d.]+)", rx.stdout)
    rec = {"file": fn, "prog": cert["prog"], "B": cert["B"], "npass": len(cert["passes"]), "nports": len(cert["ports"]),
           "identities": int(ident.group(1)) if ident else None, "roots": roots, "valid": valid,
           "o5_surplus": float(surplus.group(1)) if surplus else None,
           "o5_exact_ok": exact_ok, "o5_exact_surplus_lb": float(ex_m.group(1)) if ex_m else None,
           "indep_ok": indep_ok}
    out.append(rec)
    print(fn, rec["prog"], "valid" if valid else "INVALID", "roots", [(x["id"], x["type"], x["N"][:3], x["a"], x["b"]) for x in roots],
          "O5 exact", "OK" if exact_ok else "FAIL", flush=True)
json.dump(out, open(os.path.join(CD, "certs_roots.json"), "w"), indent=1)
print("written", len(out), "records")
