"""collect_cert_roots.py -- run verify_cert.py, indep_check.py and verify_cert_exact.py on every certificate and record,
per file, the program, the base level, the sizes, the designated root pass families with their type, lengths and exact
recurrence, and the results of the three checkers.  Output: results/unbounded_certs/certs_roots.json (used by
notebook_src/programs_md.wl to generate PROGRAMS.md so that every row refers to one certificate file).

Input files are the certificates results/unbounded_certs/cert_*.json and cert2_*.json; the output file certs_roots.json
and summary.json are never read as certificates.  The script stops with a non-zero exit code if a checker fails to run,
if a file does not have the structure of a certificate, or if the root families reported by verify_cert.py are not
exactly the roots designated in the certificate.

usage: python collect_cert_roots.py
"""
import json, os, re, subprocess, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
CD = os.path.join(HERE, "results", "unbounded_certs")
OUT = os.path.join(CD, "certs_roots.json")
ROOT_RX = re.compile(r"^\s+root (\d+) type \((-?\d+), (-?\d+)\) parity (\d): .*?N = \[([^\]]*)\]; "
                     r"O6 exact recurrence N\(m\+2\) = (-?\d+) N\(m\) ([+-]) (-?\d+) \(N\((\d+)\) = (\d+)\)", re.M)
# the constant term is printed by verify_cert.py as "+ 66" or "+ -6": the sign and the (possibly signed) number are both read


def run(script, path):
    r = subprocess.run([sys.executable, os.path.join(HERE, script), path], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("%s failed on %s (exit code %d)\n%s\n%s" % (script, os.path.basename(path), r.returncode, r.stdout[-2000:], r.stderr[-2000:]))
    return r.stdout


files = sorted(p for p in glob.glob(os.path.join(CD, "cert_*.json")) + glob.glob(os.path.join(CD, "cert2_*.json")))
if not files:
    sys.exit("no certificate files found in " + CD)
out = []
for path in files:
    fn = os.path.basename(path)
    cert = json.load(open(path))
    if not (isinstance(cert, dict) and all(k in cert for k in ("prog", "B", "ports", "passes", "roots"))):
        sys.exit("%s does not have the structure of a certificate (keys prog, B, ports, passes, roots)" % fn)
    txt = run("verify_cert.py", path)
    valid = "RESULT: CERTIFICATE VALID" in txt
    ident = re.search(r"O2 port identities checked: (\d+)", txt)
    surplus = re.search(r"smallest margin surplus ([-\d.]+)", txt)
    roots = []
    for m in ROOT_RX.finditer(txt):
        b = int(m.group(8))
        if m.group(7) == "-":
            b = -b
        roots.append({"id": int(m.group(1)), "type": [int(m.group(2)), int(m.group(3))], "parity": int(m.group(4)),
                      "N": [int(x) for x in m.group(5).split(",")], "a": int(m.group(6)), "b": b, "m0": int(m.group(9))})
    if sorted(r["id"] for r in roots) != sorted(cert["roots"]):
        sys.exit("%s designates the roots %s but verify_cert.py reported the root families %s"
                 % (fn, cert["roots"], [r["id"] for r in roots]))
    ri = run("indep_check.py", path)
    summary_line = ri.split("SUMMARY")[-1].splitlines()[0] if "SUMMARY" in ri else ""
    indep_ok = ("SUMMARY" in ri) and ("failures: 0" in ri) and ("False" not in summary_line)
    rx = run("verify_cert_exact.py", path)
    exact_ok = "O5 EXACT OK" in rx
    ex_m = re.search(r"smallest margin surplus = ([-\d.]+)", rx)
    rec = {"file": fn, "prog": cert["prog"], "B": cert["B"], "npass": len(cert["passes"]), "nports": len(cert["ports"]),
           "identities": int(ident.group(1)) if ident else None, "roots": roots, "valid": valid,
           "o5_surplus": float(surplus.group(1)) if surplus else None,
           "o5_exact_ok": exact_ok, "o5_exact_surplus_lb": float(ex_m.group(1)) if ex_m else None,
           "indep_ok": indep_ok}
    out.append(rec)
    print(fn, rec["prog"], "valid" if valid else "INVALID", "indep", "OK" if indep_ok else "FAIL", "O5 exact", "OK" if exact_ok else "FAIL",
          "roots", [(x["id"], x["type"], x["N"][:3], x["a"], x["b"]) for x in roots], flush=True)
json.dump(out, open(OUT, "w"), indent=1)
nroots = sum(len(e["roots"]) for e in out)
print("written", len(out), "records for", len(files), "certificate files;", nroots, "designated roots;",
      "all valid:", all(e["valid"] for e in out), "all indep OK:", all(e["indep_ok"] for e in out), "all O5 exact OK:", all(e["o5_exact_ok"] for e in out))
