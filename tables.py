"""markdown tables for the report from the survey JSON files."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))


def mirror(p):
    return p.translate(str.maketrans("RL", "LR"))


def rot_to_start_S(p):
    """rotate the cyclic word so that it begins with the longest run of S (ties: lexicographic)."""
    n = len(p)
    cands = [p[i:] + p[:i] for i in range(n)]
    def lead(w):
        k = 0
        while k < len(w) and w[k] == "S": k += 1
        return k
    return min(cands, key=lambda w: (-lead(w), w))


def table43():
    rows = {N: {r["prog"]: r for r in json.load(open(os.path.join(HERE, f"surveyall_{N}_8.json")))} for N in (46, 160, 300)}
    out = ["| 命令列 (R形) | 閉軌道率 % (46/160/300) | 類数 (46/160/300) | 最大周期 (300) | 最大半径 (46/160/300) | 判定 |", "|---|---|---|---|---|---|"]
    for p in rows[46]:
        c = [100 * rows[N][p]["closed"] / rows[N][p]["starts"] for N in (46, 160, 300)]
        k = [rows[N][p]["nclasses"] for N in (46, 160, 300)]
        mr = [rows[N][p]["maxradius"] for N in (46, 160, 300)]
        mp = rows[300][p]["maxperiod"]
        pm = rot_to_start_S(mirror(p))
        if mr[2] < 25 and abs(mr[2] - mr[1]) < 0.1 and c[2] > 95:
            verdict = "有界 (一様)"
        elif c[2] < 30 and abs(c[2] - c[1]) < 12:
            verdict = "拡散 (ほぼ全て脱出)"
        elif mr[2] > 100:
            verdict = "階層 (窓に比例)"
        else:
            verdict = "混合"
        out.append(f"| ({pm})* | {c[0]:.1f} / {c[1]:.1f} / {c[2]:.1f} | {k[0]} / {k[1]} / {k[2]} | {mp} | {mr[0]:.1f} / {mr[1]:.1f} / {mr[2]:.1f} | {verdict} |")
    return "\n".join(out)


def table_exits(maxlen=2):
    rows = json.load(open(os.path.join(HERE, "surveyexits_300.json")))
    extra = []
    try:
        for line in open(os.path.join(HERE, "surveyexits_single_300.log")):
            pass
    except Exception:
        pass
    out = ["| 出口列 | 命令列 | 閉軌道率 % | 類数 | 最大周期 | 最大半径 | 判定 |", "|---|---|---|---|---|---|---|"]
    for r in rows:
        if len(r["exits"]) > maxlen:
            continue
        c = 100 * r["closed"] / r["starts"]
        mr = r["maxradius"]
        verdict = "有界候補" if (mr < 25 and c > 95) else ("拡散" if c < 30 else ("階層" if mr > 100 else "混合"))
        out.append(f"| {tuple(r['exits'])} | ({r['prog']})* | {c:.1f} | {r['nclasses']} | {r['maxperiod']} | {mr:.1f} | {verdict} |")
    return "\n".join(out)


def bounded_candidates():
    """all surveyed programs (any survey file) with window-independent small radius and >95% closed at N=300"""
    seen = {}
    for f in os.listdir(HERE):
        if f.startswith("survey_300_") and f.endswith(".json"):
            r = json.load(open(os.path.join(HERE, f)))
            cls = r["classes"]
            mr = max((x["radius"] for x in cls), default=0); c = 100 * r["closed"] / (4 * r["nStarts"])
            if mr < 25 and c > 95:
                seen[r["prog"]] = (c, len(cls), max(x["period"] for x in cls), mr)
    return seen


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "43":
        print(table43())
    elif what == "exits":
        print(table_exits(int(sys.argv[2]) if len(sys.argv) > 2 else 2))
    elif what == "bounded":
        for p, v in sorted(bounded_candidates().items(), key=lambda kv: (len(kv[0]), kv[0])):
            print(f"{p:36s} closed={v[0]:.1f}% classes={v[1]} maxP={v[2]} maxR={v[3]:.1f}")


def table_proofs():
    """table of all programs verified by prove_many.py (proof_many.json) + the four main ones."""
    import json
    from survey_exits import prog_from_exits
    pm = json.load(open(os.path.join(HERE, "proof_many.json")))
    ex = {}
    try:
        for r in json.load(open(os.path.join(HERE, "surveyexits_300.json"))):
            ex[r["prog"]] = tuple(r["exits"])
    except Exception:
        pass
    # reconstruct exit sequences for pure-R programs not in the exits survey
    def exits_of(p):
        if p in ex: return str(ex[p])
        if "L" in p: return "—"
        w = p
        if "SS" not in w: return "(SR)^k"
        k = (w + w).find("SSR")
        if k < 0: return "?"
        w = (w + w)[k + 3:k + 3 + len(w)]  # start right after a block-ending R
        seq = []; i = 0
        while i < len(w):
            j = i
            while w[j:j + 2] == "SR": j += 2
            b = (j - i) // 2
            m = 0
            while j + m < len(w) and w[j + m] == "S": m += 1
            if m < 2 or j + m >= len(w) or w[j + m] != "R": return "?"
            seq.append(3 + b); seq.extend([2] * (m - 2))
            i = j + m + 1
        n = len(seq)
        seq = min(tuple(seq[t:] + seq[:t]) for t in range(n))
        return str(seq)
    out = ["| 命令列 | 出口列 | n | 開始状態数 | 脱出 | 周期の集合 |", "|---|---|---|---|---|---|"]
    for p, d in sorted(pm.items(), key=lambda kv: (len(kv[0]), kv[0])):
        if "L" in p and mirror(p) in pm:
            continue  # show the R form only
        per = sorted(int(k) for k in d["spectrum"])
        out.append(f"| ({rot_to_start_S(p)})* | {exits_of(p)} | {d['n']} | {d['starts']:,} | {d['escape']} | {{{', '.join(map(str, per))}}} |" if d["ok"]
                   else f"| ({rot_to_start_S(p)})* | {exits_of(p)} | {d['n']} | {d['starts']:,} | {d['escape']} | **FAILED** |")
    return "\n".join(out)
