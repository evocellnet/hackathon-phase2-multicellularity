#!/usr/bin/env python
"""Apply the pre-registered W5 decision rule (results/r2/w5/decision_operationalisation.txt) to the R1 outputs of
r1_compare.py; writes results/r2/w5/r1/verdict.txt and results/r2/w5/r1/DONE.json (md5 of all r1 outputs)."""
import hashlib, json, os, sys, time
import numpy as np
import pandas as pd

OUT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/w5/r1"
DEEP = ["C03", "C06", "C07"]


def diff_ci(a1, n1, a0, n0):
    p1, p0 = a1 / n1, a0 / n0
    se = np.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    d = p1 - p0
    return d, d - 1.96 * se, d + 1.96 * se


def main():
    for f in ("verdict.txt", "DONE.json"):
        if os.path.exists(os.path.join(OUT, f)):
            sys.exit(f"refusing to overwrite {f}")
    S = pd.read_csv(os.path.join(OUT, "xval_agreement_summary.tsv"), sep="\t")
    ctl = {v: pd.read_csv(os.path.join(OUT, f"controls_{v}.tsv"), sep="\t").set_index("control")
           for v in ("default", "3di", "3di_guard")}
    L = []
    base = ctl["default"]
    for v in ("3di", "3di_guard"):
        c = ctl[v]
        n_deep = int((c.loc[DEEP, "recovered"] == "yes").sum())
        n_deep0 = int((base.loc[DEEP, "recovered"] == "yes").sum())
        lost = [k for k in base.index if base.loc[k, "recovered"] == "yes" and c.loc[k, "recovered"] != "yes"]
        comp = "; ".join(f"{k}: {base.loc[k, 'components']} -> {c.loc[k, 'components']}" for k in DEEP)
        crit1 = n_deep > n_deep0 and not lost
        L.append(f"[{v}] (1) deep-node controls recovered {n_deep}/3 vs default {n_deep0}/3; lost elsewhere: {lost or 'none'}"
                 f"\n      components {comp}\n      -> improved: {crit1}")
        crit2 = True
        for evt in ("gain", "dup"):
            for direc in ("ancestra->eukaryoma", "eukaryoma->ancestra"):
                r1 = S[(S.variant == v) & (S.event == evt) & (S.direction == direc) & (S.group == "both_stems")].iloc[0]
                r0 = S[(S.variant == "default") & (S.event == evt) & (S.direction == direc) & (S.group == "both_stems")].iloc[0]
                d, lo, hi = diff_ci(r1["same"], r1["n_mapped"], r0["same"], r0["n_mapped"])
                prim = direc == "ancestra->eukaryoma"
                if prim and r1["frac_same"] < r0["frac_same"]:
                    crit2 = False
                L.append(f"      {'PRIMARY' if prim else 'secondary'} {evt:4s} {direc}: {r1['same']}/{r1['n_mapped']} = "
                         f"{r1['frac_same']:.4f} vs default {r0['same']}/{r0['n_mapped']} = {r0['frac_same']:.4f}; "
                         f"difference {d * 100:+.2f} pp [95% CI {lo * 100:+.2f}, {hi * 100:+.2f}]")
        L.append(f"      -> (2) EukaryOMA agreement not lowered: {crit2}")
        L.append(f"      VERDICT {v}: {'REPLACES default' if crit1 and crit2 else 'does NOT replace the default'}")
    txt = (f"W5 task A verdict (generated {time.strftime('%Y-%m-%d %H:%M')} by scripts/r2/w5/r1_verdict.py)\n"
           + "\n".join(L) + "\n")
    open(os.path.join(OUT, "verdict.txt"), "w").write(txt)
    print(txt)
    md5 = {f: hashlib.md5(open(os.path.join(OUT, f), "rb").read()).hexdigest() for f in sorted(os.listdir(OUT))}
    json.dump(dict(status="done", finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"), md5=md5),
              open(os.path.join(OUT, "DONE.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
