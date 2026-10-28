#!/usr/bin/env python
"""Holm adjustment across the primary family (PREREG.md: H1/H2, H4, H5, H9, H12).

    python scripts/analysis/holm_primary.py [--registry results/analysis_registry.tsv] [--force]

Inputs: results/primary/euk/primary_tests.tsv (rows set == 'primary' and analysis == 'primary', written by
primary_tests.py; H12 counts only if its row says tier == 'primary', i.e. >= 2 loss-expected stems remained)
and results/h9/h9_tests.tsv (row label == 'primary').
Output: results/primary/euk/holm.tsv; one registry row per test (holm_p filled).
Refuses to run (exit 2) unless every member of the family is present, unless --force (then the table is
labelled PROVISIONAL and nothing is written to the registry).
"""
from __future__ import annotations

import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import branchstats as bs  # noqa: E402

PRIMARY = os.path.join(ROOT, "results", "primary", "euk", "primary_tests.tsv")
H9 = os.path.join(ROOT, "results", "h9", "h9_tests.tsv")
OUT = os.path.join(ROOT, "results", "primary", "euk", "holm.tsv")
FAMILY = ("H1H2", "H4", "H5", "H9", "H12")


def collect(primary_path: str = PRIMARY, h9_path: str = H9) -> pd.DataFrame:
    """One row per primary test with its p (and why a test is missing)."""
    rows = []
    if os.path.exists(primary_path):
        t = pd.read_csv(primary_path, sep="\t")
        t = t[(t["set"] == "primary") & (t["analysis"] == "primary")]
        for _, r in t.iterrows():
            rows.append(dict(hypothesis=r["hypothesis"], p=float(r["p"]), p_mcse=r.get("p_mcse", np.nan),
                             tier=r.get("tier", "primary"), source=os.path.relpath(primary_path, ROOT),
                             n=r.get("n"), k=r.get("k"), scheme=r.get("scheme")))
    if os.path.exists(h9_path):
        h = pd.read_csv(h9_path, sep="\t")
        h = h[h["label"] == "primary"]
        if len(h):
            r = h.iloc[0]
            rows.append(dict(hypothesis="H9", p=float(r["p"]), p_mcse=r["p_mcse"], tier="primary",
                             source=os.path.relpath(h9_path, ROOT), n=r["n"], k=r["k"], scheme=r["scheme"]))
    return pd.DataFrame(rows, columns=["hypothesis", "p", "p_mcse", "tier", "source", "n", "k", "scheme"])


def compute(primary_path: str = PRIMARY, h9_path: str = H9) -> tuple[pd.DataFrame, list]:
    """(table with holm_p over the tests whose tier is primary, list of missing family members)."""
    t = collect(primary_path, h9_path)
    have = set(t["hypothesis"])
    demoted = set(t.loc[t["tier"] != "primary", "hypothesis"])
    missing = [h for h in FAMILY if h not in have]
    fam = t[t["tier"] == "primary"].copy()
    fam["holm_p"] = bs.holm(fam["p"].to_numpy())
    out = t.merge(fam[["hypothesis", "holm_p"]], on="hypothesis", how="left")
    out["family_size"] = len(fam)
    out["note"] = ["demoted to exploratory (not in Holm family)" if h in demoted else "" for h in out["hypothesis"]]
    return out, missing


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--registry", default=os.path.join(ROOT, "results", "analysis_registry.tsv"))
    ap.add_argument("--primary", default=PRIMARY)
    ap.add_argument("--h9", default=H9)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--force", action="store_true", help="write a PROVISIONAL table even if tests are missing")
    a = ap.parse_args(argv)
    out, missing = compute(a.primary, a.h9)
    if missing and not a.force:
        print(f"Holm not computed: missing primary tests {missing}", file=sys.stderr)
        return 2
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    if missing:
        out["note"] = out["note"] + f" PROVISIONAL: missing {missing}"
    out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(out.to_string(index=False))
    if not missing:
        for _, r in out[out["holm_p"].notna()].iterrows():
            bs.append_registry(a.registry, dict(
                hypothesis=r["hypothesis"], statistic="Holm-adjusted primary p (holm_primary.py)", n=r["n"],
                k=r["k"], p=r["p"], p_mcse=r["p_mcse"], holm_p=r["holm_p"], tier="primary", scheme=r["scheme"],
                input_md5s=[x for x in (a.primary, a.h9) if os.path.exists(x)],
                notes=f"Holm across {int(r['family_size'])} primary tests ({', '.join(out.loc[out['holm_p'].notna(), 'hypothesis'])}); source {r['source']}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
