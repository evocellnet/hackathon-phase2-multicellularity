#!/usr/bin/env python
"""Export the per-branch H1/H2 statistic (NB-GLM residual) for every branch, for Figure 3.

results/primary/euk/branch_values.tsv holds the sequence-novel counts but not the H1/H2 residual
s_b = log((y_b + 0.5)/(yhat_b + 0.5)) of the control branches. This script recomputes it with exactly the
code path of scripts/analysis/primary_tests.py (mcell.analysis.load_setup, h1h2_counts on
results/gainladder/euk/branch_gain_ladder.tsv, mcell.branchstats.count_residual_statistic; no randomness), then
checks it against the committed tables before writing:
  * stem s_b and yhat equal resid_s / yhat in results/primary/euk/per_transition_H1H2.tsv;
  * the mid-rank percentile u of each primary stem among its matched controls
    (results/primary/euk/matched_sets.tsv) equals the committed u.
It aborts if either check fails, so the exported table can only ever restate the committed analysis.

Usage: python scripts/report/export_h1h2_residuals.py  -> results/report/h1h2_branch_residuals.tsv
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402

PR = os.path.join(ROOT, "results", "primary", "euk")
OUT = os.path.join(ROOT, "results", "report", "h1h2_branch_residuals.tsv")


def main():
    p = an.run_paths("euk")
    _, _, _, bt = an.load_setup("euk", p)
    y = an.h1h2_counts(p["ladder"], bt)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rs = bs.count_residual_statistic(y, bt)
    tab = pd.DataFrame({"seq_novel": y, "yhat": rs["yhat"], "resid_s": rs["s"], "in_U": bt["in_U"]})
    tab.index.name = "branch"

    pt = pd.read_csv(os.path.join(PR, "per_transition_H1H2.tsv"), sep="\t")
    ms = pd.read_csv(os.path.join(PR, "matched_sets.tsv"), sep="\t")
    ms = ms[ms["hypothesis"] == "H1H2"]
    for _, r in pt.iterrows():
        b = r["stem"]
        for col in ("resid_s", "yhat"):
            if not np.isclose(tab.at[b, col], r[col], rtol=1e-4, atol=1e-4):
                sys.exit(f"mismatch {r['transition_id']} {col}: {tab.at[b, col]} vs committed {r[col]}")
        ctrl = ms[(ms["set"] == r["set"]) & (ms["transition_id"] == r["transition_id"])]["control"]
        u = bs.u_score(tab.at[b, "resid_s"], tab.loc[ctrl, "resid_s"].to_numpy())
        if not np.isclose(u, r["u"], atol=1e-4):
            sys.exit(f"u mismatch {r['transition_id']}: {u} vs committed {r['u']}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tab.to_csv(OUT, sep="\t", float_format="%.6g")
    print(f"wrote {OUT}: {len(tab)} branches; stems and u match per_transition_H1H2.tsv")


if __name__ == "__main__":
    main()
