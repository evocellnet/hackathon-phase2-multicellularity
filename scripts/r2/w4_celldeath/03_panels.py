#!/usr/bin/env python
"""W4 cell-death follow-up, step 3 (descriptive, exploratory): are cell-death families gained (sequence-novel) or
duplicated on the simple->complex event branches of the 8 dense panels?

Input : results/r2/panels/<R>/w3_families.tsv.gz (all laddered branches; family Pfam set = union over up to 3
        representatives, transposon domains removed; architecture = longest annotated representative),
        results/r2/panels/<R>/event_table.tsv, branch_table.tsv. Primary panels only (R4, not R4_curated).
For every event branch and kind (seq_novel, dup): number of families, number with any Pfam domain, number carrying
>= 1 of the 22 cell-death families. Background = the panel's control branches (internal controls + tip controls),
pooled. Test (exploratory): one-sided Fisher exact on cell-death vs other Pfam-annotated families, event vs
background; BH over all event x kind rows of the simple->complex events (kind nested or origin_complex).
Outputs (results/r2/w4_celldeath/):
  panels_event_celldeath.tsv     event x kind counts, background fraction, odds ratio, p, q
  panels_celldeath_families.tsv  every cell-death family on an event branch (family id, n_genes, Pfam, architecture)
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = os.path.join(P, "results/r2/w4_celldeath")
PANELS = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
CD = pd.read_csv(os.path.join(P, "data/r2_sets/cell_death.tsv"), sep="\t")
CDACC = set(CD.pfam)


def out(name):
    path = os.path.join(OUT, name)
    if os.path.exists(path):
        sys.exit(f"REFUSE: {path} exists")
    return path


def main():
    rows, fam_rows = [], []
    for R in PANELS:
        d = os.path.join(P, "results/r2/panels", R)
        w = pd.read_csv(os.path.join(d, "w3_families.tsv.gz"), sep="\t", dtype=str).fillna("")
        w["accs"] = w.pfam.map(lambda s: set(x for x in s.split(";") if x))
        w["annot"] = w.accs.map(bool)
        w["cd"] = w.accs.map(lambda s: bool(s & CDACC))
        et = pd.read_csv(os.path.join(d, "event_table.tsv"), sep="\t")
        bt = pd.read_csv(os.path.join(d, "branch_table.tsv"), sep="\t")
        ctrl = set(bt.loc[(bt.control == True) | (bt.tip_control == True), "branch_child"])  # noqa: E712
        for kind in ("seq_novel", "dup"):
            bg = w[(w.kind == kind) & w.branch_child.isin(ctrl)]
            bg_cd, bg_ann = int(bg.cd.sum()), int(bg.annot.sum())
            for e in et.itertuples():
                x = w[(w.kind == kind) & (w.branch_child == e.M)]
                n_cd, n_ann = int(x.cd.sum()), int(x.annot.sum())
                if n_ann > 0 and bg_ann > 0:
                    orr, p = fisher_exact([[n_cd, n_ann - n_cd], [bg_cd, bg_ann - bg_cd]], alternative="greater")
                else:
                    orr, p = np.nan, np.nan
                rows.append(dict(panel=R, event_id=e.event_id, event_kind=e.kind, branch=e.M,
                                 from_state=e.from_state, to_state=e.to_state, kind=kind, n_families=len(x),
                                 n_pfam_annotated=n_ann, n_celldeath=n_cd,
                                 frac_celldeath_of_annotated=n_cd / n_ann if n_ann else np.nan,
                                 n_control_branches=len(ctrl), bg_annotated=bg_ann, bg_celldeath=bg_cd,
                                 bg_frac=bg_cd / bg_ann if bg_ann else np.nan, odds_ratio=orr, p_one_sided=p,
                                 genes_in_celldeath_families=int(x.loc[x.cd, "n_genes"].astype(int).sum())))
                for f in x[x.cd].itertuples():
                    fam_rows.append(dict(panel=R, event_id=e.event_id, event_kind=e.kind, kind=kind,
                                         family_id=f.family_id, n_genes=int(f.n_genes), rep_species=f.rep_species,
                                         celldeath_pfam=";".join(sorted(f.accs & CDACC)), pfam_names=f.pfam_names,
                                         architecture=f.architecture, go=f.go))
    t = pd.DataFrame(rows)
    sc = t.event_kind.isin(["nested", "origin_complex"]) & t.p_one_sided.notna()
    t["q_bh_simple_to_complex"] = np.nan
    t.loc[sc, "q_bh_simple_to_complex"] = multipletests(t.loc[sc, "p_one_sided"], method="fdr_bh")[1]
    t.to_csv(out("panels_event_celldeath.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(fam_rows).to_csv(out("panels_celldeath_families.tsv"), sep="\t", index=False)
    print(t[sc].drop(columns=["branch", "from_state", "to_state", "n_control_branches"]).to_string(index=False))


if __name__ == "__main__":
    main()
