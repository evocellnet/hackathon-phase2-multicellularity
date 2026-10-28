#!/usr/bin/env python
"""H9 per-unit relative-rate statistics on real data (DESIGN.md §4.6 H9) -- READY, NOT YET RUN.

Guarded: refuses to run unless PREREG.md exists and ``--confirm-prereg-frozen``
is given (the PI decides when). Computes, from the gene-wise fits of
scripts/phylo/gene_trees.py:

* stems: (P, M, S) of transitions in ``--sets`` (default primary) via
  mcell.rates.stem_triples;
* controls: nodes of mcell.branchstats.h9_control_branches (U, two child clades
  of ≥ 2 tips) in fixed orientation (mcell.rates.control_triples);
* per unit × gene r_g table and per-unit median r with 2,000 gene bootstraps
  (seed 11): results/h9/h9_gene_table[_genome_only].tsv, h9_units[...].tsv.

``--genome-only`` drops transcriptome taxa from every gene tree first (M6 rerun).
The relabelling test itself is done with branchstats.h9_values (random control
orientation, seed 11) + branchstats.relabel_test, by the PI's analysis script.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from mcell import rates  # noqa: E402
from mcell.phylo import read_tree  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="euk")
    ap.add_argument("--sets", default="primary", help="comma-separated analysis_set values for stems")
    ap.add_argument("--genome-only", action="store_true")
    ap.add_argument("--exclude-mc-in-controls", action="store_true",
                    help="sensitivity: drop multicellular tips from control units (9/10 euk H9 controls contain some)")
    ap.add_argument("--confirm-prereg-frozen", action="store_true")
    a = ap.parse_args()
    if not (a.confirm_prereg_frozen and os.path.exists(os.path.join(ROOT, "PREREG.md"))):
        sys.exit("H9 on real data is blocked until PREREG.md is frozen and --confirm-prereg-frozen is given.")
    from mcell import branchstats as bs
    guide = read_tree(os.path.join(ROOT, "data", f"guide_tree_{a.set}.nwk"))
    tr = pd.read_csv(os.path.join(ROOT, "data", "transitions.tsv"), sep="\t", dtype=str)
    sp = pd.read_csv(os.path.join(ROOT, "data", "species_final.tsv" if a.set == "euk" else "species_final_prok.tsv"),
                     sep="\t", dtype=str)
    gts = rates.load_gene_trees(os.path.join(ROOT, "runs", "phylo", a.set, "genetrees"))
    tag = ""
    if a.genome_only:
        tx = set(sp.loc[sp["data_type_final"] == "transcriptome", "code"])
        gts = {g: t for g, t in ((g, rates.drop_tips(t, tx)) for g, t in gts.items()) if t is not None}
        tag = "_genome_only"
    stems = rates.stem_triples(tr, guide, sets=a.sets.split(","), kingdom=a.set)
    bt = bs.build_branch_table(guide, sp, tr[tr["kingdom"] == a.set])  # U membership does not use lengths
    ctrl = rates.control_triples(guide, bs.h9_control_branches(guide, bt))
    excl = None
    if a.exclude_mc_in_controls:
        mc = set(sp.loc[sp["mc_state"].isin(["simple", "complex"]), "code"])
        excl = {c[0]: mc for c in ctrl}
        tag += "_ctrl_unicellular_tips"
    tab = rates.h9_gene_table(gts, guide, stems + ctrl, unit_exclude=excl)
    tab["kind"] = tab["unit"].map(lambda u: "stem" if u.startswith("T_") else "control")
    out = os.path.join(ROOT, "results", "h9")
    os.makedirs(out, exist_ok=True)
    tab.to_csv(os.path.join(out, f"h9_gene_table{tag}.tsv"), sep="\t", index=False, float_format="%.6g")
    summ = rates.h9_summary(tab)
    summ.to_csv(os.path.join(out, f"h9_units{tag}.tsv"), sep="\t", float_format="%.6g")
    print(summ.to_string())


if __name__ == "__main__":
    main()
