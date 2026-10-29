#!/usr/bin/env python
"""Errata item 1: where does 'depth SMD 1.77' come from, and what is the correct number?

SMD (mcell.branchstats.match_controls) = (stem value - mean of its matched controls) / SD over the
control pool U; the 'ALL' row = (mean over the 7 stems - mean over all 7 x 20 matched controls, with
multiplicity) / SD_U. 'depth' = number of edges from the root of the guide tree to the stem's lower
node. Recomputes the depth SMD under four settings:
  A  round-1 analysis setting: current transitions (D7 path stems), real branch lengths
     -> must equal results/primary/euk/smd.tsv (ALL = 2.06)
  B  path stems, placeholder branch lengths 1.0 (setting of the frozen power table)
  C  pre-D7 transitions (git 8b77079: fungal stem = Blastocladiomycota_Amastigomycota, volvocine
     stem = TGV_clade), placeholder branch lengths (setting of decision D6, t = 1.6 h)
  D  unmatched: 7 stems vs the whole pool U (no matching), path stems
Output: results/round1_errata/smd_depth_provenance.tsv
"""
import io
import os
import subprocess
import sys
import warnings

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from mcell import branchstats as bs  # noqa: E402


def depth_smd(tr, bl_path):
    traits = pd.read_csv(os.path.join(ROOT, "data/species_traits.tsv"), sep="\t")
    tree = bs.load_tree(os.path.join(ROOT, "data/guide_tree_euk.nwk"))
    bl = bs.branch_lengths_from_tree(tree, bl_path)
    bt = bs.build_branch_table(tree, traits, tr[tr.kingdom == "euk"], branch_lengths=bl)
    stems = bs.stem_branches(bt)
    m = bs.match_controls(bt, stems)
    d = m.smd[m.smd.covariate == "depth"].copy()
    pool = list(bt.index[bt.in_U])
    x = bt.loc[pool, "depth"]
    xs = bt.loc[list(stems.values), "depth"].mean()
    unmatched = dict(transition_id="ALL_unmatched_vs_U", stem="", covariate="depth", stem_value=xs,
                     control_mean=x.mean(), smd=(xs - x.mean()) / x.std(ddof=1))
    return d, unmatched, m.covariates


def main():
    warnings.simplefilter("ignore")
    cur = pd.read_csv(os.path.join(ROOT, "data/transitions.tsv"), sep="\t")
    old = pd.read_csv(io.StringIO(subprocess.check_output(
        ["git", "-C", ROOT, "show", "8b77079:data/transitions.tsv"], text=True)), sep="\t")
    real = os.path.join(ROOT, "data/species_tree_bl_euk.nwk")
    out = []
    for lab, tr, blp, desc in [
        ("A", cur, real, "round-1 analysis: path stems (D7), real branch lengths [= results/primary/euk/smd.tsv]"),
        ("B", cur, None, "path stems, placeholder branch lengths 1.0 (frozen power-table setting)"),
        ("C", old, None, "pre-D7 stems (git 8b77079), placeholder branch lengths (decision D6, t=1.6 h)"),
    ]:
        d, um, cov = depth_smd(tr, blp)
        d["setting"], d["description"], d["matching_covariates"] = lab, desc, ",".join(cov)
        out.append(d)
        if lab == "A":
            um.update(setting="D", description="no matching: 7 stems vs all 30 pool branches (SD_U)",
                      matching_covariates="none")
            out.append(pd.DataFrame([um]))
    df = pd.concat(out, ignore_index=True)
    ref = pd.read_csv(os.path.join(ROOT, "results/primary/euk/smd.tsv"), sep="\t")
    ref = ref[(ref.hypothesis == "H1H2") & (ref.set == "primary") & (ref.covariate == "depth")].set_index("transition_id")
    a = df[df.setting == "A"].set_index("transition_id")
    assert (abs(a.smd - ref.loc[a.index, "smd"]) < 1e-3).all(), "setting A does not reproduce smd.tsv"
    df.to_csv(os.path.join(ROOT, "results/round1_errata/smd_depth_provenance.tsv"), sep="\t", index=False,
              float_format="%.4g")
    print(df[["setting", "transition_id", "stem", "stem_value", "control_mean", "smd"]].to_string(index=False))


if __name__ == "__main__":
    main()
