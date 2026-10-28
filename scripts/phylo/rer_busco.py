#!/usr/bin/env python
"""EXPLORATORY gene-level relative evolutionary rates (RER) on the 250 BUSCO markers (DESIGN §4.5).

mcell.rates.relative_rates on the per-gene IQ-TREE -te fits (runs/phylo/euk/genetrees), min 15 species
per gene (round-1 review M8). Branch sets on the guide tree (branchstats.build_branch_table classes):

* foreground = the 7 primary transitions: for each, the clade rooted at the first stem edge below P
  (M, or the intermediate node X of a D7 path stem); inside it, edges whose class is stem (of a
  primary transition), nested_stem or within_multicellular;
* excluded (NaN) = other edges inside those clades (reversals, e.g. yeasts; unicellular Allomyces under
  the fungal path node; the exploratory Tetrabaena stem) and every edge outside them that is not
  within_unicellular (exploratory multicellular lineages, grey zones, other stems);
* background = within_unicellular edges outside the foreground clades. Root-child edges are dropped.

Statistic: Kendall tau_b(RER, foreground) per gene; permulation null = random non-overlapping clade
placements matched on tip count (1,000, seed 11; excluded edges stay excluded); two-sided p; BH-FDR over
genes with >= 2 transitions represented. Consistency: per gene and transition, the mean RER over the
foreground edges of that clade; a gene "shifts consistently" if >= 4 transitions have clade-mean RER
of the same sign with |mean| >= 0.5 and none >= 0.5 in the opposite direction (RER is a standardised
residual). The expected number of such genes
is obtained from the same permuted clade placements.
Outputs: results/rer/{rer_matrix.tsv, rer_genes.tsv, rer_clade_means.tsv, rer_consistency_null.tsv,
summary.md}; registry rows (tier exploratory).
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "phylo"))
from mcell import rates, branchstats as bs  # noqa: E402
import h9_rates as H  # noqa: E402

OUT = os.path.join(ROOT, "results", "rer")
SEED, N_PERM, MIN_SP, THR, MIN_T = 11, 1000, 15, 0.5, 4


def busco_descriptions():
    d = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "data", "busco_runs", "*", "run_eukaryota_odb10", "full_table.tsv")))[:40]:
        for line in open(f):
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) >= 7 and x[0] not in d and x[6]:
                d[x[0]] = x[6]
    return d


def branch_sets(tr, tree, bt):
    prim = tr[(tr.kingdom == "euk") & (tr.analysis_set == "primary")]
    roots, tid_of = [], {}
    for _, r in prim.iterrows():
        path = [x for x in str(r["stem_path"]).split(">") if x and x != "nan"]
        root = path[1] if len(path) > 1 else r["M"]
        roots.append(root)
        tid_of[root] = r["transition_id"]
    primary_ids = set(prim.transition_id)
    in_clade = {}
    for root in roots:
        for n in rates.foreground_nodes(tree, [root]):
            in_clade[n] = root
    fg, excl = set(), set()
    for b, row in bt.iterrows():
        cls = row["branch_class"]
        if b in in_clade:
            tids = set(str(row["transition_id"]).split(";")) if row["transition_id"] else set()
            if cls in ("within_multicellular", "nested_stem") or (cls == "stem" and tids & primary_ids):
                fg.add(b)
            else:
                excl.add(b)
        elif cls != "within_unicellular":
            excl.add(b)
    return roots, tid_of, fg, excl


def consistent(cm: np.ndarray, thr=THR, min_t=MIN_T):
    """>= min_t transitions with clade-mean RER beyond +thr (or -thr) and none beyond the opposite threshold."""
    up = np.nansum(cm >= thr, axis=1)
    dn = np.nansum(cm <= -thr, axis=1)
    return ((up >= min_t) & (dn == 0)) | ((dn >= min_t) & (up == 0)), up, dn


def main():
    warnings.simplefilter("ignore")
    os.makedirs(OUT, exist_ok=True)
    tr, traits, tree, bt, gts = H.setup("euk")
    roots, tid_of, fg, excl = branch_sets(tr, tree, bt)
    rer = rates.relative_rates(gts, tree, min_sp=MIN_SP)
    rer = rer.dropna(how="all")
    rer.loc[:, [c for c in rer.columns if c in excl]] = np.nan
    rer.to_csv(os.path.join(OUT, "rer_matrix.tsv"), sep="\t", float_format="%.4g")
    # observed foreground = clades minus excluded (excluded columns are NaN in every gene)
    res = rates.rer_permulation_test(rer, tree, roots, n_perm=N_PERM, seed=SEED)
    cm = rates.clade_mean_rer(rer, tree, roots)
    cm.columns = [tid_of[c] for c in cm.columns]
    n_trans = cm.notna().sum(1)
    res["n_transitions"] = n_trans
    res.loc[n_trans < 2, ["p_perm"]] = np.nan
    res["q_bh"] = rates.bh(res["p_perm"].to_numpy())
    is_c, up, dn = consistent(cm.to_numpy())
    res["n_trans_up"] = up
    res["n_trans_down"] = dn
    res["consistent_ge4"] = is_c
    desc = busco_descriptions()
    res.insert(0, "description", [desc.get(g, "") for g in res.index])
    res = res.sort_values("p_perm")
    res.to_csv(os.path.join(OUT, "rer_genes.tsv"), sep="\t", float_format="%.5g")
    cm.to_csv(os.path.join(OUT, "rer_clade_means.tsv"), sep="\t", float_format="%.4g")
    # null for the consistency count: same permuted clade placements
    perm_roots = rates.permuted_clade_roots(tree, roots, n_perm=1000, seed=SEED)
    null = []
    for pr in perm_roots:
        c0, _, _ = consistent(rates.clade_mean_rer(rer, tree, pr).to_numpy())
        null.append(int(c0.sum()))
    null = np.array(null)
    obs_c = int(is_c.sum())
    pd.DataFrame({"n_consistent_genes": null}).to_csv(os.path.join(OUT, "rer_consistency_null.tsv"), sep="\t",
                                                       index=False)
    p_c = (1 + (null >= obs_c).sum()) / (1 + len(null))
    tested = res["p_perm"].notna()
    meta = dict(finished=dt.datetime.now().isoformat(timespec="seconds"), n_genes_rer=len(rer),
                n_genes_tested=int(tested.sum()), n_fg_edges=len(fg & set(rer.columns)),
                n_excluded_edges=len(excl & set(rer.columns)),
                n_bg_edges=int(len(rer.columns) - len(fg & set(rer.columns)) - len(excl & set(rer.columns))),
                clade_roots=tid_of, n_perm=res.attrs.get("n_perm"), seed=SEED, min_sp=MIN_SP,
                n_p05=int((res["p_perm"] < 0.05).sum()), n_q10=int((res["q_bh"] < 0.10).sum()),
                n_q05=int((res["q_bh"] < 0.05).sum()), expected_p05=0.05 * int(tested.sum()),
                consistency=dict(threshold=THR, min_transitions=MIN_T, observed=obs_c,
                                 null_mean=float(null.mean()), null_q95=float(np.quantile(null, 0.95)),
                                 p=float(p_c), n_null=len(null)))
    json.dump(meta, open(os.path.join(OUT, "rer_run.json"), "w"), indent=1)
    bs.append_registry(H.REG, dict(
        hypothesis="Q5 gene-level RER (BUSCO 250)", statistic="Kendall tau_b(RER, fg=primary stems+mc crowns), "
        "permulation p per gene; estimate = n genes q<0.10", n=int(tested.sum()), k="", estimate=meta["n_q10"],
        ci_low="NA", ci_high="NA", p="NA", p_mcse="NA", holm_p="NA", tier="exploratory", seed=SEED, n_draws=N_PERM,
        scheme="clade-placement permulation", input_md5s=[os.path.join(ROOT, "runs/phylo/euk/genetrees/branch_lengths.tsv")],
        notes=f"label=exploratory; BH over genes; n_p05={meta['n_p05']} (expected {meta['expected_p05']:.1f})"))
    bs.append_registry(H.REG, dict(
        hypothesis="Q5 gene-level RER consistency (>=4 transitions)", statistic="n genes with >=4 transitions of "
        "same-sign clade-mean RER |.|>=0.5, none opposite; null = permuted clade placements", n=int(tested.sum()), k=len(roots),
        estimate=obs_c, ci_low="NA", ci_high="NA", p=float(p_c), p_mcse="NA", holm_p="NA", tier="exploratory",
        seed=SEED, n_draws=len(null), scheme="clade-placement permulation",
        input_md5s=[os.path.join(ROOT, "runs/phylo/euk/genetrees/branch_lengths.tsv")],
        notes=f"label=exploratory; null mean {null.mean():.2f}, q95 {np.quantile(null, 0.95):.0f}"))
    print(json.dumps(meta, indent=1))
    print(res.head(15).to_string())


if __name__ == "__main__":
    main()
