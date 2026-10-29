#!/usr/bin/env python
"""Descriptive re-computation of the round-1 H1/H2 gene-content statistic with EukaryOMA counts.

Round 1 (scripts/analysis/primary_tests.py, mcell.branchstats): y_b = sequence-novel gains on branch
b (path-summed for path stems); NB-GLM y ~ log_bl + log_n_desc + log_n_sister + frac_tx_desc fitted
on the control pool U (leave-one-out for pool branches); residual s_b = log((y+.5)/(yhat+.5));
u_b = rank of s_b among its k matched controls; T = mean u over the primary stems.

Here, on the 125-species pruned tree, the same machinery (bs.count_residual_statistic, bs.u_score)
is applied to y = supported HOG gains (gain events of families with >= 2 species) from
(a) Ancestra projected onto the pruned tree and (b) EukaryOMA. Covariates are recomputed on the
pruned tree (tips below / in the sister clades, transcriptome fraction; branch length = sum of the
round-1 branch lengths of the original edges the pruned branch covers). Pool = round-1 controls that
map exactly; matched sets = round-1 H1/H2 matched sets restricted to those controls. No p-values:
the species set, the statistic (HOG gains, not sequence-novel gains) and the pool differ from the
pre-registered test, so this is a direction/size comparison only.
Output: results/eukaryoma_xval/h1h2_descriptive.tsv (+ printed T values)
"""
import os, sys, warnings

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import branchstats as bs  # noqa: E402
import xval_lib as xl  # noqa: E402

OUT = os.path.join(ROOT, "results/eukaryoma_xval")


def main():
    orig = xl.load_tree(os.path.join(ROOT, "runs/euk_fast/stages/input/guide_tree.nwk"))
    pr = xl.load_tree(os.path.join(ROOT, "runs/eukaryoma_xval/tree_pruned.nwk"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bl = bs.branch_lengths_from_tree(bs.load_tree(os.path.join(ROOT, "data/guide_tree_euk.nwk")),
                                         os.path.join(ROOT, "data/species_tree_bl_euk.nwk"))
    sp = pd.read_csv(os.path.join(ROOT, "data/species_final.tsv"), sep="\t").set_index("code")
    is_tx = (sp["data_type_final"].astype(str).str.lower() == "transcriptome").to_dict()
    C = pd.read_csv(os.path.join(OUT, "counts_per_branch.tsv"), sep="\t", keep_default_na=False)
    bm = pd.read_csv(os.path.join(OUT, "branch_mapping.tsv"), sep="\t", keep_default_na=False)
    tr = pd.read_csv(os.path.join(ROOT, "data/transitions.tsv"), sep="\t", keep_default_na=False)
    ms = pd.read_csv(os.path.join(ROOT, "results/primary/euk/matched_sets.tsv"), sep="\t")
    ms = ms[ms["hypothesis"] == "H1H2"]
    r1 = pd.read_csv(os.path.join(ROOT, "results/primary/euk/per_transition_H1H2.tsv"), sep="\t").set_index("transition_id")
    pn = {n.name: n for n in pr.traverse()}
    leaves = {n.name: n.get_leaf_names() for n in pr.traverse()}

    def cov(top, bottom):
        nd = pn[bottom]
        # sister tips: other children of the path's top-most pruned node below `top`
        x = nd
        while x.up.name != top:
            x = x.up
        sis = [l for s in x.up.children if s is not x for l in leaves[s.name]]
        desc = leaves[bottom]
        covered = xl.edges_on_path(orig, bottom, top)
        L = sum(bl.get(c, np.nan) for _, c in covered)
        ftx = np.mean([is_tx.get(s, False) for s in desc])
        return dict(log_bl=np.log(L) if L > 0 else np.nan, log_n_desc=np.log(len(desc)),
                    log_n_sister=np.log(len(sis)) if sis else np.nan, frac_tx_desc=ftx,
                    n_desc=len(desc), n_sister=len(sis), bl=L)

    rows = []
    br = C[C["method"] == "ancestra"].drop_duplicates(["role", "branch_id"])
    for _, b in br.iterrows():
        d = cov(b["top"], b["bottom"])
        rows.append(dict(key=f"{b['role']}:{b['branch_id']}", branch_id=b["branch_id"], role=b["role"],
                         status=b["status"], top=b["top"], bottom=b["bottom"],
                         bottom_is_leaf=b["bottom"] not in leaves or len(leaves[b["bottom"]]) == 1, **d))
    bt = pd.DataFrame(rows).set_index("key")
    bt["in_U"] = (bt["role"] == "control") & (bt["status"] == "exact")
    out = []
    for meth in ("ancestra", "eukaryoma"):
        y = C[C["method"] == meth].assign(key=lambda d: d["role"] + ":" + d["branch_id"]).set_index("key")["gains_supported"]
        y = y.reindex(bt.index).astype(float)
        sub = bt[(bt["role"] != "other")].copy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rs = bs.count_residual_statistic(y.reindex(sub.index), sub)
        s = rs["s"]
        for tset in ("primary", "secondary"):
            for tid in tr.loc[(tr["kingdom"] == "euk") & (tr["analysis_set"] == tset), "transition_id"]:
                k = f"stem:{tid}"
                if k not in bt.index:
                    out.append(dict(method=meth, set=tset, transition_id=tid, status="missing"))
                    continue
                ctrl = [f"control:{c}" for c in ms.loc[(ms["set"] == tset) & (ms["transition_id"] == tid), "control"]]
                ctrl = [c for c in ctrl if c in bt.index and bt.loc[c, "in_U"]]
                ok = not bt.loc[k, "bottom_is_leaf"]
                u = bs.u_score(s[k], s[ctrl].dropna().to_numpy()) if ok and np.isfinite(s[k]) else np.nan
                out.append(dict(method=meth, set=tset, transition_id=tid, status=bt.loc[k, "status"],
                                pruned_stem=f"{bt.loc[k, 'top']}>{bt.loc[k, 'bottom']}",
                                bottom_is_leaf=bt.loc[k, "bottom_is_leaf"], y=y[k], yhat=round(rs.loc[k, "yhat"], 2),
                                s=round(s[k], 4), k_controls=len(ctrl), u=u,
                                round1_seqnovel_y=r1.loc[tid, "y"] if tid in r1.index else np.nan,
                                round1_s=r1.loc[tid, "resid_s"] if tid in r1.index else np.nan,
                                round1_u=r1.loc[tid, "u"] if tid in r1.index else np.nan,
                                glm=f"{rs.attrs['method']} covs={rs.attrs['covariates']} n_pool={rs.attrs['n_pool']}"))
    O = pd.DataFrame(out)
    O.to_csv(os.path.join(OUT, "h1h2_descriptive.tsv"), sep="\t", index=False)
    bt.to_csv(os.path.join(OUT, "h1h2_branch_covariates.tsv"), sep="\t")
    print(O.drop(columns=["glm"]).to_string())
    summ = []
    for (meth, tset), g in O.groupby(["method", "set"]):
        for lab, gg in (("exact stems", g[g["status"] == "exact"]),
                        ("exact + approx (internal bottom)", g[g["status"].isin(["exact", "approx"])])):
            v = gg["u"].dropna()
            r1u = gg.loc[v.index, "round1_u"]
            summ.append(dict(method=meth, set=tset, stems=lab, n=len(v), T_mean_u=round(v.mean(), 3) if len(v) else np.nan,
                             round1_T_same_stems=round(r1u.mean(), 3) if len(v) else np.nan,
                             transitions=",".join(gg.loc[v.index, "transition_id"])))
    S = pd.DataFrame(summ)
    S.to_csv(os.path.join(OUT, "h1h2_descriptive_T.tsv"), sep="\t", index=False)
    print(S.to_string())


if __name__ == "__main__":
    main()
