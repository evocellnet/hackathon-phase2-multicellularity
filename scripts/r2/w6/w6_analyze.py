#!/usr/bin/env python
"""W6 (exploratory) analysis: 3Di branch length per AA branch length, multicellular vs unicellular branches.

Input: results/r2/w6/<R>/branch_lengths.tsv.gz (w6_prep.py fit), panel branch tables (results/r2/panels/<R>/branch_table.tsv),
transitions (data/dense/<R>/transitions.tsv).

Per gene g and branch b (only pruned edges that equal one guide edge; root-child edges excluded; AA length >= 0.005):
    y = log(L3Di / LAA)
Saturation / non-linearity control: within each stratum (tips, internal branches) the AA lengths of the reference branches
(unicellular: both ends 'uni', not an event branch) are cut into deciles; residual e = y - median(y of reference branches
in the same decile). Sensitivity: residual from a linear fit of log L3Di on log LAA among reference branches.
Branch value = median of e over genes; 95 % CIs by bootstrapping genes (2,000, seed 11).
(a) stems: origin / origin_complex / nested event branches vs matched controls (internal stems -> panel CONTROL branches;
    tip stems -> unicellular non-event tips). d = stem - median(controls); one-sided permutation p = (1+#ctrl>=stem)/(n+1).
    Pooled: sign test across events (one-sided, d > 0), Stouffer on p.
(b) crowns: per gene, mean e on branches inside multicellular crowns minus mean e on unicellular branches, within stratum,
    averaged over strata present; panel effect = median over genes; one-sided Wilcoxon signed-rank over genes;
    branch-level Mann-Whitney (branches as units) as the phylogenetically more conservative version.
Nulls: the same pipeline with (i) shuffled 3Di (column-wise permutation among taxa) and (ii) AA split-half
(ratio of two random column halves of the same AA alignment, both LG+G4). Sensitivity: 3Di under Poisson+FO+G4.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = os.path.join(ROOT, "results", "r2", "w6")
SEED = 11
NB = 2000
MIN_AA = 0.005
EVENT_KINDS = {"origin", "origin_complex", "nested"}
PIPES = {"real": ("bl_aa", "bl_di"), "shuffled3di": ("bl_aa", "bl_shuf"), "aa_splithalf": ("bl_aaH1", "bl_aaH2"),
         "real_poisson3di": ("bl_aa", "bl_dipois"), "sim_constratio": ("bl_aa", "bl_sim")}
PSEUDO = 0.5   # pseudo-count (substitutions) added to both expected counts L * n_cols


def classes(R):
    bt = pd.read_csv(os.path.join(ROOT, "results", "r2", "panels", R, "branch_table.tsv"), sep="\t").set_index("branch_child")
    tr = pd.read_csv(os.path.join(ROOT, "data", "dense", R.split("_")[0], "transitions.tsv"), sep="\t")
    kind = dict(zip(tr["event_id"], tr["kind"]))
    sp = pd.read_csv(os.path.join(ROOT, "data", "dense", R.split("_")[0], "species.tsv"), sep="\t")
    tx = set(sp.loc[~sp["data_type"].astype(str).str.lower().eq("genome"), "code"])
    rows = []
    for b, r in bt.iterrows():
        ev = [e for e in str(r["events"]).split(";") if e and e != "nan"]
        stem = [e for e in ev if kind.get(e) in EVENT_KINDS]
        rev = [e for e in ev if kind.get(e) == "reversal"]
        if r["root_child"]:
            c = "root_child"
        elif stem:
            c = "stem"
        elif rev:
            c = "reversal_stem"
        elif r["inside_crown"]:
            c = "crown"
        elif r["state_child"] == "uni" and r["state_parent"] == "uni" and not r["is_event"]:
            c = "uni"
        else:
            c = "other"
        rows.append(dict(node=b, cls=c, event=";".join(stem or rev), is_leaf=bool(r["is_leaf"]), control=bool(r["control"]),
                         depth=int(r["depth"]), transcriptome_tip=bool(r["is_leaf"] and b in tx)))
    return pd.DataFrame(rows).set_index("node")


def residuals(df, xa, xd):
    """Add y and residual e (decile-bin primary; linear sensitivity) per stratum."""
    df = df.copy()
    df["x"] = np.log(df[xa])
    # log ratio of expected substitution counts with a 0.5 pseudo-count (short branches: 3Di length often ~0)
    df["y"] = np.log((df[xd] * df["n_cols"] + PSEUDO) / (df[xa] * df["n_cols"] + PSEUDO))
    df["e"] = np.nan
    df["e_lin"] = np.nan
    for leaf, sub in df.groupby("is_leaf"):
        ref = sub[sub["cls"] == "uni"]
        edges = np.quantile(ref["x"], np.linspace(0, 1, 11))[1:-1]
        bref = np.digitize(ref["x"], edges)
        med = pd.Series(ref["y"].to_numpy()).groupby(bref).median()
        b = np.digitize(sub["x"], edges)
        df.loc[sub.index, "e"] = sub["y"].to_numpy() - med.reindex(b).to_numpy()
        sl, ic = np.polyfit(ref["x"], ref["y"], 1)
        df.loc[sub.index, "e_lin"] = sub["y"] - (ic + sl * sub["x"])
    return df


def branch_medians(df, col, genes=None):
    d = df if genes is None else df.merge(pd.DataFrame({"gene": genes}), on="gene")
    return d.groupby("node")[col].median()


def boot(df, col, fn, nb=NB, seed=SEED):
    """fn(df) -> scalar or dict; bootstrap over genes (resampled with replacement, duplicates kept)."""
    rng = np.random.default_rng(seed)
    genes = df["gene"].unique()
    by = {g: d for g, d in df.groupby("gene")}
    out = []
    for _ in range(nb):
        s = rng.choice(genes, len(genes))
        d = pd.concat([by[g].assign(gene=f"{g}#{i}") for i, g in enumerate(s)])
        out.append(fn(d))
    return out


def stem_test(df, cl, col):
    """(a) per event stem vs matched controls."""
    bm = branch_medians(df, col)
    ctrl_int = [n for n in cl.index[cl["control"]] if n in bm.index]
    ctrl_tip = [n for n in cl.index[(cl["cls"] == "uni") & cl["is_leaf"]] if n in bm.index]
    rows = []
    for n in cl.index[cl["cls"] == "stem"]:
        if n not in bm.index:
            continue
        cset = ctrl_tip if cl.loc[n, "is_leaf"] else ctrl_int
        if len(cset) < 2:
            continue
        c = bm[cset].to_numpy()
        s = bm[n]
        rows.append(dict(event=cl.loc[n, "event"], branch=n, tip=cl.loc[n, "is_leaf"], n_genes=int((df["node"] == n).sum()),
                         stem=s, ctrl_median=float(np.median(c)), d=s - float(np.median(c)), n_ctrl=len(c),
                         p_one_sided=(1 + int((c >= s).sum())) / (len(c) + 1)))
    return pd.DataFrame(rows), ctrl_int, ctrl_tip


def crown_gene_effect(df, col):
    """(b) per-gene crown-minus-uni contrast averaged over strata."""
    d = df[df["cls"].isin(["crown", "uni"])]
    m = d.groupby(["gene", "is_leaf", "cls"])[col].mean().unstack("cls")
    m = m.dropna()
    diff = (m["crown"] - m["uni"]).groupby(level="gene").mean()
    return diff


def analyse_panel(R, pipe, col="e", drop_tx=False):
    bl = pd.read_csv(os.path.join(OUT, R, "branch_lengths.tsv.gz"), sep="\t")
    cl = classes(R)
    xa, xd = PIPES[pipe]
    sim = os.path.join(OUT, R, "sim_lengths.tsv.gz")
    if os.path.exists(sim):
        bl = bl.merge(pd.read_csv(sim, sep="\t"), on=["gene", "node"], how="left")
    nc = pd.read_csv(os.path.join(OUT, R, "work", "alignments.tsv"), sep="\t")[["gene", "n_cols"]]
    bl = bl.merge(nc, on="gene")
    bl = bl[bl["single_edge"] & ~bl["root_child"]]
    n0 = len(bl)
    # >= 1 expected AA substitution on the branch in this gene
    bl = bl[(bl[xa] >= MIN_AA) & (bl[xa] * bl["n_cols"] >= 1)].join(cl, on="node")
    if drop_tx:
        bl = bl[~bl["transcriptome_tip"]]
    df = residuals(bl, xa, xd)
    st, ci_, ct_ = stem_test(df, cl, col)
    ge = crown_gene_effect(df, col)
    bm = branch_medians(df, col)
    crown_b = bm[[n for n in cl.index[cl["cls"] == "crown"] if n in bm.index]]
    uni_b = bm[[n for n in cl.index[cl["cls"] == "uni"] if n in bm.index]]
    res = dict(panel=R, pipeline=pipe, residual=col, drop_transcriptome_tips=drop_tx, n_gene_branches=len(df),
               n_dropped_short=n0 - len(bl), n_genes=int(df["gene"].nunique()),
               spearman_log3di_logaa=float(stats.spearmanr(df[xa], df[xd])[0]),
               spearman_branch_medians=float(stats.spearmanr(df.groupby("node")[xa].median(), df.groupby("node")[xd].median())[0]),
               b_effect=float(ge.median()), b_n_genes=int(len(ge)),
               b_wilcoxon_p_one_sided=float(stats.wilcoxon(ge, alternative="greater")[1]) if len(ge) > 5 else np.nan,
               b_frac_genes_positive=float((ge > 0).mean()),
               b_branch_level_diff=float(crown_b.median() - uni_b.median()), b_n_crown_branches=len(crown_b),
               b_n_uni_branches=len(uni_b),
               b_branch_mw_p_one_sided=float(stats.mannwhitneyu(crown_b, uni_b, alternative="greater")[1]))
    return df, st, res


def branch_pooled(R, pipe):
    """Sensitivity: per-branch ratio of summed lengths over genes, z_b = log(sum L3Di / sum LAA); AA-length control by a
    linear fit z ~ log(mean LAA) + leaf on the unicellular reference branches; residual compared as in (a)/(b)."""
    bl = pd.read_csv(os.path.join(OUT, R, "branch_lengths.tsv.gz"), sep="\t")
    sim = os.path.join(OUT, R, "sim_lengths.tsv.gz")
    if os.path.exists(sim):
        bl = bl.merge(pd.read_csv(sim, sep="\t"), on=["gene", "node"], how="left")
    xa, xd = PIPES[pipe]
    bl = bl[bl["single_edge"] & ~bl["root_child"]]
    cl = classes_cache[R]
    g = bl.groupby("node").agg(sa=(xa, "sum"), sd=(xd, "sum"), ma=(xa, "mean")).join(cl)
    g["z"] = np.log(g["sd"] / g["sa"])
    ref = g[g["cls"] == "uni"]
    X = np.c_[np.ones(len(ref)), np.log(ref["ma"]), ref["is_leaf"].astype(float)]
    beta = np.linalg.lstsq(X, ref["z"], rcond=None)[0]
    g["r"] = g["z"] - np.c_[np.ones(len(g)), np.log(g["ma"]), g["is_leaf"].astype(float)] @ beta
    rows = []
    ci_ = g.index[g["control"]]
    ct_ = g.index[(g["cls"] == "uni") & g["is_leaf"]]
    for n in g.index[g["cls"] == "stem"]:
        c = g.loc[ct_ if g.loc[n, "is_leaf"] else ci_, "r"].to_numpy()
        rows.append(dict(panel=R, pipeline=pipe, event=g.loc[n, "event"], branch=n, d=g.loc[n, "r"] - np.median(c),
                         p_one_sided=(1 + int((c >= g.loc[n, "r"]).sum())) / (len(c) + 1)))
    cr, un = g.loc[g["cls"] == "crown", "r"], g.loc[g["cls"] == "uni", "r"]
    crown = dict(panel=R, pipeline=pipe, crown_minus_uni=float(cr.median() - un.median()),
                 mw_p_one_sided=float(stats.mannwhitneyu(cr, un, alternative="greater")[1]), n_crown=len(cr), n_uni=len(un))
    return pd.DataFrame(rows), crown


def run_all(panels):
    os.makedirs(OUT, exist_ok=True)
    summ, stems = [], []
    gene_eff = []
    for R in panels:
        for pipe in PIPES:
            for col, dtx in (("e", False), ("e_lin", False), ("e", True)):
                if pipe != "real" and (col != "e" or dtx):
                    continue
                df, st, res = analyse_panel(R, pipe, col, dtx)
                tag = "primary" if (col == "e" and not dtx) else ("linear_residual" if col == "e_lin" else "no_transcriptome_tips")
                res["analysis"] = tag
                # bootstrap CIs over genes: (b) gene-median effect; (a) per-event d
                ge = crown_gene_effect(df, col)
                rng = np.random.default_rng(SEED)
                bs = [np.median(rng.choice(ge.to_numpy(), len(ge))) for _ in range(NB)]
                res["b_effect_lo"], res["b_effect_hi"] = np.quantile(bs, [0.025, 0.975])
                if tag == "primary":
                    gene_eff.append(pd.DataFrame({"panel": R, "pipeline": pipe, "gene": ge.index, "crown_minus_uni": ge.to_numpy()}))
                if len(st):
                    names = list(st["branch"])

                    def fn(d):
                        s2, _, _ = stem_test(d, classes_cache[R], col)
                        return s2.set_index("branch")["d"].reindex(names).to_numpy()
                    B = np.array(boot(df, col, fn, nb=500 if pipe == "real" and tag == "primary" else 200))
                    st["d_lo"], st["d_hi"] = np.nanquantile(B, 0.025, axis=0), np.nanquantile(B, 0.975, axis=0)
                    st["panel"], st["pipeline"], st["analysis"] = R, pipe, tag
                    stems.append(st)
                summ.append(res)
                print(R, pipe, tag, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()
                                     if k.startswith("b_") or k.startswith("spear")}, file=sys.stderr, flush=True)
            if pipe == "real":
                df, _, _ = analyse_panel(R, "real")
                df.to_csv(os.path.join(OUT, R, "gene_branch_residuals.tsv.gz"), sep="\t", index=False)
    bp_st, bp_cr = [], []
    for R in panels:
        for pipe in PIPES:
            x, y = branch_pooled(R, pipe)
            bp_st.append(x)
            bp_cr.append(y)
    bp_st = pd.concat(bp_st, ignore_index=True)
    bp_st.to_csv(os.path.join(OUT, "branch_pooled_stem_tests.tsv"), sep="\t", index=False)
    bpc = pd.DataFrame(bp_cr)
    for pipe, gg in bp_st.groupby("pipeline"):
        k, npos = len(gg), int((gg["d"] > 0).sum())
        bpc = pd.concat([bpc, pd.DataFrame([dict(panel="pooled_stems", pipeline=pipe, crown_minus_uni=float(gg["d"].median()),
                                                 mw_p_one_sided=float(stats.binomtest(npos, k, 0.5, alternative="greater").pvalue),
                                                 n_crown=npos, n_uni=k)])])
    bpc.to_csv(os.path.join(OUT, "branch_pooled_tests.tsv"), sep="\t", index=False)
    S = pd.DataFrame(summ)
    ST = pd.concat(stems, ignore_index=True)
    GE = pd.concat(gene_eff, ignore_index=True)
    # pooled (a): sign test over events and Stouffer, per pipeline/analysis
    pooled = []
    for (pipe, tag), g in ST.groupby(["pipeline", "analysis"]):
        k, npos = len(g), int((g["d"] > 0).sum())
        z = stats.norm.isf(g["p_one_sided"].clip(upper=1 - 1e-9)).sum() / np.sqrt(k)
        rng = np.random.default_rng(SEED)
        dv = g["d"].to_numpy()
        bs = [np.median(rng.choice(dv, k)) for _ in range(NB)]
        pooled.append(dict(test="(a) stems vs matched controls", pipeline=pipe, analysis=tag, k_events=k, n_positive=npos,
                           sign_test_p_one_sided=float(stats.binomtest(npos, k, 0.5, alternative="greater").pvalue),
                           stouffer_p_one_sided=float(stats.norm.sf(z)), median_d=float(np.median(dv)),
                           median_d_lo=float(np.quantile(bs, 0.025)), median_d_hi=float(np.quantile(bs, 0.975)),
                           ratio_scale=float(np.exp(np.median(dv)))))
    for (pipe, tag), g in S.groupby(["pipeline", "analysis"]):
        k, npos = len(g), int((g["b_effect"] > 0).sum())
        pooled.append(dict(test="(b) crown vs unicellular branches", pipeline=pipe, analysis=tag, k_events=k, n_positive=npos,
                           sign_test_p_one_sided=float(stats.binomtest(npos, k, 0.5, alternative="greater").pvalue),
                           median_d=float(g["b_effect"].mean()), ratio_scale=float(np.exp(g["b_effect"].mean()))))
    S.to_csv(os.path.join(OUT, "panel_tests.tsv"), sep="\t", index=False)
    ST.to_csv(os.path.join(OUT, "stem_tests.tsv"), sep="\t", index=False)
    GE.to_csv(os.path.join(OUT, "gene_crown_effects.tsv.gz"), sep="\t", index=False)
    # calibrated contrasts: real minus constant-ratio simulation (same genes, same branches, same binning), which removes
    # estimation artefacts of the pipeline (regression to the mean from binning on the noisy AA length, saturation, gaps)
    cal = []
    w = GE.pivot_table(index=["panel", "gene"], columns="pipeline", values="crown_minus_uni")
    for R, g in w.groupby(level="panel"):
        dv = (g["real"] - g["sim_constratio"]).dropna().to_numpy()
        rng = np.random.default_rng(SEED)
        bs = [np.median(rng.choice(dv, len(dv))) for _ in range(NB)]
        cal.append(dict(test="(b) crown vs unicellular, real minus constant-ratio null (per gene)", panel=R, k=len(dv),
                        effect=float(np.median(dv)), lo=float(np.quantile(bs, 0.025)), hi=float(np.quantile(bs, 0.975)),
                        frac_positive=float((dv > 0).mean()), p_one_sided=float(stats.wilcoxon(dv, alternative="greater")[1])))
    sp = ST[ST["analysis"] == "primary"].pivot_table(index=["panel", "event"], columns="pipeline", values="d")
    dd = (sp["real"] - sp["sim_constratio"]).dropna()
    rng = np.random.default_rng(SEED)
    bs = [np.median(rng.choice(dd.to_numpy(), len(dd))) for _ in range(NB)]
    cal.append(dict(test="(a) stems, d(real) minus d(constant-ratio null), sign test over events", panel="pooled", k=len(dd),
                    effect=float(dd.median()), lo=float(np.quantile(bs, 0.025)), hi=float(np.quantile(bs, 0.975)),
                    frac_positive=float((dd > 0).mean()),
                    p_one_sided=float(stats.binomtest(int((dd > 0).sum()), len(dd), 0.5, alternative="greater").pvalue)))
    for R, g in dd.groupby(level="panel"):
        cal.append(dict(test="(a) stems, d(real) minus d(constant-ratio null), sign test over events", panel=R, k=len(g),
                        effect=float(g.median()), frac_positive=float((g > 0).mean()),
                        p_one_sided=float(stats.binomtest(int((g > 0).sum()), len(g), 0.5, alternative="greater").pvalue)))
    pd.DataFrame(cal).to_csv(os.path.join(OUT, "calibrated_tests.tsv"), sep="\t", index=False)
    pd.DataFrame(pooled).to_csv(os.path.join(OUT, "pooled_tests.tsv"), sep="\t", index=False)
    calibration(panels)


def calibration(panels):
    """Per-species ProstT5 accuracy (W5, vs AFDB-derived 3Di) against the species' tip residual."""
    acc = pd.read_csv(os.path.join(ROOT, "results", "r2", "w5_3di_validation", "per_protein_accuracy.tsv"), sep="\t")
    acc["species"] = acc["protein"].str.split("_").str[0]
    sa = acc.groupby("species").agg(acc=("acc_p", "mean"), n=("acc_p", "size"))
    rows = []
    for R in panels:
        df = pd.read_csv(os.path.join(OUT, R, "gene_branch_residuals.tsv.gz"), sep="\t")
        tips = df[df["is_leaf"]].groupby("node")["e"].median()
        for s, v in tips.items():
            if s in sa.index and sa.loc[s, "n"] >= 20:
                rows.append(dict(panel=R, species=s, tip_residual=v, prostt5_accuracy=sa.loc[s, "acc"], n_acc=int(sa.loc[s, "n"])))
    c = pd.DataFrame(rows)
    c.to_csv(os.path.join(OUT, "calibration_species_accuracy.tsv"), sep="\t", index=False)
    if len(c) >= 4:
        r, p = stats.spearmanr(c["prostt5_accuracy"], c["tip_residual"])
        json.dump(dict(n_species=len(c), spearman=float(r), p=float(p)), open(os.path.join(OUT, "calibration_summary.json"), "w"), indent=1)


classes_cache = {}

if __name__ == "__main__":
    panels = sys.argv[1:] or ["R1", "R2"]
    for R in panels:
        classes_cache[R] = classes(R)
    run_all(panels)
