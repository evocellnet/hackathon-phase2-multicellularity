#!/usr/bin/env python
"""Orthology-free Pfam-count version of H4/H5 (DESIGN §4.4 + §4.6: "Expansions (H4/H5) must agree with the
orthology-free Pfam-count version to be called robust"). Pre-registered robustness analysis; not in the Holm family.

    python scripts/robustness/h4h5_pfam_counts.py [--no-registry]

Per-branch statistic (b = u -> v; path stems P > ... > M use v = M and u = P):
    child  = species below v (M for path stems)
    sister = species below u that are not below the first edge of the stem (single edges: u's other children;
             path stems: P's other children; equals data/transitions.tsv `sister_species` for every stem)
    x_C(S) = mean over species in S of the number of proteins carrying >= 1 class-C domain (frozen
             data/toolkit_pfam.tsv; Pfam-A GA hits from data/annotations/pfam_prereg/euk), x_all(S) = mean proteome size
    s_b    = log((x_C(child) + 0.5) / (x_C(sister) + 0.5)) - log((x_all(child) + 0.5) / (x_all(sister) + 0.5))
C = TF-DBD (no C2H2) for H4, adhesion/ECM U signalling (LRR conditional rule) for H5 -- the same protein-level class
flags that define the HOG classes of the primary analysis (mcell.analysis.protein_flags). No HOGs are used.

Test: the primary matched sets (results/primary/euk/matched_sets.tsv, per hypothesis and set), relabelling
scheme 'permute', T = mean u, one-sided (greater), 10,000 draws, seed 11 (mcell.branchstats.relabel_test).
Effect: s_stem - mean s(N_i) per transition, REML pooled (z and HKSJ CI).
Exploratory sensitivities (same sets/test): clade-balanced (equal-split) species means; genome-only species.

Outputs: results/robustness/h4h5/pfam_counts_{branch_values,per_transition,tests,species_counts}.tsv
"""
from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "prereg"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402

OUT = os.path.join(ROOT, "results", "robustness", "h4h5")
CLASSES = {"H4": "tf_dbd", "H5": "adh_sig"}


def species_counts(p) -> pd.DataFrame:
    """Per species: n proteins (proteome as used in the run), n TF-DBD proteins, n adh/ECM U signalling proteins."""
    genes = pd.read_parquet(os.path.join(p["cache"], "genes.parquet"), columns=["gene", "species"])
    spec = an.load_class_spec(p["toolkit"], p["h12_sets"])
    hits = an.read_pfam_dir(p["pfam"], sorted(genes["species"].astype(str).unique()), pfams=spec.all_pfams())
    if hits.attrs["missing"]:
        raise SystemExit(f"species without Pfam file: {hits.attrs['missing']}")
    hits = hits[hits["protein"].isin(set(genes["gene"]))]
    fl = an.protein_flags(hits, spec)
    sp_of = genes.set_index("gene")["species"].astype(str)
    out = genes.groupby(genes["species"].astype(str), observed=True).size().rename("n_proteins").to_frame()
    for h, c in CLASSES.items():
        prots = fl.index[fl[c].to_numpy()]
        out[f"n_{c}"] = sp_of.reindex(prots).value_counts().reindex(out.index).fillna(0).astype(int)
    # cross-check against the cached gene-level class flags (built from the same files)
    gc = pd.read_parquet(os.path.join(p["cache"], "gene_classes.parquet"), columns=["species", "tf_dbd", "adh_sig"])
    chk = gc.groupby(gc["species"].astype(str), observed=True)[["tf_dbd", "adh_sig"]].sum()
    for c in ("tf_dbd", "adh_sig"):
        if not (chk[c].reindex(out.index).astype(int) == out[f"n_{c}"]).all():
            raise SystemExit(f"class counts differ from runs/analysis cache for {c}")
    return out


def clade_sets(tree, bt):
    """branch -> (child species, sister species) as described in the module doc."""
    tips = {n.name: frozenset(n.get_leaf_names()) for n in tree.traverse()}
    out = {}
    for b in bt.index:
        if bt.at[b, "is_terminal"] or not bt.at[b, "parent"]:
            continue
        edges = bs.path_edges(bt, b)
        first = edges[0]
        u = bt.at[first, "parent"]
        out[b] = (tips[b], tips[u] - tips[first], u)
    return out


def balanced_weights(tree, allowed: set, root_name: str, exclude: frozenset = frozenset()) -> dict:
    """Equal-split weights over the species of ``allowed`` below ``root_name`` (minus ``exclude``)."""
    node = tree & root_name
    w = {}

    def rec(n, mass):
        if n.is_leaf():
            if n.name in allowed and n.name not in exclude:
                w[n.name] = w.get(n.name, 0) + mass
            return
        kids = [c for c in n.children if any(x in allowed and x not in exclude for x in c.get_leaf_names())]
        for c in kids:
            rec(c, mass / len(kids))
    rec(node, 1.0)
    return w


def stat_values(cnt: pd.DataFrame, sets: dict, tree, mode: str, genome: set) -> pd.DataFrame:
    rows = []
    for b, (ch, si, u) in sets.items():
        r = dict(branch=b)
        if mode == "genome_only":
            ch2, si2 = [x for x in ch if x in genome and x in cnt.index], [x for x in si if x in genome and x in cnt.index]
        else:
            ch2, si2 = [x for x in ch if x in cnt.index], [x for x in si if x in cnt.index]
        r.update(n_sp_child=len(ch2), n_sp_sister=len(si2))
        if not ch2 or not si2:
            for h in CLASSES:
                r[f"{h}_s"] = np.nan
            rows.append(r)
            continue
        if mode == "balanced":
            # child: equal split from v; sister: equal split from u with the child subtree removed
            wc = balanced_weights(tree, set(ch2), b)
            ws = balanced_weights(tree, set(si2), u, exclude=frozenset(ch))
            wc = pd.Series(wc)
            ws = pd.Series(ws)
            mc = cnt.loc[wc.index].mul(wc / wc.sum(), axis=0).sum()
            ms = cnt.loc[ws.index].mul(ws / ws.sum(), axis=0).sum()
        else:
            mc = cnt.loc[ch2].mean()
            ms = cnt.loc[si2].mean()
        lr_all = np.log((mc["n_proteins"] + 0.5) / (ms["n_proteins"] + 0.5))
        r.update(mean_proteins_child=mc["n_proteins"], mean_proteins_sister=ms["n_proteins"], lr_all=lr_all)
        for h, c in CLASSES.items():
            lr = np.log((mc[f"n_{c}"] + 0.5) / (ms[f"n_{c}"] + 0.5))
            r[f"{h}_mean_child"] = mc[f"n_{c}"]
            r[f"{h}_mean_sister"] = ms[f"n_{c}"]
            r[f"{h}_lr_class"] = lr
            r[f"{h}_s"] = lr - lr_all
        rows.append(r)
    return pd.DataFrame(rows).set_index("branch")


def run_test(values: pd.Series, sets: dict, stems: dict, n_draws: int, seed: int):
    r = bs.relabel_test(values, sets, stems=stems, alternative="greater", n_draws=n_draws, seed=seed,
                        scheme="permute")
    eff = bs.effect_table(values, sets, stems=stems, scale="difference")
    eff["u"] = r.u.reindex(eff.index)
    ok = eff.dropna(subset=["est"])
    pooled = bs.random_effects(ok["est"], ok["var"])
    hk = bs.random_effects(ok["est"], ok["var"], ci_method="hksj")
    return r, eff, pooled, hk


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-draws", type=int, default=bs.N_DRAWS)
    ap.add_argument("--seed", type=int, default=bs.SEED)
    ap.add_argument("--registry", default=os.path.join(ROOT, "results", "analysis_registry.tsv"))
    ap.add_argument("--no-registry", action="store_true")
    a = ap.parse_args(argv)
    import check_prereg
    if not check_prereg.verify(strict=False):
        sys.exit("PREREG hash mismatch")
    os.makedirs(OUT, exist_ok=True)
    p = an.run_paths("euk")
    tr, traits, tree, bt = an.load_setup("euk", p)
    cnt = species_counts(p)
    cnt.to_csv(os.path.join(OUT, "pfam_counts_species_counts.tsv"), sep="\t")
    genome = set(traits.loc[(traits["kingdom"] == "euk") & (traits["data_type"] == "genome"), "code"])
    sets_cs = clade_sets(tree, bt)
    # stems' sister sets must equal the frozen sister_species
    ms = pd.read_csv(os.path.join(ROOT, "results", "primary", "euk", "matched_sets.tsv"), sep="\t")
    for _, r in tr[tr["analysis_set"].isin(["primary", "secondary"])].iterrows():
        stem = ms.loc[ms["transition_id"] == r["transition_id"], "stem"]
        if stem.empty:
            continue
        b = stem.iloc[0]
        frozen = set(str(r["sister_species"]).split(","))
        if set(sets_cs[b][1]) != frozen or set(sets_cs[b][0]) != set(str(r["mc_species"]).split(",")):
            sys.exit(f"{r['transition_id']}: child/sister sets differ from data/transitions.tsv")
    modes = {"species_mean": "robustness_prereg (DESIGN 4.4/4.6 orthology-free Pfam-count H4/H5); main definition",
             "balanced": "sensitivity of the robustness check: clade-balanced (equal-split) species means",
             "genome_only": "sensitivity of the robustness check: genome-only species"}
    V = {m: stat_values(cnt, sets_cs, tree, m, genome) for m in modes}
    bv = bt[["parent", "depth", "n_desc", "n_sister", "branch_class", "transition_id", "role", "stem_path",
             "in_U"]].join(V["species_mean"], how="inner")
    for m in ("balanced", "genome_only"):
        bv = bv.join(V[m][["H4_s", "H5_s"]].add_suffix(f"_{m}"))
    bv.to_csv(os.path.join(OUT, "pfam_counts_branch_values.tsv"), sep="\t", float_format="%.6g")
    prim = pd.read_csv(os.path.join(ROOT, "results", "primary", "euk", "primary_tests.tsv"), sep="\t")
    prim = prim[prim["analysis"] == "primary"].set_index(["hypothesis", "set"])
    inputs = [os.path.join(ROOT, "data", "toolkit_pfam.tsv"), os.path.join(ROOT, "data", "transitions.tsv"),
              os.path.join(ROOT, "results", "primary", "euk", "matched_sets.tsv"),
              os.path.join(OUT, "pfam_counts_species_counts.tsv")]
    rows, pts = [], []
    for h in ("H4", "H5"):
        for setname in ("primary", "secondary"):
            g = ms[(ms["hypothesis"] == h) & (ms["set"] == setname)]
            sets = {t: list(x.sort_values("rank")["control"]) for t, x in g.groupby("transition_id")}
            stems = {t: x["stem"].iloc[0] for t, x in g.groupby("transition_id")}
            for m, label in modes.items():
                vals = V[m][f"{h}_s"]
                r, eff, pooled, hk = run_test(vals, sets, stems, a.n_draws, a.seed)
                hog_p = prim.at[(h, setname), "p"]
                hog_est = prim.at[(h, setname), "pooled_est"]
                agree = ((r.p < 0.05) == (hog_p < 0.05)) and (np.sign(pooled["mu"]) == np.sign(hog_est))
                tier = "secondary" if setname == "secondary" else "exploratory"
                rows.append(dict(hypothesis=h, set=setname, variant=m, n=r.n, k=int(np.median(list(r.k.values()))),
                                 T_obs=r.T_obs, p=r.p, p_mcse=r.p_mcse, pooled_est=pooled["mu"],
                                 pooled_ci_low=pooled["ci_low"], pooled_ci_high=pooled["ci_high"],
                                 pooled_tau2=pooled["tau2"], hksj_ci_low=hk["ci_low"], hksj_ci_high=hk["ci_high"],
                                 hog_T_obs=prim.at[(h, setname), "T_obs"], hog_p=hog_p, hog_pooled_est=hog_est,
                                 agree_with_hog=bool(agree), dropped=json.dumps(r.dropped), tier=tier, label=label))
                e = eff.copy()
                e.insert(0, "variant", m)
                e.insert(0, "set", setname)
                e.insert(0, "hypothesis", h)
                pts.append(e.reset_index())
                if not a.no_registry:
                    bs.append_registry(a.registry, dict(
                        hypothesis=h,
                        statistic=("orthology-free Pfam-count: s_b = log((mean class-protein count child+.5)/(sister+.5)) "
                                   "- same log ratio for all proteins; class " + ("TF-DBD (no C2H2)" if h == "H4" else
                                                                                 "adhesion/ECM U signalling")
                                   + f"; T = mean u; one-sided; variant={m}"),
                        n=r.n, k=int(np.median(list(r.k.values()))), estimate=pooled["mu"], ci_low=pooled["ci_low"],
                        ci_high=pooled["ci_high"], p=r.p, p_mcse=r.p_mcse, holm_p=np.nan, tier=tier,
                        input_md5s=inputs, seed=a.seed, n_draws=a.n_draws, scheme="permute",
                        notes=(f"label={label}; set={setname}; T_obs={r.T_obs:.4f}; matched sets = primary {h} sets; "
                               f"HOG-based {h}: p={hog_p:.4g}, pooled={hog_est:.3g}; agreement(same alpha=0.05 "
                               f"conclusion and pooled sign)={agree}; not in Holm family; "
                               "script=scripts/robustness/h4h5_pfam_counts.py")))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "pfam_counts_tests.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.concat(pts, ignore_index=True).to_csv(os.path.join(OUT, "pfam_counts_per_transition.tsv"), sep="\t",
                                              index=False, float_format="%.6g")
    with pd.option_context("display.width", 250):
        print(T.drop(columns=["dropped", "label"]).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
