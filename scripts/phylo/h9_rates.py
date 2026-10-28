#!/usr/bin/env python
"""Primary H9 test (DESIGN.md §4.6 H9 + v3.1 amendments; PREREG.md, tag prereg-v1).

Aborts unless scripts/prereg/check_prereg.py verifies every frozen input.

Primary (euk): stems = analysis_set 'primary' (branchstats.stem_branches; incl. the D7 path stems);
controls = branchstats.h9_control_branches (pool U of build_branch_table restricted to nodes with
two child clades of >= 2 tips), matched with branchstats.match_controls(pool=...) exactly as in
scripts/prereg_power_table.py (k = min(20, |pool|)). Per gene r_g = log(L̄(P→M)/L̄(P→S))
(mcell.rates, M6 inclusion rules), r_i = median over genes, 2,000 gene bootstraps (seed 11);
control r in random orientation (branchstats.h9_values, seed 11); null = 'permute', 10,000
draws, seed 11, two-sided.

Sensitivities: (S1) literal v3 'spec' scheme; (S2) genome-only taxa (transcriptome tips removed
from every gene tree, all units recomputed); (S3) control r computed on unicellular tips only
(multicellular tips inside the control child clades removed). Descriptive: secondary nested
transitions (r_i, CI, percentile vs the H9 controls) and prokaryote transitions (r_i, CI).
Every test is appended to results/analysis_registry.tsv.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "prereg"))
from mcell import rates, branchstats as bs  # noqa: E402
from mcell.phylo import read_tree  # noqa: E402

OUT = os.path.join(ROOT, "results", "h9")
REG = os.environ.get("H9_REGISTRY", os.path.join(ROOT, "results", "analysis_registry.tsv"))
SEED, N_DRAWS, N_BOOT = 11, 10_000, 2_000


def verify_prereg():
    import check_prereg
    if not check_prereg.verify():
        sys.exit("PREREG hash mismatch: aborting (scripts/prereg/check_prereg.py)")


def setup(kingdom):
    tr = pd.read_csv(os.path.join(ROOT, "data", "transitions.tsv"), sep="\t")
    traits = pd.read_csv(os.path.join(ROOT, "data", "species_traits.tsv"), sep="\t")
    tree = bs.load_tree(os.path.join(ROOT, "data", f"guide_tree_{kingdom}.nwk"))
    bl = bs.branch_lengths_from_tree(tree, os.path.join(ROOT, "data", f"species_tree_bl_{kingdom}.nwk"))
    bt = bs.build_branch_table(tree, traits, tr[tr.kingdom == kingdom], branch_lengths=bl)
    gts = rates.load_gene_trees(os.path.join(ROOT, "runs", "phylo", kingdom, "genetrees"))
    return tr, traits, tree, bt, gts


def r_table(gts, tree, units, exclude=None):
    tab = rates.h9_gene_table(gts, tree, units, unit_exclude=exclude)
    return tab, rates.h9_summary(tab, n_boot=N_BOOT, seed=SEED)


def run_test(summ, stems, ctrl_nodes, m9, scheme="permute"):
    r_st = {t: summ.loc[t, "r"] for t in stems.index}
    r_c = {c: summ.loc[c, "r"] for c in ctrl_nodes}
    vals = bs.h9_values(r_st, dict(stems), r_c, orient=True, seed=SEED)
    res = bs.relabel_test(vals, m9, alternative="two-sided", n_draws=N_DRAWS, seed=SEED, scheme=scheme)
    eff = bs.effect_table(vals, m9, scale="difference")
    ok = eff.dropna(subset=["est"])
    re = bs.random_effects(ok["est"], ok["var"])
    return vals, res, eff, re


def reg(row, kingdom="euk"):
    base = dict(hypothesis="H9", seed=SEED, n_draws=N_DRAWS,
                input_md5s=[os.path.join(ROOT, p) for p in ("data/transitions.tsv", f"data/guide_tree_{kingdom}.nwk",
                                                            f"data/species_tree_bl_{kingdom}.nwk",
                                                            f"runs/phylo/{kingdom}/genetrees/branch_lengths.tsv")])
    base.update(row)
    bs.append_registry(REG, base)


def main():
    warnings.simplefilter("ignore")
    verify_prereg()
    os.makedirs(OUT, exist_ok=True)
    tr, traits, tree, bt, gts = setup("euk")
    stems = bs.stem_branches(bt)                                  # primary stems -> stem branch (= M)
    prim = sorted(tr.loc[(tr.kingdom == "euk") & (tr.analysis_set == "primary"), "transition_id"])
    assert sorted(stems.index) == prim, (sorted(stems.index), prim)
    pool = bs.h9_control_branches(tree, bt)
    assert all(bt.loc[c, "in_U"] for c in pool)
    m9 = bs.match_controls(bt, stems, pool=pool)
    st_units = rates.stem_triples(tr, tree, sets=["primary"], kingdom="euk")
    assert {u[0]: u[2] for u in st_units} == dict(stems), "stem M nodes differ from branchstats stems"
    ctrl_units = rates.control_triples(tree, pool)
    assert [u[0] for u in ctrl_units] == list(pool) or set(u[0] for u in ctrl_units) == set(pool)
    mc = set(traits.loc[traits["mc_state"].isin(["simple", "complex"]), "code"])
    tips = {n.name: set(n.get_leaf_names()) for n in tree.traverse()}
    ctrl_info = pd.DataFrame([dict(control=v, child1=c1, child2=c2[0], n_tips1=len(tips[c1]), n_tips2=len(tips[c2[0]]),
                                   n_mc_tips1=len(tips[c1] & mc), n_mc_tips2=len(tips[c2[0]] & mc))
                              for v, _, c1, c2 in ctrl_units]).set_index("control")
    n_pure = int(((ctrl_info.n_mc_tips1 == 0) & (ctrl_info.n_mc_tips2 == 0)).sum())

    # ---------------- primary
    tab, summ = r_table(gts, tree, st_units + ctrl_units)
    tab.to_csv(os.path.join(OUT, "h9_gene_table.tsv"), sep="\t", index=False, float_format="%.6g")
    vals, res, eff, re = run_test(summ, stems, pool, m9)
    rs = res.summary()

    # ---------------- sensitivities
    s1 = bs.relabel_test(vals, m9, alternative="two-sided", n_draws=N_DRAWS, seed=SEED, scheme="spec").summary()
    tx = set(traits.loc[traits["data_type"] == "transcriptome", "code"])
    gts_g = {g: t for g, t in ((g, rates.drop_tips(t, tx)) for g, t in gts.items()) if t is not None}
    tab_g, summ_g = r_table(gts_g, tree, st_units + ctrl_units)
    tab_g.to_csv(os.path.join(OUT, "h9_gene_table_genome_only.tsv"), sep="\t", index=False, float_format="%.6g")
    _, res_g, eff_g, re_g = run_test(summ_g, stems, pool, m9)
    tab_u, summ_u = r_table(gts, tree, ctrl_units, exclude={c: mc for c in pool})
    summ_u = pd.concat([summ.loc[list(stems.index)], summ_u])
    _, res_u, eff_u, re_u = run_test(summ_u, stems, pool, m9)

    # ---------------- per-transition table
    rows = []
    for t in stems.index:
        rows.append(dict(set="primary", transition_id=t, P=summ.loc[t, "P"], M=summ.loc[t, "M"], S=summ.loc[t, "S"],
                         n_genes=summ.loc[t, "n_genes"], r=summ.loc[t, "r"], ci_lo=summ.loc[t, "ci_lo"],
                         ci_hi=summ.loc[t, "ci_hi"], rate_ratio=summ.loc[t, "rate_ratio"],
                         frac_genes_pos=summ.loc[t, "frac_pos"], u_percentile=res.u.get(t, np.nan),
                         diff_vs_controls=eff.loc[t, "est"], ctrl_pi_low=eff.loc[t, "pi_low"],
                         ctrl_pi_high=eff.loc[t, "pi_high"],
                         r_genome_only=summ_g.loc[t, "r"], ci_lo_genome_only=summ_g.loc[t, "ci_lo"],
                         ci_hi_genome_only=summ_g.loc[t, "ci_hi"], n_genes_genome_only=summ_g.loc[t, "n_genes"]))
    # secondary nested (descriptive): percentile vs the oriented H9 controls
    sec_units = rates.stem_triples(tr, tree, sets=["secondary"], kingdom="euk")
    tab_s, summ_s = r_table(gts, tree, sec_units)
    cvals = vals.reindex(pool).to_numpy(float)
    for t, P, M, S in sec_units:
        rows.append(dict(set="secondary_nested", transition_id=t, P=P, M=M, S=";".join(S),
                         n_genes=summ_s.loc[t, "n_genes"], r=summ_s.loc[t, "r"], ci_lo=summ_s.loc[t, "ci_lo"],
                         ci_hi=summ_s.loc[t, "ci_hi"], rate_ratio=summ_s.loc[t, "rate_ratio"],
                         frac_genes_pos=summ_s.loc[t, "frac_pos"],
                         u_percentile=bs.u_score(summ_s.loc[t, "r"], cvals),
                         diff_vs_controls=summ_s.loc[t, "r"] - np.nanmean(cvals)))
    # prokaryotes (descriptive only, |U| = 2)
    trp, _, treep, _, gtsp = setup("prok")
    pk_units = rates.stem_triples(trp, treep, sets=["secondary"], kingdom="prok")
    tab_p, summ_p = r_table(gtsp, treep, pk_units)
    tab_p.to_csv(os.path.join(OUT, "h9_gene_table_prok.tsv"), sep="\t", index=False, float_format="%.6g")
    for t, P, M, S in pk_units:
        rows.append(dict(set="prok_descriptive", transition_id=t, P=P, M=M, S=";".join(S),
                         n_genes=summ_p.loc[t, "n_genes"], r=summ_p.loc[t, "r"], ci_lo=summ_p.loc[t, "ci_lo"],
                         ci_hi=summ_p.loc[t, "ci_hi"], rate_ratio=summ_p.loc[t, "rate_ratio"],
                         frac_genes_pos=summ_p.loc[t, "frac_pos"]))
    per = pd.DataFrame(rows)
    per.to_csv(os.path.join(OUT, "h9_transitions.tsv"), sep="\t", index=False, float_format="%.5g")
    ctrl = ctrl_info.join(summ[["n_genes", "r", "ci_lo", "ci_hi"]]).join(
        pd.Series({c: vals[c] / summ.loc[c, "r"] if summ.loc[c, "r"] else np.nan for c in pool}, name="orientation_sign"))
    ctrl["r_oriented"] = vals.reindex(ctrl.index)
    ctrl["r_unicellular_tips_only"] = summ_u.reindex(ctrl.index)["r"]
    ctrl["n_genes_unicellular_tips_only"] = summ_u.reindex(ctrl.index)["n_genes"]
    ctrl["r_genome_only"] = summ_g.reindex(ctrl.index)["r"]
    ctrl.to_csv(os.path.join(OUT, "h9_controls.tsv"), sep="\t", float_format="%.5g")

    tests = [("primary", "permute", rs, re, "primary"),
             ("S1 literal v3 'spec' scheme", "spec", s1, re, "sensitivity"),
             ("S2 genome-only taxa", "permute", res_g.summary(), re_g, "sensitivity"),
             ("S3 control r on unicellular tips only", "permute", res_u.summary(), re_u, "sensitivity")]
    trows = []
    for name, scheme, r_, re_, label in tests:
        trows.append(dict(analysis=name, label=label, scheme=scheme, n=r_["n"], k=r_["k"], T_obs=r_["T_obs"], p=r_["p"],
                          p_mcse=r_["p_mcse"], min_attainable_p=r_["min_attainable_p"],
                          pooled_diff=re_["mu"], pooled_ci_low=re_["ci_low"], pooled_ci_high=re_["ci_high"],
                          pooled_tau2=re_.get("tau2", np.nan)))
        reg(dict(statistic="median per-gene log rate ratio r_i (P->M vs P->S); relabelling T = mean u; "
                           "estimate = REML pooled r_i - mean control r",
                 n=r_["n"], k=r_["k"], estimate=float(re_["mu"]), ci_low=float(re_["ci_low"]),
                 ci_high=float(re_["ci_high"]), p=float(r_["p"]), p_mcse=float(r_["p_mcse"]), holm_p="NA",
                 tier="primary" if label == "primary" else "exploratory", scheme=scheme,
                 notes=f"label={label}; {name}; " + ("NON-INFORMATIVE: all N_i identical (k=|pool|), spec null cannot reach |T_obs-1/2|; " if scheme == "spec" else "") + f"T_obs={r_['T_obs']:.4f}; min_attainable_p={r_['min_attainable_p']:.3g}; "
                       f"two-sided; pool={len(pool)} H9 controls; Holm across 5 primary tests applied by PI"))
    tests_df = pd.DataFrame(trows)
    tests_df.to_csv(os.path.join(OUT, "h9_tests.tsv"), sep="\t", index=False, float_format="%.5g")
    for _, r in per[per.set != "primary"].iterrows():
        reg(dict(statistic="median per-gene log rate ratio r_i (descriptive, no test)", n=1, k="",
                 estimate=float(r["r"]) if pd.notna(r["r"]) else float("nan"),
                 ci_low=float(r["ci_lo"]) if pd.notna(r["ci_lo"]) else float("nan"),
                 ci_high=float(r["ci_hi"]) if pd.notna(r["ci_hi"]) else float("nan"), p="NA", p_mcse="NA",
                 holm_p="NA", tier="secondary" if r["set"] == "secondary_nested" else "exploratory",
                 scheme="none", notes=f"label=descriptive; {r['set']}; {r['transition_id']}; n_genes={r['n_genes']}"),
            kingdom="prok" if r["set"] == "prok_descriptive" else "euk")
    meta = dict(finished=dt.datetime.now().isoformat(timespec="seconds"), pool=list(pool), n_pool=len(pool),
                n_pool_without_mc_tips=n_pure, U_size=int(bt["in_U"].sum()), k=m9.k, n_genes_total=len(gts),
                n_genes_genome_only=len(gts_g), seed=SEED, n_draws=N_DRAWS, n_boot=N_BOOT,
                primary=dict(T_obs=rs["T_obs"], p=rs["p"], p_mcse=rs["p_mcse"], min_attainable_p=rs["min_attainable_p"]))
    json.dump(meta, open(os.path.join(OUT, "h9_run.json"), "w"), indent=1, default=float)
    print(per.to_string())
    print(tests_df.to_string())
    print(ctrl.to_string())
    print(json.dumps(meta, indent=1, default=float))


if __name__ == "__main__":
    main()
