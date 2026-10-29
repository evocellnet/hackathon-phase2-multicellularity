#!/usr/bin/env python
"""Errata item 4: the H9 rate test over ALL control orientations, by exact enumeration.

What 'orientation' means. A transition's statistic is r = median_g log(L(P->M tips)/L(P->S tips)),
signed so that r > 0 means the multicellular side evolved faster. A control node has no
multicellular side, so its r is computed for its two child clades in an arbitrary fixed order
(mcell.rates.control_triples: alphabetical) and then multiplied by a random sign
(branchstats.h9_values, seed 11). With 10 controls there are 2^10 = 1,024 orientations; the
pre-registered analysis used one of them (seed 11), and the p-value moved from 0.19 to 0.61 across
the five pre-registered seeds.

The test (branchstats.relabel_test, scheme 'permute', two-sided). u_i = (#controls below r_i +
0.5 #ties + 0.5)/(k+1) for every stem i against all k = 10 controls (every stem's matched set is
the whole H9 pool, k = |pool| = 10); T = mean_i u_i; null = permute the 17 values over the 7 stem
and 10 control positions. Because every stem is compared with the same 10 controls, T depends only
on WHICH 7 values sit in stem positions, so the permutation null has exactly C(17,7) = 19,448
equally likely values. We enumerate them all (no Monte Carlo error) for every orientation.

Outputs (results/round1_errata/h9/):
  per_orientation.tsv      one row per orientation (signs, T_obs, exact p) for the primary analysis
  summary.tsv              per analysis: p distribution over orientations and orientation-free tests
  r_values.tsv             r recomputed from the gene trees (checked against results/h9 tables)
  seed_check.tsv           exact p at the 5 pre-registered orientation seeds vs the MC p in results/h9
Orientation-free summaries:
  (a) mean p over orientations (2 x mean p is a valid p-value under any dependence, Ruschendorf 1982);
  (b) T_bar = mean over orientations of T; its null = the same average applied to every one of the
      19,448 relabellings (T_bar(pi) = mean_o T(o, pi)), i.e. the test statistic averaged over the
      pre-registered randomisation with its permutation null computed the same way. Exact.
Analyses: primary; S2 genome-only; S3 unicellular tips only; Dictyostelia without the stem;
Dictyostelia with S = Physarum only and S = Ceratiomyxa only (both transcriptomes in our data).
Usage: python scripts/errata/h9_all_orientations.py
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "phylo"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "prereg"))
from mcell import rates, branchstats as bs  # noqa: E402
import h9_rates as H  # noqa: E402  (setup, r_table: identical code path to the round-1 analysis)

OUT = os.path.join(ROOT, "results", "round1_errata", "h9")
SEEDS = [11, 23, 37, 41, 53]


def exact_all_orientations(a: np.ndarray, c: np.ndarray, label: str):
    """a: stem r (n), c: control r in the fixed orientation (k). Returns per-orientation frame and
    summary dict. Integer arithmetic: 2(k+1)n T = sum_i (2 #less + #eq + 1)."""
    a, c = np.asarray(a, float), np.asarray(c, float)
    n, k = len(a), len(c)
    N = n + k
    combos = np.array(list(itertools.combinations(range(N), n)), dtype=np.int16)   # stem positions
    mask = np.zeros((len(combos), N), bool)
    mask[np.arange(len(combos))[:, None], combos] = True
    ctrl = np.array([np.flatnonzero(~m) for m in mask], dtype=np.int16)
    ident = int(np.flatnonzero((combos == np.arange(k, N)).all(1))[0])             # observed labelling
    signs = np.array(list(itertools.product([1.0, -1.0], repeat=k)))               # 2^k orientations
    half = (k + 1) * n                                                             # integer value of T = 1/2
    S = np.empty((len(signs), len(combos)), dtype=np.int32)
    for j, sg in enumerate(signs):
        v = np.concatenate([sg * c, a])
        X, C = v[combos], v[ctrl]                                                  # (M, n), (M, k)
        less = (C[:, None, :] < X[:, :, None]).sum((1, 2))
        eq = (C[:, None, :] == X[:, :, None]).sum((1, 2))
        S[j] = 2 * less + eq + n
    D = np.abs(S - half)                                                           # |T - 1/2| * 2(k+1)n
    d_obs = D[:, ident]
    p = (D >= d_obs[:, None]).mean(1)
    T_obs = S[:, ident] / (2.0 * half)
    Sbar = S.mean(0)
    dbar = np.abs(Sbar - half)
    p_bar = float((dbar >= dbar[ident] - 1e-9).mean())
    per = pd.DataFrame(signs.astype(int), columns=[f"sign_c{i}" for i in range(k)])
    per.insert(0, "analysis", label)
    per["T_obs"], per["p_exact"] = T_obs, p
    summ = dict(analysis=label, n_stems=n, k_controls=k, n_orientations=len(signs), n_relabellings=len(combos),
                p_min=p.min(), p_q05=np.quantile(p, .05), p_q25=np.quantile(p, .25), p_median=np.median(p),
                p_q75=np.quantile(p, .75), p_q95=np.quantile(p, .95), p_max=p.max(), p_mean=p.mean(),
                p_2xmean_valid_bound=min(1.0, 2 * p.mean()), frac_orient_p_le_0_05=(p <= 0.05).mean(),
                frac_orient_p_le_0_10=(p <= 0.10).mean(),
                T_obs_min=T_obs.min(), T_obs_max=T_obs.max(), T_bar_obs=Sbar[ident] / (2.0 * half),
                p_T_bar_exact=p_bar)
    return per, summ, (signs, S, ident, half)


def main():
    warnings.simplefilter("ignore")
    os.makedirs(OUT, exist_ok=True)
    tr, traits, tree, bt, gts = H.setup("euk")
    stems = bs.stem_branches(bt)
    pool = bs.h9_control_branches(tree, bt)
    m9 = bs.match_controls(bt, stems, pool=pool)
    assert all(sorted(m9.sets[t]) == sorted(pool) for t in stems.index), "H9 matched sets != whole pool"
    st_units = rates.stem_triples(tr, tree, sets=["primary"], kingdom="euk")
    ctrl_units = rates.control_triples(tree, pool)
    keys = sorted(pool)                                   # h9_values orders controls by sorted name
    tids = list(stems.index)
    mc = set(traits.loc[traits["mc_state"].isin(["simple", "complex"]), "code"])
    tx = set(traits.loc[traits["data_type"] == "transcriptome", "code"])

    # --- r values, recomputed with the round-1 code path
    tab, summ = H.r_table(gts, tree, st_units + ctrl_units)
    gts_g = {g: t for g, t in ((g, rates.drop_tips(t, tx)) for g, t in gts.items()) if t is not None}
    _, summ_g = H.r_table(gts_g, tree, st_units + ctrl_units)
    _, summ_u = H.r_table(gts, tree, ctrl_units, exclude={c: mc for c in pool})
    dic = [u for u in st_units if u[0] == "T_dictyostelia"][0]
    alt = [("T_dicty_S_Physpol", dic[1], dic[2], dic[3]), ("T_dicty_S_Cerafru", dic[1], dic[2], dic[3])]
    _, summ_d = H.r_table(gts, tree, alt, exclude={"T_dicty_S_Physpol": {"Cerafru"}, "T_dicty_S_Cerafru": {"Physpol"}})

    old_t = pd.read_csv(os.path.join(ROOT, "results/h9/h9_transitions.tsv"), sep="\t").set_index("transition_id")
    old_c = pd.read_csv(os.path.join(ROOT, "results/h9/h9_controls.tsv"), sep="\t").set_index("control")
    rv = []
    for u in tids:
        rv.append(dict(unit=u, kind="stem", r=summ.loc[u, "r"], r_round1=old_t.loc[u, "r"],
                       r_genome_only=summ_g.loc[u, "r"], n_genes=summ.loc[u, "n_genes"]))
    for u in keys:
        rv.append(dict(unit=u, kind="control", r=summ.loc[u, "r"], r_round1=old_c.loc[u, "r"],
                       r_genome_only=summ_g.loc[u, "r"], r_unicellular_tips_only=summ_u.loc[u, "r"],
                       n_genes=summ.loc[u, "n_genes"],
                       n_mc_tips_child1=old_c.loc[u, "n_mc_tips1"], n_mc_tips_child2=old_c.loc[u, "n_mc_tips2"]))
    for u in summ_d.index:
        rv.append(dict(unit=u, kind="stem_variant", r=summ_d.loc[u, "r"], n_genes=summ_d.loc[u, "n_genes"],
                       ci_lo=summ_d.loc[u, "ci_lo"], ci_hi=summ_d.loc[u, "ci_hi"]))
    rv = pd.DataFrame(rv)
    rv.to_csv(os.path.join(OUT, "r_values.tsv"), sep="\t", index=False, float_format="%.6g")
    chk = rv.dropna(subset=["r_round1"])
    maxdiff = float((chk.r - chk.r_round1).abs().max())
    assert maxdiff < 1e-4, f"recomputed r differ from round-1 tables by {maxdiff}"

    a = np.array([summ.loc[t, "r"] for t in tids])
    c = np.array([summ.loc[x, "r"] for x in keys])

    # --- seed check: exact p at the pre-registered orientations vs the MC p of round 1
    per, s_prim, (signs, S, ident, half) = exact_all_orientations(a, c, "primary")
    stab = pd.read_csv(os.path.join(ROOT, "results/h9/h9_seed_stability.tsv"), sep="\t")
    stab = stab[stab.varied == "orientation_seed"].set_index("orientation_seed")
    rows = []
    for sd in SEEDS:
        sg = np.random.default_rng(sd).choice([-1.0, 1.0], size=len(keys))
        j = int(np.flatnonzero((signs == sg).all(1))[0])
        res = bs.relabel_test(bs.h9_values(dict(zip(tids, a)), dict(stems), dict(zip(keys, c)), orient=True, seed=sd),
                              m9, alternative="two-sided", n_draws=10_000, seed=11)
        rows.append(dict(orientation_seed=sd, orientation_index=j, T_obs_exact=per.loc[j, "T_obs"],
                         T_obs_relabel_test=res.T_obs, p_exact=per.loc[j, "p_exact"], p_mc_now=res.p,
                         p_mc_round1=stab.loc[sd, "p"] if sd in stab.index else np.nan))
    seedchk = pd.DataFrame(rows)
    seedchk.to_csv(os.path.join(OUT, "seed_check.tsv"), sep="\t", index=False, float_format="%.5g")
    assert np.allclose(seedchk.T_obs_exact, seedchk.T_obs_relabel_test)

    # pooled effect (r_i - mean oriented control r, REML random effects) for every orientation
    pe = []
    for j, sg in enumerate(signs):
        vals = bs.h9_values(dict(zip(tids, a)), dict(stems), dict(zip(keys, sg * c)), orient=False)
        eff = bs.effect_table(vals, m9, scale="difference").dropna(subset=["est"])
        re = bs.random_effects(eff["est"], eff["var"])
        oc = sg * c                     # shared-control CI as in results/h9/pooled_ci_shared_controls.tsv
        d_sh, se_sh = a.mean() - oc.mean(), np.sqrt(np.var(oc, ddof=1) * (1 / len(oc) + 1 / len(a)))
        pe.append(dict(orientation_index=j, pooled_diff=re["mu"], ci_low=re["ci_low"], ci_high=re["ci_high"],
                       shared_rr=np.exp(d_sh), shared_rr_ci_low=np.exp(d_sh - 1.959964 * se_sh),
                       shared_rr_ci_high=np.exp(d_sh + 1.959964 * se_sh),
                       rate_ratio=np.exp(re["mu"]), rr_ci_low=np.exp(re["ci_low"]), rr_ci_high=np.exp(re["ci_high"])))
    pe = pd.DataFrame(pe)
    pe.to_csv(os.path.join(OUT, "pooled_effect_per_orientation.tsv"), sep="\t", index=False, float_format="%.5g")
    pe_s = dict(rr_min=pe.rate_ratio.min(), rr_median=pe.rate_ratio.median(), rr_max=pe.rate_ratio.max(),
                rr_ci_low_min=pe.rr_ci_low.min(), rr_ci_high_max=pe.rr_ci_high.max(),
                frac_ci_inside_0_67_1_5=float(((pe.rr_ci_low > 1 / 1.5) & (pe.rr_ci_high < 1.5)).mean()),
                frac_ci_excludes_1=float(((pe.rr_ci_low > 1) | (pe.rr_ci_high < 1)).mean()),
                shared_ci_low_min=pe.shared_rr_ci_low.min(), shared_ci_high_max=pe.shared_rr_ci_high.max(),
                frac_shared_ci_inside_0_67_1_5=float(((pe.shared_rr_ci_low > 1 / 1.5) & (pe.shared_rr_ci_high < 1.5)).mean()),
                seed11=pe.loc[seedchk.loc[0, "orientation_index"]].to_dict(),
                orientation_free_rr=float(np.exp(a.mean())),
                note="orientation-averaged control mean is exactly 0, so the orientation-free pooled difference is "
                     "the mean stem r (equal weights: every stem has the same control set and variance)")
    print(json.dumps(pe_s, indent=1, default=float))

    summaries = [s_prim]
    pers = [per]

    def run(label, a_, c_):
        p_, s_, _ = exact_all_orientations(a_, c_, label)
        summaries.append(s_)
        pers.append(p_)

    # S2 genome-only: stems / controls without a value drop out (as in round 1)
    ag = {t: summ_g.loc[t, "r"] for t in tids if pd.notna(summ_g.loc[t, "r"])}
    cg = [summ_g.loc[x, "r"] for x in keys if pd.notna(summ_g.loc[x, "r"])]
    run("S2_genome_only", np.array(list(ag.values())), np.array(cg))
    cu = [summ_u.loc[x, "r"] for x in keys if pd.notna(summ_u.loc[x, "r"])]
    run("S3_unicellular_tips_only", a, np.array(cu))
    keep = [i for i, t in enumerate(tids) if t != "T_dictyostelia"]
    run("without_Dictyostelia", a[keep], c)
    idd = tids.index("T_dictyostelia")
    for u in summ_d.index:
        a2 = a.copy()
        a2[idd] = summ_d.loc[u, "r"]
        run(f"Dictyostelia_{u.split('_S_')[1]}_only_as_sister", a2, c)
    # controls whose two child clades contain no multicellular species
    pure = [x for x in keys if old_c.loc[x, "n_mc_tips1"] == 0 and old_c.loc[x, "n_mc_tips2"] == 0]

    sm = pd.DataFrame(summaries)
    sm.to_csv(os.path.join(OUT, "summary.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.concat(pers).to_csv(os.path.join(OUT, "per_orientation.tsv"), sep="\t", index=False, float_format="%.5g")
    meta = dict(controls_sorted=keys, stems=tids, n_controls=len(keys),
                n_controls_with_mc_species_in_a_child_clade=len(keys) - len(pure),
                controls_without_mc_species=pure, max_abs_diff_r_vs_round1=maxdiff,
                pooled_effect_over_orientations=pe_s)
    json.dump(meta, open(os.path.join(OUT, "run.json"), "w"), indent=1, default=float)
    pd.set_option("display.width", 250)
    print(seedchk.to_string())
    print(sm.T.to_string())
    print(rv.to_string())
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
