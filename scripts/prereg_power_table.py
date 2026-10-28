#!/usr/bin/env python
"""Minimum attainable p and simulated MDE per primary test (DESIGN v3.1, PREREG).

Real eukaryote guide tree + data/transitions.tsv (analysis_set == 'primary'),
control pool U and matched sets from mcell.branchstats, relabelling scheme
'permute' (D6). For each test a planted multiplicative effect on the stems is
simulated under a stated variance assumption; power = share of replicates with
p <= alpha; MDE = smallest effect with 80 % power (log-linear interpolation).

Output: results/prereg/power_table.tsv (+ power_curves.tsv with every point).
Usage: python scripts/prereg_power_table.py [--n-sim 200] [--n-draws 1000] [--jobs 20]
"""
import argparse
import math
import os
import sys
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mcell import branchstats as bs  # noqa: E402

SEED = 11
G = {}


def setup():
    warnings.simplefilter("ignore")
    tr = pd.read_csv(os.path.join(ROOT, "data/transitions.tsv"), sep="\t")
    traits = pd.read_csv(os.path.join(ROOT, "data/species_traits.tsv"), sep="\t")
    tree = bs.load_tree(os.path.join(ROOT, "data/guide_tree_euk.nwk"))
    blp = os.path.join(ROOT, "data/species_tree_bl_euk.nwk")
    bl = bs.branch_lengths_from_tree(tree, blp if os.path.exists(blp) else None)
    bt = bs.build_branch_table(tree, traits, tr[tr.kingdom == "euk"], branch_lengths=bl)
    stems = bs.stem_branches(bt)
    m = bs.match_controls(bt, stems)
    h9pool = bs.h9_control_branches(tree, bt)
    m9 = bs.match_controls(bt, stems, pool=h9pool)
    exp = pd.read_csv(os.path.join(ROOT, "data/h12_expected_direction_frozen.tsv"), sep="\t")
    exp = dict(zip(exp.transition_id, exp.expected))
    h12_t = [t for t in stems.index if exp.get(t) in ("loss", "no_loss")]
    signs = bs.h12_signs({t: exp[t] for t in h12_t})
    m12 = bs.MatchResult(stems={t: stems[t] for t in h12_t}, sets={t: m.sets[t] for t in h12_t},
                         distances={}, covariates=m.covariates, dropped_covariates={}, smd=m.smd,
                         k=m.k, pool=m.pool)
    G.update(bt=bt, stems=stems, m=m, m9=m9, m12=m12, signs=signs, bl_placeholder=not os.path.exists(blp))
    return G


def _init():
    setup()


def _nb(rng, mu, a):
    r = 1.0 / a
    return rng.negative_binomial(r, r / (r + np.asarray(mu, float)))


def sim_one(args):
    test, param, effect, seed, n_draws = args
    rng = np.random.default_rng(seed)
    bt, stems, m = G["bt"], G["stems"], G["m"]
    S = list(stems.values)
    if test == "H1/H2":
        mu = pd.Series(20.0, index=bt.index)
        mu[S] *= effect
        y = pd.Series(_nb(rng, mu, param), index=bt.index, dtype=float)
        s = bs.count_residual_statistic(y, bt)["s"]
        return bs.relabel_test(s, m, alternative="two-sided", n_draws=n_draws, seed=seed).p
    if test in ("H4", "H5"):
        pi_u = 0.04 if test == "H4" else 0.12
        d = 10 + _nb(rng, np.full(len(bt), 30.0), 0.5)
        lp = math.log(pi_u / (1 - pi_u)) + rng.normal(0, 0.5, len(bt))
        lp = pd.Series(lp, index=bt.index)
        lp[S] += math.log(effect)
        a = rng.binomial(d, 1 / (1 + np.exp(-lp.to_numpy())))
        s = pd.Series(np.log((a + 0.5) / (d + 1 - a - 0.5)) - math.log(pi_u / (1 - pi_u)), index=bt.index)
        mm = bs.match_controls(bt, stems, extra=pd.Series(np.log(d), index=bt.index, name="log_d"))
        return bs.relabel_test(s, mm, alternative="greater", n_draws=n_draws, seed=seed).p
    if test == "H9":
        r = pd.Series(rng.normal(0, param, len(bt)), index=bt.index)
        r[S] += math.log(effect)
        return bs.relabel_test(r, G["m9"], alternative="two-sided", n_draws=n_draws, seed=seed).p
    if test == "H12":
        m12, signs = G["m12"], G["signs"]
        p = pd.Series(rng.beta(0.5, 12, len(bt)), index=bt.index)
        for t, b in m12.stems.items():
            if signs[t] > 0:
                p[b] = min(1.0, p[b] * effect)
        f = pd.Series(rng.binomial(40, p.to_numpy()) / 40.0, index=bt.index)
        return bs.relabel_test(f, m12, alternative="greater", n_draws=n_draws, seed=seed, signs=signs).p
    raise ValueError(test)


SPECS = [
    ("H1/H2", "nb_dispersion", 0.3, (1.25, 1.5, 2, 2.5, 3, 4, 6),
     "s_b = residual log((y+0.5)/(yhat+0.5)) from NB-GLM on U (LOO for controls); y ~ NB(mean 20, alpha) "
     "on every branch, effect multiplies the stem mean; two-sided"),
    ("H1/H2", "nb_dispersion", 1.0, (1.25, 1.5, 2, 2.5, 3, 4, 6, 8),
     "as above with NB alpha = 1.0"),
    ("H4", "logit_sd", 0.5, (1.25, 1.5, 2, 3, 4, 6),
     "d = 10 + NB(mean 30, alpha 0.5); a ~ Bin(d, p_b), logit p_b = logit(0.04) + N(0, 0.5^2); effect = odds "
     "ratio on stems; matching adds log d; one-sided; all U branches assumed to pass d >= 10"),
    ("H5", "logit_sd", 0.5, (1.25, 1.5, 2, 3, 4, 6),
     "as H4 with pi_u = 0.12 (adhesion/ECM + signalling)"),
    ("H9", "sd_log_rate_ratio", 0.15, (1.05, 1.1, 1.2, 1.3, 1.5, 2),
     "r ~ N(0, 0.15^2) on controls (pool = U nodes with two child clades >= 2 tips), stems N(log effect, "
     "0.15^2); effect = rate ratio; two-sided"),
    ("H9", "sd_log_rate_ratio", 0.30, (1.1, 1.2, 1.3, 1.5, 2, 3),
     "as above with sd 0.30"),
    ("H12", "beta_loss_prob", "Beta(0.5,12)", (2, 4, 8, 16, 40),
     "40 module HOGs at u; loss prob p_b ~ Beta(0.5, 12) (mean 0.04); effect multiplies p_b on loss-expected "
     "stems only (frozen directions; T_rhodophyta not_applicable -> excluded); s = w f; one-sided"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sim", type=int, default=200)
    ap.add_argument("--n-draws", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--out", default=os.path.join(ROOT, "results/prereg"))
    a = ap.parse_args()
    setup()
    bt, m = G["bt"], G["m"]
    os.makedirs(a.out, exist_ok=True)
    jobs = []
    for si, (test, pname, pval, effects, _) in enumerate(SPECS):
        for e in (1.0,) + tuple(effects):
            for r in range(a.n_sim):
                jobs.append((si, (test, pval if not isinstance(pval, str) else None, e,
                                  SEED * 1_000_003 + si * 100_003 + int(e * 1000) * 1009 + r, a.n_draws)))
    with Pool(a.jobs, initializer=_init) as P:
        ps = P.map(sim_one, [j[1] for j in jobs], chunksize=8)
    res = pd.DataFrame([dict(spec=j[0], effect=j[1][2], p=p) for j, p in zip(jobs, ps)])
    curves, rows = [], []
    for si, (test, pname, pval, effects, assumption) in enumerate(SPECS):
        sub = res[res.spec == si]
        pw = sub.groupby("effect")["p"].agg(power_0_01=lambda x: np.mean(x <= 0.01),
                                            power_0_05=lambda x: np.mean(x <= 0.05)).reset_index()
        pw = pw.rename(columns={"power_0_01": "power_0.01", "power_0_05": "power_0.05"})
        for _, r in pw.iterrows():
            curves.append(dict(test=test, param=pname, param_value=pval, effect=r["effect"],
                               power_0_01=r["power_0.01"], power_0_05=r["power_0.05"], n_sim=a.n_sim))
        mde = bs.mde(pw[pw.effect >= 1.0])
        mt = G["m9"] if test == "H9" else (G["m12"] if test == "H12" else m)
        n_ctrl = len({c for t in mt.stems for c in mt.sets[t]})
        G_perm = n_ctrl + len(mt.stems)
        rows.append(dict(
            test=test, scheme="permute", n_stems=len(mt.stems), k=mt.k, U_size=int(bt["in_U"].sum()),
            pool_size=len(mt.pool), n_controls_used=n_ctrl,
            min_attainable_p=max(1 / 10_001, 1 / (math.perm(G_perm, len(mt.stems)) + 1)),
            variance_param=pname, variance_value=pval,
            mde_alpha_0_01=mde.get("alpha_0.01"), mde_alpha_0_05=mde.get("alpha_0.05"),
            power_at_1_5x_alpha_0_05=float(pw.loc[np.isclose(pw.effect, 1.5), "power_0.05"].iloc[0])
            if np.isclose(pw.effect, 1.5).any() else np.nan,
            size_at_effect_1_alpha_0_05=float(pw.loc[pw.effect == 1.0, "power_0.05"].iloc[0]),
            n_sim=a.n_sim, n_draws=a.n_draws, seed=SEED,
            branch_lengths="placeholder_1.0" if G["bl_placeholder"] else "species_tree_bl_euk.nwk",
            assumption=assumption, effects_grid=",".join(str(e) for e in effects)))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out, "power_table.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(curves).to_csv(os.path.join(a.out, "power_curves.tsv"), sep="\t", index=False,
                                float_format="%.4g")
    print(out[["test", "variance_value", "n_stems", "k", "pool_size", "min_attainable_p",
               "mde_alpha_0_01", "mde_alpha_0_05", "size_at_effect_1_alpha_0_05"]].to_string())


if __name__ == "__main__":
    main()
