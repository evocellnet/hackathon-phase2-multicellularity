#!/usr/bin/env python
"""Null calibration of the relabelling schemes (DESIGN v3.1, D6).

Part A (real tree): the 7 primary eukaryote stems and their matched sets N_i
from the frozen guide tree; under H0 every branch value is iid N(0, 1); share
of replicates with p <= alpha for scheme 'spec' (v3: draws without replacement)
and 'permute' (v3.1 primary).
Part B (synthetic covariate shift, H1/H2 residual statistic): |U| = 30, 7
stems shifted by 0 / 1 / 2.5 pool SDs on log branch length and log n_desc, NB
counts (alpha 0.3) depending on both; permutation test on the residuals vs the
parametric bootstrap (residual_bootstrap_test).

Part C (real tree, H1/H2 residual statistic): NB(mean 20, alpha 0.3) counts on
every branch, no covariate effect, permutation vs parametric bootstrap.
Bootstrap rows use 99 resamples, so their alpha = 0.01 row has no resolution.

Output: results/validation_branchstats/null_calibration.tsv
"""
import argparse
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
ALPHAS = (0.01, 0.05, 0.10)
G = {}


def _init():
    warnings.simplefilter("ignore")
    tr = pd.read_csv(os.path.join(ROOT, "data/transitions.tsv"), sep="\t")
    traits = pd.read_csv(os.path.join(ROOT, "data/species_traits.tsv"), sep="\t")
    tree = bs.load_tree(os.path.join(ROOT, "data/guide_tree_euk.nwk"))
    bt = bs.build_branch_table(tree, traits, tr[tr.kingdom == "euk"])
    G["bt"] = bt
    G["m"] = bs.match_controls(bt)
    G["ids"] = sorted({c for s in G["m"].sets.values() for c in s}) + list(G["m"].stems.values())


def real_one(args):
    scheme, alt, seed, n_draws = args
    rng = np.random.default_rng(seed)
    v = pd.Series(rng.normal(size=len(G["ids"])), index=G["ids"])
    return bs.relabel_test(v, G["m"], alternative=alt, n_draws=n_draws, seed=seed + 1, scheme=scheme).p


def real_resid_one(args):
    method, seed, n_draws = args
    warnings.simplefilter("ignore")
    rng = np.random.default_rng(seed)
    bt = G["bt"]
    y = pd.Series(rng.negative_binomial(1 / 0.3, (1 / 0.3) / (1 / 0.3 + 20.0), len(bt)), index=bt.index,
                  dtype=float)
    if method == "permute_on_residual":
        return bs.relabel_test(bs.count_residual_statistic(y, bt)["s"], G["m"], n_draws=n_draws,
                               seed=seed + 1).p
    return bs.residual_bootstrap_test(y, bt, G["m"], n_boot=n_draws, seed=seed + 1)["p"]


def shift_data(shift, seed):
    rng = np.random.default_rng(seed)
    n_pool, n_st = 30, 7
    pool = [f"u{i}" for i in range(n_pool)]
    stems = {f"T{i}": f"s{i}" for i in range(n_st)}
    idx = pool + list(stems.values())
    isst = np.r_[np.zeros(n_pool, bool), np.ones(n_st, bool)]
    bt = pd.DataFrame(index=idx)
    bt["log_bl"] = rng.normal(0, 0.4, len(idx)) + isst * shift * 0.4
    bt["log_n_desc"] = rng.normal(1.5, 0.5, len(idx)) + isst * shift * 0.5
    bt["log_n_sister"] = rng.normal(2, 0.5, len(idx))
    bt["frac_tx_desc"] = rng.uniform(0, 0.3, len(idx))
    bt["in_U"] = ~isst
    mu = np.exp(1.5 + 0.8 * bt.log_bl + 0.5 * bt.log_n_desc)
    m = bs.match_controls(bt, stems, covariates=list(bs.RESID_COVARIATES), k=20)
    y = pd.Series(rng.negative_binomial(1 / 0.3, (1 / 0.3) / (1 / 0.3 + mu)), index=idx, dtype=float)
    return bt, m, y


def shift_one(args):
    method, shift, seed, n_draws = args
    warnings.simplefilter("ignore")
    bt, m, y = shift_data(shift, seed)
    if method == "permute_on_residual":
        s = bs.count_residual_statistic(y, bt)["s"]
        return bs.relabel_test(s, m, n_draws=n_draws, seed=seed + 1).p
    if method == "permute_on_raw_counts":
        return bs.relabel_test(y, m, n_draws=n_draws, seed=seed + 1).p
    return bs.residual_bootstrap_test(y, bt, m, n_boot=n_draws, seed=seed + 1)["p"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--shift-reps", type=int, default=400)
    ap.add_argument("--boot-reps", type=int, default=300)
    ap.add_argument("--jobs", type=int, default=20)
    ap.add_argument("--only-boot-shifts", default="",
                    help="comma-separated shifts: run only these bootstrap shift rows and append them")
    ap.add_argument("--n-boot", type=int, default=99)
    a = ap.parse_args()
    if a.only_boot_shifts:
        return boot_only(a)
    out = os.path.join(ROOT, "results/validation_branchstats")
    os.makedirs(out, exist_ok=True)
    rows = []
    with Pool(a.jobs, initializer=_init) as P:
        for scheme in ("spec", "permute"):
            for alt in ("two-sided", "greater"):
                seeds = [SEED * 100_000 + r for r in range(a.reps)]
                ps = np.array(P.map(real_one, [(scheme, alt, s, 2000) for s in seeds], chunksize=10))
                for al in ALPHAS:
                    rows.append(dict(setting="real_euk_7_primary_stems_|U|=30_k=20", statistic="iid_normal",
                                     scheme=scheme, alternative=alt, alpha=al,
                                     rejection_rate=float(np.mean(ps <= al)),
                                     se=float(np.sqrt(al * (1 - al) / len(ps))), n_reps=len(ps),
                                     n_draws=2000, seed=SEED))
                print(scheme, alt, [round(np.mean(ps <= al), 4) for al in ALPHAS], flush=True)
        for method, nrep, nd in (("permute_on_residual", a.shift_reps, 500),
                                 ("parametric_bootstrap_on_residual", a.boot_reps, 99)):
            seeds = [SEED * 300_000 + r for r in range(nrep)]
            ps = np.array(P.map(real_resid_one, [(method, s, nd) for s in seeds], chunksize=2))
            for al in ALPHAS:
                rows.append(dict(setting="real_euk_7_primary_stems_|U|=30_k=20_bl_placeholder",
                                 statistic="H1H2_NB_counts_mean20_alpha0.3", scheme=method,
                                 alternative="two-sided", alpha=al, rejection_rate=float(np.mean(ps <= al)),
                                 se=float(np.sqrt(al * (1 - al) / len(ps))), n_reps=len(ps), n_draws=nd,
                                 seed=SEED))
            print("real", method, [round(np.mean(ps <= al), 4) for al in ALPHAS], flush=True)
        for method, nrep, nd in (("permute_on_raw_counts", a.shift_reps, 500),
                                 ("permute_on_residual", a.shift_reps, 500),
                                 ("parametric_bootstrap_on_residual", a.boot_reps, 99)):
            for shift in (0.0, 1.0, 2.5):
                seeds = [SEED * 200_000 + int(shift * 10) * 10_000 + r for r in range(nrep)]
                ps = np.array(P.map(shift_one, [(method, shift, s, nd) for s in seeds], chunksize=2))
                for al in ALPHAS:
                    rows.append(dict(setting=f"synthetic_shift_{shift:g}SD_|U|=30_n=7_k=20",
                                     statistic="H1H2_NB_counts_alpha0.3", scheme=method,
                                     alternative="two-sided", alpha=al,
                                     rejection_rate=float(np.mean(ps <= al)),
                                     se=float(np.sqrt(al * (1 - al) / len(ps))), n_reps=len(ps),
                                     n_draws=nd, seed=SEED))
                print(method, shift, [round(np.mean(ps <= al), 4) for al in ALPHAS], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(out, "null_calibration.tsv"), sep="\t", index=False,
                              float_format="%.4g")


def boot_only(a):
    """Append parametric-bootstrap shift rows (used when the node was oversubscribed)."""
    out = os.path.join(ROOT, "results/validation_branchstats/null_calibration.tsv")
    rows = []
    with Pool(a.jobs) as P:
        for shift in [float(x) for x in a.only_boot_shifts.split(",")]:
            seeds = [SEED * 200_000 + int(shift * 10) * 10_000 + r for r in range(a.boot_reps)]
            ps = np.array(P.map(shift_one, [("parametric_bootstrap_on_residual", shift, s, a.n_boot)
                                            for s in seeds], chunksize=1))
            for al in ALPHAS:
                rows.append(dict(setting=f"synthetic_shift_{shift:g}SD_|U|=30_n=7_k=20",
                                 statistic="H1H2_NB_counts_alpha0.3", scheme="parametric_bootstrap_on_residual",
                                 alternative="two-sided", alpha=al, rejection_rate=float(np.mean(ps <= al)),
                                 se=float(np.sqrt(al * (1 - al) / len(ps))), n_reps=len(ps), n_draws=a.n_boot,
                                 seed=SEED))
            print("parametric_bootstrap_on_residual", shift, [round(np.mean(ps <= al), 4) for al in ALPHAS],
                  flush=True)
            pd.DataFrame(rows[-len(ALPHAS):]).to_csv(out, sep="\t", index=False, header=False, mode="a",
                                                    float_format="%.4g")


if __name__ == "__main__":
    main()
