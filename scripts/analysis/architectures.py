#!/usr/bin/env python
"""EXPLORATORY (H11, DESIGN Q4). New ordered Pfam architectures on transition stems vs control branches.

    python scripts/analysis/architectures.py --run {euk,prok} [--transitions default] [--allow-partial]

Architecture of a protein = Pfam accessions (full Pfam, gathering cutoffs) ordered by alignment start after
greedy overlap resolution by score, tandem repeats collapsed (A-A-B -> A-B), as in mcell.gainladder.
An architecture is *new on branch p -> c* when c is the MRCA of the annotated species carrying it and it
occurs in >= 2 species (so it is present in >= 2 descendant species of c and absent from every annotated species
outside c). Per branch (D7 path stems summed over their edges): new_arch_all, new_arch_multi (>= 2 distinct
domains; the H11 statistic), new_arch_recombined (multi-domain and every component domain occurs outside the clade).

Stems vs controls (mcell.branchstats machinery; exploratory): s_b = log((y+.5)/(yhat+.5)) with yhat from the
NB-GLM on the control pool (count_residual_statistic, as H1/H2 in v3.1), relabelling scheme 'permute', one-sided
(greater), 10,000 draws, seed 11; also the raw-count statistic; effect = log O/E vs the matched controls, REML pooled.
Controls = U; if |U| < 10 (prokaryotes) a relaxed pool is used and the test is labelled dry-run/descriptive.
Requires full Pfam for >= 95 % of species unless --allow-partial (then coverage columns must be read).

Outputs: results/exploratory/<run>/architectures/: branch_new_architectures.tsv, per_transition.tsv,
tests.tsv, stem_architectures.tsv.gz (new architectures on each stem), summary.json; registry rows.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402
from mcell import explore as ex  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[architectures {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def relabel(name, vals, y, bt, stems, pool, n_draws, seed):
    """(RelabelResult, effects, pooled) via an.run_relabel_suite; on degenerate data (e.g. zero pooled SE)
    the relabel test and effect table are run directly and pooling is skipped when it fails."""
    try:
        su = an.run_relabel_suite(name, vals, bt, stems, "greater", effect_values=y, scale="log_ratio",
                                  n_draws=n_draws, seed=seed, spec_sensitivity=False, pool=pool)
        return su.relabel, su.effects, su.pooled or {}
    except ZeroDivisionError:
        st = {t: b for t, b in stems.items() if np.isfinite(vals.get(b, np.nan))}
        pl = [b for b in pool if np.isfinite(vals.get(b, np.nan))]
        m = bs.match_controls(bt, st, pool=pl)
        r = bs.relabel_test(vals, m, alternative="greater", n_draws=n_draws, seed=seed)
        eff = bs.effect_table(y, m, scale="log_ratio")
        eff["u"] = r.u.reindex(eff.index)
        ok = eff.dropna(subset=["est"]) if "est" in eff.columns else eff.iloc[0:0]
        try:
            pooled = bs.random_effects(ok["est"], ok["var"]) if len(ok) else {}
        except ZeroDivisionError:
            pooled = {}
        return r, eff, pooled


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--transitions", default="default")
    ap.add_argument("--n-draws", type=int, default=bs.N_DRAWS)
    ap.add_argument("--seed", type=int, default=bs.SEED)
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--min-coverage", type=float, default=0.95)
    ap.add_argument("--no-registry", action="store_true")
    ap.add_argument("--skip-prereg-check", action="store_true")
    a = ap.parse_args(argv)
    if not a.skip_prereg_check and not ex.prereg_ok(ROOT):
        sys.exit("PREREG hash mismatch: aborting")
    ctx = ex.load_context(a.run, a.transitions, need_cache=False)
    od = ex.out_dir(a.run, "architectures")
    tabdir = os.path.join(ROOT, "runs", "analysis", a.run, "pfam_full")
    man = ex.build_pfam_tables(ctx.paths["pfam_full"], ctx.species, tabdir, log=log)
    ann = man["annotated"]
    cov = len(ann) / len(ctx.species)
    if cov < a.min_coverage and not a.allow_partial:
        sys.exit(f"full Pfam for {len(ann)}/{len(ctx.species)} species (< {a.min_coverage:.0%}); "
                 "rerun when annotation is complete or pass --allow-partial")
    arch = pd.read_parquet(os.path.join(tabdir, "arch.parquet"))
    log(f"{len(arch):,} annotated proteins, {arch['arch'].nunique():,} architectures, {len(ann)} species")
    org = ex.architecture_origins(arch, ctx.tree, ann)
    br = ex.branch_architecture_counts(org, ctx.bt, ctx.tree, ann)
    bt = ctx.bt
    tab = bt[["parent", "depth", "n_desc", "n_sister", "is_terminal", "branch_class", "transition_id", "stem_path",
              "in_U"]].join(br)
    tab["in_pool"] = tab.index.isin(ctx.pool)
    stems = ctx.stems
    rev = {b: t for t, b in stems.items()}
    tab["selected_transition"] = [rev.get(b, "") for b in tab.index]
    tab.to_csv(os.path.join(od, "branch_new_architectures.tsv"), sep="\t", index_label="branch")
    # architectures new on the stems
    srows = []
    for t, b in stems.items():
        edges = bs.path_edges(bt, b)
        s = org[org["origin"].isin(edges)].assign(transition=t)
        srows.append(s)
    if srows:
        pd.concat(srows).to_csv(os.path.join(od, "stem_architectures.tsv.gz"), sep="\t", index=False)

    tests, pts = [], []
    dry = ctx.pool_label != "U"
    for stat in ("new_arch_multi", "new_arch_recombined", "new_arch_all"):
        y = tab[stat].astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = bs.count_residual_statistic(y, bt, pool=ctx.pool)
        s = res["s"]
        at = res.attrs
        fits = (f"method={at.get('method')}; alpha={at.get('alpha')}; covariates={at.get('covariates')}; "
                f"n_fallback={at.get('n_fallback')}")
        for label, vals in (("residual", s), ("raw", y)):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    r, effects, pooled = relabel(f"H11 {stat} {label}", vals, y, bt, stems, ctx.pool, a.n_draws, a.seed)
            except Exception as e:
                log(f"{stat}/{label}: relabel failed: {e}")
                continue
            tests.append(dict(statistic=stat, version=label, n=r.n, k=int(np.median(list(r.k.values()))),
                              T_obs=r.T_obs, p=r.p, p_mcse=r.p_mcse, min_attainable_p=r.min_attainable_p,
                              pooled_log_OE=pooled.get("mu"), ci_low=pooled.get("ci_low"), ci_high=pooled.get("ci_high"),
                              pool=ctx.pool_label, glm=fits if label == "residual" else ""))
            if label == "residual":
                e = effects.copy()
                e["statistic"] = stat
                e["y_stem"] = [y.get(stems[t], np.nan) for t in e.index]
                e["s_stem"] = [s.get(stems[t], np.nan) for t in e.index]
                pts.append(e.reset_index())
    T = pd.DataFrame(tests)
    if len(T):
        T["q_bh"] = ex.bh(T["p"].to_numpy())
    T.to_csv(os.path.join(od, "tests.tsv"), sep="\t", index=False)
    PT = pd.concat(pts, ignore_index=True) if pts else pd.DataFrame()
    PT.to_csv(os.path.join(od, "per_transition.tsv"), sep="\t", index=False)
    if not a.no_registry:
        for _, r in T.iterrows():
            bs.append_registry(ctx.paths["registry"], ex.registry_row(
                "H11 new domain architectures",
                f"{r['statistic']} per branch ({r['version']}: "
                + ("NB-GLM residual log((y+.5)/(yhat+.5))" if r["version"] == "residual" else "raw count")
                + "); relabelling T = mean u, one-sided greater; estimate = pooled log O/E",
                n=int(r["n"]), k=int(r["k"]), estimate=r["pooled_log_OE"], ci=(r["ci_low"], r["ci_high"]), p=r["p"],
                seed=a.seed, n_draws=a.n_draws, scheme="permute",
                notes=f"run={a.run}; T_obs={r['T_obs']:.4g}; q_BH(within H11 rows)={r['q_bh']:.3g}; pool={ctx.pool_label}; "
                      f"pfam_annotated={len(ann)}/{len(ctx.species)}"
                      + ("; DRY RUN prok, relaxed pool: not interpretable" if dry else ""),
                inputs=[os.path.join(tabdir, "manifest.json"), ctx.paths["transitions"]]))
    summ = dict(run=a.run, pool=ctx.pool_label, n_pool=len(ctx.pool), n_annotated=len(ann), n_species=len(ctx.species),
                n_architectures=int(len(org)), n_multi=int(org["multi"].sum()),
                n_new_multi_any_branch=int(((org["origin"] != "") & org["multi"]).sum()),
                stems={t: dict(branch=b, **{c: float(tab.at[b, c]) for c in ("new_arch_all", "new_arch_multi",
                                                                           "new_arch_recombined")})
                       for t, b in stems.items()},
                pool_mean={c: float(tab.loc[ctx.pool, c].mean()) for c in ("new_arch_all", "new_arch_multi",
                                                                         "new_arch_recombined")},
                tests=T.round(5).to_dict("records"))
    ex.write_summary(od, summ)
    log(f"done -> {os.path.relpath(od, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
