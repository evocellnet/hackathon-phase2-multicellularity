#!/usr/bin/env python
"""Primary tests H1/H2, H4, H5, H12 (DESIGN §4.6 + v3.1; PREREG.md incl. Dev-1) and the prokaryote descriptive run.

    python scripts/analysis/primary_tests.py --run euk  [--n-boot 9999] [--n-boot-secondary 1999]
    python scripts/analysis/primary_tests.py --run prok            # descriptive only (|U| = 2)

Aborts unless scripts/prereg/check_prereg.verify() passes. Needs runs/analysis/<run>/ (build_branch_events.py) and,
for H1/H2, results/gainladder/<run>/branch_gain_ladder.tsv (aborts if missing unless --allow-missing-ladder).

euk, per test (primary stems = analysis_set 'primary' incl. D7 path stems; secondary = nested transitions, same
statistics, not in Holm):
  H1/H2  s_b = log((y+.5)/(yhat+.5)), y = sequence-novel gains (path-summed), yhat from the NB-GLM on U (D6).
         PRIMARY p = parametric bootstrap (Dev-1, residual_bootstrap_test, two-sided); permute and v3 'spec'
         relabelling p are sensitivities. Effect: log O/E of raw counts vs the matched controls, REML pooled.
  H4/H5  compositional logit statistic (d >= 10), C = TF-DBD (no C2H2) / adhesion-ECM U signalling HOGs (>= 50 %
         of members), WGD-flagged branches excluded (mcell.analysis.wgd_flags; sensitivity with them included);
         matching adds log d; one-sided (greater), permute.
  H12    signed flagellar-module loss f_b (4 modules), w from data/h12_expected_direction_frozen.tsv
         (not_applicable excluded); one-sided, permute. Primary only if >= 2 loss-expected stems keep a value.
All: k = min(20, |U valid|) Mahalanobis-nearest controls, 10,000 draws, seed 11; per-transition table (estimate,
5-95 % control range, percentile u), REML pooled effect (z and HKSJ CI), leave-one-transition-out and
leave-one-supergroup-out. Holm across H1/H2, H4, H5, H9 (results/h9/h9_tests.tsv), H12 when all exist
(scripts/analysis/holm_primary.py). Every test row is appended to results/analysis_registry.tsv.

Outputs: results/primary/euk/ (primary_tests.tsv, per_transition_<test>.tsv, branch_values.tsv, matched_sets.tsv,
smd.tsv, loo.tsv, wgd_flags.tsv, h12_module_coverage.tsv, summary.md, DONE.json) or results/descriptive/prok/.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "prereg"))
sys.path.insert(0, os.path.join(ROOT, "scripts", "analysis"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[primary_tests {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def verify_prereg():
    import check_prereg
    if not check_prereg.verify(strict=False):
        sys.exit("PREREG hash mismatch: aborting (scripts/prereg/check_prereg.py)")
    log("PREREG hashes verified")


# =============================================================================
# shared value computation
# =============================================================================

def compute_values(cache, bt, p, ladder_ok):
    """All per-branch statistics (indexed by branch = child name)."""
    out = {}
    wgd = an.wgd_flags(cache, bt, p["qc_flags"])
    out["wgd"] = wgd
    excl = set(wgd.index[wgd["wgd_flag"]])
    if ladder_ok:
        y = an.h1h2_counts(p["ladder"], bt)
        out["y_seq"] = y
    for name, col in (("H4", "tf_dbd"), ("H5", "adh_sig")):
        C = cache.class_hogs(col)
        out[name] = bs.compositional_statistic(cache.events, cache.hogs, C, bt, min_d=10, exclude=excl)
        out[name + "_with_wgd"] = bs.compositional_statistic(cache.events, cache.hogs, C, bt, min_d=10)
    mh, cov = an.h12_module_hogs(cache)
    out["h12_cov"] = cov
    if len(mh):
        out["H12"] = bs.module_loss_statistic(cache.events, mh, {m: 0 for m in an.H12_MODULES}, bt, hogs=cache.hogs)
    else:
        out["H12"] = pd.DataFrame(index=bt.index, data={"f": np.nan, "n_modules": 0})
    # raw counts per branch (path-summed): HOG gains, dups (HOGs), losses
    br = cache.branches.set_index("branch_child")
    raw = pd.DataFrame(index=bt.index)
    for c in ("gains", "losses", "dups", "dups_at_origin", "gains_ge2_child_clades"):
        raw[c] = bs.path_sum(br[c].reindex(bt.index).astype(float), bt)
    raw["n_parent"] = br["n_parent"].reindex(bt.index)
    raw["n_child"] = br["n_child"].reindex(bt.index)
    for col in ("gain", "dup", "loss"):
        e = cache.events[(cache.events["event"] == col) & ~((cache.events["event"] == "dup") & cache.events["at_origin"])]
        raw[f"{col}_hogs"] = bs.path_sum(e.groupby("branch_child")["hog_id"].nunique().reindex(bt.index)
                                         .fillna(0).astype(float), bt)
    out["raw"] = raw
    return out


def branch_values_table(bt, V):
    cols = ["parent", "depth", "n_desc", "n_sister", "is_terminal", "branch_length", "branch_class",
            "transition_id", "role", "stem_path", "primary", "in_U", "U_exclusion"]
    t = bt[cols].copy()
    t = t.join(V["raw"])
    if "y_seq" in V:
        t["seq_novel"] = V["y_seq"]
    for k in ("H4", "H5"):
        t[f"{k}_d"] = V[k]["d"]
        t[f"{k}_a"] = V[k]["a"]
        t[f"{k}_pi_u"] = V[k]["pi_u"]
        t[f"{k}_s"] = V[k]["s"]
        t[f"{k}_s_with_wgd"] = V[k + "_with_wgd"]["s"]
        t[f"{k}_reason"] = V[k]["reason"]
    t["H12_f"] = V["H12"]["f"]
    t["H12_n_modules"] = V["H12"]["n_modules"]
    t = t.join(V["wgd"][["dup_rate", "ratio_to_median", "wgd_flag", "flag_dup_excess", "flag_qc"]])
    return t


# =============================================================================
# euk
# =============================================================================

def suite_rows(hyp, setname, suite, analysis, scheme_p=None, extra_notes="", tier=None, statistic=""):
    r = suite.relabel
    rows = []
    pooled = suite.pooled or {}
    ph = suite.pooled_hksj or {}
    base = dict(hypothesis=hyp, set=setname, statistic=statistic, n=r.n, k=int(np.median(list(r.k.values()))),
                k_min=min(r.k.values()), T_obs=r.T_obs, pooled_est=pooled.get("mu"), pooled_ci_low=pooled.get("ci_low"),
                pooled_ci_high=pooled.get("ci_high"), pooled_tau2=pooled.get("tau2"), pooled_hksj_ci_low=ph.get("ci_low"),
                pooled_hksj_ci_high=ph.get("ci_high"), min_attainable_p=r.min_attainable_p,
                excluded_stems=";".join(f"{t}:{w}" for t, w in suite.excluded_stems.items()),
                stems=";".join(f"{t}={b}" for t, b in suite.stems.items()))
    for analysis_name, scheme, p, mcse, T, label_tier in scheme_p:
        d = dict(base)
        d.update(analysis=analysis_name, scheme=scheme, p=p, p_mcse=mcse, tier=label_tier, notes=extra_notes)
        if T is not None:
            d["T_obs"] = T
        rows.append(d)
    return rows


def per_transition(suite, extra: pd.DataFrame, setname, hyp):
    e = suite.effects.copy()
    e.insert(0, "set", setname)
    e.insert(0, "hypothesis", hyp)
    e = e.join(extra, on="stem", rsuffix="_branch")
    for t, why in suite.excluded_stems.items():
        e.loc[t, ["hypothesis", "set", "stem"]] = [hyp, setname, ""]
        e.loc[t, "excluded"] = why
    return e.reset_index().rename(columns={"index": "transition_id"})


def run_euk(a, p, cache, tr, traits, bt):
    out_dir = a.out or os.path.join(ROOT, "results", "primary", "euk")
    os.makedirs(out_dir, exist_ok=True)
    ladder_ok = os.path.exists(p["ladder"])
    if not ladder_ok and not a.allow_missing_ladder:
        sys.exit(f"gain ladder missing ({p['ladder']}); H1/H2 needs it. Wait for scripts/gain_ladder.py "
                 "or pass --allow-missing-ladder to run H4/H5/H12 only.")
    V = compute_values(cache, bt, p, ladder_ok)
    bv = branch_values_table(bt, V)
    bv.to_csv(os.path.join(out_dir, "branch_values.tsv"), sep="\t", float_format="%.6g")
    V["wgd"].to_csv(os.path.join(out_dir, "wgd_flags.tsv"), sep="\t", float_format="%.6g")
    V["h12_cov"].to_csv(os.path.join(out_dir, "h12_module_coverage.tsv"), sep="\t", index=False)
    species = pd.read_csv(p["species"], sep="\t")
    stems_p = bs.stem_branches(bt)
    prim = sorted(tr.loc[tr["analysis_set"] == "primary", "transition_id"])
    assert sorted(stems_p.index) == prim, (sorted(stems_p.index), prim)
    allst = bs.stem_branches(bt, primary_only=False)
    sec = sorted(tr.loc[tr["analysis_set"] == "secondary", "transition_id"])
    stems_s = allst.reindex(sec).dropna()
    groups = an.supergroups(tr, species)
    sets = {"primary": dict(stems_p), "secondary": dict(stems_s)}
    log(f"primary stems: {dict(stems_p)}; secondary: {dict(stems_s)}; |U| = {int(bt['in_U'].sum())}")
    wanted = set(a.tests.split(",")) if a.tests != "all" else {"H1H2", "H4", "H5", "H12"}
    rows, pts, smd, msets, loos = [], [], [], [], []
    kw = dict(n_draws=a.n_draws, seed=a.seed)

    def collect(suite, hyp, setname):
        s = suite.match.smd.copy()
        s.insert(0, "set", setname)
        s.insert(0, "hypothesis", hyp)
        smd.append(s)
        ms = an.matched_sets_table(suite.match)
        ms.insert(0, "set", setname)
        ms.insert(0, "hypothesis", hyp)
        msets.append(ms)
        lo = suite.loo.copy()
        lo.insert(0, "set", setname)
        lo.insert(0, "hypothesis", hyp)
        loos.append(lo)

    # ---------------- H1/H2
    if "H1H2" in wanted and ladder_ok:
        y = V["y_seq"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rs = bs.count_residual_statistic(y, bt)
        s = rs["s"]
        for setname, st in sets.items():
            if not st:
                continue
            log(f"H1/H2 {setname}: relabelling")
            su = an.run_relabel_suite("H1H2", s, bt, st, "two-sided", effect_values=y, scale="log_ratio",
                                      groups=groups, **kw)
            nb = a.n_boot if setname == "primary" else a.n_boot_secondary
            log(f"H1/H2 {setname}: parametric bootstrap ({nb} replicates)")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                boot = bs.residual_bootstrap_test(y, bt, su.match, alternative="two-sided", n_boot=nb, seed=a.seed)
            tier_p = "primary" if setname == "primary" else "secondary"
            note = (f"GLM {rs.attrs['method']} alpha={rs.attrs['alpha']:.3g} covariates={rs.attrs['covariates']} "
                    f"fallbacks={rs.attrs['n_fallback']}; effect = log O/E raw counts")
            sp = [("primary", "parametric_bootstrap", boot["p"], boot["p_mcse"], boot["T_obs"], tier_p),
                  ("sensitivity: permute residuals", "permute", su.relabel.p, su.relabel.p_mcse, None,
                   "exploratory" if setname == "primary" else "secondary")]
            if su.relabel_spec is not None:
                sp.append(("sensitivity: v3 spec scheme", "spec", su.relabel_spec.p, su.relabel_spec.p_mcse, None,
                           "exploratory" if setname == "primary" else "secondary"))
            rows += suite_rows("H1H2", setname, su, "", sp, note, statistic="sequence-novel gains, NB-GLM residual "
                               "s=log((y+.5)/(yhat+.5)); T=mean u; two-sided")
            for r in rows[-len(sp):]:
                r["n_boot"] = nb if r["scheme"] == "parametric_bootstrap" else np.nan
                if r["scheme"] == "parametric_bootstrap":
                    r["min_attainable_p"] = 1.0 / (nb + 1)
            ex = pd.DataFrame({"seq_novel": y, "yhat": rs["yhat"], "resid_s": s, "hog_gains": V["raw"]["gains"],
                               "n_edges": bt["stem_path"].map(lambda x: len(x.split(";")) if x else 1)})
            pts.append(per_transition(su, ex, setname, "H1H2"))
            collect(su, "H1H2", setname)
    # ---------------- H4 / H5
    for hyp in ("H4", "H5"):
        if hyp not in wanted:
            continue
        for setname, st in sets.items():
            if not st:
                continue
            log(f"{hyp} {setname}")
            comp = V[hyp]
            su = an.run_relabel_suite(hyp, comp["s"], bt, st, "greater", scale="difference",
                                      extra=comp[["log_d"]], groups=groups, **kw)
            tier_p = "primary" if setname == "primary" else "secondary"
            sp = [("primary", "permute", su.relabel.p, su.relabel.p_mcse, None, tier_p)]
            if su.relabel_spec is not None:
                sp.append(("sensitivity: v3 spec scheme", "spec", su.relabel_spec.p, su.relabel_spec.p_mcse, None,
                           "exploratory" if setname == "primary" else "secondary"))
            excl = [b for b in V["wgd"].index[V["wgd"]["wgd_flag"]] if bt.at[b, "in_U"] or b in st.values()]
            rows += suite_rows(hyp, setname, su, "", sp,
                               f"WGD-flagged branches excluded ({len(excl)} in U/stems: {';'.join(excl)}); "
                               "effect = difference in adjusted logit (log odds ratio)",
                               statistic=f"compositional logit share of stem-duplicated HOGs in "
                                         f"{'TF-DBD' if hyp == 'H4' else 'adhesion/ECM U signalling'} (d>=10); one-sided")
            # sensitivity: WGD branches included
            cw = V[hyp + "_with_wgd"]
            try:
                su2 = an.run_relabel_suite(hyp, cw["s"], bt, st, "greater", scale="difference",
                                           extra=cw[["log_d"]], spec_sensitivity=False, **kw)
                rows += suite_rows(hyp, setname, su2, "", [("sensitivity: WGD-flagged branches included", "permute",
                                                            su2.relabel.p, su2.relabel.p_mcse, None,
                                                            "exploratory" if setname == "primary" else "secondary")],
                                   "", statistic=rows[-1]["statistic"])
            except ValueError as e:
                log(f"{hyp} with WGD: {e}")
            ex = comp[["d", "a", "n_u", "c_u", "pi_u", "pi_smoothed", "s", "reason", "n_edges"]].rename(
                columns={"s": "stat_s"})
            ex["wgd_flag"] = V["wgd"]["wgd_flag"]
            pts.append(per_transition(su, ex, setname, hyp))
            collect(su, hyp, setname)
    # ---------------- H12
    if "H12" in wanted:
        d = pd.read_csv(p["h12_dir"], sep="\t").set_index("transition_id")["expected"]
        appl = {t: b for t, b in sets["primary"].items() if str(d.get(t, "not_applicable")) != "not_applicable"}
        signs = bs.h12_signs({t: d[t] for t in appl})
        f = V["H12"]["f"]
        loss_exp = [t for t in appl if signs[t] > 0]
        remain = [t for t in loss_exp if np.isfinite(f.get(appl[t], np.nan))]
        tier = "primary" if len(remain) >= 2 else "exploratory"
        log(f"H12: loss-expected stems with a value: {remain} -> tier {tier}")
        su = an.run_relabel_suite("H12", f, bt, appl, "greater", scale="difference", signs=signs, groups=groups, **kw)
        sp = [("primary", "permute", su.relabel.p, su.relabel.p_mcse, None, tier)]
        if su.relabel_spec is not None:
            sp.append(("sensitivity: v3 spec scheme", "spec", su.relabel_spec.p, su.relabel_spec.p_mcse, None,
                       "exploratory"))
        rows += suite_rows("H12", "primary", su, "", sp,
                           f"signs {signs}; not_applicable excluded: "
                           f"{[t for t in sets['primary'] if t not in appl]}; loss-expected stems with value: {remain}"
                           + ("" if tier == "primary" else "; DEMOTED to exploratory (<2 loss-expected stems)"),
                           statistic="signed flagellar-module loss w*f_b (4 modules); one-sided")
        ex = V["H12"].drop(columns=["parent"], errors="ignore").rename(columns={"f": "f_b"})
        ex["expected"] = pd.Series({b: d.get(t) for t, b in appl.items()})
        pt = per_transition(su, ex, "primary", "H12")
        pt["w"] = pt["transition_id"].map(signs)
        pts.append(pt)
        collect(su, "H12", "primary")
        # secondary: no frozen direction -> descriptive f only
        for t, b in sets["secondary"].items():
            pts.append(pd.DataFrame([dict(transition_id=t, hypothesis="H12", set="secondary", stem=b,
                                          y=f.get(b, np.nan), excluded="no frozen direction: descriptive f_b only")]))
    tests = pd.DataFrame(rows)
    # ---------------- Holm
    tests.to_csv(os.path.join(out_dir, "primary_tests.tsv"), sep="\t", index=False, float_format="%.6g")
    import holm_primary
    holm_tab, missing = holm_primary.compute(os.path.join(out_dir, "primary_tests.tsv"), holm_primary.H9)
    tests["holm_p"] = np.nan
    if not missing:
        hp = holm_tab.set_index("hypothesis")["holm_p"]
        sel = (tests["set"] == "primary") & (tests["analysis"] == "primary")
        tests.loc[sel, "holm_p"] = tests.loc[sel, "hypothesis"].map(hp)
        holm_tab.to_csv(os.path.join(out_dir, "holm.tsv"), sep="\t", index=False, float_format="%.6g")
        log(f"Holm over {holm_tab['holm_p'].notna().sum()} tests")
    else:
        log(f"Holm deferred: missing {missing} (run scripts/analysis/holm_primary.py later)")
    tests.to_csv(os.path.join(out_dir, "primary_tests.tsv"), sep="\t", index=False, float_format="%.6g")
    pt = pd.concat(pts, ignore_index=True)
    pt.to_csv(os.path.join(out_dir, "per_transition.tsv"), sep="\t", index=False, float_format="%.6g")
    for h, g in pt.groupby("hypothesis"):
        g.dropna(axis=1, how="all").to_csv(os.path.join(out_dir, f"per_transition_{h}.tsv"), sep="\t", index=False,
                                           float_format="%.6g")
    pd.concat(smd, ignore_index=True).to_csv(os.path.join(out_dir, "smd.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.concat(msets, ignore_index=True).to_csv(os.path.join(out_dir, "matched_sets.tsv"), sep="\t", index=False,
                                               float_format="%.4g")
    pd.concat(loos, ignore_index=True).to_csv(os.path.join(out_dir, "loo.tsv"), sep="\t", index=False, float_format="%.6g")
    # ---------------- registry
    inputs = [x for x in (p["transitions"], p["guide"], p["bl"], p["ladder"], os.path.join(p["cache"], "run_info.json"),
                          p["toolkit"], p["h12_sets"], p["h12_dir"]) if os.path.exists(x)]
    for _, r in tests.iterrows():
        bs.append_registry(a.registry, dict(
            hypothesis=r["hypothesis"], statistic=r["statistic"], n=r["n"], k=r["k"], estimate=r["pooled_est"],
            ci_low=r["pooled_ci_low"], ci_high=r["pooled_ci_high"], p=r["p"], p_mcse=r["p_mcse"],
            holm_p=r["holm_p"], tier=r["tier"], input_md5s=inputs, seed=a.seed,
            n_draws=r.get("n_boot") if r["scheme"] == "parametric_bootstrap" else a.n_draws, scheme=r["scheme"],
            notes=f"set={r['set']}; analysis={r['analysis'] or 'primary'}; T_obs={r['T_obs']:.4f}; "
                  f"min_attainable_p={r['min_attainable_p']:.3g}; {r['notes']}"))
    write_summary_euk(out_dir, tests, pt, V, bt)
    json.dump(dict(status="done", finished=dt.datetime.now().isoformat(timespec="seconds"), n_tests=len(tests),
                   holm_missing=missing, ladder=ladder_ok, wall_s=round(time.time() - T0)),
              open(os.path.join(out_dir, "DONE.json"), "w"), indent=1)
    log(f"wrote {out_dir}")


def write_summary_euk(out_dir, tests, pt, V, bt):
    L = ["# Primary tests (eukaryotes)", "",
         f"Generated by `scripts/analysis/primary_tests.py` on {dt.date.today()}; numbers only (no interpretation).",
         "Spec: DESIGN.md §4.6 + v3.1, PREREG.md (Dev-1: H1/H2 primary p = parametric bootstrap).", "",
         "## Tests", "",
         "| hypothesis | set | analysis | scheme | n | k | T_obs | p | Holm p | pooled est [95% CI] | tier |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tests.iterrows():
        L.append(f"| {r['hypothesis']} | {r['set']} | {r['analysis']} | {r['scheme']} | {r['n']} | {r['k']} | "
                 f"{r['T_obs']:.3f} | {an.fmt_p(r['p'])} | {an.fmt_p(r.get('holm_p'))} | "
                 f"{r['pooled_est']:.3g} [{r['pooled_ci_low']:.3g}, {r['pooled_ci_high']:.3g}] | {r['tier']} |")
    L += ["", "## Per transition", "",
          "| hypothesis | set | transition | stem | y | expected | est | control 5-95% | u | note |", "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in pt.iterrows():
        L.append(f"| {r['hypothesis']} | {r['set']} | {r['transition_id']} | {r.get('stem', '')} | {r.get('y', np.nan):.4g} | "
                 f"{r.get('expected', np.nan):.4g} | {r.get('est', np.nan):.3g} | "
                 f"[{r.get('pi_low', np.nan):.3g}, {r.get('pi_high', np.nan):.3g}] | {r.get('u', np.nan):.3f} | "
                 f"{r.get('excluded', '') if isinstance(r.get('excluded', ''), str) else ''} |")
    w = V["wgd"]
    L += ["", f"WGD/duplication-excess flags (internal branches): {', '.join(w.index[w['wgd_flag'] & ~bt['is_terminal']]) or 'none'}; "
          f"median internal duplication rate {w['median_rate_internal'].iloc[0]:.4g}.", "",
          "Files: primary_tests.tsv, per_transition*.tsv, loo.tsv, smd.tsv, matched_sets.tsv, branch_values.tsv, "
          "wgd_flags.tsv, h12_module_coverage.tsv, holm.tsv (when H9 exists)."]
    open(os.path.join(out_dir, "summary.md"), "w").write("\n".join(L) + "\n")


# =============================================================================
# prok (descriptive)
# =============================================================================

def prok_stems(tr, bt):
    """(transition_id, analysis_set, stem branch, terminal) for every prokaryote transition."""
    rows = []
    for _, r in tr.iterrows():
        Ms = [m.strip() for m in str(r["M"]).split(",") if m.strip() in bt.index]
        for m in Ms:
            rows.append(dict(transition_id=r["transition_id"], analysis_set=r["analysis_set"], label=r["label"],
                             from_state=r["from_state"], to_state=r["to_state"], stem=m,
                             stem_path=bt.at[m, "stem_path"] or m, terminal=bool(bt.at[m, "is_terminal"]),
                             M_list=len(Ms) > 1))
    return pd.DataFrame(rows)


def run_prok(a, p, cache, tr, traits, bt):
    out_dir = a.out or os.path.join(ROOT, "results", "descriptive", "prok")
    os.makedirs(out_dir, exist_ok=True)
    ladder_ok = os.path.exists(p["ladder"])
    V = compute_values(cache, bt, p, ladder_ok)
    bv = branch_values_table(bt, V)
    bv.to_csv(os.path.join(out_dir, "branch_values.tsv"), sep="\t", float_format="%.6g")
    st = prok_stems(tr, bt)
    U = list(bt.index[bt["in_U"]])
    bce = cache.branch_class_events.set_index("branch_child")
    rows = []
    for _, r in st.iterrows():
        b = r["stem"]
        edges = bs.path_edges(bt, b)
        d = r.to_dict()
        raw = V["raw"].loc[b]
        d.update(n_edges=len(edges), P=bt.at[edges[0], "parent"], n_hogs_at_P=raw["n_parent"] if len(edges) == 1 else
                 cache.branches.set_index("branch_child").at[edges[0], "n_parent"], n_hogs_at_M=raw["n_child"],
                 hog_gains=raw["gains"], gains_ge2_child_clades=raw["gains_ge2_child_clades"],
                 dup_events=raw["dups"], dup_hogs=raw["dup_hogs"], losses=raw["losses"],
                 branch_length=bt.at[b, "branch_length"], n_desc=bt.at[b, "n_desc"])
        d["seq_novel"] = V["y_seq"].get(b, np.nan) if ladder_ok else np.nan
        for c in an.CLASS_COLS:
            for ev in ("gain", "loss"):
                d[f"{ev}_{c}"] = float(bce.reindex(edges)[f"{ev}_{c}"].fillna(0).sum())
        for hyp, cls in (("H4", "tf_dbd"), ("H5", "adh_sig")):
            cmp_ = V[hyp + "_with_wgd"].loc[b]
            d[f"{hyp}_d"], d[f"{hyp}_a"], d[f"{hyp}_pi_u"] = cmp_["d"], cmp_["a"], cmp_["pi_u"]
            d[f"{hyp}_s"] = cmp_["s"]
            d[f"{hyp}_reason"] = cmp_["reason"]
        d["wgd_flag"] = bool(V["wgd"].at[b, "wgd_flag"])
        rows.append(d)
    T = pd.DataFrame(rows)
    # per-class duplication composition (union over path edges, lifted to P)
    for c in an.CLASS_COLS:
        comp = bs.compositional_statistic(cache.events, cache.hogs, cache.class_hogs(c), bt, min_d=0)
        T[f"dup_{c}"] = T["stem"].map(comp["a"])
    if ladder_ok:
        yU = V["y_seq"].reindex(U)
        T["seq_novel_U_mean"] = float(yU.mean())
        T["seq_novel_logOE_vs_U"] = np.log((T["seq_novel"] + 0.5) / (yU.mean() + 0.5))
    T.to_csv(os.path.join(out_dir, "per_transition.tsv"), sep="\t", index=False, float_format="%.6g")
    ctrl = bv.loc[U]
    ctrl.to_csv(os.path.join(out_dir, "U_controls.tsv"), sep="\t", float_format="%.6g")
    V["wgd"].to_csv(os.path.join(out_dir, "wgd_flags.tsv"), sep="\t", float_format="%.6g")
    # registry: one descriptive row per statistic
    inputs = [x for x in (p["transitions"], p["guide"], p["bl"], os.path.join(p["cache"], "run_info.json"),
                          p["ladder"]) if os.path.exists(x)]
    sec = T[(T["analysis_set"] == "secondary") & ~T["terminal"]]
    for hyp, stat, col in (("H1H2", "sequence-novel gains per stem (descriptive)", "seq_novel"),
                           ("H4", "compositional logit share TF-DBD of stem-duplicated HOGs (descriptive)", "H4_s"),
                           ("H5", "compositional logit share adhesion/ECM U signalling (descriptive)", "H5_s")):
        if col not in sec.columns or sec[col].isna().all():
            continue
        bs.append_registry(a.registry, dict(
            hypothesis=hyp, statistic=stat, n=int(sec[col].notna().sum()), k=len(U),
            estimate=float(sec[col].median()), p="NA", p_mcse="NA", holm_p="NA", tier="exploratory", scheme="none",
            seed=a.seed, input_md5s=inputs,
            notes="prok descriptive only (|U|=2, no relabelling test); estimate = median over secondary prokaryote stems; "
                  + "; ".join(f"{t}={v:.4g}" for t, v in zip(sec["transition_id"], sec[col]))))
    write_summary_prok(out_dir, T, ctrl, V, ladder_ok, cache)
    json.dump(dict(status="done", finished=dt.datetime.now().isoformat(timespec="seconds"), ladder=ladder_ok,
                   n_transitions=int(T["transition_id"].nunique())), open(os.path.join(out_dir, "DONE.json"), "w"), indent=1)
    log(f"wrote {out_dir}")
    return T


def _f(x, nd=0):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.{nd}f}" if isinstance(x, (float, np.floating)) else str(x)


def write_summary_prok(out_dir, T, ctrl, V, ladder_ok, cache):
    pc = os.path.join(ROOT, "results", "positive_controls", "prok", "controls.tsv")
    L = ["# Prokaryote run: descriptive per-transition summary", "",
         f"Generated by `scripts/analysis/primary_tests.py --run prok` on {dt.date.today()} from `runs/prok_fast` "
         f"(Ancestra fast; {cache.info['stats']['n_genes']:,} genes, {cache.info['stats']['n_species']} species, "
         f"{len(cache.hogs):,} ancestral genes at internal nodes). Descriptive only: |U| = {len(ctrl)} control branches, "
         "so no relabelling test (DESIGN v3.1). Numbers only; no interpretation.", "",
         "Definitions: HOG gains / losses / duplications from `mcell.hogs` (stems = branch into M; path stems summed over edges). "
         "dup HOGs = HOGs at the parent with >= 1 non-origin duplication. Class counts: HOGs with >= 50% of members carrying a "
         "pre-registered class domain (data/toolkit_pfam.tsv; pre-registered Pfam subset). "
         + ("Sequence-novel gains from results/gainladder/prok/branch_gain_ladder.tsv." if ladder_ok else
            "Sequence-novel gains: gain ladder not yet available (column NA; re-run when results/gainladder/prok/ exists)."),
         "", "## Gene-content change on each transition stem", "",
         "| transition | set | stem (edges) | HOGs at P | HOG gains | gains in >=2 child clades | seq-novel gains | dup HOGs | dup events | losses |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in T.iterrows():
        L.append(f"| {r['transition_id']}{' (terminal)' if r['terminal'] else ''} | {r['analysis_set']} | {r['stem']} ({r['n_edges']}) | "
                 f"{_f(r['n_hogs_at_P'])} | {_f(r['hog_gains'])} | {_f(r['gains_ge2_child_clades'])} | {_f(r['seq_novel'])} | "
                 f"{_f(r['dup_hogs'])} | {_f(r['dup_events'])} | {_f(r['losses'])} |")
    L += ["", "## Class composition on each stem (HOG counts: gained / duplicated / lost)", "",
          "| transition | stem | " + " | ".join(an.CLASS_COLS) + " |", "|---|---|" + "---|" * len(an.CLASS_COLS)]
    for _, r in T.iterrows():
        L.append(f"| {r['transition_id']} | {r['stem']} | " + " | ".join(
            f"{_f(r[f'gain_{c}'])}/{_f(r[f'dup_{c}'])}/{_f(r[f'loss_{c}'])}" for c in an.CLASS_COLS) + " |")
    L += ["", "## H4/H5 compositional statistic (descriptive; WGD-flagged branches not excluded here)", "",
          "| transition | stem | d | a (TF-DBD) | pi_u (TF) | s H4 | a (adh/ECM U sig.) | pi_u | s H5 |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in T.iterrows():
        L.append(f"| {r['transition_id']} | {r['stem']} | {_f(r['H4_d'])} | {_f(r['H4_a'])} | {_f(r['H4_pi_u'], 4)} | "
                 f"{_f(r['H4_s'], 3) if r['H4_reason'] == '' else r['H4_reason']} | {_f(r['H5_a'])} | {_f(r['H5_pi_u'], 4)} | "
                 f"{_f(r['H5_s'], 3) if r['H5_reason'] == '' else r['H5_reason']} |")
    z = T[(T["H4_a"] == 0) & (T["H4_reason"] == "")]
    if len(z):
        L += ["", f"Note on the frozen H4 statistic: {len(z)} stems have a = 0 TF-DBD duplicated HOGs but s > 0 whenever "
              "pi_u < 0.5/(d+1) (continuity correction); here pi_u(TF-DBD) <= "
              f"{T['H4_pi_u'].max():.2g}, i.e. pre-registered TF-DBD domains are nearly absent from prokaryote HOGs."]
    L += ["", f"Control pool U ({len(ctrl)} branches): " + "; ".join(
        f"{b}: gains {_f(r['gains'])}, seq-novel {_f(r.get('seq_novel', np.nan))}, dup HOGs {_f(r['dup_hogs'])}, losses {_f(r['losses'])}, "
        f"H4 s {_f(r['H4_s_with_wgd'], 3)}, H5 s {_f(r['H5_s_with_wgd'], 3)}" for b, r in ctrl.iterrows()) + ".",
        "H12 (flagellar modules) is not applicable to prokaryotes: "
        f"{int(sum(cache.hog_classes[c].sum() for c in cache.hog_classes.columns if c.startswith('mod_')))} module HOGs in this run.",
        f"Duplication-excess flags (>= 3x median internal rate; internal branches): "
        f"{', '.join(V['wgd'].index[V['wgd']['flag_dup_excess']]) or 'none'}.", ""]
    L += ["## Positive controls P1/P2", ""]
    if os.path.exists(pc):
        c = pd.read_csv(pc, sep="\t")
        for _, r in c.iterrows():
            L.append(f"- **{r['control']}** ({r['expected']}): recovered = **{r['recovered']}**. {r['observed']}")
    else:
        L.append("Not yet scored (run scripts/analysis/score_positive_controls.py --run prok).")
    L += ["", "Files: per_transition.tsv, branch_values.tsv, U_controls.tsv, wgd_flags.tsv; positive controls in "
          "results/positive_controls/prok/."]
    open(os.path.join(out_dir, "summary.md"), "w").write("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--tests", default="all", help="comma list of H1H2,H4,H5,H12 (euk)")
    ap.add_argument("--n-draws", type=int, default=bs.N_DRAWS)
    ap.add_argument("--n-boot", type=int, default=9999, help="H1/H2 parametric bootstrap replicates (primary)")
    ap.add_argument("--n-boot-secondary", type=int, default=1999)
    ap.add_argument("--seed", type=int, default=bs.SEED)
    ap.add_argument("--registry", default=os.path.join(ROOT, "results", "analysis_registry.tsv"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--allow-missing-ladder", action="store_true")
    a = ap.parse_args(argv)
    verify_prereg()
    p = an.run_paths(a.run)
    cache = an.load_cache(a.run, p)
    tr, traits, tree, bt = an.load_setup(a.run, p)
    names = set(cache.tree.names) - {cache.tree.names[0]}
    if names != set(bt.index):
        sys.exit("branch table and HOG tree differ")
    if a.run == "euk":
        run_euk(a, p, cache, tr, traits, bt)
    else:
        run_prok(a, p, cache, tr, traits, bt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
