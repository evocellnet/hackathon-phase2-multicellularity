#!/usr/bin/env python
"""POST HOC sensitivity (exploratory_posthoc): duplication re-placement for H4/H5.

Motivation: results/positive_controls/euk/failure_diagnosis.md (C07): Ancestra places duplications one branch too
deep when the sister lineage lacks the family (52 % of duplication products on the Choanozoa stem are metazoan-only),
paid for by compensating losses in the sister lineage.

Rule. Every Ancestra duplication node (OrthoXML paralogGroup) is re-assigned to the most recent species-tree node
at which all of its daughter lineages can still be ancestral, i.e. the MRCA of the species of its members (LCA
reconciliation / "latest-possible placement"). Equivalently: if all daughter HOGs are restricted to one child clade
c of the node Ancestra chose, the duplication moves down to c's stem (recursively), and the k copies that Ancestra
had to carry through the sister lineage of c, and lose there, become one lineage (k - 1 compensating losses are
dropped per sister branch). Implementation: mcell.hogs.parse_hogs with the paralogGroup ``taxonId`` ignored, which
is the parser's documented fallback (duplication at the MRCA of its members); everything else (ortholog levels,
consistency repair, conservation law) is unchanged, so losses and gains are recomputed consistently.

    python scripts/robustness/dup_replacement.py validate     # simulations with known duplication branches
    python scripts/robustness/dup_replacement.py euk [--no-registry]

Outputs: results/robustness/h4h5/replacement_validation_{metrics,branches}.tsv;
         results/robustness/h4h5/replacement_{branch_values,per_transition,tests}.tsv; cache in runs/robustness_h4h5/.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "prereg"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402
from mcell import hogs as mh  # noqa: E402

OUT = os.path.join(ROOT, "results", "robustness", "h4h5")
WORK = os.path.join(ROOT, "runs", "robustness_h4h5")
T0 = time.time()


def log(m):
    print(f"[dup_replacement {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


# =============================================================================
# re-placement parser
# =============================================================================

_ORIG_FAMILY = mh._Builder.family


def _family_mrca_dups(self, F, fam_index):
    for i, k in enumerate(F.kind):
        if k == mh.PARA:
            F.taxid[i] = None
    return _ORIG_FAMILY(self, F, fam_index)


def parse(src, newick=None, replace=False, progress=None):
    """parse_hogs; with ``replace`` paralogGroups are placed at the MRCA of their members."""
    if replace:
        mh._Builder.family = _family_mrca_dups
    try:
        return mh.parse_hogs(src, newick=newick, progress=progress)
    finally:
        mh._Builder.family = _ORIG_FAMILY


# =============================================================================
# validation on simulations
# =============================================================================

SIMS = [  # (name, truth orthoxml, ancestra orthoxml, species tree, branch classes or None)
    *[(f"{s}/{pr}", f"runs/validation_hogparser/{s}/truth/hogs.orthoxml",
       f"runs/validation_hogparser/{s}_ancestra_{pr}/hogs.orthoxml", f"runs/validation_hogparser/{s}/species_tree.nwk",
       None)
      for s in ("S1_12sp_default_hgt", "S2_16sp_highturnover_wgd", "S3_24sp_default") for pr in ("fast", "accurate")],
    *[(f"calibration_tree_{s}/fast", f"runs/calibration_tree/{s}/truth/hogs.orthoxml",
       f"runs/calibration_tree/{s}_ancestra_fast/hogs.orthoxml", f"runs/calibration_tree/{s}/species_tree.nwk", s)
      for s in ("default", "high_turnover")],
]


def _metrics(t, p):
    t, p = np.asarray(t, float), np.asarray(p, float)
    ok = len(t) > 2 and t.std() > 0 and p.std() > 0
    return dict(n_branches=len(t), true_total=int(t.sum()), inferred_total=int(p.sum()),
                ratio=p.sum() / t.sum() if t.sum() else np.nan,
                pearson=pearsonr(t, p)[0] if ok else np.nan, spearman=spearmanr(t, p)[0] if ok else np.nan,
                mae=float(np.abs(t - p).mean()))


def validate(a):
    os.makedirs(OUT, exist_ok=True)
    rows, comp = [], []
    calib = pd.read_csv(os.path.join(ROOT, "results", "calibration_tree", "branch_comparison.tsv"), sep="\t")
    vh = pd.read_csv(os.path.join(ROOT, "results", "validation_hogparser", "branch_comparison.tsv"), sep="\t")
    for name, tx, ax, nwk, cset in SIMS:
        tx, ax, nwk = (os.path.join(ROOT, x) for x in (tx, ax, nwk))
        truth = parse(tx, newick=nwk)
        orig = parse(ax, newick=nwk)
        repl = parse(ax, newick=nwk, replace=True)
        for r in (truth, orig, repl):
            if not r.check_conservation().empty:
                raise SystemExit(f"{name}: conservation violated")
        tb, ob, rb = (x.branches.set_index("branch_child") for x in (truth, orig, repl))
        # cross-check truth against the earlier validated truth tables
        if cset is None:
            s, pr = name.split("/")
            ref = vh[(vh["simulation"] == s) & (vh["preset"] == pr)].set_index("branch_child")["dups_total"]
        else:
            ref = calib[calib["setting"] == cset].set_index("branch")["true_dups"]
        tt = (tb["dups"] + tb["dups_at_origin"]).reindex(ref.index)
        if not (tt == ref).all():
            raise SystemExit(f"{name}: truth differs from the validated truth table")
        if not (ob[["dups", "dups_at_origin", "losses", "gains"]] == orig.branches.set_index("branch_child")[
                ["dups", "dups_at_origin", "losses", "gains"]]).all().all():
            raise SystemExit("bug")
        root = truth.tree.names[0]
        internal = [b for b in tb.index if b != root and not tb.at[b, "child_is_leaf"]]
        allb = [b for b in tb.index if b != root]
        scopes = {"internal_branches": internal, "all_branches": allb}
        if cset is not None:
            cl = calib[calib["setting"] == cset].set_index("branch")["branch_class"]
            scopes["stem_edges+U"] = [b for b in internal if cl.get(b) in ("primary_stem_edge", "control_U")]
        for q, f in (("dups_total", lambda d: d["dups"] + d["dups_at_origin"]), ("dups_non_origin", lambda d: d["dups"]),
                     ("losses", lambda d: d["losses"]), ("gains", lambda d: d["gains"])):
            tv, ov, rv = f(tb), f(ob), f(rb)
            for scope, bb in scopes.items():
                for method, v in (("ancestra", ov), ("replaced", rv)):
                    rows.append(dict(simulation=name, quantity=q, scope=scope, method=method,
                                     **_metrics(tv.reindex(bb), v.reindex(bb))))
        c = pd.DataFrame({"true_dups_total": tb["dups"] + tb["dups_at_origin"],
                          "ancestra_dups_total": ob["dups"] + ob["dups_at_origin"],
                          "replaced_dups_total": rb["dups"] + rb["dups_at_origin"],
                          "true_losses": tb["losses"], "ancestra_losses": ob["losses"], "replaced_losses": rb["losses"],
                          "child_is_leaf": tb["child_is_leaf"]}).drop(index=root)
        c.insert(0, "simulation", name)
        comp.append(c.reset_index())
        log(f"validated {name}")
    M = pd.DataFrame(rows)
    M.to_csv(os.path.join(OUT, "replacement_validation_metrics.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.concat(comp).to_csv(os.path.join(OUT, "replacement_validation_branches.tsv"), sep="\t", index=False)
    w = M.pivot_table(index=["simulation", "quantity", "scope"], columns="method", values=["pearson", "mae", "ratio"])
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(w.round(3).to_string())
    return 0


# =============================================================================
# euk: re-placed H4/H5
# =============================================================================

class _CacheShim:
    def __init__(self, events, branches):
        self.events, self.branches = events, branches


def build_replaced(p):
    os.makedirs(WORK, exist_ok=True)
    done = os.path.join(WORK, "DONE.json")
    if os.path.exists(done):
        log("re-placed cache exists")
    else:
        info = json.load(open(os.path.join(p["cache"], "run_info.json")))
        res = parse(info["orthoxml"], newick=p["tree"], replace=True,
                    progress=lambda n: log(f"{n:,} families") if n % 100000 == 0 else None)
        bad = res.check_conservation()
        if not bad.empty:
            raise SystemExit("conservation violated")
        log(f"parsed: {res.stats}")
        gcl = pd.read_parquet(os.path.join(p["cache"], "gene_classes.parquet"))
        if not (gcl["gene"].to_numpy() == res.genes["gene"].to_numpy()).all():
            raise SystemExit("gene order differs from the cache")
        gflags = gcl[list(an.CLASS_COLS)].astype(bool)
        hc = an.hog_class_table(res.hogs, res.genes, gflags, res.tree)
        res.hogs[["hog_id", "level", "level_idx", "family_id", "parent_hog_id", "parent_idx", "n_genes"]].to_parquet(
            os.path.join(WORK, "hogs.parquet"), index=False)
        hc[["hog_id", "level", "n_genes", "in_tf_dbd", "in_adh_sig"]].to_parquet(
            os.path.join(WORK, "hog_classes.parquet"), index=False)
        ev = res.events[res.events["event"].isin(["dup", "loss"])].drop(columns=["rep_genes"])
        ev.to_parquet(os.path.join(WORK, "events_dup_loss.parquet"), index=False)
        res.branches.to_csv(os.path.join(WORK, "branch_summary.tsv"), sep="\t", index=False)
        json.dump(dict(stats=res.stats, finished=time.strftime("%Y-%m-%dT%H:%M:%S")), open(done, "w"), indent=1,
                  default=str)
    hogs = pd.read_parquet(os.path.join(WORK, "hogs.parquet"))
    hc = pd.read_parquet(os.path.join(WORK, "hog_classes.parquet"))
    ev = pd.read_parquet(os.path.join(WORK, "events_dup_loss.parquet"))
    br = pd.read_csv(os.path.join(WORK, "branch_summary.tsv"), sep="\t", keep_default_na=False, na_values=[])
    return hogs, hc, ev, br


def suite_row(su, hyp, setname, analysis, tier, label, notes=""):
    r = su.relabel
    return dict(hypothesis=hyp, set=setname, analysis=analysis, n=r.n, k=int(np.median(list(r.k.values()))),
                T_obs=r.T_obs, p=r.p, p_mcse=r.p_mcse, pooled_est=su.pooled.get("mu"),
                pooled_ci_low=su.pooled.get("ci_low"), pooled_ci_high=su.pooled.get("ci_high"),
                hksj_ci_low=su.pooled_hksj.get("ci_low"), hksj_ci_high=su.pooled_hksj.get("ci_high"),
                tier=tier, label=label, excluded=json.dumps(su.excluded_stems), notes=notes)


def run_euk(a):
    import check_prereg
    if not check_prereg.verify(strict=False):
        sys.exit("PREREG hash mismatch")
    os.makedirs(OUT, exist_ok=True)
    p = an.run_paths("euk")
    tr, traits, tree, bt = an.load_setup("euk", p)
    hogs, hc, ev, br = build_replaced(p)
    # original (primary) per-branch values for comparison
    prim_bv = pd.read_csv(os.path.join(ROOT, "results", "primary", "euk", "branch_values.tsv"), sep="\t", index_col=0)
    orig_ev = pd.read_parquet(os.path.join(p["cache"], "events.parquet"),
                              columns=["branch_child", "event", "at_origin", "hog_id", "n_events"])
    shim = _CacheShim(ev, br)
    wgd = an.wgd_flags(shim, bt, p["qc_flags"])
    excl = set(wgd.index[wgd["wgd_flag"]])
    log(f"WGD flags (re-placed): {sorted(excl)}")
    V = {}
    for h, col in (("H4", "in_tf_dbd"), ("H5", "in_adh_sig")):
        C = set(hc.loc[hc[col].to_numpy(), "hog_id"])
        V[h] = bs.compositional_statistic(ev, hogs, C, bt, min_d=10, exclude=excl)
    # per-branch table: d, a, s before / after; total dup events before / after; losses before / after
    bvt = bt[["parent", "branch_class", "transition_id", "role", "stem_path", "in_U"]].copy()
    for h in ("H4", "H5"):
        for c in ("d", "a", "s"):
            bvt[f"{h}_{c}_orig"] = prim_bv[f"{h}_{c}"].reindex(bvt.index)
            bvt[f"{h}_{c}_repl"] = V[h][c].reindex(bvt.index)
        bvt[f"{h}_reason_repl"] = V[h]["reason"].reindex(bvt.index)
    o = orig_ev[(orig_ev["event"] == "dup") & ~orig_ev["at_origin"]].groupby("branch_child")["n_events"].sum()
    r = ev[(ev["event"] == "dup") & ~ev["at_origin"]].groupby("branch_child")["n_events"].sum()
    bvt["dup_events_orig"] = bs.path_sum(o.reindex(bt.index).fillna(0), bt)
    bvt["dup_events_repl"] = bs.path_sum(r.reindex(bt.index).fillna(0), bt)
    ol = orig_ev[orig_ev["event"] == "loss"].groupby("branch_child").size()
    rl = ev[ev["event"] == "loss"].groupby("branch_child").size()
    bvt["losses_orig"] = bs.path_sum(ol.reindex(bt.index).fillna(0), bt)
    bvt["losses_repl"] = bs.path_sum(rl.reindex(bt.index).fillna(0), bt)
    bvt["wgd_flag_repl"] = wgd["wgd_flag"].reindex(bvt.index)
    bvt.to_csv(os.path.join(OUT, "replacement_branch_values.tsv"), sep="\t", float_format="%.6g")
    del orig_ev
    # ---- tests as in the primary analysis
    species = pd.read_csv(p["species"], sep="\t")
    groups = an.supergroups(tr, species)
    stems_p = bs.stem_branches(bt)
    allst = bs.stem_branches(bt, primary_only=False)
    sec = sorted(tr.loc[tr["analysis_set"] == "secondary", "transition_id"])
    sets = {"primary": dict(stems_p), "secondary": dict(allst.reindex(sec).dropna())}
    ms = pd.read_csv(os.path.join(ROOT, "results", "primary", "euk", "matched_sets.tsv"), sep="\t")
    prim = pd.read_csv(os.path.join(ROOT, "results", "primary", "euk", "primary_tests.tsv"), sep="\t")
    prim = prim[prim["analysis"] == "primary"].set_index(["hypothesis", "set"])
    rows, pts, msets = [], [], []
    kw = dict(n_draws=a.n_draws, seed=a.seed)
    for h in ("H4", "H5"):
        comp = V[h]
        for setname, st in sets.items():
            log(f"{h} {setname}")
            su = an.run_relabel_suite(h, comp["s"], bt, st, "greater", scale="difference", extra=comp[["log_d"]],
                                      groups=groups, spec_sensitivity=False, **kw)
            mt = an.matched_sets_table(su.match)
            mt.insert(0, "set", setname)
            mt.insert(0, "hypothesis", h)
            msets.append(mt)
            tier = "exploratory_posthoc"
            lbl = "POST HOC sensitivity: duplication re-placement (MRCA of paralogGroup members)"
            hog_p, hog_est = prim.at[(h, setname), "p"], prim.at[(h, setname), "pooled_est"]
            rows.append(suite_row(su, h, setname, "re-placed; re-matched (as primary pipeline)", tier, lbl,
                                  f"orig p={hog_p:.4g}, pooled={hog_est:.3g}"))
            e = su.effects.copy()
            e.insert(0, "analysis", "re-matched")
            e.insert(0, "set", setname)
            e.insert(0, "hypothesis", h)
            pts.append(e.reset_index())
            # same matched sets as the primary analysis
            g = ms[(ms["hypothesis"] == h) & (ms["set"] == setname)]
            fs = {t: list(x.sort_values("rank")["control"]) for t, x in g.groupby("transition_id")}
            fst = {t: x["stem"].iloc[0] for t, x in g.groupby("transition_id")}
            rr = bs.relabel_test(comp["s"], fs, stems=fst, alternative="greater", scheme="permute", **kw)
            ef = bs.effect_table(comp["s"], fs, stems=fst, scale="difference")
            ef["u"] = rr.u.reindex(ef.index)
            ok = ef.dropna(subset=["est"])
            pool = bs.random_effects(ok["est"], ok["var"])
            hk = bs.random_effects(ok["est"], ok["var"], ci_method="hksj")
            rows.append(dict(hypothesis=h, set=setname, analysis="re-placed; primary matched sets fixed", n=rr.n,
                             k=int(np.median(list(rr.k.values()))), T_obs=rr.T_obs, p=rr.p, p_mcse=rr.p_mcse,
                             pooled_est=pool["mu"], pooled_ci_low=pool["ci_low"], pooled_ci_high=pool["ci_high"],
                             hksj_ci_low=hk["ci_low"], hksj_ci_high=hk["ci_high"], tier=tier, label=lbl,
                             excluded=json.dumps(rr.dropped), notes=f"orig p={hog_p:.4g}, pooled={hog_est:.3g}"))
            ef.insert(0, "analysis", "primary sets fixed")
            ef.insert(0, "set", setname)
            ef.insert(0, "hypothesis", h)
            pts.append(ef.reset_index())
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(OUT, "replacement_tests.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.concat(pts, ignore_index=True).to_csv(os.path.join(OUT, "replacement_per_transition.tsv"), sep="\t",
                                              index=False, float_format="%.6g")
    pd.concat(msets, ignore_index=True).to_csv(os.path.join(OUT, "replacement_matched_sets.tsv"), sep="\t",
                                               index=False, float_format="%.4g")
    if not a.no_registry:
        inputs = [os.path.join(p["cache"], "run_info.json"), os.path.join(WORK, "DONE.json"),
                  os.path.join(ROOT, "data", "toolkit_pfam.tsv"), os.path.join(ROOT, "data", "transitions.tsv"),
                  os.path.join(OUT, "replacement_validation_metrics.tsv")]
        for _, r in T.iterrows():
            append_posthoc(a.registry, dict(
                hypothesis=r["hypothesis"],
                statistic=("compositional logit share of stem-duplicated HOGs in "
                           + ("TF-DBD" if r["hypothesis"] == "H4" else "adhesion/ECM U signalling")
                           + " (d>=10), duplications re-placed at the MRCA of their members; one-sided"),
                n=r["n"], k=r["k"], estimate=r["pooled_est"], ci_low=r["pooled_ci_low"], ci_high=r["pooled_ci_high"],
                p=r["p"], p_mcse=r["p_mcse"], holm_p=np.nan, input_md5s=inputs, seed=a.seed, n_draws=a.n_draws,
                scheme="permute",
                notes=(f"label=exploratory_posthoc (POST HOC, motivated by C07 failure diagnosis; not pre-registered); "
                       f"set={r['set']}; analysis={r['analysis']}; T_obs={r['T_obs']:.4f}; {r['notes']}; "
                       "script=scripts/robustness/dup_replacement.py")))
    with pd.option_context("display.width", 250):
        print(T.drop(columns=["label", "excluded"]).to_string())
    return 0


def append_posthoc(path, row):
    """bs.append_registry only accepts tiers primary/exploratory/secondary: validate and format through it on a
    scratch file, then append the formatted line with tier 'exploratory_posthoc' to the real registry."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        tmp = os.path.join(d, "r.tsv")
        r = bs.append_registry(tmp, dict(row, tier="exploratory", git_commit=bs.git_commit(path)))
        line = open(tmp).read().splitlines()[1].split("\t")
    line[bs.REGISTRY_COLUMNS.index("tier")] = "exploratory_posthoc"
    with open(path) as fh:
        if fh.readline().rstrip("\n").split("\t") != bs.REGISTRY_COLUMNS:
            raise ValueError("registry header differs")
    with open(path, "a") as fh:
        fh.write("\t".join(line) + "\n")
    return r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["validate", "euk"])
    ap.add_argument("--n-draws", type=int, default=bs.N_DRAWS)
    ap.add_argument("--seed", type=int, default=bs.SEED)
    ap.add_argument("--registry", default=os.path.join(ROOT, "results", "analysis_registry.tsv"))
    ap.add_argument("--no-registry", action="store_true")
    a = ap.parse_args(argv)
    warnings.simplefilter("ignore")
    return validate(a) if a.what == "validate" else run_euk(a)


if __name__ == "__main__":
    sys.exit(main())
