#!/usr/bin/env python
"""EXPLORATORY. Pfam-wide view of transitions (DESIGN Q2, §3.3, §4.4 orthology-free Pfam-count version).

    python scripts/analysis/pfam_wide.py --run {euk,prok} [--transitions default] [--n-perm 10000] [--seed 11]

(a) Per species x Pfam domain copy numbers (proteins carrying the domain) from the full Pfam annotation
    data/annotations/pfam/<run>/<code>.pfam.tsv.gz (species without a finished file are listed and left out,
    never counted as zeros). Cached in runs/analysis/<run>/pfam_full/.
(b) Per transition (M clade vs sister species; mcell.explore.domain_contrast): expanded domains (one-sided
    Mann-Whitney on per-species counts, exact by enumeration when small, BH within transition, q <= 0.05; or
    log2((mean_M+1)/(mean_S+1)) >= 1 with presence in >= 50 % of M), gained (>= 50 % of M, absent from every
    annotated species outside M), lost (absent from all M, present in >= 50 % of S), contracted.
    HOG view (needs runs/analysis/<run>/): per domain, HOGs gained / duplicated / lost on the stem edges that
    carry the domain in >= 50 % of their members; flag 'hog_expanded' = gain+dup HOGs >= 2 and percentile
    u >= 0.95 among the matched controls.
(c) Convergence: k = number of transitions in which a domain is expanded-or-gained (also: gained only, lost,
    HOG-expanded) vs the same statistic on pseudo-transitions (child clades of the control pool U: M = tips(v),
    S = tips(u) - tips(v)); null draws one matched control per transition (distinct), 10,000 permutations,
    seed 11 (mcell.explore.convergence_test). If |U| < 10 (prokaryotes) a relaxed pool is used and labelled.

Outputs: results/exploratory/<run>/pfam_wide/: domain_contrasts.tsv.gz (transitions), control_flags.tsv.gz,
convergence_<flag>.tsv + convergence_global.tsv, heatmap_status.tsv (domain x transition: 2 gained, 1 expanded,
-1 contracted, -2 lost), heatmap_log2fc.tsv, hog_domain_events.tsv.gz, heatmap_hog.tsv, toolkit_by_transition.tsv,
coverage.tsv, summary.json. Registry rows (tier exploratory) for every convergence test.
Single-threaded; run with nice -n 19.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402
from mcell import explore as ex  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[pfam_wide {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def pfam_names(pfam_tab_dir):
    n = pd.read_csv(os.path.join(pfam_tab_dir, "names.tsv"), sep="\t")
    names = dict(zip(n["pfam_acc"], n["pfam_name"]))
    clan = {}
    cf = os.path.join(ROOT, "data", "raw", "pfam", "Pfam-A.clans.tsv.gz")
    if os.path.exists(cf):
        c = pd.read_csv(cf, sep="\t", header=None, names=["pfam_acc", "clan", "clan_name", "pfam_name", "desc"])
        clan = dict(zip(c["pfam_acc"], c["clan"].fillna("")))
        for a, nm in zip(c["pfam_acc"], c["pfam_name"]):
            names.setdefault(a, nm)
    return names, clan


def toolkit_map():
    tk = pd.read_csv(os.path.join(ROOT, "data", "toolkit_pfam.tsv"), sep="\t")
    tk = tk[tk["role"].isin(["member", "conditional_member"])]
    return dict(zip(tk["pfam_acc"].str.split(".").str[0], tk["category"]))


def hog_view(ctx, gdom, contrasts, row_of):
    out = {}
    c = ctx.cache
    for k, con in enumerate(contrasts):
        out[con.cid] = ex.hog_domain_events(c.events, c.hogs, c.genes, gdom, con.edges, row_of)
        if (k + 1) % 10 == 0:
            log(f"HOG view {k + 1}/{len(contrasts)}")
    return out


def hog_flags(Y: pd.DataFrame, ref_sets: dict, min_count: int = 2, u_min: float = 0.95) -> pd.DataFrame:
    """Y: domain x branch (gain+dup HOG counts). Flag for column b: Y >= min_count and u(Y_b vs Y over
    ref_sets[b]) >= u_min (mid-rank percentile as in §4.6)."""
    F = pd.DataFrame(False, index=Y.index, columns=list(ref_sets))
    for b, ref in ref_sets.items():
        ref = [r for r in ref if r in Y.columns]
        if not ref or b not in Y.columns:
            continue
        y = Y[b].to_numpy()[:, None]
        R = Y[ref].to_numpy()
        u = ((R < y).sum(1) + 0.5 * (R == y).sum(1) + 0.5) / (R.shape[1] + 1)
        F[b] = (Y[b].to_numpy() >= min_count) & (u >= u_min)
    return F


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--transitions", default="default")
    ap.add_argument("--n-perm", type=int, default=ex.N_PERM)
    ap.add_argument("--seed", type=int, default=ex.SEED)
    ap.add_argument("--no-hog", action="store_true", help="skip the HOG view (no cache needed)")
    ap.add_argument("--no-registry", action="store_true")
    ap.add_argument("--skip-prereg-check", action="store_true")
    a = ap.parse_args(argv)
    if not a.skip_prereg_check and not ex.prereg_ok(ROOT):
        sys.exit("PREREG hash mismatch: aborting")
    ctx = ex.load_context(a.run, a.transitions, need_cache=not a.no_hog)
    od = ex.out_dir(a.run, "pfam_wide")
    tabdir = os.path.join(ROOT, "runs", "analysis", a.run, "pfam_full")
    man = ex.build_pfam_tables(ctx.paths["pfam_full"], ctx.species, tabdir, log=log)
    ann = man["annotated"]
    log(f"Pfam: {len(ann)} annotated species, {len(man['missing'])} missing")
    if len(ann) < 4:
        sys.exit("fewer than 4 species with full Pfam annotation; nothing to do")
    counts = pd.read_parquet(os.path.join(tabdir, "counts.parquet"))
    X = ex.count_matrix(counts, ann)
    names, clan = pfam_names(tabdir)
    tk = toolkit_map()
    X.to_csv(os.path.join(od, "pfam_counts_species_x_domain.tsv.gz"), sep="\t")

    # ---- coverage
    cov = []
    for con in ctx.trans + ctx.ctrls:
        r = con.restrict(ann)
        cov.append(dict(contrast=con.cid, kind=con.kind, branch=con.branch, n_M=len(con.M), n_M_annot=len(r.M),
                        n_S=len(con.S), n_S_annot=len(r.S), n_out=len(con.outside), n_out_annot=len(r.outside),
                        evaluable=ex.evaluable(con, ann)))
    cov = pd.DataFrame(cov)
    ev = dict(zip(cov["contrast"], cov["evaluable"]))
    n_ev_t = sum(ev[c.cid] for c in ctx.trans)
    n_ev_c = sum(ev[c.cid] for c in ctx.ctrls)
    log(f"evaluable contrasts: {n_ev_t}/{len(ctx.trans)} transitions, {n_ev_c}/{len(ctx.ctrls)} controls")
    cov.to_csv(os.path.join(od, "coverage.tsv"), sep="\t", index=False)

    # ---- (b) per-contrast domain tables
    tabs = {con.cid: ex.domain_contrast(X, con) for con in ctx.trans}
    ctab = {con.branch: ex.domain_contrast(X, con) for con in ctx.ctrls}
    log(f"domain contrasts: {len(tabs)} transitions, {len(ctab)} controls")
    D = pd.concat(tabs.values(), ignore_index=True)
    D["pfam_name"] = D["pfam_acc"].map(names)
    D["clan"] = D["pfam_acc"].map(clan)
    D["toolkit"] = D["pfam_acc"].map(tk).fillna("")
    D.to_csv(os.path.join(od, "domain_contrasts.tsv.gz"), sep="\t", index=False)

    def flagmat(tables, fn):
        cols = {k: t.set_index("pfam_acc").pipe(fn) for k, t in tables.items()}
        return pd.DataFrame(cols).fillna(False).astype(bool)

    flag_fns = {
        "exp_or_gain": lambda t: t["status"].isin(["expanded", "gained"]),
        "gained": lambda t: t["status"] == "gained",
        "expanded": lambda t: t["status"] == "expanded",
        "lost": lambda t: t["status"] == "lost",
    }
    tabs_ev = {k: v for k, v in tabs.items() if ev[k]}
    ctab_ev = {k: v for k, v in ctab.items() if ev[f"ctrl:{k}"]}
    if len(tabs_ev) < 2 or len(ctab_ev) < 2:
        sys.exit(f"too few evaluable contrasts for convergence ({len(tabs_ev)} transitions, {len(ctab_ev)} controls); "
                 "wait for the Pfam annotation")
    obs_flags = {k: flagmat(tabs_ev, f) for k, f in flag_fns.items()}
    ctrl_flags = {k: flagmat(ctab_ev, f) for k, f in flag_fns.items()}
    cf = pd.concat({k: v.stack().rename("flag") for k, v in ctrl_flags.items()}, names=["flag_type"]).reset_index()
    cf.columns = ["flag_type", "pfam_acc", "control", "flag"]
    cf[cf["flag"]].to_csv(os.path.join(od, "control_flags.tsv.gz"), sep="\t", index=False)

    # ---- heatmaps
    sm = ex.status_matrix(tabs)
    sm = sm[(sm != 0).any(axis=1)]
    sm.insert(0, "pfam_name", sm.index.map(names))
    sm.to_csv(os.path.join(od, "heatmap_status.tsv"), sep="\t", index_label="pfam_acc")
    lf = pd.DataFrame({k: t.set_index("pfam_acc")["log2fc"] for k, t in tabs.items()}).reindex(sm.index)
    lf.insert(0, "pfam_name", lf.index.map(names))
    lf.to_csv(os.path.join(od, "heatmap_log2fc.tsv"), sep="\t", index_label="pfam_acc")

    # ---- (c) convergence
    sets = {t: ctx.match.sets[t] for t in tabs_ev}
    glob_rows, conv = [], {}
    for k in ("exp_or_gain", "gained", "lost"):
        per, glob, info = ex.convergence_test(obs_flags[k], ctrl_flags[k], sets, n_perm=a.n_perm, seed=a.seed)
        per["pfam_name"] = per["pfam_acc"].map(names)
        per["clan"] = per["pfam_acc"].map(clan)
        per["toolkit"] = per["pfam_acc"].map(tk).fillna("")
        per.to_csv(os.path.join(od, f"convergence_{k}.tsv"), sep="\t", index=False)
        glob.insert(0, "flag", k)
        glob["n_controls"] = info["n_controls"]
        glob_rows.append(glob)
        conv[k] = (per, glob, info)
        log(f"convergence {k}: {info}")

    # ---- HOG view
    hog_summary = {}
    if not a.no_hog:
        c = ctx.cache
        pdm = pd.read_parquet(os.path.join(tabdir, "protein_domains.parquet"))
        gdom = ex.gene_domain_codes(c.genes, pdm)
        row_of = pd.Series(np.arange(len(c.hogs)), index=c.hogs["hog_id"].to_numpy())
        hv = hog_view(ctx, gdom, ctx.trans, row_of)
        hvc = hog_view(ctx, gdom, ctx.ctrls, row_of)
        H = pd.concat({k: v for k, v in hv.items()}, names=["transition"]).reset_index()
        H["pfam_name"] = H["pfam_acc"].map(names)
        H.to_csv(os.path.join(od, "hog_domain_events.tsv.gz"), sep="\t", index=False)
        # gain + dup HOG counts: domain x branch (transitions by stem branch, controls by branch)
        Y = pd.DataFrame({**{con.branch: (hv[con.cid]["gain"] + hv[con.cid]["dup"]) for con in ctx.trans},
                          **{con.branch: (hvc[con.cid]["gain"] + hvc[con.cid]["dup"]) for con in ctx.ctrls}}).fillna(0)
        ref_t = {con.branch: ctx.match.sets[con.cid] for con in ctx.trans}
        ref_c = ex.control_self_sets(ctx.bt, ctx.pool)
        Ft = hog_flags(Y, ref_t)
        Ft.columns = [con.cid for con in ctx.trans]
        Fc = hog_flags(Y, ref_c)
        Ft = Ft[[t for t in Ft.columns if ev[t]]]
        Fc = Fc[[b for b in Fc.columns if ev[f"ctrl:{b}"]]]
        per, glob, info = ex.convergence_test(Ft, Fc, sets, n_perm=a.n_perm, seed=a.seed)
        per["pfam_name"] = per["pfam_acc"].map(names)
        per["toolkit"] = per["pfam_acc"].map(tk).fillna("")
        per.to_csv(os.path.join(od, "convergence_hog_expanded.tsv"), sep="\t", index=False)
        glob.insert(0, "flag", "hog_expanded")
        glob["n_controls"] = info["n_controls"]
        glob_rows.append(glob)
        conv["hog_expanded"] = (per, glob, info)
        hm = pd.DataFrame({con.cid: Y[con.branch] for con in ctx.trans})
        hm = hm[(hm > 0).any(axis=1)]
        hm.insert(0, "pfam_name", hm.index.map(names))
        hm.to_csv(os.path.join(od, "heatmap_hog.tsv"), sep="\t", index_label="pfam_acc")
        # concordance orthology-free vs HOG view, per transition
        for con in [x for x in ctx.trans if ev[x.cid]]:
            of = set(obs_flags["exp_or_gain"].index[obs_flags["exp_or_gain"][con.cid]])
            hf = set(Ft.index[Ft[con.cid]])
            hog_summary[con.cid] = dict(n_pfam_exp_or_gain=len(of), n_hog_expanded=len(hf), n_both=len(of & hf))
        log("HOG view done")
    G = pd.concat(glob_rows, ignore_index=True)
    G.to_csv(os.path.join(od, "convergence_global.tsv"), sep="\t", index=False)

    # ---- toolkit summary per transition (robustness check for H4/H5: orthology-free Pfam-count version)
    rows = []
    for t, tab in tabs.items():
        tt = tab.assign(cat=tab["pfam_acc"].map(tk)).dropna(subset=["cat"])
        for cat, g in tt.groupby("cat"):
            rows.append(dict(transition=t, category=cat, n_domains_present=int((g["frac_M"] > 0).sum()),
                             n_expanded=int((g["status"] == "expanded").sum()),
                             n_gained=int((g["status"] == "gained").sum()),
                             n_contracted=int((g["status"] == "contracted").sum()),
                             n_lost=int((g["status"] == "lost").sum()),
                             expanded=";".join(g.loc[g["status"].isin(["expanded", "gained"]), "pfam_acc"]
                                               .map(lambda x: names.get(x, x)))))
    pd.DataFrame(rows).to_csv(os.path.join(od, "toolkit_by_transition.tsv"), sep="\t", index=False)

    # ---- per-transition counts and summary
    pt = []
    for t, tab in tabs.items():
        vc = tab["status"].value_counts()
        r = cov.set_index("contrast").loc[t]
        pt.append(dict(transition=t, evaluable=bool(r["evaluable"]), n_M_annot=int(r["n_M_annot"]), n_S_annot=int(r["n_S_annot"]),
                       n_out_annot=int(r["n_out_annot"]), n_domains=int(len(tab)),
                       n_tested_mwu=int(tab["p_mwu"].notna().sum()),
                       n_expanded=int(vc.get("expanded", 0)), n_expanded_by_test=int(tab["expanded_by_test"].sum()),
                       n_gained=int(vc.get("gained", 0)), n_contracted=int(vc.get("contracted", 0)),
                       n_lost=int(vc.get("lost", 0)), **hog_summary.get(t, {})))
    pt = pd.DataFrame(pt)
    pt.to_csv(os.path.join(od, "per_transition.tsv"), sep="\t", index=False)
    cm = pd.DataFrame([dict(control=k, n_expanded=int((v["status"] == "expanded").sum()),
                            n_gained=int((v["status"] == "gained").sum()), n_lost=int((v["status"] == "lost").sum()))
                       for k, v in ctab.items()])
    cm.to_csv(os.path.join(od, "per_control.tsv"), sep="\t", index=False)

    # ---- registry
    inputs = [os.path.join(tabdir, "manifest.json"), ctx.paths["transitions"]]
    if not a.no_registry:
        for _, g in G[G["scheme"] == "permute"].iterrows():
            gm = G[(G["scheme"] == "matched") & (G["flag"] == g["flag"]) & (G["k_min"] == g["k_min"])]
            pm = f"; p_matched_sensitivity={gm['p_perm'].iloc[0]:.3g}" if len(gm) else ""
            bs.append_registry(ctx.paths["registry"], ex.registry_row(
                f"Q2 Pfam-wide convergence ({g['flag']})",
                f"n Pfam domains flagged '{g['flag']}' in >= {g['k_min']} transitions vs pseudo-transitions "
                f"(control-pool child clades), permute null; estimate = observed n, null mean in notes",
                n=n_ev_t, k=int(g["n_controls"]), estimate=float(g["n_domains_obs"]), p=float(g["p_perm"]),
                seed=a.seed, n_draws=a.n_perm, scheme="permute",
                notes=f"run={a.run}; null_mean={g['null_mean']:.3g}; null_q95={g['null_q95']:.3g}{pm}; pool={ctx.pool_label}; "
                      f"evaluable transitions={n_ev_t}/{len(ctx.trans)}, controls={n_ev_c}/{len(ctx.ctrls)}; "
                      f"pfam_annotated={len(ann)}/{len(ctx.species)}" + ("; DRY RUN prok" if a.run == "prok" else ""),
                inputs=inputs))
    summ = dict(run=a.run, transitions=ctx.tids, pool=ctx.pool_label, n_pool=len(ctx.pool),
                n_evaluable_transitions=n_ev_t, n_evaluable_controls=n_ev_c,
                n_species=len(ctx.species), n_pfam_annotated=len(ann), pfam_missing=man["missing"],
                n_domains_total=int(X.shape[1]), per_transition=pt.to_dict("records"),
                convergence_global=G.to_dict("records"),
                top_convergent=conv["exp_or_gain"][0].head(15)[["pfam_acc", "pfam_name", "k_obs", "k_null_mean",
                                                                 "p_perm", "q", "transitions"]].to_dict("records"),
                top_lost=conv["lost"][0].head(10)[["pfam_acc", "pfam_name", "k_obs", "k_null_mean", "p_perm", "q",
                                                   "transitions"]].to_dict("records"),
                control_flag_means=cm[["n_expanded", "n_gained", "n_lost"]].mean().round(2).to_dict())
    ex.write_summary(od, summ)
    log(f"done -> {os.path.relpath(od, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
