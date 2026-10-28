#!/usr/bin/env python
"""EXPLORATORY / DESCRIPTIVE. Co-option index C_i (DESIGN §4.3, §4.6 'Co-option index (descriptive)').

    python scripts/analysis/cooption.py --run {euk,prok} [--transitions default]

For each transition i and every control pseudo-transition (child clade of the control pool, as in
pfam_wide.py), F_i = families with log2((mean_M + 1)/(mean_S + 1)) >= 1 present in >= 50 % of M species;
C_i = fraction of F_i present at P (>= 1 genome or >= 2 transcriptome species carrying it among S U outgroups,
i.e. every species outside M); C_i^parent = the same among species outside P (present at P's parent).
Reported: C_i, median C over the matched controls N_i, C_i - median(C_{N_i}), percentile u (§4.6), and
T = mean_i u. No p-value (descriptive by design: presence in S implies presence at P).
Versions:
  pfam_expanded      all Pfam domains (full Pfam, data/annotations/pfam/<run>; annotated species only)
  toolkit_expanded   toolkit Pfam domains (data/toolkit_pfam.tsv members) with the expansion rule
  toolkit_present    toolkit Pfam domains present in >= 50 % of M (no expansion rule; §4.6 'toolkit-based version')
  hog_toolkit        Ancestra families (top-level HOGs, family_id) whose genes carry a toolkit domain in >= 50 %;
                     presence at P also when the family's gain on a branch below P was split-repaired
                     (results/gainladder/<run>/gains_classified.tsv.gz; without it: raw presence, noted)
The toolkit versions use the pre-registered Pfam subset (cache pfam_hits) when full Pfam is missing.

Outputs: results/exploratory/<run>/cooption/: cooption_index.tsv (version x contrast), cooption_families.tsv.gz
(families in F_i with presence flags), per_transition.tsv (vs matched controls), summary.json; one registry row per
version (estimate = mean_i [C_i - median C_{N_i}], p = NA).
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

from mcell import branchstats as bs  # noqa: E402
from mcell import explore as ex  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[cooption {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def toolkit_accs():
    tk = pd.read_csv(os.path.join(ROOT, "data", "toolkit_pfam.tsv"), sep="\t")
    tk = tk[tk["role"].isin(["member", "conditional_member"])]
    return set(tk["pfam_acc"].str.split(".").str[0])


def family_matrix(cache, toolkit_cols=ex.TOOLKIT_CLASSES):
    """species x family_id copy numbers for toolkit families (>= 50 % of genes carry a toolkit domain)."""
    gc = cache.gene_classes
    tkg = gc[[c for c in toolkit_cols if c in gc.columns]].any(axis=1).to_numpy()
    g = cache.genes[["gene", "species", "family_id"]].copy()
    g["tk"] = pd.Series(tkg, index=gc["gene"]).reindex(g["gene"]).fillna(False).to_numpy()
    fam = g.groupby("family_id")["tk"].agg(["mean", "size"])
    keep = set(fam.index[(fam["mean"] >= 0.5) & (fam["size"] >= 2)])
    g = g[g["family_id"].isin(keep)]
    X = g.groupby(["species", "family_id"], observed=True).size().unstack(fill_value=0)
    X.index = X.index.astype(str)
    return X


def repaired_families(path, tree, P):
    if not os.path.exists(path):
        return None
    g = pd.read_csv(path, sep="\t", usecols=lambda c: c in ("family_id", "branch_child", "split_repaired"))
    g = g[g["split_repaired"].astype(bool)]
    below = {n.name for n in (tree & P).traverse()} if P in {n.name for n in tree.traverse()} else set()
    return set(g.loc[g["branch_child"].isin(below), "family_id"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--transitions", default="default")
    ap.add_argument("--no-registry", action="store_true")
    ap.add_argument("--skip-prereg-check", action="store_true")
    a = ap.parse_args(argv)
    if not a.skip_prereg_check and not ex.prereg_ok(ROOT):
        sys.exit("PREREG hash mismatch: aborting")
    ctx = ex.load_context(a.run, a.transitions, need_cache=True)
    od = ex.out_dir(a.run, "cooption")
    tabdir = os.path.join(ROOT, "runs", "analysis", a.run, "pfam_full")
    man = ex.build_pfam_tables(ctx.paths["pfam_full"], ctx.species, tabdir, log=log)
    ann = man["annotated"]
    tk = toolkit_accs()
    versions = {}
    notes = {}
    if len(ann) >= 4:
        X = ex.count_matrix(pd.read_parquet(os.path.join(tabdir, "counts.parquet")), ann)
        versions["pfam_expanded"] = (X, None, "expanded")
        versions["toolkit_expanded"] = (X, tk, "expanded")
        versions["toolkit_present"] = (X, tk, "present")
        notes["pfam"] = f"full Pfam, {len(ann)}/{len(ctx.species)} species annotated"
    else:
        h = ctx.cache.hits
        c = h[["protein", "pfam_acc", "species"]].drop_duplicates().groupby(["species", "pfam_acc"]).size() \
            .rename("n_proteins").reset_index()
        Xp = ex.count_matrix(c, ctx.species)
        versions["toolkit_expanded"] = (Xp, tk, "expanded")
        versions["toolkit_present"] = (Xp, tk, "present")
        notes["pfam"] = "full Pfam missing: toolkit versions on the pre-registered subset; pfam_expanded skipped"
    Xf = family_matrix(ctx.cache)
    versions["hog_toolkit"] = (Xf, None, "expanded")
    ladder = ctx.paths["ladder_gains"]
    notes["hog"] = "split repair from gain ladder" if os.path.exists(ladder) else "gain ladder missing: raw presence"
    log(f"versions: {list(versions)}; {notes}")

    rows, fam_rows = [], []
    for v, (X, fams, rule) in versions.items():
        for con in ctx.trans + ctx.ctrls:
            extra = ()
            if v == "hog_toolkit":
                rp = repaired_families(ladder, ctx.tree, con.P)
                extra = rp or ()
            r = ex.coption_index(X, con, ctx.dtype, families=fams, rule=rule, extra_present=extra)
            F = r.pop("F")
            ok = ex.evaluable(con, X.index) if v != "hog_toolkit" else True
            if not ok:
                r.update(C=np.nan, C_Pparent=np.nan)
                F = []
            rows.append(dict(version=v, kind=con.kind, branch=con.branch, evaluable=ok, **r))
            if con.kind == "transition" and F:
                c2 = con.restrict(X.index)
                pres = ex.present_rule(X[F], c2.outside, ctx.dtype)
                presP = ex.present_rule(X[F], c2.outside_P, ctx.dtype)
                for f in F:
                    fam_rows.append(dict(version=v, transition=con.cid, family=f, present_at_P=bool(pres[f]) or f in set(extra),
                                         present_at_P_parent=bool(presP[f]),
                                         mean_M=float(X.loc[c2.M, f].mean()), mean_S=float(X.loc[c2.S, f].mean())
                                         if c2.S else np.nan))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(od, "cooption_index.tsv"), sep="\t", index=False)
    pd.DataFrame(fam_rows).to_csv(os.path.join(od, "cooption_families.tsv.gz"), sep="\t", index=False)

    # vs matched controls
    out = []
    for v in versions:
        rv = R[R["version"] == v]
        cval = rv[rv["kind"] == "control"].set_index("branch")["C"]
        cvalP = rv[rv["kind"] == "control"].set_index("branch")["C_Pparent"]
        for _, r in rv[rv["kind"] == "transition"].iterrows():
            N = [x for x in ctx.match.sets[r["contrast"]] if np.isfinite(cval.get(x, np.nan))]
            cn = cval.reindex(N).to_numpy(float)
            med = float(np.median(cn)) if len(cn) else np.nan
            u = bs.u_score(r["C"], cn) if (len(cn) and np.isfinite(r["C"])) else np.nan
            cnp = cvalP.reindex(N).dropna().to_numpy(float)
            out.append(dict(version=v, transition=r["contrast"], n_F=r["n_F"], n_P=r["n_P"], C=r["C"],
                            C_Pparent=r["C_Pparent"], k_controls=len(cn), median_C_controls=med,
                            C_minus_median=r["C"] - med if np.isfinite(med) else np.nan,
                            median_C_Pparent_controls=float(np.median(cnp)) if len(cnp) else np.nan,
                            u=u, n_M=r["n_M"], n_S=r["n_S"]))
    PT = pd.DataFrame(out)
    PT.to_csv(os.path.join(od, "per_transition.tsv"), sep="\t", index=False)
    summ = dict(run=a.run, pool=ctx.pool_label, notes=notes, per_transition=PT.round(4).to_dict("records"), versions={})
    for v, g in PT.groupby("version"):
        d = g["C_minus_median"].dropna()
        uu = g["u"].dropna()
        summ["versions"][v] = dict(n=int(len(d)), mean_C=float(g["C"].mean()), mean_diff=float(d.mean()) if len(d) else None,
                                   T_mean_u=float(uu.mean()) if len(uu) else None)
        if not a.no_registry and len(d):
            bs.append_registry(ctx.paths["registry"], ex.registry_row(
                f"Q3 co-option index ({v})",
                "C_i = fraction of M-expanded (or M-present) families present at P; estimate = mean_i [C_i - median C(N_i)]",
                n=int(len(d)), k=int(g["k_controls"].median()), estimate=float(d.mean()),
                ci=(float(d.min()), float(d.max())),
                notes=f"run={a.run}; descriptive, no p; T=mean u={uu.mean():.3g}; ci columns = min/max over transitions; "
                      f"pool={ctx.pool_label}; {notes.get('pfam') if v != 'hog_toolkit' else notes['hog']}"
                      + ("; DRY RUN prok" if a.run == "prok" else ""),
                inputs=[ctx.paths["transitions"]]))
    ex.write_summary(od, summ)
    log(f"done -> {os.path.relpath(od, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
