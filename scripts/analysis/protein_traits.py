#!/usr/bin/env python
"""EXPLORATORY. Protein length and disorder, multicellular vs sister species (DESIGN §4.5 'Protein traits';
H9 companion, H10).

    python scripts/analysis/protein_traits.py --run {euk,prok} [--transitions default]

Traits: protein length (from data/proteomes/<run>/<code>.fa; cached in runs/analysis/<run>/explore/) and
metapredict disorder (frac_disordered; data/annotations/disorder/<run>/<code>.disorder.tsv.gz, finished species
only -- missing species are dropped from the disorder tests and counted).
Per transition (M clade vs sister species S; HOGs = ancestral genes at P from runs/analysis/<run>/):
  length_sc     HOGs present in >= 50 % of M and of S species and single-copy in every M U S species that has
                them; per HOG d = log2(median length M / median length S); Wilcoxon signed-rank (two-sided) +
                paired Hodges-Lehmann with 95 % CI.
  disorder_sc   the same on species with disorder data; d = median frac_disordered M - S.
  tf_disorder_sc / tf_disorder_all   H10: TF-DBD HOGs at P (>= 50 % of members carry a frozen TF-DBD domain,
                C2H2 excluded); single-copy, and any copy number (per-side median over all member genes).
  whole_length, whole_disorder (+ _genome: genome-only species)   secondary, confounded by data type: per-species
                median trait, M vs S species, Mann-Whitney (two-sided) + two-sample Hodges-Lehmann.
  tf_disorder_excess (+ _genome)   per species median disorder of TF-DBD proteins minus the proteome median.
For the paired single-copy tests the same statistic is computed on every control pseudo-transition (control pool
child clades); u = percentile of the transition's HL among its matched controls (§4.6 mid-rank).
BH across all tests of this script (q). Registry: one row per test with a p-value.

Outputs: results/exploratory/<run>/protein_traits/: tests.tsv, hog_pairs.tsv.gz, per_species.tsv,
controls.tsv, summary.json.
"""
from __future__ import annotations

import argparse
import glob
import json
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
    print(f"[protein_traits {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def protein_lengths(run, species, cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    f = os.path.join(cache_dir, "protein_lengths.parquet")
    if os.path.exists(f):
        L = pd.read_parquet(f)
        if set(species) <= set(L["species"].astype(str).unique()):
            return L
    rows = []
    for sp in species:
        fa = os.path.join(ROOT, "data", "proteomes", run, f"{sp}.fa")
        if not os.path.exists(fa):
            continue
        name, n = None, 0
        names, lens = [], []
        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        names.append(name)
                        lens.append(n)
                    name, n = line[1:].split()[0], 0
                else:
                    n += len(line.strip())
        if name is not None:
            names.append(name)
            lens.append(n)
        rows.append(pd.DataFrame({"gene": names, "species": sp, "length": lens}))
    L = pd.concat(rows, ignore_index=True)
    L.to_parquet(f, index=False)
    return L


def disorder_table(run, species):
    d = os.path.join(ROOT, "data", "annotations", "disorder", run)
    rows, have = [], []
    for sp in species:
        f = os.path.join(d, f"{sp}.disorder.tsv.gz")
        if os.path.exists(f) and os.path.exists(os.path.join(d, f"{sp}.DONE.json")):
            t = pd.read_csv(f, sep="\t", usecols=["protein", "frac_disordered"])
            t["species"] = sp
            rows.append(t)
            have.append(sp)
    if not rows:
        return pd.DataFrame(columns=["protein", "frac_disordered", "species"]), have
    return pd.concat(rows, ignore_index=True), have


class Members:
    """gene membership at a level, cached per level."""

    def __init__(self, cache):
        self.c = cache
        self._m = {}
        self.genes = cache.genes[["gene", "species"]].copy()
        self.genes["species"] = self.genes["species"].astype(str)

    def at(self, level):
        if level not in self._m:
            li = self.c.tree.index[level]
            r = ex.gene_hog_at(self.c.hogs, self.c.genes, li)
            k = r >= 0
            self._m[level] = pd.DataFrame({"hog_row": r[k], "gene": self.genes["gene"].to_numpy()[k],
                                           "species": self.genes["species"].to_numpy()[k]})
        return self._m[level]


def paired(memb, M, S, trait, kind, hog_rows=None, single_copy=True):
    """per-HOG differences and the paired test for one contrast."""
    m = memb if hog_rows is None else memb[memb["hog_row"].isin(hog_rows)]
    sc = ex.single_copy_hogs(m, M, S, single_copy=single_copy)
    if sc.empty:
        return pd.DataFrame(columns=["hog_row", "M", "S", "d"]), dict(n=0, p=np.nan, hl=np.nan, hl_lo=np.nan,
                                                                     hl_hi=np.nan, median_diff=np.nan)
    med = ex.per_hog_side_medians(sc, trait).dropna()
    if kind == "log2ratio":
        med = med[(med["M"] > 0) & (med["S"] > 0)]
        d = np.log2(med["M"] / med["S"])
    else:
        d = med["M"] - med["S"]
    med["d"] = d
    return med.reset_index(), ex.paired_test(d.to_numpy())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--transitions", default="default")
    ap.add_argument("--no-controls", action="store_true")
    ap.add_argument("--no-registry", action="store_true")
    ap.add_argument("--skip-prereg-check", action="store_true")
    a = ap.parse_args(argv)
    if not a.skip_prereg_check and not ex.prereg_ok(ROOT):
        sys.exit("PREREG hash mismatch: aborting")
    ctx = ex.load_context(a.run, a.transitions, need_cache=True)
    c = ctx.cache
    od = ex.out_dir(a.run, "protein_traits")
    sp_all = ctx.species
    L = protein_lengths(a.run, sp_all, os.path.join(ROOT, "runs", "analysis", a.run, "explore"))
    length = L.set_index("gene")["length"].astype(float)
    D, have_dis = disorder_table(a.run, sp_all)
    dis = D.set_index("protein")["frac_disordered"].astype(float)
    log(f"lengths for {L['species'].nunique()} species; disorder for {len(have_dis)}/{len(sp_all)}")
    mem = Members(c)
    tree_names = set(c.tree.names)
    tf_rows = set(np.flatnonzero(c.hog_classes["in_tf_dbd"].to_numpy()))
    gc = c.gene_classes
    tf_genes = set(gc.loc[gc["tf_dbd"].to_numpy(), "gene"])
    genome = {s for s in sp_all if ctx.dtype.get(s, "genome") == "genome"}

    # per-species whole-proteome table
    ps = L.groupby("species")["length"].median().rename("median_length").to_frame()
    ps["n_proteins"] = L.groupby("species").size()
    if len(D):
        ps["median_disorder"] = D.groupby("species")["frac_disordered"].median()
        Dt = D[D["protein"].isin(tf_genes)]
        ps["n_tf"] = Dt.groupby("species").size()
        ps["median_tf_disorder"] = Dt.groupby("species")["frac_disordered"].median()
        ps["tf_disorder_excess"] = ps["median_tf_disorder"] - ps["median_disorder"]
    else:
        for col in ("median_disorder", "n_tf", "median_tf_disorder", "tf_disorder_excess"):
            ps[col] = np.nan
    ps.index = ps.index.astype(str)
    ps["data_type"] = [ctx.dtype.get(s, "genome") for s in ps.index]
    ps.to_csv(os.path.join(od, "per_species.tsv"), sep="\t", index_label="species")

    def sc_tests(con, with_tf):
        if con.P not in tree_names:
            log(f"{con.cid}: P={con.P} not in the Ancestra tree; skipped")
            return [], []
        memb = mem.at(con.P)
        out, pairs = [], []
        specs = [("length_sc", length, "log2ratio", con.M, con.S, None, True)]
        Md = [s for s in con.M if s in set(have_dis)]
        Sd = [s for s in con.S if s in set(have_dis)]
        if Md and Sd:
            specs.append(("disorder_sc", dis, "difference", Md, Sd, None, True))
            if with_tf:
                specs.append(("tf_disorder_sc", dis, "difference", Md, Sd, tf_rows, True))
                specs.append(("tf_disorder_all", dis, "difference", Md, Sd, tf_rows, False))
        if with_tf:
            specs.append(("tf_length_all", length, "log2ratio", con.M, con.S, tf_rows, False))
        for name, tr_, kind, M, S, rows, scp in specs:
            pr, t = paired(memb, M, S, tr_, kind, rows, scp)
            out.append(dict(contrast=con.cid, test=name, subset="all", n_M_species=len(M), n_S_species=len(S),
                            n_M_missing_trait=len(con.M) - len(M), n_S_missing_trait=len(con.S) - len(S), **t))
            if len(pr):
                pairs.append(pr.assign(contrast=con.cid, test=name))
        return out, pairs

    rows, pairs = [], []
    for con in ctx.trans:
        o, p = sc_tests(con, True)
        rows += o
        pairs += p
        # whole-proteome and per-species TF excess
        for col, test in (("median_length", "whole_length"), ("median_disorder", "whole_disorder"),
                          ("tf_disorder_excess", "tf_disorder_excess")):
            for subset, keep in (("all", None), ("genome", genome)):
                M = [s for s in con.M if (keep is None or s in keep)]
                S = [s for s in con.S if (keep is None or s in keep)]
                x = ps[col].reindex(M).to_numpy(float)
                y = ps[col].reindex(S).to_numpy(float)
                if col == "median_length":
                    x, y = np.log2(x), np.log2(y)
                t = ex.two_sample_test(x, y)
                rows.append(dict(contrast=con.cid, test=test, subset=subset, n=t["n_M"] + t["n_S"],
                                 n_M_species=t["n_M"], n_S_species=t["n_S"], hl=t["hl"], hl_lo=t["hl_lo"],
                                 hl_hi=t["hl_hi"], p=t["p"], median_diff=t["median_M"] - t["median_S"]))
        log(f"{con.cid} done")
    T = pd.DataFrame(rows)

    # controls for the paired single-copy statistics
    C = pd.DataFrame()
    if not a.no_controls:
        crow = []
        for k, con in enumerate(ctx.ctrls):
            o, _ = sc_tests(con, False)
            crow += o
        C = pd.DataFrame(crow)
        C["branch"] = C["contrast"].str.replace("ctrl:", "", regex=False)
        C.to_csv(os.path.join(od, "controls.tsv"), sep="\t", index=False)
        T["u_vs_controls"] = np.nan
        T["median_hl_controls"] = np.nan
        T["k_controls"] = 0
        for i, r in T.iterrows():
            if r["test"] not in ("length_sc", "disorder_sc") or not np.isfinite(r["hl"]):
                continue
            N = ctx.match.sets[r["contrast"]]
            cv = C[(C["test"] == r["test"]) & C["branch"].isin(N)]["hl"].dropna().to_numpy(float)
            if len(cv):
                T.at[i, "u_vs_controls"] = bs.u_score(r["hl"], cv)
                T.at[i, "median_hl_controls"] = float(np.median(cv))
                T.at[i, "k_controls"] = len(cv)
    T["q_bh"] = ex.bh(T["p"].to_numpy())
    T = T.rename(columns={"contrast": "transition"})
    T.to_csv(os.path.join(od, "tests.tsv"), sep="\t", index=False)
    if pairs:
        pd.concat(pairs, ignore_index=True).to_csv(os.path.join(od, "hog_pairs.tsv.gz"), sep="\t", index=False)

    # across transitions (descriptive)
    acr = []
    for (test, subset), g in T.groupby(["test", "subset"]):
        h = g["hl"].dropna()
        acr.append(dict(test=test, subset=subset, n_transitions=int(len(h)), n_positive=int((h > 0).sum()),
                        median_hl=float(h.median()) if len(h) else np.nan))
    A = pd.DataFrame(acr)
    A.to_csv(os.path.join(od, "across_transitions.tsv"), sep="\t", index=False)

    if not a.no_registry:
        hyp = {"length_sc": "Q5 protein length (paired single-copy HOGs)",
               "disorder_sc": "Q5 protein disorder (paired single-copy HOGs)",
               "tf_disorder_sc": "H10 TF-DBD disorder (paired single-copy HOGs)",
               "tf_disorder_all": "H10 TF-DBD disorder (paired TF HOGs, any copy number)",
               "tf_length_all": "H10 TF-DBD length (paired TF HOGs, any copy number)",
               "whole_length": "Q5 whole-proteome median length (secondary, data-type confounded)",
               "whole_disorder": "Q5 whole-proteome median disorder (secondary, data-type confounded)",
               "tf_disorder_excess": "H10 TF disorder minus proteome disorder per species (secondary)"}
        for _, r in T[T["p"].notna()].iterrows():
            unit = "log2 ratio" if "length" in r["test"] else "difference in frac_disordered"
            paired_ = r["test"].endswith("_sc") or r["test"].endswith("_all")
            bs.append_registry(ctx.paths["registry"], ex.registry_row(
                hyp[r["test"]] + ("" if r["subset"] == "all" else f" [{r['subset']}-only]"),
                ("Wilcoxon signed-rank over HOGs + paired Hodges-Lehmann" if paired_ else
                 "Mann-Whitney over species + two-sample Hodges-Lehmann") + f" ({unit}, M minus S); {r['transition']}",
                n=int(r["n"]), estimate=r["hl"], ci=(r["hl_lo"], r["hl_hi"]), p=r["p"],
                notes=f"run={a.run}; transition={r['transition']}; q_BH(within protein_traits)={r['q_bh']:.3g}"
                      + (f"; u_vs_controls={r['u_vs_controls']:.3g}" if np.isfinite(r.get("u_vs_controls", np.nan)) else "")
                      + f"; disorder species {len(have_dis)}/{len(sp_all)}"
                      + ("; DRY RUN prok" if a.run == "prok" else ""),
                inputs=[ctx.paths["transitions"]]))
    summ = dict(run=a.run, pool=ctx.pool_label, n_species=len(sp_all), n_disorder_species=len(have_dis),
                tests=T.round(5).replace({np.nan: None}).to_dict("records"),
                across=A.round(5).replace({np.nan: None}).to_dict("records"))
    ex.write_summary(od, summ)
    log(f"done -> {os.path.relpath(od, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
