#!/usr/bin/env python
"""Validation 2b: how well do Ancestra HOGs recover per-branch gene-content change?

For each fixed simulation (``ancestra.simulate``, seeds below):

1. simulate and write proteomes (deterministic; re-simulated in memory for the truth);
2. truth = per-branch events from the true gene trees (:mod:`mcell.simtruth`); it is
   cross-checked against :func:`mcell.hogs.parse_hogs` on the true OrthoXML (must be equal);
3. ``ancestra run --preset {fast,accurate}`` (skipped if the output exists);
4. inferred = :func:`mcell.hogs.parse_hogs` on Ancestra's ``hogs.orthoxml``;
5. metrics per event type across branches (Pearson r, Spearman rho, totals, MAE), exact
   recovery of the ancestral gene sets at internal nodes, and per-gene family-origin accuracy.

Writes ``results/validation_hogparser/{summary.md, branch_comparison.tsv, metrics.tsv}``.
Large intermediate files go to ``runs/validation_hogparser/`` (git-ignored).

    python scripts/validate_hogparser_sim.py [--threads 4] [--presets fast,accurate]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import io  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

from ancestra.orthoxml import write_orthoxml  # noqa: E402
from ancestra.simulate import SimParams, simulate, write_simulation  # noqa: E402
from mcell.hogs import parse_hogs  # noqa: E402
from mcell.simtruth import true_ancestral_genes, true_branch_events  # noqa: E402

# fixed list of simulations (seeds are part of the protocol; do not change after results exist)
SIMS = {
    "S1_12sp_default_hgt": (11, dict(n_species=12, n_genes=1000, hgt_rate=0.01)),
    "S2_16sp_highturnover_wgd": (23, dict(n_species=16, n_genes=1500, dup_rate=0.3, loss_rate=0.3,
                                          wgd_branch="random")),
    "S3_24sp_default": (7, dict(n_species=24, n_genes=1000)),
}
EVENTS = ["gains", "losses", "dups_total", "n_child"]


def run_ancestra(sim_dir: str, out: str, preset: str, threads: int) -> float:
    if os.path.exists(os.path.join(out, "hogs.orthoxml")):
        return float("nan")
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), NUMBA_NUM_THREADS=str(threads),
               PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""))
    t0 = time.time()
    with open(out + ".log", "w") as log:
        subprocess.run([os.path.join(os.path.dirname(sys.executable), "ancestra"), "run",
                        "--proteomes", os.path.join(sim_dir, "proteomes"),
                        "--tree", os.path.join(sim_dir, "species_tree.nwk"), "--out", out,
                        "--threads", str(threads), "--preset", preset],
                       check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
    return time.time() - t0


def branch_table(df: pd.DataFrame) -> pd.DataFrame:
    d = df.set_index("branch_child").copy()
    d["dups_total"] = d["dups"] + d["dups_at_origin"]
    return d[["branch_parent", "child_is_leaf", "gains", "losses", "dups_total", "n_child",
              "leaf_specific_genes", "gains_ge2_child_clades"]]


def corr_metrics(t: pd.Series, p: pd.Series) -> dict:
    t, p = t.astype(float), p.astype(float)
    ok = len(t) > 2 and t.std() > 0 and p.std() > 0
    return dict(n_branches=len(t), true_total=int(t.sum()), inferred_total=int(p.sum()),
                ratio=round(p.sum() / t.sum(), 3) if t.sum() else np.nan,
                pearson=round(pearsonr(t, p)[0], 3) if ok else np.nan,
                spearman=round(spearmanr(t, p)[0], 3) if ok else np.nan,
                mae=round(float((t - p).abs().mean()), 2))


def origin_of_genes(res) -> pd.Series:
    """gene -> species-tree node where its family originates (gain branch child)."""
    g = res.events[res.events.event == "gain"].drop_duplicates("family_id")
    fam2origin = dict(zip(g["family_id"], g["branch_child"]))
    return res.genes.set_index("gene")["family_id"].map(fam2origin)


def origin_accuracy(truth_res, inf_res) -> dict:
    T = truth_res.tree
    to, io_ = origin_of_genes(truth_res), origin_of_genes(inf_res).reindex(origin_of_genes(truth_res).index)
    exact = younger = older = other = 0
    for a, b in zip(to.to_numpy(), io_.to_numpy()):
        if a == b:
            exact += 1
        else:
            ia, ib = T.index[a], T.index[b]
            if T.is_anc(ia, ib):
                younger += 1
            elif T.is_anc(ib, ia):
                older += 1
            else:
                other += 1
    n = len(to)
    return dict(genes=n, origin_exact=round(exact / n, 4), origin_too_young=round(younger / n, 4),
                origin_too_old=round(older / n, 4), origin_other=round(other / n, 4))


def hog_set_recovery(truth_sets: dict, inf_res) -> dict:
    """Exact recovery of ancestral gene sets at internal nodes, pooled over nodes."""
    tp = nt = ni = 0
    for i in inf_res.tree.internal:
        lv = inf_res.tree.names[i]
        m = inf_res.members(lv)
        inf = {frozenset(x[1]) for x in m.groupby("hog_id")["gene"]}
        tr = truth_sets.get(lv, set())
        tp += len(inf & tr)
        nt += len(tr)
        ni += len(inf)
    return dict(hog_sets_true=nt, hog_sets_inferred=ni, hog_set_recall=round(tp / nt, 4),
                hog_set_precision=round(tp / ni, 4))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--presets", default="fast,accurate")
    ap.add_argument("--work", default=os.path.join(ROOT, "runs", "validation_hogparser"))
    ap.add_argument("--results", default=os.path.join(ROOT, "results", "validation_hogparser"))
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.results, exist_ok=True)
    presets = [p for p in a.presets.split(",") if p]
    rows, comp, runtimes = [], [], {}
    for name, (seed, kw) in SIMS.items():
        res = simulate(SimParams(**kw), seed)
        sim_dir = os.path.join(a.work, name)
        if not os.path.exists(os.path.join(sim_dir, "species_tree.nwk")):
            write_simulation(res, sim_dir)
        truth_df = true_branch_events(res.gene_trees, res.species_tree)
        truth_sets = true_ancestral_genes(res.gene_trees, res.species_tree)
        buf = io.StringIO()
        write_orthoxml(res.orthoxml, buf)
        truth_res = parse_hogs(io.BytesIO(buf.getvalue().encode()))
        cols = ["n_child", "gains", "losses", "dups", "dups_at_origin"]
        same = (truth_df.set_index("branch_child")[cols]
                == truth_res.branches.set_index("branch_child")[cols].loc[truth_df.branch_child]).all().all()
        if not same:
            raise SystemExit(f"{name}: parser on true HOGs != gene-tree truth (validation 2a failed)")
        tb = branch_table(truth_df)
        root = truth_res.tree.names[0]
        for preset in presets:
            out = os.path.join(a.work, f"{name}_ancestra_{preset}")
            dt = run_ancestra(sim_dir, out, preset, a.threads)
            if dt == dt:
                runtimes[f"{name}/{preset}"] = round(dt, 1)
            inf = parse_hogs(os.path.join(out, "hogs.orthoxml"), newick=os.path.join(sim_dir, "species_tree.nwk"))
            assert inf.check_conservation().empty
            ib = branch_table(inf.branches).loc[tb.index]
            c = tb.join(ib.drop(columns=["branch_parent", "child_is_leaf"]), rsuffix="_inferred")
            c.insert(0, "preset", preset)
            c.insert(0, "simulation", name)
            comp.append(c.reset_index())
            nonstem = tb.index != root
            internal = nonstem & ~tb["child_is_leaf"].to_numpy()
            leafb = tb["child_is_leaf"].to_numpy()
            for ev in EVENTS:
                for scope, mask in (("all_branches", nonstem), ("internal_branches", internal),
                                    ("leaf_branches", leafb)):
                    if ev == "n_child" and scope != "internal_branches":
                        continue
                    r = corr_metrics(tb.loc[mask, ev], ib.loc[mask, ev])
                    rows.append(dict(simulation=name, preset=preset, quantity=ev, scope=scope, **r))
            # root: number of ancestral genes at the root (no branch above)
            rows.append(dict(simulation=name, preset=preset, quantity="n_root_genes", scope="root",
                             n_branches=1, true_total=int(tb.loc[root, "n_child"]),
                             inferred_total=int(ib.loc[root, "n_child"]),
                             ratio=round(ib.loc[root, "n_child"] / tb.loc[root, "n_child"], 3)))
            extra = dict(simulation=name, preset=preset, **origin_accuracy(truth_res, inf),
                         **hog_set_recovery(truth_sets, inf), stats=json.dumps(inf.stats))
            rows.append(dict(simulation=name, preset=preset, quantity="gene_level", scope="genes",
                             **{k: v for k, v in extra.items() if k not in ("simulation", "preset")}))
            print(f"[validate] {name} {preset} done", file=sys.stderr, flush=True)
    met = pd.DataFrame(rows)
    met.to_csv(os.path.join(a.results, "metrics.tsv"), sep="\t", index=False)
    pd.concat(comp).to_csv(os.path.join(a.results, "branch_comparison.tsv"), sep="\t", index=False)
    write_summary(met, a.results, runtimes, presets)
    return 0


def write_summary(met: pd.DataFrame, out: str, runtimes: dict, presets: list[str]) -> None:
    L = ["# Validation of the HOG -> per-branch event parser (mcell.hogs)", "",
         "Generated by `scripts/validate_hogparser_sim.py` (fixed seeds; numbers below are copied from "
         "`metrics.tsv`).", "",
         "## (a) Parser on the TRUE HOGs vs. events from the true gene trees", "",
         "Exact agreement (all branches, all counts, and all ancestral gene sets) on the 3 simulations "
         "below and on the 4 configurations in `tests/test_hogs_simulation.py` (incl. HGT, WGD, a "
         "multifurcating tree, 40 species). Gains are *observable* gains: a family that survives only in a "
         "subclade originates at the MRCA of its survivors (the true root-genome origin is unobservable).", "",
         "## (b) Ancestra-inferred HOGs vs. truth", "",
         "Per-branch counts compared across branches (root stem excluded). `ratio` = inferred total / true "
         "total. `dups_total` = duplications incl. those at a family's origin. `n_child` = number of "
         "ancestral genes at each internal node.", ""]
    for sim in met.simulation.unique():
        L.append(f"### {sim}")
        L.append("")
        L.append("| preset | quantity | scope | true total | inferred total | ratio | Pearson r | Spearman | MAE |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        sub = met[(met.simulation == sim) & ~met.quantity.isin(["gene_level"])]
        for _, r in sub.iterrows():
            if r.scope == "leaf_branches" and r.quantity != "gains":
                pass
            L.append(f"| {r.preset} | {r.quantity} | {r.scope} | {r.true_total} | {r.inferred_total} | "
                     f"{r.ratio} | {r.get('pearson', '')} | {r.get('spearman', '')} | {r.get('mae', '')} |")
        L.append("")
    g = met[met.quantity == "gene_level"]
    L += ["### Ancestral gene sets and family origins", "",
          "`hog_set_recall/precision`: fraction of true / inferred ancestral genes at internal nodes whose "
          "extant gene set is reproduced exactly. `origin_*`: per extant gene, is its family's inferred "
          "origin node the true one, a descendant of it (too young = spurious gain + missed losses), or an "
          "ancestor (too old).", "",
          "| simulation | preset | genes | origin exact | too young | too old | other | HOG-set recall | "
          "HOG-set precision |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in g.iterrows():
        L.append(f"| {r.simulation} | {r.preset} | {int(r.genes)} | {r.origin_exact} | {r.origin_too_young} | "
                 f"{r.origin_too_old} | {r.origin_other} | {r.hog_set_recall} | {r.hog_set_precision} |")
    L += ["", f"Ancestra wall times (s, this run; empty if reused): {json.dumps(runtimes)}", ""]
    with open(os.path.join(out, "summary.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
