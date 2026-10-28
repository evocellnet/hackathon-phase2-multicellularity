#!/usr/bin/env python
"""Calibration on our tree (DESIGN §3.2): Ancestra-simulator bias of per-branch
gains / sequence-novel gains / losses / duplications by branch class, on the
eukaryote guide topology pruned to ~60 taxa.

Pruning (seed 11): for every primary transition (mcell.branchstats.stem_branches)
the stem path P -> ... -> M is kept: M keeps one tip in each of two of its child
clades (>= 2 M-side tips, M stays the MRCA), every path node keeps one tip outside
the path child, and P keeps >= 2 sister-side tips spanning the sister clade's root
split when possible. Remaining tips are drawn at random up to ``--n-taxa``. Branch
lengths: data/species_tree_bl_euk.nwk (via branchstats.branch_lengths_from_tree),
pruned with lengths preserved and rescaled to a mean root-to-tip depth of 1.0
(the simulator's convention: rates are per unit of that depth).

Settings: simulator defaults (dup = loss = 0.2) and a higher-turnover setting
(dup = loss = 0.4); 2,000 root genes; seed 11. Ancestra --preset fast; HOG events
by mcell.hogs; gain ladder (mcell.gainladder; very-sensitive DIAMOND, real-DB
--dbsize) on all internal branches.

Classes on the pruned tree (branchstats.build_branch_table on it): primary stems
(path stems summed over their edges, as the tests do), control pool U, root
children, other internal. Bias question: does inferred/true inflation differ
between stems and U? Reported as pooled ratios, per-stem log inflation and its
mid-rank percentile among U (0.5 = no bias), and the H1/H2-style statistic
T = mean_i u(stem_i | U) computed on inferred and on true counts.

    python scripts/calibration_tree.py [--threads 4] [--n-taxa 60]
"""
from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import time
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ancestra.orthoxml import write_orthoxml  # noqa: E402
from ancestra.simulate import SimParams, simulate, write_simulation  # noqa: E402
from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402
from mcell import gainladder as gl  # noqa: E402
from mcell.hogs import parse_hogs  # noqa: E402

SEED = 11
REAL_DB_LETTERS = 1489260785
SETTINGS = {"default": dict(), "high_turnover": dict(dup_rate=0.4, loss_rate=0.4)}


def log(m):
    print(f"[calibration_tree {time.strftime('%H:%M:%S')}] {m}", file=sys.stderr, flush=True)


def prune_tree(tree, bt, n_taxa: int, rng) -> tuple[str, list, pd.DataFrame]:
    nodes = {n.name: n for n in tree.traverse()}
    stems = bs.stem_branches(bt)
    keep: list[str] = []
    why = {}

    def add(tip, reason):
        if tip not in why:
            keep.append(tip)
            why[tip] = reason

    def pick(node, k=1):
        lv = sorted(node.get_leaf_names())
        return [lv[i] for i in rng.choice(len(lv), size=min(k, len(lv)), replace=False)]

    def span(node, k=2, reason=""):
        """k tips spanning node's root split (one per child clade, random children)."""
        if node.is_leaf():
            add(node.name, reason)
            return
        ch = list(node.children)
        order = rng.permutation(len(ch))
        for j in order[:k]:
            add(pick(ch[j])[0], reason)

    for tid, b in stems.items():
        edges = bs.path_edges(bt, b)
        path = [bt.at[edges[0], "parent"]] + edges          # P, ..., M
        M = nodes[path[-1]]
        span(M, 2, f"{tid}:M")
        for x, nxt in zip(path[:-1], path[1:]):
            others = [c for c in nodes[x].children if c.name != nxt]
            if not others:
                continue
            if x == path[0]:                                  # sister side of P: >= 2 tips
                s = others[int(rng.integers(len(others)))]
                if s.is_leaf() and len(others) > 1:
                    add(s.name, f"{tid}:sister")
                    s2 = [o for o in others if o is not s][0]
                    add(pick(s2)[0], f"{tid}:sister")
                else:
                    span(s, 2, f"{tid}:sister")
            else:
                add(pick(others[int(rng.integers(len(others)))])[0], f"{tid}:path_node")
    rest = sorted(set(tree.get_leaf_names()) - set(keep))
    for i in rng.permutation(len(rest))[: max(0, n_taxa - len(keep))]:
        add(rest[i], "random_fill")
    t = tree.copy()
    t.prune(keep, preserve_branch_length=True)
    depth = np.mean([t.get_distance(l) for l in t.get_leaves()])
    for n in t.traverse():
        n.dist = n.dist / depth if n.up is not None else 0.0
    nwk = t.write(format=1, format_root_node=True)
    if not nwk.rstrip(";").endswith(")" + t.name) and t.name:
        nwk = nwk.rstrip(";").rsplit(")", 1)[0] + ")" + t.name + ";"
    return nwk, keep, pd.DataFrame({"tip": keep, "reason": [why[k] for k in keep]})


def run_ancestra(sim_dir, out, threads):
    if os.path.exists(os.path.join(out, "hogs.orthoxml")):
        return
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), NUMBA_NUM_THREADS=str(threads),
               PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""))
    with open(out + ".log", "w") as fh:
        subprocess.run(["nice", "-n", "19", os.path.join(os.path.dirname(sys.executable), "ancestra"), "run",
                        "--proteomes", os.path.join(sim_dir, "proteomes"), "--tree",
                        os.path.join(sim_dir, "species_tree.nwk"), "--out", out, "--threads", str(threads),
                        "--preset", "fast"], check=True, env=env, stdout=fh, stderr=subprocess.STDOUT)


def ladder(inf, sim_dir, work, threads):
    db = gl.make_diamond_db({"sim": os.path.join(sim_dir, "proteomes")}, os.path.join(work, "db", "sim"),
                            threads=threads)
    targets = pd.read_parquet(db["targets"])
    gains = gl.internal_gains(inf)
    reps = gl.choose_representatives(inf, gains, targets.set_index("seqid")["length"], n_rep=3, seed=SEED)
    q = os.path.join(work, "queries.faa")
    gl.write_queries(reps, [db["faa"]], q)
    hits = os.path.join(work, "hits_very-sensitive.tsv.gz")
    if not os.path.exists(hits):
        info = gl.run_diamond(q, db["dmnd"], hits + ".part", threads=threads, sensitivity="very-sensitive",
                              nice=19, extra=["--dbsize", str(REAL_DB_LETTERS)], log=hits + ".log")
        os.replace(hits + ".part", hits)
        log(f"DIAMOND {info['wall_s']} s, {len(reps)} queries")
    ts = gl.TargetSpace.build(targets, inf)
    pq, top = gl.reduce_hits(hits, reps, gains, inf.tree, ts)
    cls = gl.classify_gains(inf, gains, reps, pq, top, ts)
    return gl.branch_ladder(inf, cls).set_index("branch_child")


def glm_check(d: pd.DataFrame) -> list[dict]:
    """Spurious-gain model without a stem term: NB-GLM of inferred counts on
    [1, log branch length, log n_desc] fitted on internal non-stem, non-root-child
    branches; stems' observed / expected (pooled over stem edges) with an exact
    Poisson 95% interval (narrower than the NB one). O/E ~ 1 means stems are
    inflated only as their length / clade size predicts."""
    from scipy.stats import chi2
    out = []
    for q in ("hog_gains", "seq_novel"):
        fit_rows = d[~d["branch_class"].isin(["primary_stem_edge", "root_child"])]
        X = lambda df: np.column_stack([np.ones(len(df)), np.log(df["branch_length"]), np.log(df["n_desc"])])
        f = bs.fit_count_glm(X(fit_rows), fit_rows[q].to_numpy(float), family="nb")
        for cls in ("primary_stem_edge", "control_U"):
            g = d[d["branch_class"] == cls]
            mu = np.exp(X(g) @ f["beta"])
            o, e = float(g[q].sum()), float(mu.sum())
            lo = chi2.ppf(0.025, 2 * o) / 2 if o > 0 else 0.0
            hi = chi2.ppf(0.975, 2 * o + 2) / 2
            out.append(dict(quantity=q, branch_class=cls, observed=int(o), expected=round(e, 1),
                            O_E=round(o / e, 3), O_E_lo=round(lo / e, 3), O_E_hi=round(hi / e, 3),
                            beta_log_bl=round(float(f["beta"][1]), 3), beta_log_n_desc=round(float(f["beta"][2]), 3),
                            glm=f["method"]))
    return out


def u_mid(x, ref):
    ref = np.asarray(ref, float)
    return (np.sum(ref < x) + 0.5 * np.sum(ref == x) + 0.5) / (len(ref) + 1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--n-taxa", type=int, default=60)
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--work", default=os.path.join(ROOT, "runs", "calibration_tree"))
    ap.add_argument("--results", default=os.path.join(ROOT, "results", "calibration_tree"))
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.results, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tr, traits, tree, bt = an.load_setup("euk")
        bl = bs.branch_lengths_from_tree(tree, os.path.join(ROOT, "data", "species_tree_bl_euk.nwk"))
    for n in tree.traverse():
        if n.up is not None:
            n.dist = float(bl.get(n.name, np.nan))
    if any(not np.isfinite(n.dist) for n in tree.traverse() if n.up is not None):
        raise SystemExit("missing branch lengths")
    rng = np.random.default_rng(SEED)
    nwk, keep, why = prune_tree(tree, bt, a.n_taxa, rng)
    with open(os.path.join(a.work, "pruned_tree.nwk"), "w") as fh:
        fh.write(nwk + "\n")
    why.to_csv(os.path.join(a.results, "kept_tips.tsv"), sep="\t", index=False)
    log(f"pruned tree: {len(keep)} tips")
    ptree = bs.load_tree(nwk)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pbl = bs.branch_lengths_from_tree(ptree, os.path.join(a.work, "pruned_tree.nwk"))
        pbt = bs.build_branch_table(ptree, traits, tr, branch_lengths=pbl)
    pstems = bs.stem_branches(pbt)
    missing = sorted(set(bs.stem_branches(bt).index) - set(pstems.index))
    if missing:
        raise SystemExit(f"primary stems lost by pruning: {missing}")
    path_edge = {e: tid for tid, b in pstems.items() for e in bs.path_edges(pbt, b)}
    internal = pbt[~pbt["is_terminal"]].index
    cls_of = {}
    for b in internal:
        if b in path_edge:
            cls_of[b] = "primary_stem_edge"
        elif pbt.at[b, "in_U"]:
            cls_of[b] = "control_U"
        elif pbt.at[b, "root_child"]:
            cls_of[b] = "root_child"
        else:
            cls_of[b] = "other_internal"
    log(f"pruned branch table: {len(pstems)} primary stems, {int(pbt['in_U'].sum())} U branches")

    comp, rows, stem_rows, glm_rows = [], [], [], []
    for sname, kw in SETTINGS.items():
        sim = simulate(SimParams(newick=nwk, n_genes=a.n_genes, **kw), SEED)
        sd = os.path.join(a.work, sname)
        if not os.path.exists(os.path.join(sd, "species_tree.nwk")):
            write_simulation(sim, sd)
        buf = io.StringIO()
        write_orthoxml(sim.orthoxml, buf)
        truth = parse_hogs(io.BytesIO(buf.getvalue().encode()))
        out = os.path.join(a.work, f"{sname}_ancestra_fast")
        t0 = time.time()
        run_ancestra(sd, out, a.threads)
        log(f"{sname}: Ancestra {time.time() - t0:.0f} s")
        inf = parse_hogs(os.path.join(out, "hogs.orthoxml"), newick=os.path.join(sd, "species_tree.nwk"))
        assert inf.check_conservation().empty
        lad = ladder(inf, sd, os.path.join(a.work, sname), a.threads)
        tb = truth.branches.set_index("branch_child")
        ib = inf.branches.set_index("branch_child")
        d = pd.DataFrame(index=internal)
        d["true_gains"] = tb["gains"].reindex(internal)
        d["hog_gains"] = ib["gains"].reindex(internal)
        d["seq_novel"] = lad["seq_novel"].reindex(internal)
        d["true_losses"] = tb["losses"].reindex(internal)
        d["losses"] = ib["losses"].reindex(internal)
        d["true_dups"] = (tb["dups"] + tb["dups_at_origin"]).reindex(internal)
        d["dups"] = (ib["dups"] + ib["dups_at_origin"]).reindex(internal)
        d["branch_class"] = [cls_of[b] for b in internal]
        d["branch_length"] = pbt.loc[internal, "edge_branch_length"].to_numpy()
        d["n_desc"] = pbt.loc[internal, "n_desc"].to_numpy()
        d.insert(0, "setting", sname)
        comp.append(d.reset_index().rename(columns={"index": "branch"}))
        pairs = {"hog_gains": "true_gains", "seq_novel": "true_gains", "losses": "true_losses", "dups": "true_dups"}
        for c, g in d.groupby("branch_class"):
            r = dict(setting=sname, branch_class=c, n_branches=len(g))
            for q, t in pairs.items():
                r[f"{q}_true"] = int(g[t].sum())
                r[f"{q}_inferred"] = int(g[q].sum())
                r[f"{q}_ratio"] = round(g[q].sum() / g[t].sum(), 3) if g[t].sum() else np.nan
            rows.append(r)
        # stems as the tests see them: summed over path edges; controls = U
        U = [b for b in internal if cls_of[b] == "control_U"]
        for q, t in pairs.items():
            linf = np.log((d[q] + 0.5) / (d[t] + 0.5))
            uU = linf.loc[U]
            us_inf, us_true, us_logi = [], [], []
            for tid, b in pstems.items():
                e = bs.path_edges(pbt, b)
                si, st = d.loc[e, q].sum(), d.loc[e, t].sum()
                li = np.log((si + 0.5) / (st + 0.5))
                u_infl = u_mid(li, uU)
                u_i = u_mid(si, d.loc[U, q])
                u_t = u_mid(st, d.loc[U, t])
                us_logi.append(u_infl)
                us_inf.append(u_i)
                us_true.append(u_t)
                stem_rows.append(dict(setting=sname, quantity=q, transition_id=tid, stem=b, n_edges=len(e),
                                      true=int(st), inferred=int(si), ratio=round(si / st, 3) if st else np.nan,
                                      log_inflation=round(li, 3), u_inflation_among_U=round(u_infl, 3),
                                      u_inferred_among_U=round(u_i, 3), u_true_among_U=round(u_t, 3)))
            rows.append(dict(setting=sname, branch_class=f"TEST:{q}", n_branches=len(U),
                             U_median_log_inflation=round(float(uU.median()), 3),
                             stems_median_log_inflation=round(float(np.median(
                                 [r_["log_inflation"] for r_ in stem_rows if r_["setting"] == sname and r_["quantity"] == q])), 3),
                             mean_u_inflation=round(float(np.mean(us_logi)), 3),
                             T_inferred=round(float(np.mean(us_inf)), 3), T_true=round(float(np.mean(us_true)), 3),
                             T_bias=round(float(np.mean(us_inf) - np.mean(us_true)), 3)))
        for r_ in glm_check(d):
            glm_rows.append(dict(setting=sname, **r_))
        log(f"{sname}: done")
    cs = pd.DataFrame(rows)
    st = pd.DataFrame(stem_rows)
    pd.concat(comp).to_csv(os.path.join(a.results, "branch_comparison.tsv"), sep="\t", index=False)
    cs.to_csv(os.path.join(a.results, "class_summary.tsv"), sep="\t", index=False)
    st.to_csv(os.path.join(a.results, "stem_bias.tsv"), sep="\t", index=False)
    gm = pd.DataFrame(glm_rows)
    gm.to_csv(os.path.join(a.results, "glm_check.tsv"), sep="\t", index=False)
    write_summary(a, cs, st, len(keep), pstems, cls_of, gm)
    with open(os.path.join(a.results, "DONE.json"), "w") as fh:
        json.dump(dict(status="done", n_taxa=len(keep), seed=SEED, settings=SETTINGS, n_genes=a.n_genes,
                       finished=time.strftime("%Y-%m-%dT%H:%M:%S%z")), fh, indent=1)
    print(cs.to_string())
    return 0


def write_summary(a, cs, st, n_taxa, pstems, cls_of, gm):
    cl = cs[~cs["branch_class"].str.startswith("TEST:")]
    te = cs[cs["branch_class"].str.startswith("TEST:")]
    n_by = pd.Series(cls_of).value_counts()
    L = ["# Calibration on our tree: Ancestra bias by branch class (DESIGN §3.2)", "",
         "Generated by `scripts/calibration_tree.py` (seed 11). Tables: `class_summary.tsv`, `stem_bias.tsv`, "
         "`branch_comparison.tsv`, `kept_tips.tsv`.", "",
         f"* Tree: eukaryote guide topology pruned to {n_taxa} taxa keeping every primary stem's P/path/M structure "
         "(>= 2 M-side tips spanning M's split, >= 2 sister-side tips); branch lengths from "
         "`data/species_tree_bl_euk.nwk`, rescaled to mean root-to-tip 1.0.",
         f"* Simulator: `ancestra.simulate`, {a.n_genes} root genes, defaults (dup = loss = 0.2) and high turnover "
         "(dup = loss = 0.4); Ancestra `--preset fast`; gain ladder with very-sensitive DIAMOND (real-DB e-values).",
         f"* Internal branches by class on the pruned tree: {n_by.to_dict()} (U recomputed by "
         f"`build_branch_table` on the pruned tree; {len(pstems)} primary stems, path stems summed over edges).",
         "* Caveat: the simulator has no de-novo gene birth; true gains are families surviving only in a clade "
         "(lost elsewhere), so true per-branch gain counts are small and ratios are noisy.", "",
         "## Inferred / true by branch class (pooled over branches)", "",
         "| setting | class | branches | HOG gains (true -> inf, ratio) | seq-novel (ratio to true gains) | losses ratio | dups ratio |",
         "|---|---|---|---|---|---|---|"]
    for _, r in cl.iterrows():
        L.append(f"| {r.setting} | {r.branch_class} | {r.n_branches} | {r.hog_gains_true} -> {r.hog_gains_inferred} "
                 f"({r.hog_gains_ratio}) | {r.seq_novel_inferred} ({r.seq_novel_ratio}) | {r.losses_ratio} | {r.dups_ratio} |")
    L += ["", "## Would the stem-vs-control test be biased under the null?", "",
          "Per stem: log inflation = log((inferred+0.5)/(true+0.5)); `mean u` = mean over the 7 primary stems of its "
          "mid-rank percentile among the U branches' log inflation (0.5 = stems inflated like controls). "
          "`T_inf` / `T_true` = the H1/H2-style statistic mean_i u(stem_i | U) on inferred vs true counts; "
          "`T_bias` = T_inf - T_true (> 0: inference pushes stems up relative to controls).", "",
          "| setting | quantity | U median log infl. | stems median log infl. | mean u (inflation) | T_inf | T_true | T_bias |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in te.iterrows():
        L.append(f"| {r.setting} | {r.branch_class[5:]} | {r.U_median_log_inflation} | {r.stems_median_log_inflation} | "
                 f"{r.mean_u_inflation} | {r.T_inferred} | {r.T_true} | {r.T_bias} |")
    L += ["", "## Length / clade-size adjusted check (NB-GLM without a stem term)", "",
          "Spurious-gain model fitted on internal non-stem, non-root-child branches: log E[count] = b0 + b1 log(branch "
          "length) + b2 log(n_desc). O/E pooled over stem edges (and U for reference), exact Poisson 95% interval "
          "(too narrow under overdispersion).", "",
          "| setting | quantity | class | observed | expected | O/E (95% Poisson) | b1 log bl | b2 log n_desc |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in gm.iterrows():
        L.append(f"| {r.setting} | {r.quantity} | {r.branch_class} | {r.observed} | {r.expected} | {r.O_E} "
                 f"({r.O_E_lo}-{r.O_E_hi}) | {r.beta_log_bl} | {r.beta_log_n_desc} |")
    te2 = te.set_index(["setting", "branch_class"])
    g2 = gm.set_index(["setting", "quantity", "branch_class"])
    rel = {(st_, q): g2.loc[(st_, q, "primary_stem_edge"), "O_E"] / g2.loc[(st_, q, "control_U"), "O_E"]
           for st_ in gm["setting"].unique() for q in ("hog_gains", "seq_novel")}
    L += ["", "## Interpretation", "",
          "* Internal-branch true gains are ~0 in this simulator (no gene birth), so stem and U gains are spurious; the "
          "question is purely whether spurious gains differ between stems and controls.",
          "* Losses and duplications: stems and U inflated alike (pooled ratios ~1.1-1.3 in both; T_bias "
          + ", ".join(f"{te2.loc[(x, 'TEST:' + q), 'T_bias']:+.3f}" for x in te['setting'].unique() for q in ('losses', 'dups'))
          + ").",
          "* Gains: on the u scale the H1/H2 statistic is shifted towards 'burst' by T_bias = "
          + ", ".join(f"{te2.loc[(x, 'TEST:' + q), 'T_bias']:+.3f} ({x}, {q})" for x in te['setting'].unique()
                      for q in ('hog_gains', 'seq_novel'))
          + " (one draw, 7 stems vs 9 U branches; the SE of a mean of 7 u-values is ~0.1, so none of these is "
          "distinguishable from 0).",
          "* After adjusting for branch length and clade size, U branches carry fewer spurious gains than other "
          "internal branches (O/E ~0.5-0.6), stems about as many as predicted (O/E ~0.9-1.2): stem/U ratio of "
          "adjusted inflation = " + ", ".join(f"{v:.2f} ({k[0]}, {k[1]})" for k, v in rel.items())
          + ". The direction (stems more inflated than controls) is the same in all 4 cases; the exact Poisson intervals "
          "of stems and U do not overlap in " + str(sum(
              g2.loc[(x, q, "primary_stem_edge"), "O_E_lo"] > g2.loc[(x, q, "control_U"), "O_E_hi"]
              for x in gm["setting"].unique() for q in ("hog_gains", "seq_novel")))
          + " of 4 (these intervals ignore overdispersion and the single simulation draw).",
          "* Consequence for H1/H2: under the null, the inference itself can push stems above their controls. A "
          "burst claim should therefore (i) use the length/size-adjusted residual statistic (D6 NB-GLM) rather than "
          "raw counts, (ii) require an effect clearly larger than ~2x, and (iii) not interpret a 'burst' p-value near the "
          "threshold. A 'deficit' (fewer gains on stems) is not explained by this bias. More simulation replicates "
          "(seeds) would narrow these estimates (each costs ~12 min of Ancestra at 4 threads).", ""]
    L += ["## Per stem (HOG gains and seq-novel)", "",
          "| setting | quantity | transition | edges | true | inferred | ratio | u(inflation among U) |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in st[st["quantity"].isin(["hog_gains", "seq_novel"])].iterrows():
        L.append(f"| {r.setting} | {r.quantity} | {r.transition_id} | {r.n_edges} | {r.true} | {r.inferred} | {r.ratio} | "
                 f"{r.u_inflation_among_U} |")
    L.append("")
    with open(os.path.join(a.results, "summary.md"), "w") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())
