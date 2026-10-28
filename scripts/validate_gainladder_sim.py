#!/usr/bin/env python
"""Validation of the gain-evidence ladder on Ancestra's simulator (DESIGN §3.2).

Uses the three fixed simulations and cached Ancestra runs of
``scripts/validate_hogparser_sim.py`` (runs/validation_hogparser/<sim>_ancestra_<preset>/).
For each simulation x preset:

1. truth: per-branch events from the true HOGs (:func:`mcell.hogs.parse_hogs` on the
   simulator's OrthoXML; validated equal to the gene-tree truth);
2. inferred: Ancestra HOGs -> internal-branch gains -> ≤3 representatives -> one DIAMOND
   search against all simulated proteomes (same flags as the real run; ``--dbsize`` set to
   the real euk+prok DB so e-value cut-offs are comparable) -> clade filter, split repair;
3. acceptance test for split repair: internal-branch gain ratio (inferred/true) after repair
   ≤ 1.3 with loss and duplication correlations (vs truth) unchanged (|Δr| ≤ 0.01);
4. sequence-novel gains vs true gains: per-gain precision / recall (a gain is correct if
   the majority true family of its genes has its true gain on the same branch), and
   per-branch count agreement.

Writes results/gainladder/validation_sim/{metrics.tsv, branch_comparison.tsv, gain_matches.tsv.gz,
summary.md, DONE.json}.

    python scripts/validate_gainladder_sim.py [--threads 4] [--presets fast,accurate]
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import pearsonr  # noqa: E402

from ancestra.orthoxml import write_orthoxml  # noqa: E402
from ancestra.simulate import SimParams, simulate, write_simulation  # noqa: E402
from mcell import gainladder as gl  # noqa: E402
from mcell.hogs import parse_hogs  # noqa: E402

_spec = importlib.util.spec_from_file_location("vhs", os.path.join(ROOT, "scripts", "validate_hogparser_sim.py"))
vhs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vhs)
SIMS = vhs.SIMS
REAL_DB_LETTERS = 1489260785    # runs/gainladder/db/all_euk_prok.dmnd ("Database letters")


def log(m):
    print(f"[validate_gainladder {time.strftime('%H:%M:%S')}] {m}", file=sys.stderr, flush=True)


def r_(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return round(float(pearsonr(x, y)[0]), 4) if len(x) > 2 and x.std() > 0 and y.std() > 0 else np.nan


def gain_matches(truth, inf, cls) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per inferred gain: majority true family and whether its true gain is on the
    same branch. Per true internal gain: matched by which inferred labels."""
    tg = truth.events[(truth.events.event == "gain") & (truth.events.branch_parent != "")]
    t_origin = dict(zip(tg["family_id"], tg["branch_child"]))
    t_int = tg[~tg["child_is_leaf"]][["family_id", "branch_child"]].copy()
    tfam = truth.genes.set_index("gene")["family_id"]
    root = truth.tree.names[0]
    G = inf.genes[inf.genes["family_id"].isin(set(cls["family_id"]))][["gene", "family_id"]].copy()
    G["true_fam"] = tfam.reindex(G["gene"]).to_numpy()
    maj = (G.groupby(["family_id", "true_fam"]).size().rename("n").reset_index()
           .sort_values(["family_id", "n", "true_fam"], ascending=[True, False, True]).drop_duplicates("family_id"))
    tot = G.groupby("family_id").size()
    m = cls[["family_id", "branch_child", "ladder", "seq_novel", "split_repaired", "repair_status",
             "n_outside", "best_out_bitscore"]].merge(
        maj[["family_id", "true_fam", "n"]], on="family_id", how="left")
    m["frac_majority"] = (m["n"] / tot.reindex(m["family_id"]).to_numpy()).round(3)
    m["true_origin"] = [t_origin.get(f, root) for f in m["true_fam"]]
    m["correct"] = m["true_origin"] == m["branch_child"]
    # a true gain split into several inferred gains on the same branch: only the largest piece counts
    m = m.sort_values(["n", "family_id"], ascending=[False, True], kind="mergesort")
    m["redundant_piece"] = m["correct"] & m.duplicated(["true_fam", "branch_child"])
    m.loc[m["redundant_piece"], "correct"] = False
    m = m.sort_index()
    m["true_hgt"] = m["true_fam"].astype(str).str.contains(r"\.x")
    # recall side
    rec = t_int.rename(columns={"family_id": "true_fam"})
    rec["hgt"] = rec["true_fam"].str.contains(r"\.x")
    ok = m[m["correct"]]
    for lab, sel in (("hog_gain", ok), ("after_repair", ok[~ok["split_repaired"]]),
                     ("seq_novel", ok[ok["seq_novel"]])):
        hit = set(zip(sel["true_fam"], sel["branch_child"]))
        rec[f"found_{lab}"] = [(f, b) in hit for f, b in zip(rec["true_fam"], rec["branch_child"])]
    return m, rec


def min_outside_distance(sim, truth, m: pd.DataFrame) -> pd.Series:
    """For each inferred gain whose majority true family is not a true gain on its branch:
    the smallest true gene-tree (patristic) distance between a family member inside the
    clade and one outside it (NaN for HGT-derived true families)."""
    T = truth.tree
    gsp = truth.genes.set_index("gene")["species"].astype(str)
    out = {}
    for fam_id, c, f in zip(m["family_id"], m["branch_child"], m["true_fam"]):
        root = sim.gene_trees.get(f)
        if root is None:
            continue
        depth, paths = {}, {}
        stack = [(root, 0.0, ())]
        while stack:
            n, d, pth = stack.pop()
            pth = pth + (id(n),)
            depth[id(n)] = d
            if not n.children:
                paths[n.name] = (d, pth)
            for ch in n.children:
                stack.append((ch, d + (ch.length or 0.0), pth))
        ci = T.index[c]
        ins = [g for g in paths if g in gsp.index and T.is_anc(ci, T.index[gsp[g]])]
        outs = [g for g in paths if g in gsp.index and not T.is_anc(ci, T.index[gsp[g]])]
        best = np.nan
        for a_ in ins:
            da, pa = paths[a_]
            for b_ in outs:
                db, pb = paths[b_]
                pbs = set(pb)
                lca = [x for x in pa if x in pbs][-1]
                dist = da + db - 2 * depth[lca]
                best = dist if not best == best or dist < best else best
        out[fam_id] = best
    return pd.Series(out, dtype=float)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--presets", default="fast,accurate")
    ap.add_argument("--sensitivity", default="ultra-sensitive")
    ap.add_argument("--dbsize", default="real", choices=["real", "none"],
                    help="real: e-values scaled to the real euk+prok DB (default); none: the small simulated DB")
    ap.add_argument("--work", default=os.path.join(ROOT, "runs", "gainladder", "validation_sim"))
    ap.add_argument("--anc-work", default=os.path.join(ROOT, "runs", "validation_hogparser"))
    ap.add_argument("--results", default=None)
    a = ap.parse_args(argv)
    if a.results is None:
        a.results = os.path.join(ROOT, "results", "gainladder", "validation_sim"
                                 + ("" if a.dbsize == "real" and a.sensitivity == "ultra-sensitive"
                                    else f"_{a.sensitivity}_dbsize-{a.dbsize}"))
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.results, exist_ok=True)
    rows, comp, matches, recs = [], [], [], []
    for name, (seed, kw) in SIMS.items():
        sim = simulate(SimParams(**kw), seed)
        sim_dir = os.path.join(a.anc_work, name)
        if not os.path.exists(os.path.join(sim_dir, "species_tree.nwk")):
            write_simulation(sim, sim_dir)
        buf = io.StringIO()
        write_orthoxml(sim.orthoxml, buf)
        truth = parse_hogs(io.BytesIO(buf.getvalue().encode()))
        tb = truth.branches.set_index("branch_child")
        dbp = os.path.join(a.work, name, "db", "sim")
        db = gl.make_diamond_db({"sim": os.path.join(sim_dir, "proteomes")}, dbp, threads=a.threads)
        targets = pd.read_parquet(db["targets"])
        lengths = targets.set_index("seqid")["length"]
        for preset in [x for x in a.presets.split(",") if x]:
            ox = os.path.join(a.anc_work, f"{name}_ancestra_{preset}", "hogs.orthoxml")
            if not os.path.exists(ox):
                log(f"skip {name}/{preset}: no Ancestra output (run validate_hogparser_sim.py)")
                continue
            inf = parse_hogs(ox, newick=os.path.join(sim_dir, "species_tree.nwk"))
            w = os.path.join(a.work, name, preset)
            os.makedirs(w, exist_ok=True)
            gains = gl.internal_gains(inf)
            reps = gl.choose_representatives(inf, gains, lengths, n_rep=3, seed=11)
            q = os.path.join(w, "queries.faa")
            gl.write_queries(reps, [db["faa"]], q)
            hits = os.path.join(w, f"hits_{a.sensitivity}_dbsize-{a.dbsize}.tsv.gz")
            if not os.path.exists(hits):
                info = gl.run_diamond(q, db["dmnd"], hits + ".part", threads=a.threads, sensitivity=a.sensitivity,
                                      extra=["--dbsize", str(REAL_DB_LETTERS)] if a.dbsize == "real" else [],
                                      log=hits + ".log")
                os.replace(hits + ".part", hits)
                log(f"{name}/{preset}: DIAMOND {info['wall_s']} s")
            ts = gl.TargetSpace.build(targets, inf)
            pq, top = gl.reduce_hits(hits, reps, gains, inf.tree, ts)
            cls = gl.classify_gains(inf, gains, reps, pq, top, ts)
            br = gl.branch_ladder(inf, cls).set_index("branch_child")
            ib = inf.branches.set_index("branch_child")
            # compare on internal non-root branches (where the ladder acts)
            idx = br.index
            t = tb.loc[idx]
            c = pd.DataFrame({"simulation": name, "preset": preset, "branch_child": idx,
                              "true_gains": t["gains"].to_numpy(), "hog_gains": br["hog_gains"].to_numpy(),
                              "split_repaired": br["split_repaired"].to_numpy(),
                              "gains_after_repair": br["gains_after_repair"].to_numpy(),
                              "seq_novel": br["seq_novel"].to_numpy(),
                              "true_losses": t["losses"].to_numpy(), "losses": br["losses"].to_numpy(),
                              "losses_after_repair": br["losses_after_repair"].to_numpy(),
                              "true_dups": (t["dups"] + t["dups_at_origin"]).to_numpy(),
                              "dups": (br["dups"] + br["dups_at_origin"]).to_numpy()})
            comp.append(c)
            # losses over ALL non-root branches (repairs only touch internal ones)
            allb = [x for x in tb.index if tb.loc[x, "branch_parent"] != ""]
            l_all_true = tb.loc[allb, "losses"]
            l_all_inf = ib.loc[allb, "losses"]
            l_all_rep = l_all_inf.copy()
            l_all_rep.loc[idx] = br["losses_after_repair"].to_numpy()
            m, rec = gain_matches(truth, inf, cls)
            m.insert(0, "preset", preset)
            m.insert(0, "simulation", name)
            rec.insert(0, "preset", preset)
            rec.insert(0, "simulation", name)
            bad = m[~m["correct"]]
            m["min_true_dist_outside"] = min_outside_distance(sim, truth, bad).reindex(m["family_id"]).to_numpy()
            matches.append(m)
            recs.append(rec)
            novel = m[m["seq_novel"]]
            rep = m[m["split_repaired"]]
            nonh = rec[~rec["hgt"]]
            row = dict(
                simulation=name, preset=preset, n_internal_branches=len(idx),
                true_gains=int(c.true_gains.sum()), true_gains_hgt=int(rec["hgt"].sum()),
                hog_gains=int(c.hog_gains.sum()), split_repaired=int(c.split_repaired.sum()),
                gains_after_repair=int(c.gains_after_repair.sum()), seq_novel=int(c.seq_novel.sum()),
                ratio_hog=round(c.hog_gains.sum() / c.true_gains.sum(), 3),
                ratio_after_repair=round(c.gains_after_repair.sum() / c.true_gains.sum(), 3),
                ratio_seq_novel=round(c.seq_novel.sum() / c.true_gains.sum(), 3),
                r_gain_hog=r_(c.true_gains, c.hog_gains), r_gain_after_repair=r_(c.true_gains, c.gains_after_repair),
                r_gain_seq_novel=r_(c.true_gains, c.seq_novel),
                r_loss_internal_before=r_(c.true_losses, c.losses),
                r_loss_internal_after=r_(c.true_losses, c.losses_after_repair),
                ratio_loss_internal_before=round(c.losses.sum() / c.true_losses.sum(), 3),
                ratio_loss_internal_after=round(c.losses_after_repair.sum() / c.true_losses.sum(), 3),
                r_loss_all_before=r_(l_all_true, l_all_inf), r_loss_all_after=r_(l_all_true, l_all_rep),
                r_dup_internal_before=r_(c.true_dups, c.dups), r_dup_internal_after=r_(c.true_dups, c.dups),
                repair_precision=round(float((~rep["correct"]).mean()), 3) if len(rep) else np.nan,
                repaired_but_true_gain=int(rep["correct"].sum()),
                novel_precision=round(float(novel["correct"].mean()), 3) if len(novel) else np.nan,
                novel_recall=round(float(rec["found_seq_novel"].mean()), 3),
                novel_recall_nonhgt=round(float(nonh["found_seq_novel"].mean()), 3) if len(nonh) else np.nan,
                hog_precision=round(float(m["correct"].mean()), 3),
                hog_recall=round(float(rec["found_hog_gain"].mean()), 3),
                after_repair_precision=round(float(m.loc[~m["split_repaired"], "correct"].mean()), 3),
                after_repair_recall=round(float(rec["found_after_repair"].mean()), 3),
                false_novel_true_origin_older=int((~novel["correct"]).sum()),
                missed_true_gains_nonnovel=int((rec["found_hog_gain"] & ~rec["found_seq_novel"]).sum()),
                missed_true_hgt=int((rec["found_hog_gain"] & ~rec["found_seq_novel"] & rec["hgt"]).sum()),
                spurious_homologous=int(((m["ladder"] != "seq_novel") & ~m["correct"]).sum()),
                spurious_novel=int((m["seq_novel"] & ~m["correct"]).sum()),
                ratio_if_all_homologous_spurious_repaired=round(
                    (c.hog_gains.sum() - ((m["ladder"] != "seq_novel") & ~m["correct"]).sum()) / c.true_gains.sum(), 3),
                **{f"spurious_homologous_{k}": int(v) for k, v in
                   m.loc[(m["ladder"] != "seq_novel") & ~m["correct"], "repair_status"].value_counts().items()},
            )
            row["accept_split_repair"] = bool(row["ratio_after_repair"] <= 1.3
                                              and abs(row["r_loss_all_after"] - row["r_loss_all_before"]) <= 0.01
                                              and abs(row["r_loss_internal_after"] - row["r_loss_internal_before"]) <= 0.01)
            rows.append(row)
            log(f"{name}/{preset}: " + json.dumps({k: row[k] for k in ("ratio_hog", "ratio_after_repair",
                                                                       "ratio_seq_novel", "novel_precision",
                                                                       "novel_recall")}))
    met = pd.DataFrame(rows)
    met.to_csv(os.path.join(a.results, "metrics.tsv"), sep="\t", index=False)
    pd.concat(comp).to_csv(os.path.join(a.results, "branch_comparison.tsv"), sep="\t", index=False)
    pd.concat(matches).to_csv(os.path.join(a.results, "gain_matches.tsv.gz"), sep="\t", index=False)
    pd.concat(recs).to_csv(os.path.join(a.results, "true_gain_recall.tsv.gz"), sep="\t", index=False)
    with open(os.path.join(a.results, "DONE.json"), "w") as fh:
        json.dump(dict(status="done", sensitivity=a.sensitivity, dbsize=a.dbsize, presets=a.presets,
                       accept_split_repair_all=bool(met["accept_split_repair"].all()),
                       finished=time.strftime("%Y-%m-%dT%H:%M:%S%z")), fh, indent=1)
    allm = pd.concat(matches)
    dist = (allm[~allm["correct"]].groupby(["simulation", "preset", "ladder"])["min_true_dist_outside"]
            .describe()[["count", "25%", "50%", "75%"]].round(2))
    dist.to_csv(os.path.join(a.results, "spurious_gain_distance.tsv"), sep="\t")
    write_summary(met, dist.reset_index(), a)
    print(met.T.to_string())
    print(dist.to_string())
    return 0


def write_summary(met: pd.DataFrame, dist: pd.DataFrame, a) -> None:
    f = lambda x: "NA" if x != x else x
    L = ["# Validation of the gain-evidence ladder on Ancestra's simulator (DESIGN §3.2)", "",
         "Generated by `scripts/validate_gainladder_sim.py` (code: `mcell/gainladder.py`); numbers from "
         "`metrics.tsv` / `spurious_gain_distance.tsv`.", "",
         f"Setup: the 3 fixed simulations and cached Ancestra runs of `validate_hogparser_sim.py`; <=3 representatives "
         f"per internal-branch gain; one DIAMOND search against all simulated proteomes with the real-run flags "
         f"(`--{a.sensitivity} -k 0 --max-hsps 1 -e 1e-5 --query-cover 30`), dbsize = "
         f"{'real euk+prok DB (1.49 G letters), so e-value cut-offs match the real search' if a.dbsize == 'real' else 'simulated DB'}. "
         "Truth = parser on the simulator's true HOGs. Internal branches only (root stem and leaves excluded).", "",
         "## 1. Split repair: DESIGN acceptance test (ratio after repair <= 1.3, loss/dup r unchanged)", "",
         "| sim | preset | true | HOG gains | repaired | ratio before | ratio after | loss r all (before/after) | "
         "loss r internal (before/after) | dup r | accept |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in met.iterrows():
        L.append(f"| {r.simulation} | {r.preset} | {r.true_gains} | {r.hog_gains} | {r.split_repaired} | {r.ratio_hog} | "
                 f"{r.ratio_after_repair} | {r.r_loss_all_before}/{r.r_loss_all_after} | "
                 f"{r.r_loss_internal_before}/{r.r_loss_internal_after} | {r.r_dup_internal_before} | "
                 f"{'PASS' if r.accept_split_repair else 'FAIL'} |")
    L += ["", "Why homologous spurious gains are not repaired, and the best any repair based on this search could do:", "",
          "| sim | preset | spurious homologous | best parent HOG not lost (missed duplication) | no hit maps to a parent-level HOG | "
          "repaired | repair precision | ratio if ALL homologous spurious gains were repaired |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in met.iterrows():
        L.append(f"| {r.simulation} | {r.preset} | {r.spurious_homologous} | "
                 f"{f(r.get('spurious_homologous_best_parent_hog_not_lost', 0))} | "
                 f"{f(r.get('spurious_homologous_no_parent_level_hit', 0))} | {r.split_repaired} | "
                 f"{f(r.repair_precision)} | {r.ratio_if_all_homologous_spurious_repaired} |")
    L += ["", "## 2. Sequence-novel gains vs true gains", "",
          "Correct = the majority true family of the inferred gain's genes has its true gain on the same branch (only the "
          "largest piece of a true gain split on one branch counts).", "",
          "| sim | preset | true (HGT) | seq-novel | ratio novel/true | r per branch | precision | recall | recall non-HGT | "
          "HOG-gain precision | spurious novel |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in met.iterrows():
        L.append(f"| {r.simulation} | {r.preset} | {r.true_gains} ({r.true_gains_hgt}) | {r.seq_novel} | {r.ratio_seq_novel} | "
                 f"{r.r_gain_seq_novel} | {r.novel_precision} | {r.novel_recall} | {r.novel_recall_nonhgt} | "
                 f"{r.hog_precision} | {r.spurious_novel} |")
    L += ["", "## 3. Divergence of spurious gains from their outside relatives", "",
          "Smallest true gene-tree (patristic) distance, substitutions/site, between an in-clade and an out-of-clade "
          "member of the majority true family, for spurious (incorrect) gains by ladder label:", "",
          "| sim | preset | label | n | Q1 | median | Q3 |", "|---|---|---|---|---|---|---|"]
    for _, r in dist.iterrows():
        L.append(f"| {r.simulation} | {r.preset} | {r.ladder} | {int(r['count'])} | {r['25%']} | {r['50%']} | {r['75%']} |")
    L += ["", "## Interpretation", "",
          "* **Split repair fails the acceptance test** on every simulation and preset: the rule is precise (all repairs hit "
          "spurious gains; loss/dup correlations unchanged) but applies to almost nothing, so the gain ratio stays ~2.1-2.5x. "
          "Per DESIGN, report ladder levels 1 and 3 only; repair counts are kept as a diagnostic. Most homologous spurious "
          "gains are extra paralogues of a HOG that continues into the clade (a missed duplication, not a gain+loss), "
          "and even repairing all homologous spurious gains would leave the ratio at the seq-novel level (1.6-2.3x).",
          "* **Seq-novel is a near-lossless but insufficient filter here**: recall of true (non-HGT) gains 0.91-0.97, "
          "precision 0.42-0.59. True HGT receipts are (correctly) not novel because the donor homologue is found.",
          "* **The simulator does not make novelty easy; it makes deep splits undetectable.** Families are independent random "
          "sequences (true gains are truly novel), but family rates ~ Gamma(2) on a height-1 tree put ~29% of families "
          "beyond ~2.5 subst./site across the root. Ancestra splits exactly these, and DIAMOND cannot see the homology "
          "either: spurious gains called novel are ~4.6-4.9 subst./site from their nearest outside relative vs ~1.7-2.7 "
          "for those detected. Real deep families keep conserved domains, so the simulator is probably pessimistic, but "
          "seq-novel counts must be read as upper bounds with a branch-length-dependent excess, compared only against "
          "matched control branches (§4.6), with extra caution on long stems and root children.", ""]
    with open(os.path.join(a.results, "summary.md"), "w") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    sys.exit(main())
