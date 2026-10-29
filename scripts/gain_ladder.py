#!/usr/bin/env python
"""Gain-evidence ladder (DESIGN §3.2) for one Ancestra run.

    python scripts/gain_ladder.py makedb                       # euk+prok DIAMOND DB (shared by both runs)
    python scripts/gain_ladder.py prepare  --run prok          # parse HOGs, gains, representatives, queries
    python scripts/gain_ladder.py search   --run prok [--threads 4] [--sensitivity ultra-sensitive]
    python scripts/gain_ladder.py classify --run prok          # clade filter, novelty, split repair, tables
    python scripts/gain_ladder.py arch     --run prok          # optional: Pfam architecture novelty (re-runnable)
    python scripts/gain_ladder.py all      --run prok          # prepare + search + classify (+ arch if Pfam exists)
    python scripts/gain_ladder.py pilot    --n 500             # DIAMOND timing on random queries

``--run`` names runs/<run>_fast (override with --run-dir). Intermediates go to
runs/gainladder/<run>/, results to results/gainladder/<run>/ (branch_gain_ladder.tsv,
gains_classified.tsv.gz, summary.md, DONE.json). Definitions: ``mcell/gainladder.py``.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import random
import resource
import shutil
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import __version__  # noqa: E402
from mcell import gainladder as gl  # noqa: E402
from mcell.hogs import parse_hogs  # noqa: E402

DB_PREFIX = os.path.join(ROOT, "runs", "gainladder", "db", "all_euk_prok")
PROTEOMES = {"euk": os.path.join(ROOT, "data", "proteomes", "euk"),
             "prok": os.path.join(ROOT, "data", "proteomes", "prok")}
SEED = 11


def log(msg: str) -> None:
    print(f"[gain_ladder {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def paths(a) -> dict:
    run_dir = a.run_dir
    if run_dir is None:
        cands = [os.path.join(ROOT, "runs", f"{a.run}_fast{x}") for x in ("", "_p")]
        done = [c for c in cands if os.path.exists(os.path.join(c, "DONE.json"))]
        run_dir = done[0] if done else cands[0]
    work = os.path.join(ROOT, "runs", "gainladder", a.run)
    res = os.path.join(ROOT, "results", "gainladder", a.run)
    ox = [p for p in (os.path.join(run_dir, "hogs.orthoxml.gz"), os.path.join(run_dir, "hogs.orthoxml"))
          if os.path.exists(p)]
    return dict(run_dir=run_dir, work=work, res=res, orthoxml=ox[0] if ox else None,
                tree=os.path.join(run_dir, "stages", "input", "guide_tree.nwk"),
                gains=os.path.join(work, "gains.parquet"), reps=os.path.join(work, "reps.parquet"),
                queries=os.path.join(work, "queries.faa"), hits=os.path.join(work, "hits.tsv.gz"),
                per_query=os.path.join(work, "per_query.parquet"), top=os.path.join(work, "top_hits.parquet"),
                classified=os.path.join(work, "classified.parquet"),
                search_info=os.path.join(work, "search_info.json"),
                branches=os.path.join(work, "branches.tsv"))


def parse_run(p: dict):
    if p["orthoxml"] is None:
        raise SystemExit(f"no hogs.orthoxml in {p['run_dir']} (run not finished?)")
    t0 = time.time()
    res = parse_hogs(p["orthoxml"], newick=p["tree"], progress=lambda n: log(f"{n:,} families"))
    if not res.check_conservation().empty:
        raise SystemExit("conservation law violated (parser bug)")
    log(f"parsed {res.stats['n_genes']:,} genes, {res.stats['n_families']:,} families in {time.time() - t0:.0f}s")
    return res


def targets() -> pd.DataFrame:
    return pd.read_parquet(DB_PREFIX + ".targets.parquet")


def transcriptome_species() -> set:
    out = set()
    for f in ("species_traits.tsv",):
        t = pd.read_csv(os.path.join(ROOT, "data", f), sep="\t")
        out |= set(t.loc[t["data_type"].astype(str).str.contains("transcriptome"), "code"])
    return out


# ----------------------------------------------------------------------------

def cmd_makedb(a) -> None:
    out = gl.make_diamond_db(PROTEOMES, DB_PREFIX, threads=a.threads)
    t = targets()
    log(f"DB: {len(t):,} proteins, {t['species'].nunique()} species; {out}")


def cmd_prepare(a) -> None:
    p = paths(a)
    os.makedirs(p["work"], exist_ok=True)
    res = parse_run(p)
    sel = selected_branches(p)
    gains = gl.internal_gains(res, branches=sel)
    if sel is not None:
        log(f"restricted to {len(sel)} selected branches ({p['branches']})")
    lengths = targets().set_index("seqid")["length"]
    reps = gl.choose_representatives(res, gains, lengths, n_rep=a.n_rep, seed=SEED)
    gains.to_parquet(p["gains"], index=False)
    reps.to_parquet(p["reps"], index=False)
    n = gl.write_queries(reps, [DB_PREFIX + ".faa"], p["queries"])
    log(f"{len(gains):,} internal-branch gains, {n:,} representatives -> {p['queries']}")


def cmd_search(a) -> None:
    """One DIAMOND search of all representatives, or (``--two-pass``) rank-0 representatives
    first and the remaining representatives only for gains without an outside-clade hit in
    pass 1. Novelty labels are identical by construction (one outside hit suffices); split-
    repair scores of non-novel gains then rest on one representative."""
    p = paths(a)
    if not a.two_pass:
        info = gl.run_diamond(p["queries"], DB_PREFIX + ".dmnd", p["hits"] + ".part", threads=a.threads,
                              sensitivity=a.sensitivity, nice=a.nice, log=os.path.join(p["work"], "diamond.log"))
        os.replace(p["hits"] + ".part", p["hits"])
        info.update(sensitivity=a.sensitivity, nice=a.nice, threads=a.threads, two_pass=False,
                    n_queries=sum(1 for _ in gl.read_fasta(p["queries"])))
    else:
        reps = pd.read_parquet(p["reps"])
        gains = pd.read_parquet(p["gains"])
        seqs = dict(gl.read_fasta(p["queries"]))
        passes = []
        h1 = os.path.join(p["work"], "hits_pass1.tsv.gz")
        q1 = reps[reps["rep_rank"] == 0]
        f1 = os.path.join(p["work"], "queries_pass1.faa")
        with open(f1, "w") as fh:
            fh.writelines(f">{g}\n{seqs[g]}\n" for g in q1["gene"])
        if not os.path.exists(h1):
            i1 = gl.run_diamond(f1, DB_PREFIX + ".dmnd", h1 + ".part", threads=a.threads, sensitivity=a.sensitivity, nice=a.nice,
                                log=os.path.join(p["work"], "diamond_pass1.log"))
            os.replace(h1 + ".part", h1)
            i1["n_queries"] = len(q1)
            passes.append(i1)
        res = parse_run(p)
        ts = gl.TargetSpace.build(targets(), res)
        pq1, _ = gl.reduce_hits(h1, reps, gains, res.tree, ts, top_n=1)
        hit_fams = set(pq1.loc[pq1["n_outside"] > 0, "family_id"])
        q2 = reps[(reps["rep_rank"] > 0) & ~reps["family_id"].isin(hit_fams)]
        log(f"pass 1: {len(hit_fams):,}/{len(gains):,} gains have an outside hit; pass 2: {len(q2):,} queries")
        h2 = os.path.join(p["work"], "hits_pass2.tsv.gz")
        if len(q2):
            f2 = os.path.join(p["work"], "queries_pass2.faa")
            with open(f2, "w") as fh:
                fh.writelines(f">{g}\n{seqs[g]}\n" for g in q2["gene"])
            i2 = gl.run_diamond(f2, DB_PREFIX + ".dmnd", h2 + ".part", threads=a.threads, sensitivity=a.sensitivity, nice=a.nice,
                                log=os.path.join(p["work"], "diamond_pass2.log"))
            os.replace(h2 + ".part", h2)
            i2["n_queries"] = len(q2)
            passes.append(i2)
        with open(p["hits"] + ".part", "wb") as out:      # gzip members concatenate
            for h in (h1, h2):
                if os.path.exists(h):
                    with open(h, "rb") as fh:
                        shutil.copyfileobj(fh, out)
        os.replace(p["hits"] + ".part", p["hits"])
        info = dict(cmd=passes[0]["cmd"] if passes else "", wall_s=round(sum(x["wall_s"] for x in passes), 1),
                    passes=passes, sensitivity=a.sensitivity, nice=a.nice, threads=a.threads, two_pass=True,
                    n_queries=len(q1) + len(q2), n_queries_pass1=len(q1), n_queries_pass2=len(q2),
                    searched_genes=sorted(set(q1["gene"]) | set(q2["gene"])))
        with open(os.path.join(p["work"], "searched_genes.txt"), "w") as fh:
            fh.write("\n".join(info.pop("searched_genes")) + "\n")
    with open(p["search_info"], "w") as fh:
        json.dump(info, fh, indent=1)
    log(f"search done in {info['wall_s']:.0f}s")


def branch_classes(res, kingdom: str) -> pd.DataFrame:
    """Descriptive branch classes (not the §4.6 control pool): stem of a transition
    (M = MRCA of its multicellular species in this tree), inside a multicellular clade,
    root child, clade of unicellular species only, other."""
    T = res.tree
    tr = pd.read_csv(os.path.join(ROOT, "data", "transitions.tsv"), sep="\t")
    # only transitions whose M is a clade (M_kind node/mrca_node/internal); tip-level M (scattered tips) has no stem
    tr = tr[(tr["kingdom"] == kingdom) & tr["mc_species"].notna()
            & tr["M_kind"].isin(["node", "mrca_node", "internal"])]
    traits = pd.read_csv(os.path.join(ROOT, "data", "species_traits.tsv"), sep="\t")
    mc = set(traits.loc[traits["multicellular"] == 1, "code"])
    stems: dict[str, list] = {}
    for _, r in tr.iterrows():
        sps = [T.index[s] for s in str(r["mc_species"]).split(",") if s in T.index]
        if len(sps) < 2:
            continue
        m = T.mrca(min(sps), max(T.end[s] for s in sps))
        if m == 0:
            continue
        stems.setdefault(T.names[m], []).append((r["transition_id"], r["analysis_set"]))
    rows = []
    for c in T.internal:
        if c == 0:
            continue
        nm = T.names[c]
        leaves = [T.names[s] for s in T.leaves if T.is_anc(c, s)]
        n_mc = sum(s in mc for s in leaves)
        if nm in stems:
            cls = "stem"
        elif T.parent[c] == 0:
            cls = "root_child"
        elif n_mc == 0:
            cls = "unicellular_clade"
        elif n_mc == len(leaves):
            cls = "multicellular_clade"
        else:
            cls = "mixed_clade"
        rows.append(dict(branch_child=nm, branch_class=cls,
                         transitions=";".join(t for t, _ in stems.get(nm, [])),
                         analysis_set=";".join(s for _, s in stems.get(nm, [])),
                         n_tips=len(leaves), n_mc_tips=n_mc))
    return pd.DataFrame(rows)


def cmd_classify(a) -> None:
    p = paths(a)
    os.makedirs(p["res"], exist_ok=True)
    res = parse_run(p)
    gains = pd.read_parquet(p["gains"])
    reps = pd.read_parquet(p["reps"])
    ts = gl.TargetSpace.build(targets(), res, transcriptome_species())
    t0 = time.time()
    pq, top = gl.reduce_hits(p["hits"], reps, gains, res.tree, ts, top_n=a.top_n)
    log(f"reduced hits in {time.time() - t0:.0f}s ({int(pq['n_hits'].sum()):,} hits)")
    pq.to_parquet(p["per_query"], index=False)
    top.to_parquet(p["top"], index=False)
    cls = gl.classify_gains(res, gains, reps, pq, top, ts)
    cls.to_parquet(p["classified"], index=False)
    write_outputs(a, p, res, cls)


def cmd_arch(a) -> None:
    p = paths(a)
    res = parse_run(p)
    cls = pd.read_parquet(p["classified"]).drop(columns=["arch", "arch_novel", "pfam_outside_coverage"],
                                                 errors="ignore")
    files = {}
    for d in a.pfam_dirs.split(","):
        for f in sorted(glob.glob(os.path.join(ROOT, d, "*.pfam.tsv.gz"))):
            files[os.path.basename(f).split(".pfam")[0]] = f
    t = targets()
    all_sp = sorted(t["species"].astype(str).unique())
    log(f"Pfam tables for {len(files)}/{len(all_sp)} species")
    if not files:
        raise SystemExit("no Pfam outputs found; re-run `arch` when they exist")
    arch = gl.pfam_architectures(files)
    cls = gl.architecture_novelty(res, cls, arch, files.keys(), all_sp)
    cls.to_parquet(p["classified"], index=False)
    write_outputs(a, p, res, cls, pfam_info=dict(n_species_with_pfam=len(files), n_species_total=len(all_sp)))


def selected_branches(p) -> list | None:
    if os.path.exists(p["branches"]):
        return pd.read_csv(p["branches"], sep="\t")["branch"].tolist()
    return None


def write_outputs(a, p, res, cls, pfam_info=None) -> None:
    os.makedirs(p["res"], exist_ok=True)
    sel = selected_branches(p)
    br = gl.branch_ladder(res, cls, branches=sel)
    if sel is not None:
        why = pd.read_csv(p["branches"], sep="\t").set_index("branch")["reason"]
        br["selection_reason"] = why.reindex(br["branch_child"]).to_numpy()
    bc = branch_classes(res, a.run)
    br = br.merge(bc, on="branch_child", how="left")
    br.to_csv(os.path.join(p["res"], "branch_gain_ladder.tsv"), sep="\t", index=False)
    keep = [c for c in ["family_id", "branch_parent", "branch_child", "hog_id", "n_genes", "n_species",
                        "n_child_clades_present", "n_reps", "n_hits", "n_outside", "n_outside_robust",
                        "n_outside_transcriptome", "n_outside_species_max", "n_outside_other_domain", "best_out_target",
                        "best_out_species", "best_out_domain", "best_out_pident", "best_out_bitscore",
                        "best_parent_hog", "best_parent_hog_lost", "repair_hog", "repair_score",
                        "seq_novel", "seq_novel_tx_sens", "seq_novel_genomes_only", "outside_other_domain_only", "split_repaired",
                        "repair_status", "ladder", "arch", "arch_novel", "pfam_outside_coverage"] if c in cls.columns]
    cls[keep].to_csv(os.path.join(p["res"], "gains_classified.tsv.gz"), sep="\t", index=False)
    info = json.load(open(p["search_info"])) if os.path.exists(p["search_info"]) else {}
    write_summary(a, p, br, cls, info, pfam_info)
    done = dict(status="done", run=a.run, mcell_version=__version__, orthoxml=p["orthoxml"],
                n_gains=int(len(cls)), n_split_repaired=int(cls["split_repaired"].sum()),
                n_seq_novel=int(cls["seq_novel"].sum()), search=info, pfam=pfam_info,
                finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2))
    with open(os.path.join(p["res"], "DONE.json"), "w") as fh:
        json.dump(done, fh, indent=1, default=str)
    log(f"wrote {p['res']}")


def write_summary(a, p, br, cls, info, pfam_info) -> None:
    n = len(cls)
    rep = int(cls["split_repaired"].sum())
    nov = int(cls["seq_novel"].sum())
    novt = int(cls["seq_novel_tx_sens"].sum())
    novg = int(cls["seq_novel_genomes_only"].sum())
    hom = int((cls["ladder"] == "homologous").sum())
    od = int(cls["outside_other_domain_only"].sum())
    pct = lambda x: f"{100 * x / n:.1f}%" if n else "NA"
    L = [f"# Gain-evidence ladder: {a.run} run", "",
         f"Generated by `scripts/gain_ladder.py` (mcell {__version__}); definitions in `mcell/gainladder.py` "
         "and DESIGN §3.2. Descriptive only (no hypothesis tests).", "",
         f"* HOGs: `{os.path.relpath(p['orthoxml'], ROOT)}`; tree `{os.path.relpath(p['tree'], ROOT)}`.",
         (f"* Scope: all {len(br)} internal branches (root stem and leaf branches excluded)." if "selection_reason" not in br
          else f"* Scope: {len(br)} selected branches (`{os.path.relpath(p['branches'], ROOT)}`: transition stems, path "
               "edges, grey/reversal stems, control pool U, and the edges above them; tip stems included)."),
         f"* DIAMOND: `{info.get('cmd', 'n/a').replace(ROOT + '/', '')}`; {info.get('n_queries', 'NA')} "
         f"queries, wall {info.get('wall_s', 'NA')} s at {info.get('threads', 'NA')} threads. Target DB: all "
         "eukaryote + prokaryote proteomes (identical for both runs)."
         + (f" Two-pass: {info.get('n_queries_pass1')} rank-0 representatives, then {info.get('n_queries_pass2')} "
            "further representatives of gains still without an outside hit." if info.get("two_pass") else ""),
         "* Deviation from DESIGN §3.2: `--very-sensitive` instead of `--ultra-sensitive` (timing pilot: ~3.5x "
         "cheaper; on the simulations it calls 1-3% more gains novel; see `recheck_ultra.json` for the real-data check)."
         if info.get("sensitivity") == "very-sensitive" else "", ""]
    rc = os.path.join(p["res"], "recheck_ultra.json")
    if os.path.exists(rc):
        r = json.load(open(rc))
        L += [f"* Ultra-sensitive recheck: of {r['n_sampled']} random seq-novel gains (seed 11), "
              f"{r['n_gain_outside_hit_with_ultra']} ({100 * (r['frac_gain_outside_hit_with_ultra'] or 0):.1f}%) gain an "
              "outside-clade hit with `--ultra-sensitive`.", ""]
    L += ["* **Split repair failed the DESIGN acceptance test on the simulations** (gain ratio stays 2.1-2.5x; "
          "`results/gainladder/validation_sim/summary.md`): level 2 is a diagnostic only; report levels 1 and 3. "
          "On the simulations seq-novel recall of true gains is 0.91-0.97 but precision only 0.42-0.56 "
          "(deep splits beyond detectable homology), so seq-novel counts are upper bounds.", "",
          "## Ladder totals (all selected branches, including terminal branches)", "",
         "| level | n | % of HOG gains |", "|---|---|---|",
         f"| 1. HOG-level gains | {n} | 100% |",
         f"| 2. split-repaired (gain re-merged with a HOG lost on the same branch) | {rep} | {pct(rep)} |",
         f"| homologous outside the clade, not repaired | {hom} | {pct(hom)} |",
         f"| 3. sequence-novel (no outside-clade hit) | {nov} | {pct(nov)} |",
         f"| 3'. sequence-novel ignoring transcriptome hits at >=90% id (sensitivity) | {novt} | {pct(novt)} |",
         f"| 3''. sequence-novel ignoring ALL hits to transcriptome-only species (contamination sensitivity) | {novg} | {pct(novg)} |",
         f"| non-novel with outside hits only in the other domain (HGT candidates) | {od} | {pct(od)} |"]
    if "arch_novel" in cls.columns:
        sn = cls[cls["seq_novel"]]
        L.append(f"| 4. architecture-novel (of seq-novel with a Pfam architecture: "
                 f"{int(sn['arch_novel'].notna().sum())}) | {int((sn['arch_novel'] == True).sum())} | "  # noqa: E712
                 f"{pct(int((sn['arch_novel'] == True).sum()))} |")  # noqa: E712
        if pfam_info:
            L.append(f"\nPfam tables available for {pfam_info['n_species_with_pfam']}/"
                     f"{pfam_info['n_species_total']} species; outside-clade coverage per gain: median "
                     f"{cls['pfam_outside_coverage'].median():.2f}.")
    else:
        L.append("| 4. architecture-novel | not computed (Pfam outputs missing; run `gain_ladder.py arch`) | |")
    L += ["", "## By branch class (descriptive)", "",
          "Classes from `data/transitions.tsv` + `data/species_traits.tsv` (stem = branch into the MRCA of a "
          "transition's multicellular species in this tree; NOT the §4.6 matched control pool).", "",
          "| class | branches | HOG gains | split-repaired | seq-novel | frac seq-novel (pooled) | median frac per branch |",
          "|---|---|---|---|---|---|---|"]
    for k, d in br.groupby("branch_class"):
        g = int(d["hog_gains"].sum())
        L.append(f"| {k} | {len(d)} | {g} | {int(d['split_repaired'].sum())} | {int(d['seq_novel'].sum())} | "
                 f"{d['seq_novel'].sum() / g if g else float('nan'):.3f} | {d['frac_seq_novel'].median():.3f} |")
    st = br[br["branch_class"] == "stem"]
    if len(st):
        L += ["", "### Stems", "", "| branch | transitions | analysis set | tips | HOG gains | split-repaired | "
              "seq-novel | frac |", "|---|---|---|---|---|---|---|---|"]
        for _, r in st.iterrows():
            L.append(f"| {r.branch_parent}->{r.branch_child} | {r.transitions} | {r.analysis_set} | {r.n_tips} | "
                     f"{r.hog_gains} | {r.split_repaired} | {r.seq_novel} | {r.frac_seq_novel:.3f} |")
    L += ["", "Per-branch table: `branch_gain_ladder.tsv`; per-gain labels: `gains_classified.tsv.gz`.", ""]
    with open(os.path.join(p["res"], "summary.md"), "w") as fh:
        fh.write("\n".join(L))


def cmd_recheck(a) -> None:
    """Re-search the representatives of a random sample (seed 11) of sequence-novel gains
    with --ultra-sensitive; reports the fraction that then gain an outside-clade hit
    (sensitivity loss of the main search). Writes recheck_ultra.json in the results dir."""
    p = paths(a)
    cls = pd.read_parquet(p["classified"])
    reps = pd.read_parquet(p["reps"])
    gains = pd.read_parquet(p["gains"])
    nov = sorted(cls.loc[cls["seq_novel"], "family_id"])
    rng = np.random.default_rng(SEED)
    pick = set(rng.choice(nov, size=min(a.n, len(nov)), replace=False)) if nov else set()
    r = reps[reps["family_id"].isin(pick)]
    seqs = dict(gl.read_fasta(p["queries"]))
    q = os.path.join(p["work"], "recheck_queries.faa")
    with open(q, "w") as fh:
        fh.writelines(f">{g}\n{seqs[g]}\n" for g in r["gene"])
    h = os.path.join(p["work"], "recheck_ultra.tsv.gz")
    info = gl.run_diamond(q, DB_PREFIX + ".dmnd", h, threads=a.threads, sensitivity="ultra-sensitive", nice=a.nice,
                          log=h + ".log")
    res = parse_run(p)
    ts = gl.TargetSpace.build(targets(), res, transcriptome_species())
    pq, _ = gl.reduce_hits(h, r.reset_index(drop=True), gains, res.tree, ts, top_n=1)
    fam_hit = pq.groupby("family_id")["n_outside"].sum() > 0
    out = dict(n_novel=len(nov), n_sampled=len(pick), n_queries=len(r), wall_s=info["wall_s"], threads=a.threads,
               n_gain_outside_hit_with_ultra=int(fam_hit.sum()),
               frac_gain_outside_hit_with_ultra=round(float(fam_hit.mean()), 4) if len(fam_hit) else None,
               cmd=info["cmd"].replace(ROOT + "/", ""))
    os.makedirs(p["res"], exist_ok=True)
    with open(os.path.join(p["res"], "recheck_ultra.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    log(json.dumps(out))


def branch_table(kingdom: str):
    import warnings
    from mcell import analysis as an
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return an.load_setup(kingdom)[3]


def cmd_select_branches(a) -> None:
    """Branches the tests use (mcell.branchstats.build_branch_table on the guide tree):
    every branch carrying a transition (stem, nested/grey/reversal stem, path-stem edge;
    tip stems too unless --no-tips), every control-pool branch (in_U), and the edge
    above each of those (not the root stem). Written to runs/gainladder/<run>/branches.tsv."""
    p = paths(a)
    os.makedirs(p["work"], exist_ok=True)
    bt = branch_table(a.run)
    why: dict[str, list] = {}
    for b, r in bt.iterrows():
        if r["role"] and not (a.no_tips and r["is_terminal"]):
            why.setdefault(b, []).append(f"transition:{r['role']}:{r['transition_id']}")
        if r["in_U"]:
            why.setdefault(b, []).append("control_pool_U")
    for b in list(why):
        par = bt.at[b, "parent"]
        if par in bt.index:                          # root has no row (no edge above it)
            why.setdefault(par, []).append(f"parent_edge_of:{b}")
    df = pd.DataFrame({"branch": list(why), "reason": [";".join(v) for v in why.values()]})
    df["is_terminal"] = bt.loc[df["branch"], "is_terminal"].to_numpy()
    df["n_desc"] = bt.loc[df["branch"], "n_desc"].to_numpy()
    df.sort_values("branch").to_csv(p["branches"], sep="\t", index=False)
    log(f"{len(df)} branches ({int(df['is_terminal'].sum())} terminal) -> {p['branches']}")


def cmd_plan(a) -> None:
    """Query-count / CPU estimate for the selected branches, from the prokaryote run:
    gains per internal branch = mean proteome size of the clade x exp(b0) x n_desc^b1
    (log-log fit on prok internal non-root-child branches); terminal branches =
    leaf-gain families / proteome (prok); queries = 2.0 x gains (prok two-pass ratio)."""
    import warnings
    from mcell.hogs import parse_hogs as _ph
    pp = os.path.join(ROOT, "runs", "prok_fast")
    res = _ph(os.path.join(pp, "hogs.orthoxml"), newick=os.path.join(pp, "stages", "input", "guide_tree.nwk"))
    T = res.tree
    tg = targets()
    psize = tg.groupby("species", observed=True).size()
    rows = []
    for _, r in res.branches[res.branches["branch_parent"] != ""].iterrows():
        c = T.index[r["branch_child"]]
        lv = [T.names[x] for x in T.leaves if T.is_anc(c, x)]
        rows.append(dict(leaf=r["child_is_leaf"], root_child=T.parent[c] == 0, n_desc=len(lv),
                         mprot=float(psize.reindex(lv).mean()), gains=r["gains"]))
    d = pd.DataFrame(rows)
    i = d[~d["leaf"] & ~d["root_child"]]
    b1, b0 = np.polyfit(np.log(i["n_desc"]), np.log(i["gains"] / i["mprot"]), 1)
    leaf_rate = float((d.loc[d["leaf"], "gains"] / d.loc[d["leaf"], "mprot"]).median())
    reps = pd.read_parquet(os.path.join(ROOT, "runs", "gainladder", "prok", "reps.parquet"))
    si = json.load(open(os.path.join(ROOT, "runs", "gainladder", "prok", "search_info.json")))
    q_per_gain = si["n_queries"] / reps["family_id"].nunique()
    cpu_prok = a.threads_prok * si["wall_s"] / si["n_queries"]
    sel = pd.read_csv(paths(a)["branches"], sep="\t")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from mcell import analysis as an
        tree = an.load_setup(a.run)[2]
    est = []
    for _, r in sel.iterrows():
        nd = tree.search_nodes(name=r["branch"])[0]
        lv = nd.get_leaf_names()
        mp = float(psize.reindex(lv).mean())
        g = mp * leaf_rate if r["is_terminal"] else mp * np.exp(b0) * len(lv) ** b1
        est.append(dict(branch=r["branch"], terminal=r["is_terminal"], n_desc=len(lv), mean_proteome=round(mp),
                        est_gains=round(g), reason=r["reason"]))
    e = pd.DataFrame(est)
    G = e["est_gains"].sum()
    Gi, Gt = e.loc[~e["terminal"], "est_gains"].sum(), e.loc[e["terminal"], "est_gains"].sum()
    Q = G * q_per_gain
    lo, hi = Q * cpu_prok / 3600, Q * 1.1 / 3600
    allq = None
    L = ["# Eukaryote gain-ladder plan (restricted to the branches the tests use)", "",
         f"Generated by `scripts/gain_ladder.py plan` before the eukaryote HOGs exist (PREREG Dev-2/Dev-3).", "",
         f"* Branches: {len(sel)} ({int((~sel['is_terminal']).sum())} internal, {int(sel['is_terminal'].sum())} "
         "terminal tip stems), from `mcell.branchstats.build_branch_table` on `data/guide_tree_euk.nwk` (identical "
         "to the Ancestra input tree): every transition branch in `data/transitions.tsv` (stems, nested/grey/reversal "
         "stems, path-stem edges; primary, secondary, exploratory, reversal), every control-pool branch U, and the "
         "edge above each. List: `runs/gainladder/euk/branches.tsv`.",
         f"* Model (prokaryote run): gains on an internal branch = mean proteome size of the clade x "
         f"{np.exp(b0):.4f} x n_desc^{b1:.2f} (log-log fit, {len(i)} prok internal non-root-child branches); tip stems = "
         f"{leaf_rate:.3f} x proteome (median prok leaf); {q_per_gain:.2f} queries per gain (prok two-pass); "
         f"{cpu_prok:.2f}-1.1 CPU-s per query (prok actual - random-protein pilot, very-sensitive).",
         f"* **Estimated gains: {G:,.0f}** ({Gi:,.0f} internal, {Gt:,.0f} on tip stems); **queries ~{Q:,.0f}**.",
         f"* **CPU: {lo:,.0f}-{hi:,.0f} CPU-h = {lo / 20:.1f}-{hi / 20:.1f} h at 20 threads** (nice 5). Uncertainty "
         "is at least x2 either way: eukaryote HOG splitting at the 7,500-gene family cap and deep, long stems "
         "probably raise the gain rate above the prokaryote model.",
         "* Output size: ~157 hits per pass-1 query (prok) to 640 (random eukaryote proteins): "
         f"{Q * 157 / 1e6:,.0f}-{Q * 640 / 1e6:,.0f} M rows, {Q * 157 * 16.8 / 1e9:.1f}-{Q * 640 * 16.8 / 1e9:.1f} GB gzipped.",
         "", "## Launch", "",
         "`runs/gainladder/euk/run_euk.sh` (started with nohup) waits for `runs/euk_fast/DONE.json` or "
         "`runs/euk_fast_p/DONE.json`, uses the first one present (records which; if both exist it checks that the "
         "two `hogs.orthoxml` are byte-identical), then runs `prepare` -> `search --two-pass --threads 20 --nice 5` -> "
         "`classify` (-> `arch` if Pfam tables exist) and writes `results/gainladder/euk/DONE.json`. Log: "
         "`runs/gainladder/euk/run_euk.log`.", "",
         "## Largest estimated branches", "", "| branch | terminal | n_desc | mean proteome | est. gains | reason |",
         "|---|---|---|---|---|---|"]
    for _, r in e.sort_values("est_gains", ascending=False).head(15).iterrows():
        L.append(f"| {r.branch} | {r.terminal} | {r.n_desc} | {r.mean_proteome} | {r.est_gains} | {r.reason[:80]} |")
    os.makedirs(os.path.join(ROOT, "results", "gainladder"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "gainladder", f"{a.run}_plan.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    e.to_csv(os.path.join(paths(a)["work"], "plan_estimate.tsv"), sep="\t", index=False)
    log(f"estimate: {G:,.0f} gains, {Q:,.0f} queries, {lo:.0f}-{hi:.0f} CPU-h")


def cmd_pilot(a) -> None:
    work = os.path.join(ROOT, "runs", "gainladder", "pilot")
    os.makedirs(work, exist_ok=True)
    q = os.path.join(work, f"q{a.n}.faa")
    if not os.path.exists(q):
        random.seed(SEED)
        res, k = [], 0
        for f in sorted(glob.glob(os.path.join(PROTEOMES["euk"], "*.fa"))):
            for h, s in gl.read_fasta(f):
                k += 1
                if len(res) < a.n:
                    res.append((h, s))
                elif (j := random.randrange(k)) < a.n:
                    res[j] = (h, s)
        with open(q, "w") as fh:
            fh.writelines(f">{h}\n{s}\n" for h, s in res)
    out = os.path.join(work, f"{a.sensitivity}_q{a.n}.tsv.gz")
    info = gl.run_diamond(q, DB_PREFIX + ".dmnd", out, threads=a.threads, sensitivity=a.sensitivity, nice=a.nice,
                          log=out + ".log")
    log(json.dumps(info))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["makedb", "select-branches", "plan", "prepare", "search", "classify", "arch", "all",
                                       "pilot", "recheck"])
    ap.add_argument("--run", default="prok", help="run name (euk|prok); also the kingdom in transitions.tsv")
    ap.add_argument("--run-dir", default=None)
    ap.add_argument("--db", default=None, help="DB prefix (default runs/gainladder/db/all_euk_prok; for testing)")
    ap.add_argument("--threads", type=int, default=2,
                    help="DIAMOND threads (PI cap while the eukaryote Ancestra run is active: 2)")
    ap.add_argument("--nice", type=int, default=19, help="run DIAMOND under nice -n (default 19)")
    ap.add_argument("--sensitivity", default="very-sensitive",
                    choices=["fast", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"])
    ap.add_argument("--n-rep", type=int, default=3)
    ap.add_argument("--no-tips", action="store_true", help="select-branches: leave out terminal (tip) stems")
    ap.add_argument("--threads-prok", type=int, default=2, help="plan: threads of the prok search (CPU-s/query)")
    ap.add_argument("--two-pass", action="store_true",
                    help="search rank-0 representatives first, other representatives only for still-novel gains")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--n", type=int, default=500,
                    help="pilot: number of random queries; recheck: number of novel gains sampled")
    ap.add_argument("--pfam-dirs", default=None,
                    help="comma-separated Pfam output dirs (default data/annotations/pfam/{euk,prok})")
    a = ap.parse_args(argv)
    global DB_PREFIX
    if a.db:
        DB_PREFIX = a.db
    if a.pfam_dirs is None:
        a.pfam_dirs = "data/annotations/pfam/euk,data/annotations/pfam/prok"
    if a.cmd == "makedb":
        cmd_makedb(a)
    elif a.cmd == "prepare":
        cmd_prepare(a)
    elif a.cmd == "search":
        cmd_search(a)
    elif a.cmd == "classify":
        cmd_classify(a)
    elif a.cmd == "arch":
        cmd_arch(a)
    elif a.cmd == "pilot":
        cmd_pilot(a)
    elif a.cmd == "recheck":
        cmd_recheck(a)
    elif a.cmd == "select-branches":
        cmd_select_branches(a)
    elif a.cmd == "plan":
        cmd_plan(a)
    elif a.cmd == "all":
        a.two_pass = True
        cmd_prepare(a)
        cmd_search(a)
        cmd_classify(a)
        if any(glob.glob(os.path.join(ROOT, d, "*.pfam.tsv.gz")) for d in a.pfam_dirs.split(",")):
            cmd_arch(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
