#!/usr/bin/env python
"""EXPLORATORY. Reversals to unicellularity (DESIGN §4.2 secondary/exploratory list: "genes gained at the earlier
transition lost preferentially").

    python scripts/analysis/reversals.py --run euk
    python scripts/analysis/reversals.py --run prok --pairs Synechococcales:T_cyano_filament   # pseudo-reversal demo

For every row of data/transitions.tsv with to_state 'reversal' in this kingdom (yeasts: Saccharomycotina,
Schizosaccharomyces, Saitoella, Cryptococcus, Rhodotorula; Zygnematophyceae): reversal branch = rev_parent ->
rev_child (rev_parent = first node of stem_path, else parent(M); rev_child = M, a node or a tip); preceding
transition = the nearest transition whose M is an ancestor of the reversal's M and whose role is an origin
(unicellular -> simple/complex), its D7 path edges included. --pairs <node>:<transition_id>,... adds or replaces
pairs (pseudo-reversals, labelled as such).
HOGs present at rev_parent are classified as stem-gained (an ancestor, or itself, was gained on a stem edge of the
preceding transition), stem-duplicated (an ancestor is the parent HOG of a non-origin duplication on a stem edge)
or other; lost = no descendant at rev_child. Fisher exact test (two-sided and one-sided 'greater' = preferential
loss), odds ratio with Haldane/Woolf 95 % CI, for gained, duplicated and gained-or-duplicated vs all other HOGs at
rev_parent. Within-tree control: the same test on every sister branch rev_parent -> s (lineages that kept the
multicellular state); ratio of odds ratios = OR(reversal) / geometric mean OR(sisters in the multicellular state;
all sisters if none). Losses on transcriptome-only tips are uncertain (noted).

Outputs: results/exploratory/<run>/reversals/: reversal_tests.tsv, summary.json; registry rows (gained-or-dup,
one-sided p) for every real reversal (pseudo-reversals only with --register-pseudo).
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
    print(f"[reversals {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def preceding_transition(tr, tree, bt, m_node):
    """Nearest ancestral origin transition (role 'stem') whose M is an ancestor of m_node."""
    nodes = {n.name: n for n in tree.traverse()}
    anc = []
    x = nodes[m_node].up
    while x is not None:
        anc.append(x.name)
        x = x.up
    best = None
    for _, r in tr.iterrows():
        if str(r["to_state"]) not in ("simple", "complex") or str(r["from_state"]) != "unicellular":
            continue
        Ms = ex._split(r["M"])
        if len(Ms) != 1 or Ms[0] not in anc:
            continue
        d = anc.index(Ms[0])
        if best is None or d < best[0]:
            best = (d, r["transition_id"], Ms[0])
    return (best[1], best[2]) if best else (None, None)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--pairs", default="", help="node:transition_id[,node:transition_id] (pseudo-reversals)")
    ap.add_argument("--register-pseudo", action="store_true")
    ap.add_argument("--no-registry", action="store_true")
    ap.add_argument("--skip-prereg-check", action="store_true")
    a = ap.parse_args(argv)
    if not a.skip_prereg_check and not ex.prereg_ok(ROOT):
        sys.exit("PREREG hash mismatch: aborting")
    ctx = ex.load_context(a.run, "default", need_cache=True)
    c = ctx.cache
    od = ex.out_dir(a.run, "reversals")
    tr = ctx.tr
    tree, bt = ctx.tree, ctx.bt
    nodes = {n.name: n for n in tree.traverse()}
    names = list(c.tree.names)
    is_leaf = list(c.tree.is_leaf)
    traits = ctx.traits.set_index("code") if "code" in ctx.traits.columns else None

    mcs = ctx.traits.set_index("code")["mc_state"].astype(str).to_dict() if "mc_state" in ctx.traits.columns else {}

    def tip_state(node):
        """'multicellular' if >= 50 % of the tips below have trait mc_state simple/complex, else 'unicellular'."""
        t = [x.name for x in nodes[node].get_leaves()]
        f = np.mean([mcs.get(x, "") in ("simple", "complex") for x in t])
        return "multicellular" if f >= 0.5 else "unicellular"

    jobs = []
    for _, r in tr[tr["to_state"] == "reversal"].iterrows():
        Ms = ex._split(r["M"])
        if len(Ms) != 1 or Ms[0] not in nodes:
            continue
        m = Ms[0]
        sp = r.get("stem_path")
        rp = sp.split(">")[0].strip() if isinstance(sp, str) and ">" in sp else nodes[m].up.name
        tid, pm = preceding_transition(tr, tree, bt, m)
        jobs.append(dict(reversal=r["transition_id"], rev_parent=rp, rev_child=m, preceding=tid, preceding_M=pm,
                         pseudo=False))
    for item in [x for x in a.pairs.split(",") if x]:
        node, tid = item.split(":")
        rp = nodes[node].up.name
        pm = ex._split(tr.set_index("transition_id").loc[tid, "M"])[0]
        jobs = [j for j in jobs if j["rev_child"] != node]
        jobs.append(dict(reversal=f"pseudo:{node}", rev_parent=rp, rev_child=node, preceding=tid, preceding_M=pm,
                         pseudo=True))
    if not jobs:
        log("no reversal rows for this kingdom and no --pairs: nothing to do")
        ex.write_summary(od, dict(run=a.run, n_reversals=0, note="no reversal rows in data/transitions.tsv"))
        return 0

    rows = []
    for j in jobs:
        if j["preceding"] is None:
            log(f"{j['reversal']}: no preceding origin transition found; skipped")
            continue
        edges = bs.path_edges(bt, j["preceding_M"]) if j["preceding_M"] in bt.index else [j["preceding_M"]]
        if j["rev_parent"] not in names or j["rev_child"] not in names:
            log(f"{j['reversal']}: node missing from the Ancestra tree; skipped")
            continue
        sisters = [ch.name for ch in nodes[j["rev_parent"]].children if ch.name != j["rev_child"]]
        tips = [x.name for x in nodes[j["rev_child"]].get_leaves()]
        tx_only = traits is not None and all(str(traits["data_type"].get(t, "genome")) == "transcriptome" for t in tips)
        for role, child in [("reversal", j["rev_child"])] + [("sister_control", s) for s in sisters]:
            res = ex.reversal_test(c.hogs, c.genes, c.events, names, is_leaf, edges, j["rev_parent"], child)
            for grp in ("gained", "dup", "gained_or_dup"):
                g = res[grp]
                rows.append(dict(reversal=j["reversal"], pseudo=j["pseudo"], role=role, branch_child=child,
                                 child_state=tip_state(child),
                                 rev_parent=j["rev_parent"], preceding=j["preceding"], preceding_edges=";".join(edges),
                                 group=grp, n_hogs_parent=res["n_hogs_parent"], n_lost=res["n_lost"],
                                 transcriptome_only_child=bool(tx_only) if role == "reversal" else None, **g))
        log(f"{j['reversal']}: {len(sisters)} sister control(s)")
    R = pd.DataFrame(rows)
    # ratio of odds ratios: reversal vs mean log OR of its sister controls
    R["log_or"] = np.log(R["odds_ratio"])
    sis = R[R["role"] == "sister_control"].copy()
    mc = sis.groupby("reversal")["child_state"].apply(lambda x: (x == "multicellular").any())
    sis = sis[[(st == "multicellular") or not mc.get(rv, False) for rv, st in zip(sis["reversal"], sis["child_state"])]]
    ref = sis.groupby(["reversal", "group"])["log_or"].mean().rename("sister_mean_log_or")
    R = R.join(ref, on=["reversal", "group"])
    R["sister_reference"] = R["reversal"].map(lambda r: "multicellular sisters" if mc.get(r, False) else "all sisters")
    R["ratio_of_or_vs_sisters"] = np.where(R["role"] == "reversal", np.exp(R["log_or"] - R["sister_mean_log_or"]), np.nan)
    sel = (R["role"] == "reversal")
    R.loc[sel, "q_bh"] = ex.bh(R.loc[sel, "p_greater"].to_numpy())
    R.to_csv(os.path.join(od, "reversal_tests.tsv"), sep="\t", index=False)
    if not a.no_registry:
        for _, r in R[(R["role"] == "reversal") & (R["group"] == "gained_or_dup")].iterrows():
            if r["pseudo"] and not a.register_pseudo:
                continue
            bs.append_registry(ctx.paths["registry"], ex.registry_row(
                "reversal: preferential loss of stem-gained/duplicated HOGs",
                f"Fisher exact one-sided (greater); odds ratio lost|stem-gained-or-dup vs other HOGs at "
                f"{r['rev_parent']}; reversal {r['branch_child']} after {r['preceding']}",
                n=int(r["a_lost"] + r["b_kept"] + r["c_bg_lost"] + r["d_bg_kept"]), estimate=r["odds_ratio"],
                ci=(r["or_lo"], r["or_hi"]), p=r["p_greater"],
                notes=f"run={a.run}; reversal={r['reversal']}; p_two_sided={r['p_two_sided']:.3g}; "
                      f"q_BH(within reversals)={r['q_bh']:.3g}; ratio_of_OR_vs_sisters={r['ratio_of_or_vs_sisters']:.3g}"
                      + ("; PSEUDO-REVERSAL demo" if r["pseudo"] else "")
                      + ("; child transcriptome-only: losses uncertain" if r["transcriptome_only_child"] else ""),
                inputs=[ctx.paths["transitions"]]))
    summ = dict(run=a.run, n_reversals=int(R.loc[R["role"] == "reversal", "reversal"].nunique()),
                tests=R.round(5).replace({np.nan: None}).to_dict("records"))
    ex.write_summary(od, summ)
    log(f"done -> {os.path.relpath(od, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
