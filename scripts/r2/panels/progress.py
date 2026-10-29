#!/usr/bin/env python
"""Rewrite results/r2/panels/PROGRESS.txt (interim progress; the only file that is rewritten)."""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, PANELS, ROOT  # noqa: E402

import pandas as pd  # noqa: E402

lines = [f"ROUND-2 DENSE PANELS - PROGRESS (rewritten {time.strftime('%Y-%m-%d %H:%M %Z')})", ""]
for R in PANELS:
    o = os.path.join(OUT, R)
    anc = os.path.exists(os.path.join(ROOT, "runs", "r2", R, "DONE.json"))
    steps = [(n, os.path.exists(os.path.join(o, f))) for n, f in
             (("parse", "branch_events.tsv"), ("ladder", "branch_gain_ladder.tsv"), ("branch lengths", "tree_bl.nwk"),
              ("tables", "event_table.tsv"), ("controls", "positive_controls.tsv"), ("W3", "w3_families.tsv.gz"))]
    lines.append(f"{R}: Ancestra {'done' if anc else 'running/queued'}; " + ", ".join(f"{n} {'done' if d else '-'}" for n, d in steps))
    if os.path.exists(os.path.join(o, "ladder_info.json")):
        li = json.load(open(os.path.join(o, "ladder_info.json")))
        lines.append(f"    ladder: {li['n_gains']:,} HOG gains on {li['n_branches']} branches; {li['n_seq_novel']:,} sequence-novel")
    if os.path.exists(os.path.join(o, "event_table.tsv")):
        e = pd.read_csv(os.path.join(o, "event_table.tsv"), sep="\t")
        b = pd.read_csv(os.path.join(o, "branch_table.tsv"), sep="\t")
        c = b[b["control"]]
        lines.append(f"    controls: {len(c)} internal non-crown branches; median seq-novel gains {c['seq_novel'].median():.0f} "
                     f"(range {c['seq_novel'].min():.0f}-{c['seq_novel'].max():.0f})")
        for r in e.itertuples():
            lines.append(f"    {r.event_id:32s} {r.kind:15s} M={r.M:30s} raw {r.hog_gains:6.0f} supported "
                         f"{r.supported_gains if r.supported_gains == r.supported_gains else float('nan'):6.0f} seq-novel {r.seq_novel:6.0f} dups {r.dups:5.0f}")
    if os.path.exists(os.path.join(o, "positive_controls.tsv")):
        pc = pd.read_csv(os.path.join(o, "positive_controls.tsv"), sep="\t")
        lines.append("    positive controls: " + ", ".join(f"{r.control} {r.panel_recovered} (r1 {r.round1_recovered})" for r in pc.itertuples()))
pdir = os.path.join(OUT, "pooled")
if os.path.isdir(pdir):
    for t in sorted(os.listdir(pdir)):
        f = os.path.join(pdir, t, "p1_and_meta.json")
        if os.path.exists(f):
            d = json.load(open(f))
            lines.append("")
            lines.append(f"pooled/{t}: panels {d['panels']}")
            if "p1" in d:
                lines.append(f"    P1: {d['p1']}")
                lines.append(f"    P1 meta: {d.get('p1_meta')}")
            lines.append(f"    secondary all nested: {d.get('secondary_all_nested')}")
lines.append("")
lines.append("CAVEAT (interim): tip events (species-specific orphans) exceed internal controls by 4-235x regardless of event type "
             "(reversals too); internal-branch-only and tip-vs-tip-control estimates are the interpretable ones (DEVIATIONS D-P1, D-P6).")
lines.append("Residuals, P1 and pooled estimates are computed only by scripts/r2/panels/pooled.py (interim tags are labelled).")
with open(os.path.join(OUT, "PROGRESS.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
