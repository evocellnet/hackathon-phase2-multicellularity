#!/usr/bin/env python
"""Step 2: concatenated supermatrix -> IQ-TREE LG+G4 branch lengths on the FIXED guide topology.

DESIGN.md §3.1. ``iqtree3 -s supermatrix -te guide -m LG+G4``: IQ-TREE 3.1.4
accepts multifurcating fixed trees and keeps the polytomies (checked on a toy
tree: no zero-length resolution, no constraint search needed). Only branch
lengths and model parameters are optimised; the topology is never changed.
Lengths are transferred back to the named guide tree by bipartition matching
(mcell.phylo.transfer_branch_lengths; a bifurcating root edge is split
equally). Output: data/species_tree_bl_<set>.nwk (guide internal names kept),
runs/phylo/<set>/concat/ (supermatrix, partition file, IQ-TREE files) and
runs/phylo/<set>/species_tree_bl.DONE.json.

All trimmed markers of step 1 are concatenated (missing = gaps); taxa with low
coverage are kept (the supermatrix tolerates missing data).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
os.environ["PATH"] = os.path.join(ROOT, "conda_envs", "ancestra", "bin") + os.pathsep + os.environ["PATH"]

from mcell.phylo import read_tree, transfer_branch_lengths, write_named, assert_same_topology  # noqa: E402
from scripts.phylo.build_markers import read_fasta  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["euk", "prok"], required=True)
    ap.add_argument("--threads", default="4")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    t0 = time.time()
    base = os.path.join(ROOT, "runs", "phylo", a.set)
    mdir, cdir = os.path.join(base, "markers"), os.path.join(base, "concat")
    os.makedirs(cdir, exist_ok=True)
    guide = read_tree(os.path.join(ROOT, "data", f"guide_tree_{a.set}.nwk"))
    taxa = guide.get_leaf_names()
    with open(os.path.join(mdir, "markers.tsv")) as fh:
        markers = [r["marker"] for r in csv.DictReader(fh, delimiter="\t") if int(r["trim_len"]) > 0]
    seqs = {t: [] for t in taxa}
    parts, pos = [], 1
    for m in markers:
        al = read_fasta(os.path.join(mdir, f"{m}.trim.fa"))
        L = len(next(iter(al.values())))
        for t in taxa:
            seqs[t].append(al.get(t, "-" * L))
        parts.append((m, pos, pos + L - 1))
        pos += L
    sm = os.path.join(cdir, "supermatrix.fa")
    with open(sm, "w") as fh:
        for t in taxa:
            fh.write(f">{t}\n{''.join(seqs[t])}\n")
    with open(os.path.join(cdir, "partitions.nex"), "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for m, s, e in parts:
            fh.write(f"  charset {m} = {s}-{e};\n")
        fh.write("end;\n")
    nsites = pos - 1
    occ = {t: sum(ch not in "-X?" for ch in "".join(seqs[t])) / nsites for t in taxa}
    with open(os.path.join(cdir, "taxon_occupancy.tsv"), "w") as fh:
        fh.write("taxon\tfrac_sites_present\n")
        for t in taxa:
            fh.write(f"{t}\t{occ[t]:.4f}\n")
    te = os.path.join(cdir, "guide_topology.nwk")
    with open(te, "w") as fh:
        fh.write(guide.write(format=9) + "\n")
    prefix = os.path.join(cdir, "concat")
    cmd = ["iqtree3", "-s", sm, "-te", te, "-m", "LG+G4", "-T", str(a.threads), "-seed", str(a.seed),
           "-keep-ident", "--prefix", prefix, "-redo", "-quiet"]
    print(" ".join(cmd), flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"iqtree3 failed: {r.stderr[-3000:]}")
    fitted = read_tree(prefix + ".treefile")
    lengths = transfer_branch_lengths(guide, fitted)
    out_nwk = write_named(guide, lengths)
    assert_same_topology(guide, read_tree(out_nwk))
    dst = os.path.join(ROOT, "data", f"species_tree_bl_{a.set}.nwk")
    with open(dst, "w") as fh:
        fh.write(out_nwk + "\n")
    lnl = None
    for line in open(prefix + ".iqtree"):
        if line.startswith("Log-likelihood of the tree"):
            lnl = float(line.split(":")[1].split()[0])
    done = dict(step="species_tree_bl", set=a.set, n_taxa=len(taxa), n_markers=len(markers), n_sites=nsites,
                model="LG+G4", fixed_topology=True, total_tree_length=round(sum(lengths.values()), 4),
                logL=lnl, output=os.path.relpath(dst, ROOT), threads=a.threads, seed=a.seed,
                seconds=round(time.time() - t0, 1), finished=dt.datetime.now().isoformat(timespec="seconds"))
    with open(os.path.join(base, "species_tree_bl.DONE.json"), "w") as fh:
        json.dump(done, fh, indent=1)
    print(json.dumps(done, indent=1))


if __name__ == "__main__":
    main()
