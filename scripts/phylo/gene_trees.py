#!/usr/bin/env python
"""Step 3: gene-wise branch lengths on the fixed (pruned) guide topology.

DESIGN.md §4.6 H9 / reviews/round2/statistics.md M6: per marker, IQ-TREE
``-te`` (LG+G4, 1 thread, seed 11) on the guide tree pruned to the taxa present.
Fitted once; the H9 stems, the H9 controls and the RER analysis all reuse these
fits.

Site coverage (M6 inclusion "≥ 50 % site coverage") is applied per taxon: a
taxon whose sequence covers < ``--min-cov`` (0.5) of the trimmed alignment
columns is dropped from that gene before fitting (fragments give unreliable
terminal lengths and would otherwise enter L̄). Genes with < 4 remaining taxa
are skipped.

Outputs in runs/phylo/<set>/genetrees/:
  <m>.fa / <m>.te.nwk     filtered alignment and pruned topology given to IQ-TREE
  <m>.bl.nwk              rooted pruned guide tree with guide node names + lengths
  branch_lengths.tsv      gene, node, parent, length, n_tips, root_child
  genes.tsv               per gene: n_taxa, dropped (coverage), n_sites, tree length,
                          max branch, alpha, logL, seconds
and runs/phylo/<set>/gene_trees.DONE.json.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import json
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
os.environ["PATH"] = os.path.join(ROOT, "conda_envs", "ancestra", "bin") + os.pathsep + os.environ["PATH"]

from mcell.phylo import read_tree, prune_named, transfer_branch_lengths, write_named  # noqa: E402
from scripts.phylo.build_markers import read_fasta, write_fasta  # noqa: E402


def fit_one(m, mdir, gdir, guide_nwk, min_cov, seed):
    t0 = time.time()
    guide = read_tree(guide_nwk)
    al = read_fasta(os.path.join(mdir, f"{m}.trim.fa"))
    if not al:
        return dict(gene=m, status="empty_alignment")
    L = len(next(iter(al.values())))
    cov = {k: sum(ch not in "-X?" for ch in v) / L for k, v in al.items()}
    keep = {k: v for k, v in al.items() if cov[k] >= min_cov}
    dropped = sorted(set(al) - set(keep))
    if len(keep) < 4:
        return dict(gene=m, status="lt4_taxa", n_taxa=len(keep), dropped=",".join(dropped), n_sites=L)
    pruned = prune_named(guide, keep)
    fa, te, prefix = (os.path.join(gdir, f"{m}.fa"), os.path.join(gdir, f"{m}.te.nwk"), os.path.join(gdir, m))
    bl_out = os.path.join(gdir, f"{m}.bl.nwk")
    if not os.path.exists(bl_out):
        write_fasta(fa, keep)
        with open(te, "w") as fh:
            fh.write(pruned.write(format=9) + "\n")
        r = subprocess.run(["iqtree3", "-s", fa, "-te", te, "-m", "LG+G4", "-T", "1", "-seed", str(seed),
                            "--prefix", prefix, "-redo", "-quiet"], capture_output=True, text=True)
        if r.returncode != 0:
            return dict(gene=m, status="iqtree_failed", n_taxa=len(keep), err=r.stderr[-300:].replace("\n", " "))
        lengths = transfer_branch_lengths(pruned, read_tree(prefix + ".treefile"))
        with open(bl_out + ".tmp", "w") as fh:
            fh.write(write_named(pruned, lengths) + "\n")
        os.replace(bl_out + ".tmp", bl_out)
        for ext in (".ckp.gz", ".mldist", ".model.gz", ".uniqueseq.phy", ".bionj"):
            if os.path.exists(prefix + ext):
                os.remove(prefix + ext)
    t = read_tree(bl_out)
    alpha = lnl = None
    if os.path.exists(prefix + ".iqtree"):
        txt = open(prefix + ".iqtree").read()
        mm = re.search(r"Gamma shape alpha:\s*([\d.eE+-]+)", txt)
        alpha = float(mm.group(1)) if mm else None
        mm = re.search(r"Log-likelihood of the tree:\s*([\d.eE+-]+)", txt)
        lnl = float(mm.group(1)) if mm else None
    rows = [(n.name, n.up.name, n.dist, len(n), int(n.up.up is None)) for n in t.traverse() if n.up is not None]
    return dict(gene=m, status="ok", n_taxa=len(keep), dropped=",".join(dropped), n_sites=L,
                tree_length=round(sum(r[2] for r in rows), 5), max_branch=round(max(r[2] for r in rows), 5),
                max_branch_node=max(rows, key=lambda r: r[2])[0], alpha=alpha, logL=lnl,
                seconds=round(time.time() - t0, 1), rows=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["euk", "prok"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    t0 = time.time()
    base = os.path.join(ROOT, "runs", "phylo", a.set)
    mdir, gdir = os.path.join(base, "markers"), os.path.join(base, "genetrees")
    os.makedirs(gdir, exist_ok=True)
    guide_nwk = os.path.join(ROOT, "data", f"guide_tree_{a.set}.nwk")
    with open(os.path.join(mdir, "markers.tsv")) as fh:
        markers = [r["marker"] for r in csv.DictReader(fh, delimiter="\t")]
    res = []
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(fit_one, m, mdir, gdir, guide_nwk, a.min_cov, a.seed) for m in markers]
        for i, f in enumerate(cf.as_completed(futs), 1):
            res.append(f.result())
            if i % 25 == 0:
                print(f"  {i}/{len(futs)} genes ({time.time() - t0:.0f} s)", flush=True)
    res.sort(key=lambda r: r["gene"])
    with open(os.path.join(gdir, "branch_lengths.tsv"), "w") as fh:
        fh.write("gene\tnode\tparent\tlength\tn_tips\troot_child\n")
        for r in res:
            for node, par, d, nt, rc in r.get("rows", []):
                fh.write(f"{r['gene']}\t{node}\t{par}\t{d:.8g}\t{nt}\t{rc}\n")
    cols = ["gene", "status", "n_taxa", "dropped", "n_sites", "tree_length", "max_branch", "max_branch_node",
            "alpha", "logL", "seconds"]
    with open(os.path.join(gdir, "genes.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in res:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in cols) + "\n")
    ok = [r for r in res if r["status"] == "ok"]
    done = dict(step="gene_trees", set=a.set, n_markers=len(markers), n_ok=len(ok),
                n_failed=len(res) - len(ok), min_cov=a.min_cov, model="LG+G4", jobs=a.jobs, seed=a.seed,
                cpu_seconds=round(sum(r.get("seconds", 0) for r in ok), 1),
                seconds=round(time.time() - t0, 1), finished=dt.datetime.now().isoformat(timespec="seconds"))
    with open(os.path.join(base, "gene_trees.DONE.json"), "w") as fh:
        json.dump(done, fh, indent=1)
    print(json.dumps(done, indent=1))


if __name__ == "__main__":
    main()
