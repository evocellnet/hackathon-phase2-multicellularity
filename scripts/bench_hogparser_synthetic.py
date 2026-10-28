#!/usr/bin/env python
"""Write a large synthetic OrthoXML (no sequences) to benchmark ``mcell.hogs`` at scale.

Families are generated top-down on a random Yule species tree: a gene at node X
is lost on each child branch with probability ``--loss`` and duplicated with
probability ``--dup`` (both copies then descend independently). Families start
at a random node (root with probability 0.6). The file is written with
Ancestra's OrthoXML writer, so it has Ancestra's exact layout.

    python scripts/bench_hogparser_synthetic.py --species 250 --families 15000 --out bench.orthoxml.gz
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import time

import numpy as np

from ancestra.orthoxml import GeneEntry, Group, OrthoXML, write_orthoxml
from ancestra.simulate import random_species_tree


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", type=int, default=250)
    ap.add_argument("--families", type=int, default=15000)
    ap.add_argument("--loss", type=float, default=0.08)
    ap.add_argument("--dup", type=float, default=0.03)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rng = np.random.default_rng(a.seed)
    t0 = time.time()
    tree = random_species_tree(a.species, rng)
    nodes = list(tree.preorder())
    genes: dict[str, GeneEntry] = {}
    counter = [0]

    def gene(sp):
        counter[0] += 1
        nm = f"{sp}_p{counter[0]}"
        genes[nm] = GeneEntry(nm, sp)
        return nm

    def grow(nd):
        """Lineage arriving at species node nd (after the branch above it)."""
        if not nd.children:
            return gene(nd.name)
        kids = []
        for c in nd.children:
            if rng.random() < a.loss:
                continue
            if rng.random() < a.dup:
                copies = [grow(c), grow(c)]
                kids.append(Group("paralog", copies, level=c.name))
            else:
                kids.append(grow(c))
        if not kids:
            return None
        kids = [k for k in kids if k is not None]
        if len(kids) == 1:
            return kids[0]
        return Group("ortholog", kids, level=nd.name)

    groups = []
    for _ in range(a.families):
        start = tree if rng.random() < 0.6 else nodes[int(rng.integers(len(nodes)))]
        g = grow(start)
        if g is not None:
            groups.append(g)
    ox = OrthoXML(genes=genes, groups=groups, species_tree=tree)
    print(f"generated {len(genes):,} genes, {len(groups):,} families in {time.time() - t0:.0f}s",
          file=sys.stderr)
    opener = gzip.open if a.out.endswith(".gz") else open
    with opener(a.out, "wt") as fh:
        write_orthoxml(ox, fh)
    with open(os.path.splitext(a.out)[0] + ".nwk" if not a.out.endswith(".gz")
              else a.out[:-3] + ".nwk", "w") as fh:
        fh.write(tree.to_newick() + "\n")
    print(f"written {a.out} in {time.time() - t0:.0f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
