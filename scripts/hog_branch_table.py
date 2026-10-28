#!/usr/bin/env python
"""Per-branch gene gain / loss / duplication tables from Ancestra HOGs.

Usage::

    python scripts/hog_branch_table.py RUN/hogs.orthoxml[.gz] --tree RUN/stages/input/guide_tree.nwk \\
        --out OUT [--members-levels N5,N7|all] [--copy-number-levels N5|all] [--format parquet|tsv]

Outputs in OUT/ (definitions: ``mcell/hogs.py`` module docstring):

* ``branch_summary.tsv``  per branch: n_parent, n_child, gains, losses, dups, dups_at_origin,
  leaf_specific_genes, gains_ge2_child_clades, conservation_ok (root stem has branch_parent "").
* ``hogs.parquet``        ancestral genes at internal nodes (hog_id, level, family_id,
  parent_hog_id, origin, n_genes, n_species, n_child_clades_present, rep_genes).
* ``genes.parquet``       extant gene -> species, family, ancestral gene at the parent node.
* ``events.parquet``      one row per (branch, HOG, event) for enrichment analyses.
* ``members.parquet``     (level, hog_id, gene, species) for the requested levels.
* ``copy_number.parquet`` (level, hog_id, species, n_copies) for the requested levels.
* ``run_info.json``       inputs, options, parser statistics, wall time, peak RSS.

The guide tree should be the one Ancestra used (``stages/input/guide_tree.nwk``, or the
original Newick: it is pruned/collapsed the same way). Without ``--tree`` the OrthoXML
``<taxonomy>`` is used; with both, they must agree.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd  # noqa: E402

from mcell import __version__  # noqa: E402
from mcell.hogs import parse_hogs  # noqa: E402


def _levels(arg: str | None, res) -> list[str]:
    if not arg:
        return []
    internal = [res.tree.names[i] for i in res.tree.internal]
    if arg == "all":
        return internal
    out = [x.strip() for x in arg.split(",") if x.strip()]
    bad = [x for x in out if x not in internal]
    if bad:
        raise SystemExit(f"unknown or leaf level(s): {bad}; internal levels: {internal[:20]}...")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("orthoxml", help="hogs.orthoxml or hogs.orthoxml.gz")
    ap.add_argument("--tree", help="rooted Newick guide tree used for the Ancestra run")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--members-levels", default=None,
                    help="comma-separated internal levels (or 'all') for members.parquet")
    ap.add_argument("--copy-number-levels", default=None,
                    help="comma-separated internal levels (or 'all') for copy_number.parquet")
    ap.add_argument("--n-rep", type=int, default=5, help="representative gene ids per HOG/event row")
    ap.add_argument("--format", choices=("parquet", "tsv"), default="parquet")
    a = ap.parse_args(argv)

    t0 = time.time()
    log = lambda m: print(f"[hog_branch_table {time.time() - t0:7.1f}s] {m}", file=sys.stderr, flush=True)
    res = parse_hogs(a.orthoxml, newick=a.tree, n_rep=a.n_rep,
                     progress=lambda n: log(f"{n:,} families"))
    log(f"parsed: {res.stats['n_genes']:,} genes, {res.stats['n_families']:,} families "
        f"({res.stats['n_singleton_genes']:,} ungrouped genes), {len(res.hogs):,} ancestral genes")
    bad = res.check_conservation()
    if not bad.empty:
        log(f"WARNING: conservation law violated on {len(bad)} branches (bug?)")
    paths = res.write(a.out, fmt=a.format)

    def dump(frames, name):
        if not frames:
            return
        df = pd.concat(frames, ignore_index=True)
        if a.format == "parquet":
            p = os.path.join(a.out, f"{name}.parquet")
            df.to_parquet(p, index=False)
        else:
            p = os.path.join(a.out, f"{name}.tsv.gz")
            df.to_csv(p, sep="\t", index=False)
        paths[name] = p

    dump([res.members(lv) for lv in _levels(a.members_levels, res)], "members")
    dump([res.copy_number(lv) for lv in _levels(a.copy_number_levels, res)], "copy_number")
    info = dict(mcell_version=__version__, argv=sys.argv if argv is None else argv,
                orthoxml=os.path.abspath(a.orthoxml), tree=os.path.abspath(a.tree) if a.tree else None,
                orthoxml_bytes=os.path.getsize(a.orthoxml), stats=res.stats,
                outputs={k: os.path.basename(v) for k, v in paths.items()},
                wall_s=round(time.time() - t0, 2),
                peak_rss_gb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 3))
    with open(os.path.join(a.out, "run_info.json"), "w") as fh:
        json.dump(info, fh, indent=1, default=str)
    log(f"done -> {a.out}")
    return 0 if bad.empty else 2


if __name__ == "__main__":
    sys.exit(main())
