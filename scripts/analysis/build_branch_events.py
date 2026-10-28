#!/usr/bin/env python
"""Parse an Ancestra run once and cache the analysis tables (mcell.analysis.build_cache).

    python scripts/analysis/build_branch_events.py --run {euk,prok} [--levels auto|all|N1,N2]

Writes runs/analysis/<run>/:
  branch_summary.tsv        per branch: n_parent, n_child, gains, losses, dups, dups_at_origin, ... (mcell.hogs)
  events.parquet            one row per (branch, HOG, event)
  hogs.parquet, genes.parquet  ancestral genes (all internal levels, with parent links) and extant genes
  members.parquet (level, hog_row, gene_row: rows of hogs/genes.parquet), copy_number.parquet (level, hog_row, hog_id,
                            species, n_copies) at the chosen levels
                            (auto = every P, M and path node of this kingdom's transitions + positive-control nodes)
  pfam_hits.parquet         pre-registered Pfam subset hits (data/annotations/pfam_prereg/<run>)
  gene_classes.parquet, hog_classes.parquet  class annotation (>=50% of members carry a class domain)
  branch_class_events.tsv   per branch: HOGs gained / duplicated / lost, total and per class
  run_info.json, DONE.json
Runs single-threaded; use `nice -n 19`.
"""
import argparse
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from mcell import analysis as an  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, choices=["euk", "prok"])
    ap.add_argument("--levels", default=None, help="auto (default euk), all (default prok), or comma list")
    ap.add_argument("--force", action="store_true", help="rebuild even if DONE.json exists")
    a = ap.parse_args(argv)
    p = an.run_paths(a.run)
    if os.path.exists(os.path.join(p["cache"], "DONE.json")) and not a.force:
        print(f"cache exists: {p['cache']} (use --force to rebuild)")
        return 0
    if not os.path.exists(os.path.join(p["run_dir"], "DONE.json")):
        raise SystemExit(f"Ancestra run not finished: {p['run_dir']}/DONE.json missing")
    t0 = time.time()
    log = lambda m: print(f"[build_branch_events {time.time() - t0:7.1f}s] {m}", file=sys.stderr, flush=True)  # noqa: E731
    files = an.build_cache(a.run, p, levels=a.levels or ("all" if a.run == "prok" else "auto"), log=log)
    for k, v in files.items():
        log(f"{k}: {os.path.relpath(v, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
