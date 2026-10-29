#!/usr/bin/env python
"""Errata item 2: re-run the pre-registered power table (scripts/prereg_power_table.py)
with the real phylogenomic branch lengths (data/species_tree_bl_euk.nwk, the IQ-TREE
fixed-topology BUSCO-marker branch lengths used by the analyses) instead of the
placeholder 1.0 that was in force at the PREREG freeze (the bl tree did not exist yet).

Identical simulation code, specs, effect grid, n_sim=200, n_draws=1000 and seed scheme
(SEED=11; per-replicate seed = 11*1_000_003 + spec*100_003 + int(effect*1000)*1009 + r).
Only the branch-length source differs. --bl placeholder re-runs the frozen setting with
the current code (reproducibility check against results/prereg/power_table.tsv).

Usage: python scripts/errata/power_table_real_bl.py --bl real|placeholder [--jobs 12]
Output: results/round1_errata/power_<bl>/power_table.tsv, power_curves.tsv
"""
import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import prereg_power_table as ppt  # noqa: E402

MODE = {"bl": "real"}
_orig_exists = os.path.exists


def setup():
    blp = os.path.join(ROOT, "data/species_tree_bl_euk.nwk")
    if MODE["bl"] == "placeholder":
        # make the frozen script believe the bl tree is absent -> placeholder 1.0
        os.path.exists = lambda p: False if os.path.abspath(str(p)) == blp else _orig_exists(p)
    try:
        return ppt._setup_orig()
    finally:
        os.path.exists = _orig_exists


ppt._setup_orig = ppt.setup
ppt.setup = setup

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bl", choices=["real", "placeholder"], default="real")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--n-sim", type=int, default=200)
    a = ap.parse_args()
    MODE["bl"] = a.bl
    out = os.path.join(ROOT, "results/round1_errata", f"power_{a.bl}")
    sys.argv = [sys.argv[0], "--n-sim", str(a.n_sim), "--n-draws", "1000", "--jobs", str(a.jobs), "--out", out]
    ppt.main()
