#!/usr/bin/env python
"""Freeze or verify the pre-registration (DESIGN.md §7).

  check_prereg.py --freeze   write PREREG_HASHES.tsv (sha256 of every frozen input listed in FROZEN below)
  check_prereg.py            verify that every frozen input still matches its hash; exit 1 on any mismatch

Analysis scripts that produce primary results call `verify()` first, so a primary result can never be computed from
silently changed inputs. Changes after the freeze must be recorded in PREREG.md "Deviations" and re-frozen with a
new tag (prereg-v2, ...), never by editing the hash file by hand.
"""
from __future__ import annotations

import hashlib
import os
import sys

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
FROZEN = [
    "DESIGN.md",
    "data/transitions.tsv",
    "data/guide_tree_euk.nwk",
    "data/guide_tree_prok.nwk",
    "data/positive_controls.tsv",
    "data/positive_controls_nodes.tsv",
    "data/prereg_definitions.tsv",
    "data/toolkit_pfam.tsv",
    "data/h12_sets.tsv",
    "data/h12_expected_direction_frozen.tsv",
    "data/literature_scorecard_template.tsv",
    "data/species_final.tsv",
    "data/species_final_prok.tsv",
    "results/prereg/power_table.tsv",
    "scripts/prereg/apply_eligibility.py",
]
HASHFILE = os.path.join(P, "PREREG_HASHES.tsv")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def freeze() -> None:
    with open(HASHFILE, "w") as out:
        out.write("path\tsha256\n")
        for rel in FROZEN:
            out.write(f"{rel}\t{sha256(os.path.join(P, rel))}\n")
    print(f"wrote {HASHFILE} ({len(FROZEN)} files)")


def verify(strict: bool = True) -> bool:
    ok = True
    for line in open(HASHFILE).read().splitlines()[1:]:
        rel, h = line.split("\t")
        cur = sha256(os.path.join(P, rel)) if os.path.exists(os.path.join(P, rel)) else "MISSING"
        if cur != h:
            ok = False
            print(f"PREREG MISMATCH: {rel}", file=sys.stderr)
    if not ok and strict:
        raise SystemExit("pre-registered inputs changed; record a deviation in PREREG.md and re-freeze")
    return ok


if __name__ == "__main__":
    if "--freeze" in sys.argv:
        freeze()
    else:
        print("PREREG OK" if verify(strict=False) else "PREREG CHANGED")
        sys.exit(0 if verify(strict=False) else 1)
