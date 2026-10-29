"""Shared helpers for the round-2 dense-panel analysis (W2/W3).

Paths, refuse-to-overwrite guards, panel metadata, and the placement of round-1
proteomes that are not in a panel onto the panel tree (for the novelty clade filter).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
sys.path.insert(0, ROOT)
os.environ["PATH"] = os.path.join(ROOT, "conda_envs", "ancestra", "bin") + os.pathsep + os.environ["PATH"]

import ete3  # noqa: E402
import pandas as pd  # noqa: E402

PANELS = ["R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"]
SEEDS = [11, 23, 37, 41, 53]
OUT = os.path.join(ROOT, "results", "r2", "panels")


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


EXTRA = ["R4_curated"]   # sensitivity run: R4 without the two AUGUSTUS proteomes (Eudosp, Yamauni)


def base(R: str) -> str:
    """panel whose species/transitions tables a run uses"""
    return R.split("_")[0]


def kingdom(R: str) -> str:
    return "prok" if R == "R8" else "euk"


def paths(R: str) -> dict:
    rd = os.path.join(ROOT, "runs", "r2", R)
    out = os.path.join(OUT, R)
    work = os.path.join(out, "work")
    return dict(
        R=R, run_dir=rd, done=os.path.join(rd, "DONE.json"), orthoxml=os.path.join(rd, "hogs.orthoxml"),
        tree=os.path.join(rd, "stages", "input", "guide_tree.nwk"),
        guide=os.path.join(ROOT, "data", "dense", R, "guide_tree.nwk"),
        species=os.path.join(ROOT, "data", "dense", base(R), "species.tsv"),
        transitions=os.path.join(ROOT, "data", "dense", base(R), "transitions.tsv"),
        proteomes=os.path.join(ROOT, "data", "dense", R, "proteomes"),
        r1_proteomes=os.path.join(ROOT, "data", "proteomes", kingdom(R)),
        r1_tree=os.path.join(ROOT, "data", f"guide_tree_{kingdom(R)}.nwk"),
        out=out, work=work, cache=os.path.join(work, "cache"), phylo=os.path.join(work, "phylo"),
        ladder=os.path.join(work, "ladder"), pfam=os.path.join(work, "pfam"),
    )


def guard(path: str) -> None:
    """Refuse to overwrite a final output."""
    if os.path.exists(path):
        raise SystemExit(f"refusing to overwrite existing output {path}")


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_done(dirpath: str, files: list, extra: dict | None = None, name: str = "DONE.json") -> None:
    d = dict(status="done", finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
             md5={os.path.relpath(f, dirpath): md5(f) for f in files if os.path.isfile(f)})
    if extra:
        d.update(extra)
    p = os.path.join(dirpath, name)
    guard(p)
    with open(p, "w") as fh:
        json.dump(d, fh, indent=1, default=str)


def panel_codes(R: str) -> list:
    p = paths(R)
    return sorted(f[:-3] for f in os.listdir(p["proteomes"]) if f.endswith(".fa"))


def r1_codes(R: str) -> list:
    p = paths(R)
    return sorted(f[:-3] for f in os.listdir(p["r1_proteomes"]) if f.endswith(".fa"))


def place_round1(R: str, panel_tree: ete3.Tree | None = None) -> pd.DataFrame:
    """For every round-1 proteome NOT in the panel: the smallest panel node whose clade
    contains it ('' = outside every panel clade, i.e. outside the panel root).

    Rule (fixed before any novelty output): a round-1 species s lies inside panel node c iff
    c has >= 2 panel tips that are round-1 species (set O_c) and, in the round-1 guide tree,
    s is below MRCA(O_c) and that MRCA contains no other panel/round-1 overlap species (a
    consistent clade). Nodes with < 2 overlap tips are never used, so s is placed at the
    deepest consistent node; this can only place s *higher* (more inclusive) than its true
    position when the panel node is unresolved by round-1 species, which makes s count as
    in-clade for the lower node -> fewer outside targets for that node (novelty
    counts can be slightly inflated there; reported per placement)."""
    p = paths(R)
    pt = panel_tree or ete3.Tree(p["guide"], format=1)
    rt = ete3.Tree(p["r1_tree"], format=1)
    r1 = set(rt.get_leaf_names())
    pan = set(pt.get_leaf_names())
    ov = pan & r1
    extra = sorted(r1 - pan)
    r1_leaf = {n.name: n for n in rt.get_leaves()}
    rows = []
    # candidate panel nodes (internal, >=2 overlap tips) and their round-1 MRCA
    cand = []
    for n in pt.traverse("preorder"):
        if n.is_leaf():
            continue
        O = set(n.get_leaf_names()) & ov
        if len(O) < 2:
            continue
        m = rt.get_common_ancestor([r1_leaf[x] for x in O])
        under = set(m.get_leaf_names())
        consistent = (under & ov) == O
        cand.append((n, m, under, consistent))
    for s in extra:
        best, depth = "", -1
        for n, m, under, ok in cand:
            if ok and s in under:
                d = len(n.get_ancestors())
                if d > depth:
                    best, depth = n.name, d
        rows.append(dict(species=s, panel_node=best, depth=depth))
    return pd.DataFrame(rows)
