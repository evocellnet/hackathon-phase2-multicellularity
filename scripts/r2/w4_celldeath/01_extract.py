#!/usr/bin/env python
"""W4 cell-death follow-up, step 1: per-protein ordered Pfam architectures of every protein that carries at least
one of the 22 cell-death families (data/r2_sets/cell_death.tsv), for the 187 round-1 eukaryotic proteomes.

Input : data/annotations/pfam/euk/<code>.pfam.tsv.gz (pyhmmer, Pfam-A 38.2, gathering cutoffs; same tables as W4).
Output: results/r2/w4_celldeath/celldeath_proteins.tsv.gz
        columns: protein, species, cd_accs (cell-death families on the protein), arch (ordered Pfam names,
        consecutive repeats collapsed as NAME(xN); same rule as scripts/r2/w4/01_build_matrix.py), arch_accs.
Domain order = ali_from; overlapping hits are kept (as in W4 matrix/architectures.tsv.gz).
Run: nice -n 10 conda_envs/ancestra/bin/python scripts/r2/w4_celldeath/01_extract.py   (8 processes)
"""
import os
import sys
from multiprocessing import Pool

import pandas as pd

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = os.path.join(P, "results/r2/w4_celldeath")
CD = pd.read_csv(os.path.join(P, "data/r2_sets/cell_death.tsv"), sep="\t")
CDACC = set(CD.pfam)


def collapse(names):
    res, prev, n = [], None, 0
    for d in names:
        if d == prev:
            n += 1
        else:
            if prev is not None:
                res.append(prev if n == 1 else f"{prev}(x{n})")
            prev, n = d, 1
    if prev is not None:
        res.append(prev if n == 1 else f"{prev}(x{n})")
    return " | ".join(res)   # ' | ' separator: several Pfam names contain '-'


def one(code):
    pf = pd.read_csv(f"{P}/data/annotations/pfam/euk/{code}.pfam.tsv.gz", sep="\t",
                     usecols=["protein", "pfam_acc", "pfam_name", "ali_from"])
    keep = set(pf.loc[pf.pfam_acc.isin(CDACC), "protein"])
    pf = pf[pf.protein.isin(keep)].sort_values(["protein", "ali_from"])
    rows = []
    for prot, g in pf.groupby("protein", sort=False):
        accs = list(g.pfam_acc)
        rows.append(dict(protein=prot, species=code,
                         cd_accs=";".join(sorted(set(a for a in accs if a in CDACC))),
                         arch=collapse(list(g.pfam_name)),
                         arch_accs=collapse(accs)))
    return pd.DataFrame(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    dst = os.path.join(OUT, "celldeath_proteins.tsv.gz")
    if os.path.exists(dst):
        sys.exit(f"REFUSE: {dst} exists")
    tr = pd.read_csv(os.path.join(P, "data/species_traits.tsv"), sep="\t")
    codes = list(tr[tr.kingdom == "euk"].code)
    assert len(codes) == 187, len(codes)
    with Pool(8) as pool:
        res = pool.map(one, codes)
    df = pd.concat(res, ignore_index=True)
    df.to_csv(dst, sep="\t", index=False, compression="gzip")
    print("proteins:", len(df), "species with >=1:", df.species.nunique())


if __name__ == "__main__":
    main()
