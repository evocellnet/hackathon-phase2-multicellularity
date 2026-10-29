"""Build per-species per-Pfam-family copy-number matrix (proteins with >=1 GA hit of the family),
species covariates (proteome size, residue-weighted disorder fraction) and per-protein domain
architectures for the 187 round-1 eukaryotic proteomes.

Inputs: data/annotations/pfam/euk/<code>.pfam.tsv.gz (pyhmmer, Pfam 38.2, GA cutoffs),
        data/annotations/disorder/euk/<code>.disorder.tsv.gz (metapredict, threshold 0.5),
        data/proteomes/euk/<code>.fa, data/species_traits.tsv.
Outputs (results/r2/w4/matrix/): counts.tsv.gz, species.tsv, architectures.tsv.gz
"""
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import P, out  # noqa: E402


def collapse(doms):
    """Ordered domain string with consecutive repeats collapsed: SH3,SH3,SH2 -> SH3(x2)-SH2."""
    res, prev, n = [], None, 0
    for d in doms:
        if d == prev:
            n += 1
        else:
            if prev is not None:
                res.append(prev if n == 1 else f"{prev}(x{n})")
            prev, n = d, 1
    if prev is not None:
        res.append(prev if n == 1 else f"{prev}(x{n})")
    return "-".join(res)


def one(code):
    pf = pd.read_csv(f"{P}/data/annotations/pfam/euk/{code}.pfam.tsv.gz", sep="\t",
                     usecols=["protein", "pfam_acc", "pfam_name", "ali_from"])
    cnt = pf.groupby("pfam_acc").protein.nunique()
    arch = (pf.sort_values(["protein", "ali_from"]).groupby("protein").pfam_name
            .agg(lambda s: collapse(list(s))))
    acc_sets = pf.groupby("protein").pfam_acc.agg(lambda s: ";".join(sorted(set(s))))
    a = pd.DataFrame({"protein": arch.index, "species": code, "arch": arch.values,
                      "accs": acc_sets.loc[arch.index].values})
    dis = pd.read_csv(f"{P}/data/annotations/disorder/euk/{code}.disorder.tsv.gz", sep="\t")
    nfa = 0
    with open(f"{P}/data/proteomes/euk/{code}.fa") as fh:
        for line in fh:
            if line.startswith(">"):
                nfa += 1
    dfrac = float((dis.frac_disordered * dis.length).sum() / dis.length.sum())
    cov = dict(code=code, n_proteins=nfa, n_disorder_rows=len(dis), disorder_frac=dfrac,
               n_proteins_with_pfam=int(pf.protein.nunique()), n_domain_hits=len(pf))
    return code, cnt, a, cov


def main():
    tr = pd.read_csv(f"{P}/data/species_traits.tsv", sep="\t")
    tr = tr[tr.kingdom == "euk"]
    codes = list(tr.code)
    assert len(codes) == 187
    f_counts, f_species, f_arch = out("matrix/counts.tsv.gz"), out("matrix/species.tsv"), out("matrix/architectures.tsv.gz")
    with ProcessPoolExecutor(12) as ex:
        res = list(ex.map(one, codes))
    counts = pd.DataFrame({c: cnt for c, cnt, _, _ in res}).fillna(0).astype(int).T
    counts.index.name = "code"
    counts = counts.loc[codes]
    cov = pd.DataFrame([r[3] for r in res])
    sp = tr.merge(cov, on="code")
    bad = sp[sp.n_proteins != sp.n_disorder_rows]
    if len(bad):
        print("WARNING proteome/disorder count mismatch:\n", bad[["code", "n_proteins", "n_disorder_rows"]])
    counts.to_csv(f_counts, sep="\t", compression="gzip")
    sp.to_csv(f_species, sep="\t", index=False)
    pd.concat([r[2] for r in res]).to_csv(f_arch, sep="\t", index=False, compression="gzip")
    print(counts.shape, "families; species cov written")


if __name__ == "__main__":
    main()
