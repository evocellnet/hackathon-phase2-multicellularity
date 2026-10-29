# Independent verifier: build species x Pfam protein-count matrix, proteome sizes, disorder fractions (all 187 euk).
import gzip, csv, os, sys
from multiprocessing import Pool
import pandas as pd
P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = f"{P}/results/r2/verification/w4"
tr = pd.read_csv(f"{P}/data/species_traits.tsv", sep="\t")
sp = tr.loc[tr.kingdom == "euk", "code"].tolist()

def one(code):
    fams = {}
    with gzip.open(f"{P}/data/annotations/pfam/euk/{code}.pfam.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            f = line.split("\t", 3)
            fams.setdefault(f[1].split(".")[0], set()).add(f[0])
    counts = {k: len(v) for k, v in fams.items()}
    n = 0
    with open(f"{P}/data/proteomes/euk/{code}.fa") as fh:
        for line in fh:
            if line.startswith(">"): n += 1
    L = D = 0.0; nd = 0
    with gzip.open(f"{P}/data/annotations/disorder/euk/{code}.disorder.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            l = float(f[1]); L += l; D += l * float(f[2]); nd += 1
    return code, counts, n, D / L, nd

with Pool(8) as p:
    res = p.map(one, sp)
rows = []
cnt = {}
for code, counts, n, dis, nd in res:
    rows.append((code, n, dis, nd)); cnt[code] = counts
pd.DataFrame(rows, columns=["code", "n_proteins", "disorder", "n_disorder_rows"]).to_csv(f"{OUT}/species_cov.tsv", sep="\t", index=False)
M = pd.DataFrame.from_dict(cnt, orient="index").fillna(0).astype(int)
M.index.name = "code"
M.to_csv(f"{OUT}/count_matrix.tsv.gz", sep="\t", compression="gzip")
print(M.shape)
