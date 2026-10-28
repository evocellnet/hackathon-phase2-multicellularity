#!/usr/bin/env python
"""Diagnostic: connected-component sizes of Ancestra's S0 stage-2 graph, recomputed from its own files.

Replicates ancestra.families default config (pipeline.DEFAULT_FAMILY_CONFIG): gene-vs-rep DIAMOND hits filtered
at e <= 1e-5 and >= 30% coverage of the shorter sequence, lifted to stage-1 cluster representatives, connected
components; reports how many components exceed s2_gmax = 7,500 genes and the size of the largest. Read-only.

Usage: component_sizes.py <families.tmp*/s0 dir> [cov_mode=shorter|min] [cov=0.3] [evalue=1e-5]
"""
import sys

import numpy as np
import pandas as pd
import pyarrow.csv as pacsv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

d = sys.argv[1]
MODE = sys.argv[2] if len(sys.argv) > 2 else "shorter"
COV = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
EV = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-5
print(f"rule: cov_mode={MODE} cov={COV} evalue={EV}")
cl = pd.read_csv(f"{d}/stage1_clusters.tsv", sep="\t", header=None, names=["rep", "member"], dtype=str)
reps = pd.Index(sorted(set(cl["rep"])))
gene2rep = pd.Series(cl["rep"].values, index=cl["member"].values)
rep_size = cl.groupby("rep").size()
print(f"stage-1 clusters: {len(reps):,}; genes: {len(cl):,}", flush=True)

cols = ["q", "s", "ev", "bs", "ql", "sl", "qs", "qe", "ss", "se"]
tab = pacsv.read_csv(f"{d}/reps_hits.tsv", read_options=pacsv.ReadOptions(column_names=cols, block_size=1 << 26),
                     parse_options=pacsv.ParseOptions(delimiter="\t"),
                     convert_options=pacsv.ConvertOptions(column_types={"q": "string", "s": "string"}))
df = tab.to_pandas()
print(f"hits: {len(df):,}", flush=True)
df = df[(df.q != df.s) & (df.ev <= EV)]
qcov = (np.abs(df.qe - df.qs) + 1) / df.ql
scov = (np.abs(df.se - df.ss) + 1) / df.sl
cov_use = np.where(df.ql <= df.sl, qcov, scov) if MODE == 'shorter' else np.minimum(qcov, scov)
df = df[cov_use >= COV]
ra = gene2rep.reindex(df.q.values).values
rb = gene2rep.reindex(df.s.values).values
ia = reps.get_indexer(ra)
ib = reps.get_indexer(rb)
m = (ia >= 0) & (ib >= 0) & (ia != ib)
n = len(reps)
g = coo_matrix((np.ones(m.sum(), dtype=np.int8), (ia[m], ib[m])), shape=(n, n))
nc, lab = connected_components(g, directed=False)
w = rep_size.reindex(reps).fillna(1).values
comp_genes = np.bincount(lab, weights=w)
comp_reps = np.bincount(lab)
order = np.argsort(-comp_genes)
big = comp_genes > 7500
print(f"components: {nc:,}; > 7500 genes: {big.sum()}; genes in them: {comp_genes[big].sum():,.0f} "
      f"({comp_genes[big].sum() / w.sum():.1%})")
for i in order[:10]:
    print(f"  component: {comp_genes[i]:,.0f} genes, {comp_reps[i]:,} reps")
