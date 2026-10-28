#!/usr/bin/env python3
"""EXPLORATORY (DESIGN §5). Sample sequence-novel gains with a UniProt-sourced representative for the 7 primary stems
and an equal-sized matched control sample.  seed 11.  Output results/structure/query_set.tsv"""
import os, sys, numpy as np, pandas as pd
P = os.environ.get('P', os.getcwd()); SEED = 11; NMAX = 60
tr = pd.read_csv(f'{P}/data/transitions.tsv', sep='\t'); tr = tr[tr.analysis_set == 'primary']
sp = pd.read_csv(f'{P}/data/species_final.tsv', sep='\t')
uni = set(sp.loc[sp.source_db == 'UniProt', 'code'])
g = pd.read_csv(f'{P}/results/gainladder/euk/gains_classified.tsv.gz', sep='\t',
                usecols=['family_id', 'branch_parent', 'branch_child', 'ladder', 'seq_novel', 'n_reps'])
g = g[g.ladder == 'seq_novel']
r = pd.read_parquet(f'{P}/runs/gainladder/euk/reps.parquet')
r = r[r.species.isin(uni)].sort_values(['family_id', 'rep_rank'])
first = r.groupby('family_id').head(1)[['family_id', 'gene', 'species', 'rep_rank']]
first['accession'] = first.gene.str.split('_', n=1).str[1]
g = g.merge(first, on='family_id', how='left')
g['has_uniprot_rep'] = g.gene.notna()
ms = pd.read_csv(f'{P}/results/primary/euk/matched_sets.tsv', sep='\t')
ms = ms[(ms.hypothesis == 'H4') & (ms.set == 'primary')]
rng = np.random.default_rng(SEED)
out = []; cnt = []
for _, t in tr.iterrows():
    nodes = t.stem_path.split('>')
    stem_children = nodes[1:]            # path-stem edges count as the stem
    ctrl = list(ms[ms.transition_id == t.transition_id].control)
    gs = g[g.branch_child.isin(stem_children)]
    gc = g[g.branch_child.isin(ctrl)]
    es, ec = gs[gs.has_uniprot_rep], gc[gc.has_uniprot_rep]
    k = min(NMAX, len(es), len(ec))
    ss = es.sample(n=k, random_state=int(rng.integers(1 << 31))) if k else es.iloc[:0]
    sc = ec.sample(n=k, random_state=int(rng.integers(1 << 31))) if k else ec.iloc[:0]
    cnt.append(dict(transition_id=t.transition_id, stem_edges='|'.join(stem_children), n_seqnovel_stem=len(gs),
                    n_stem_uniprot_rep=len(es), n_control_branches=len(ctrl), n_seqnovel_control=len(gc),
                    n_control_uniprot_rep=len(ec), n_sampled_each=k))
    for lab, d in (('stem', ss), ('control', sc)):
        d = d.copy(); d['transition_id'] = t.transition_id; d['set'] = lab
        out.append(d)
q = pd.concat(out)
# exclude_taxids: for CONTROL gains the NCBI LCA taxon of the gaining clade's species (hits inside the gain's own clade are not
# evidence of an older fold); for stem gains blank -> transition-level M definition (data/structure_ref_exclusion.tsv) is used.
from ete3 import Tree
tree = Tree(open(f'{P}/data/guide_tree_euk.nwk').read(), format=1)
nodes = {n.name: n for n in tree.traverse() if n.name}
par = {}
for l in open(f'{P}/data/raw/taxdump/nodes.dmp'):
    f = l.split('\t|\t'); par[int(f[0])] = int(f[1])
def lin(x):
    o = []
    while x not in o: o.append(x); x = par[x]
    return o
tx = sp.set_index('code').ncbi_taxid
def lca_taxon(clade):
    ls = [lin(int(tx[l.name])) for l in nodes[clade].iter_leaves() if l.name in tx.index]
    c = set(ls[0]).intersection(*map(set, ls[1:])) if len(ls) > 1 else set(ls[0])
    return next(x for x in ls[0] if x in c)
cl = {c: lca_taxon(c) for c in set(q.branch_child)}
q['exclude_taxids'] = [str(cl[c]) if s_ == 'control' else '' for c, s_ in zip(q.branch_child, q.set)]
q = q[['accession', 'transition_id', 'set', 'family_id', 'branch_parent', 'branch_child', 'species', 'rep_rank', 'exclude_taxids']]
q.to_csv(f'{P}/results/structure/query_set.tsv', sep='\t', index=False)
pd.DataFrame(cnt).to_csv(f'{P}/results/structure/sampling_counts.tsv', sep='\t', index=False)
print(pd.DataFrame(cnt).to_string()); print(q.groupby(['transition_id', 'set']).size().unstack())
