#!/usr/bin/env python3
"""EXPLORATORY (DESIGN §5 / Q6). Stem vs matched-control sequence-novel gains: fraction 'structurally old'
and fraction without a usable AFDB model.  Wilson CIs, per-stem Fisher, Mantel-Haenszel/CMH stratified by stem.
Usage: analyze_structure.py <queries.tsv from novel_gain_structures.py> <label> [--append-registry]"""
import os, sys, math, argparse, numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.contingency_tables import StratifiedTable
P = os.environ.get('P', os.getcwd()); sys.path.insert(0, P)
ap = argparse.ArgumentParser(); ap.add_argument('queries'); ap.add_argument('label')
ap.add_argument('--only', nargs='*', help='restrict to these transition ids (e.g. covered by the proteome DB)')
ap.add_argument('--qset', default=f'{P}/results/structure/query_set.tsv', help='sampling table (clade names, exclude_taxids)')
ap.add_argument('--exclude-ref-inside', action='store_true',
                help='proteome mode: drop queries whose gaining clade contains a reference species of its transition (hit would be a descendant, not an outside fold)')
ap.add_argument('--drop-control-euk-lca', action='store_true',
                help='swissprot sensitivity: drop controls whose own-clade LCA is Eukaryota (all eukaryote hits excluded -> only prokaryote hits can count)')
ap.add_argument('--append-registry', action='store_true'); ap.add_argument('--out', default=f'{P}/results/structure')
a = ap.parse_args()
q = pd.read_csv(a.queries, sep='\t', dtype={'accession': str})
if a.only: q = q[q.transition_id.isin(a.only)]
qs = pd.read_csv(a.qset, sep='\t', dtype=str)[['accession', 'transition_id', 'set', 'branch_child']].drop_duplicates()
q = q.merge(qs, on=['accession', 'transition_id', 'set'], how='left', suffixes=('', '_qs'))
dropped = {}
if a.exclude_ref_inside:
    from ete3 import Tree
    tree = Tree(open(f'{P}/data/guide_tree_euk.nwk').read(), format=1); nodes = {n.name: n for n in tree.traverse() if n.name}
    ref = pd.read_csv(f'{P}/data/structure_species.tsv', sep='\t').groupby('transition_id').code.apply(set)
    inside = [bool({l.name for l in nodes[c].iter_leaves()} & ref[t]) for c, t in zip(q.branch_child, q.transition_id)]
    dropped['ref_inside_own_clade'] = q[inside].groupby(['transition_id', 'set']).size().to_dict()
    q = q[[not x for x in inside]]
if a.drop_control_euk_lca:
    m = (q.set == 'control') & (pd.to_numeric(q.exclude_taxids, errors='coerce') == 2759)
    dropped['control_euk_lca'] = int(m.sum()); q = q[~m]
q['searched'] = (q.model_status == 'searched').astype(int)
q['old'] = pd.to_numeric(q.structurally_old, errors='coerce').fillna(0).astype(int) * q.searched
q['usable'] = q.searched
q['no_model'] = (q.model_status == 'no_model').astype(int)
q['low_plddt'] = (q.model_status == 'low_plddt').astype(int)
q['other_fail'] = (~q.model_status.isin(['searched', 'no_model', 'low_plddt'])).astype(int)

def wil(k, n):
    return proportion_confint(k, n, method='wilson') if n else (float('nan'), float('nan'))

def two(d, num, den):  # 2x2: stem/control x num/(den-num)
    t = []
    for s in ('stem', 'control'):
        x = d[d.set == s]; t.append([int(x[num].sum()), int(x[den].sum() - x[num].sum())])
    return np.array(t)

def block(d, num, den, name):
    rows = []; tabs = []; ids = []
    for tid, x in list(d.groupby('transition_id')) + [('POOLED', d)]:
        t = two(x, num, den)
        rec = dict(measure=name, transition_id=tid)
        for s, r in zip(('stem', 'control'), t):
            n = int(r.sum()); k = int(r[0]); lo, hi = wil(k, n)
            rec.update({f'{s}_k': k, f'{s}_n': n, f'{s}_frac': k / n if n else float('nan'), f'{s}_ci_lo': lo, f'{s}_ci_hi': hi})
        if t.sum(axis=1).min() > 0:
            rec['fisher_OR'], rec['fisher_p'] = fisher_exact(t)
        rows.append(rec)
        if tid != 'POOLED' and t.sum(axis=1).min() > 0: tabs.append(t); ids.append(tid)
    cmh = {}
    if tabs:
        st = StratifiedTable([tt.astype(float) + 0.0 for tt in tabs])
        r = st.test_null_odds(correction=True); ci = st.oddsratio_pooled_confint()
        cmh = dict(measure=name, mh_OR=st.oddsratio_pooled, mh_ci_lo=ci[0], mh_ci_hi=ci[1], cmh_stat=r.statistic, cmh_p=r.pvalue,
                   n_strata=len(tabs), breslow_day_p=(st.test_equal_odds(adjust=True).pvalue if len(tabs) > 1 else float('nan')))
    return pd.DataFrame(rows), cmh

res = []; cmhs = []
b1, c1 = block(q, 'old', 'searched', 'frac_old_of_searched'); res.append(b1); cmhs.append(c1)
# no usable model: numerator = no usable model, denominator = all sampled
q['nomodel_any'] = 1 - q.searched; q['one'] = 1
b2, c2 = block(q, 'nomodel_any', 'one', 'frac_no_usable_model'); res.append(b2); cmhs.append(c2)
# components, pooled only
extra = []
for s in ('stem', 'control'):
    x = q[q.set == s]
    extra.append(dict(set=s, n=len(x), no_model=int(x.no_model.sum()), low_plddt=int(x.low_plddt.sum()), other=int(x.other_fail.sum()),
                      usable=int(x.searched.sum()), median_plddt=float(pd.to_numeric(x.mean_plddt, errors='coerce').median())))
ex = pd.DataFrame(extra)
R = pd.concat(res); C = pd.DataFrame(cmhs)
R.to_csv(f'{a.out}/per_stem_{a.label}.tsv', sep='\t', index=False); C.to_csv(f'{a.out}/cmh_{a.label}.tsv', sep='\t', index=False)
ex.to_csv(f'{a.out}/model_availability_{a.label}.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.float_format', lambda v: f'{v:.3g}')
print('dropped:', dropped); print(R.to_string(index=False)); print(C.to_string(index=False)); print(ex.to_string(index=False))
if a.append_registry:
    from mcell.branchstats import append_registry
    reg = f'{P}/results/analysis_registry.tsv'
    for name, c in zip(('frac_old_of_searched', 'frac_no_usable_model'), cmhs):
        pooled = R[(R.measure == name) & (R.transition_id == 'POOLED')].iloc[0]
        append_registry(reg, dict(
            hypothesis=f'Q6_structure_{a.label}', tier='exploratory',
            statistic=f'EXPLORATORY: {name}; stem vs matched-control sequence-novel gains (<=60/stem, seed 11); estimate = Mantel-Haenszel common OR (stem/control) stratified by stem; p = CMH (continuity corrected); CI = MH OR CI',
            n=int(pooled.stem_n + pooled.control_n), k=int(pooled.stem_k), estimate=c['mh_OR'], ci_low=c['mh_ci_lo'], ci_high=c['mh_ci_hi'],
            p=c['cmh_p'], p_mcse='NA', holm_p='NA', seed=11, n_draws='NA', scheme='cmh',
            input_md5s=[a.queries],
            notes=f"dropped={dropped}; pooled stem frac={pooled.stem_frac:.3f} ({int(pooled.stem_k)}/{int(pooled.stem_n)}), control frac={pooled.control_frac:.3f} ({int(pooled.control_k)}/{int(pooled.control_n)}); strata={c['n_strata']}; breslow_day_p={c['breslow_day_p']:.3g}; E<=1e-3, TM(q)>=0.5; ref={a.label}; exploratory, not in Holm family"))
    print('registry appended')
