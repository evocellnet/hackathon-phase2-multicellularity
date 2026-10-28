#!/usr/bin/env python3
"""Look up species names (one per line on stdin, or args) in the merged source index; print best hit per db."""
import csv, sys, os
R = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'raw', 'source_index.tsv')
rows = list(csv.DictReader(open(R), delimiter='\t'))
def f(x):
    try: return float(x)
    except: return -1
names = sys.argv[1:] or [l.strip() for l in sys.stdin if l.strip()]
for n in names:
    hits = [r for r in rows if r['name'].lower().startswith(n.lower())]
    if not hits:
        print(f'{n}\tNONE'); continue
    best = {}
    for r in hits:
        k = r['db']; 
        score = (r['data_type'] == 'genome', f(r['busco_C']))
        if k not in best or score > best[k][0]: best[k] = (score, r)
    s = ' | '.join(f"{k}:{r['id']} {r['data_type'][:5]} B{r['busco_C']} n{r['n_proteins']} tx{r['taxid']}" for k, (sc, r) in best.items())
    print(f'{n}\t{len(hits)}\t{s}')
