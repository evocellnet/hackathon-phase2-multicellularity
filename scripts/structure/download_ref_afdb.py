#!/usr/bin/env python3
"""Download AFDB models for all proteins of the reference species in data/structure_species.tsv.
Resumable: existing files are skipped; status log appended to runs/structure/afdb_ref/download_log.tsv
(acc, code, status, version, mean_plddt, bytes). Missing entries -> missing.tsv."""
import os, sys, csv, time
sys.path.insert(0, os.path.dirname(__file__))
import afdb_fetch as A
P = os.environ.get('P', os.getcwd())
OUT = f'{P}/runs/structure/afdb_ref'
os.makedirs(OUT, exist_ok=True)
workers = int(os.environ.get('WORKERS', 8))
species = []
for r in csv.DictReader(open(f'{P}/data/structure_species.tsv'), delimiter='\t'):
    if r['code'] not in species: species.append(r['code'])
only = sys.argv[1:] or species
t0 = time.time()
# entries already known to be absent from AFDB (do not re-request on resume)
known_missing = set()
if os.path.exists(f'{OUT}/download_log.tsv'):
    for l in open(f'{OUT}/download_log.tsv'):
        f = l.rstrip('\n').split('\t')
        if len(f) > 2 and f[2] == 'missing': known_missing.add(f[0])
logf = open(f'{OUT}/download_log.tsv', 'a')
if logf.tell() == 0: logf.write('acc\tcode\tstatus\tversion\tmean_plddt\tbytes\n')
for code in only:
    accs = [l[1:].strip().split('_', 1)[1] for l in open(f'{P}/data/proteomes/euk/{code}.fa') if l[0] == '>']
    d = f'{OUT}/{code}'
    todo = [a for a in accs if a not in known_missing]
    res = A.fetch_many(todo, d, workers, on_progress=lambda k, n, c=code: print(c, k, n, f'{time.time()-t0:.0f}s', flush=True))
    for a in accs:
        if a in known_missing: continue
        st, path, ver = res[a]
        if st == 'error': continue   # not logged: retried on the next (resume) run
        pl = A.mean_plddt_from_pdb_gz(path) if path else ''
        logf.write(f'{a}\t{code}\t{st}\t{ver}\t{pl if pl=="" else round(pl,2)}\t{os.path.getsize(path) if path else 0}\n')
    logf.flush()
    c = [v[0] for v in res.values()]
    if c.count('error') > 0.2 * max(1, len(c)): sys.exit(f'{code}: >20% errors, aborting (check network / rate limit)')
    print(code, {s: c.count(s) for s in set(c)}, f'{time.time()-t0:.0f}s', flush=True)
