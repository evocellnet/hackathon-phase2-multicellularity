#!/usr/bin/env python3
"""Build Foldseek reference DBs from downloaded AFDB models (runs/structure/afdb_ref/<code>/*.pdb.gz).
Creates one DB per primary transition (union of its P-side species) at runs/structure/refdb/<transition_id>/db
and a combined DB runs/structure/refdb/all/db.  Safe to re-run (rebuilds).  Usage: build_refdb.py [--threads 4]"""
import os, sys, csv, glob, shutil, subprocess, argparse
ap = argparse.ArgumentParser(); ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--only', nargs='*', help='only these transition ids (default all + all-species DB)')
ap.add_argument('--root', default=os.environ.get('P', os.getcwd()))
a = ap.parse_args()
P = a.root
FS = os.environ.get('FOLDSEEK', f'{P}/conda_envs/ancestra/bin/foldseek')
rows = list(csv.DictReader(open(f'{P}/data/structure_species.tsv'), delimiter='\t'))
by_t = {}
for r in rows: by_t.setdefault(r['transition_id'], []).append(r['code'])
by_t['all'] = sorted({r['code'] for r in rows})
os.nice(19)
if a.only: by_t = {k: v for k, v in by_t.items() if k in a.only}
for tid, codes in by_t.items():
    dirs = [f'{P}/runs/structure/afdb_ref/{c}' for c in codes if os.path.isdir(f'{P}/runs/structure/afdb_ref/{c}')]
    out = f'{P}/runs/structure/refdb/{tid}'
    shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
    out_src = f'{out}/src'; os.makedirs(out_src)   # symlink staging: complete *.pdb.gz files only (ignores .part)
    n = 0
    for d in dirs:
        for f in glob.glob(d + '/*.pdb.gz'):
            os.symlink(f, f'{out_src}/{os.path.basename(f)}'); n += 1
    print(tid, codes, n, 'models', flush=True)
    if n == 0:
        shutil.rmtree(out, ignore_errors=True); continue
    subprocess.run(['nice', '-n', '19', FS, 'createdb', out_src, f'{out}/db', '--threads', str(min(a.threads, 4)), '-v', '2'], check=True)
    shutil.rmtree(out_src)
    open(f'{out}/species.txt', 'w').write('\n'.join(codes) + f'\nn_models\t{n}\n')
