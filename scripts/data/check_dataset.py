#!/usr/bin/env python3
"""Step 8: invariant checks on the final dataset (exit 1 on any violation).

- every selected species has exactly one FASTA <code>.fa; no extra files in data/proteomes/{euk,prok}
- protein IDs: start with '<code>_', only [A-Za-z0-9_.-], unique across ALL species
- sequences: only the 20 standard letters + X, length >= 50, no '*', no exact duplicate within a species
- md5 and protein count of every FASTA equal the values in data/species_final*.tsv
- guide-tree leaves == file stems; transitions.tsv M/P names exist in the tree
- idmaps: one row per raw protein; kept rows == FASTA IDs
"""
import gzip, hashlib, os, re, sys
import pandas as pd
from ete3 import Tree

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OK = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYX]+$')
errs = []


def err(m):
    errs.append(m)
    print('FAIL', m)


def main():
    seen_ids = set()
    tr = pd.read_csv(f'{P}/data/transitions.tsv', sep='\t', dtype=str, keep_default_na=False)
    for kingdom, f in (('euk', 'species_final.tsv'), ('prok', 'species_final_prok.tsv')):
        d = pd.read_csv(f'{P}/data/{f}', sep='\t', dtype=str)
        allf = os.listdir(f'{P}/data/proteomes/{kingdom}')
        extra = [x for x in allf if not x.endswith('.fa') or x[:-3] not in set(d['code'])]
        if extra:
            err(f'{kingdom}: extra entries in the Ancestra input directory: {extra}')
        files = {x[:-3] for x in allf if x.endswith('.fa')}
        if files != set(d['code']):
            err(f'{kingdom}: files vs table differ: {files ^ set(d["code"])}')
        t = Tree(open(f'{P}/data/guide_tree_{kingdom}.nwk').read(), format=1)
        if {l.name for l in t.iter_leaves()} != files:
            err(f'{kingdom}: tree leaves != files')
        names = {n.name for n in t.traverse()}
        for _, r in tr[tr.kingdom == kingdom].iterrows():
            if r['M'] == 'NA':
                continue
            for x in (r['M'].split(',') + [r['P']]):
                if x not in names:
                    err(f'{r.transition_id}: node {x} not in tree')
        ntot = 0
        for _, r in d.iterrows():
            c = r['code']
            path = f'{P}/data/proteomes/{kingdom}/{c}.fa'
            raw = open(path, 'rb').read()
            if hashlib.md5(raw).hexdigest() != r['md5']:
                err(f'{c}: md5 mismatch')
            ids, seqs, cur = [], [], []
            for line in raw.decode().splitlines():
                if line.startswith('>'):
                    if cur:
                        seqs.append(''.join(cur))
                    ids.append(line[1:])
                    cur = []
                else:
                    cur.append(line)
            seqs.append(''.join(cur))
            if len(ids) != int(r['n_proteins_final']):
                err(f'{c}: count {len(ids)} != table {r["n_proteins_final"]}')
            for i in ids:
                if not i.startswith(c + '_') or not re.fullmatch(r'[A-Za-z0-9_.\-]+', i):
                    err(f'{c}: bad id {i}')
                    break
                if i in seen_ids:
                    err(f'{c}: duplicate id across species {i}')
                    break
                seen_ids.add(i)
            bad = [s for s in seqs if len(s) < 50 or not OK.match(s)]
            if bad:
                err(f'{c}: {len(bad)} sequences short or with invalid letters')
            if len(set(seqs)) != len(seqs):
                err(f'{c}: exact duplicate sequences')
            m = pd.read_csv(f'{P}/data/proteomes/idmaps/{c}.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
            if set(m.loc[m.status == 'kept', 'new_id']) != set(ids):
                err(f'{c}: idmap kept rows != FASTA ids')
            if len(m) != int(r['n_proteins_raw']):
                err(f'{c}: idmap rows != n_proteins_raw')
            ntot += len(ids)
        print(f'{kingdom}: {len(d)} species, {ntot} proteins checked')
    print('OK' if not errs else f'{len(errs)} failures')
    sys.exit(1 if errs else 0)


if __name__ == '__main__':
    main()
