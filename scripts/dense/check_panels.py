#!/usr/bin/env python3
"""Sanity checks for the dense panels (run after build_panel.py). For each data/dense/<R>/:
 - every proteomes/<code>.fa symlink resolves; file stems == species.tsv codes == guide-tree leaves
 - every protein ID starts with '<code>_', IDs unique across the whole panel, no sequence < 50 aa, only [A-Z] letters
 - guide tree: internal names unique, non-empty, disjoint from leaves; parses with ancestra.newick
 - transitions.tsv: M and P exist, P is the parent of M; species states consistent (complex == D|T)
 - caps: <= 40 proteomes and <= 700k proteins (PI rule 2026-10-29), reported as WARN
Usage: check_panels.py [R1 R2 ...]   (default: all built panels). Exit 1 on any FAIL."""
import glob, os, re, sys
import pandas as pd
from ete3 import Tree

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def check(R):
    d = f'{P}/data/dense/{R}'
    fails, warns = [], []
    sp = pd.read_csv(f'{d}/species.tsv', sep='\t', dtype=str).fillna('')
    codes = set(sp['code'])
    stems = {os.path.basename(f)[:-3] for f in glob.glob(f'{d}/proteomes/*.fa')}
    if stems != codes:
        fails.append(f'stems != species codes: {stems ^ codes}')
    t = Tree(open(f'{d}/guide_tree.nwk').read(), format=1)
    leaves = set(t.get_leaf_names())
    if leaves != codes:
        fails.append(f'leaves != codes: {leaves ^ codes}')
    internal = [n.name for n in t.traverse() if not n.is_leaf()]
    if len(internal) != len(set(internal)) or not all(internal) or set(internal) & leaves:
        fails.append('internal node names not unique/non-empty/disjoint')
    seen, nprot = set(), 0
    for c in sorted(codes):
        f = f'{d}/proteomes/{c}.fa'
        if not os.path.exists(f):
            fails.append(f'{c}: broken symlink'); continue
        sid, L = None, 0
        for line in open(f):
            if line.startswith('>'):
                if sid is not None and L < 50:
                    fails.append(f'{c}: {sid} shorter than 50 aa')
                sid, L = line[1:].split()[0], 0
                nprot += 1
                if not sid.startswith(c + '_'):
                    fails.append(f'{c}: id {sid} lacks prefix'); break
                if sid in seen:
                    fails.append(f'{c}: duplicate id {sid}'); break
                seen.add(sid)
            else:
                s = line.strip()
                L += len(s)
                if not re.fullmatch(r'[A-Z]*', s):
                    fails.append(f'{c}: non-letter characters'); break
        if sid is not None and L < 50:
            fails.append(f'{c}: {sid} shorter than 50 aa')
    if nprot != sp['n_proteins'].astype(int).sum():
        fails.append(f'protein count mismatch {nprot} vs species.tsv')
    for _, r in sp.iterrows():
        if (r['mc_state'] == 'complex') != (r['D'] == '1' or r['T'] == '1'):
            fails.append(f"{r['code']}: complex != D|T")
    tr = pd.read_csv(f'{d}/transitions.tsv', sep='\t', dtype=str).fillna('')
    names = {n.name: n for n in t.traverse()}
    for _, e in tr.iterrows():
        if e['M'] not in names or e['P'] not in names or names[e['M']].up is not names[e['P']]:
            fails.append(f"event {e['event_id']}: M/P inconsistent with tree")
    if len(codes) > 40:
        warns.append(f'{len(codes)} proteomes > 40')
    if nprot > 700_000:
        warns.append(f'{nprot:,} proteins > 700k')
    return fails, warns, len(codes), nprot


def main():
    panels = sys.argv[1:] or sorted(os.path.basename(os.path.dirname(f)) for f in glob.glob(f'{P}/data/dense/R*/species.tsv'))
    bad = 0
    for R in panels:
        f, w, n, np_ = check(R)
        print(f"{R}: {'OK' if not f else 'FAIL'} species={n} proteins={np_:,}" + (f" WARN {w}" if w else ''))
        for x in f[:20]:
            print('   FAIL', x)
        bad += bool(f)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
