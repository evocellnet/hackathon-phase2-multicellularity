"""Task 2 spot-check (independent of scripts/r2/panels and mcell.gainladder).
Raw gains re-derived with mcell.hogs.parse_hogs; for 50 random gained families per branch, own representatives
(longest gene of up to 3 species, spread over child clades; tips: 3 longest genes) are searched with DIAMOND
(--very-sensitive -k 0 -e 1e-5 --query-cover 30 --max-hsps 1) against the panel DB (db.dmnd; composition checked);
novel = no hit to any species outside the clade (panel leaves below the child node + round-1 species placed at/below it)."""
import os, sys, json, subprocess, numpy as np, pandas as pd
P = '/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
sys.path.insert(0, P)
from mcell.hogs import parse_hogs
OUT = P + '/results/r2/verification/panels/spotcheck'; os.makedirs(OUT, exist_ok=True)
CHECK = {'R4': ['Volvcar', 'Chlamydomonadales'], 'R6': ['Fucophycidae', 'Ochrophyta']}
DIAMOND = P + '/conda_envs/ancestra/bin/diamond'
rng = np.random.default_rng(20261029)
summary, famrows = [], []
def read_fa(path, keep):
    seqs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                cur = line[1:].split()[0]; cur = cur if cur in keep else None
                if cur: seqs[cur] = []
            elif cur: seqs[cur].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}
for R, branches in CHECK.items():
    nwk = open(f'{P}/runs/r2/{R}/stages/input/guide_tree.nwk').read().strip()
    assert nwk == open(f'{P}/data/dense/{R}/guide_tree.nwk').read().strip()
    res = parse_hogs(f'{P}/runs/r2/{R}/hogs.orthoxml', newick=nwk)
    T, ev, genes = res.tree, res.events, res.genes
    viol = len(res.check_conservation())
    bt = pd.read_csv(f'{P}/results/r2/panels/{R}/branch_table.tsv', sep='\t').set_index('branch_child')
    gc = pd.read_csv(f'{P}/results/r2/panels/{R}/gains_classified.tsv.gz', sep='\t')
    place = pd.read_csv(f'{P}/results/r2/panels/{R}/round1_placement.tsv', sep='\t').fillna('')
    tg = pd.read_parquet(f'{P}/results/r2/panels/{R}/work/ladder/db.targets.parquet', columns=['seqid', 'species', 'source'])
    sp_of = dict(zip(tg.seqid, tg.species))
    # DB composition check: panel species + round-1 euk not in panel
    r1 = set(pd.read_csv(f'{P}/data/proteomes/summary.tsv', sep='\t').query("kingdom=='euk'").code)
    panel = set(T.names[i] for i in range(len(T.names)) if T.is_leaf[i])
    dbsp = set(tg.species)
    db_ok = dbsp == (panel | r1)
    def below(node):
        out, st = [], [T.index[node]]
        while st:
            k = st.pop(); out.append(k); st += list(T.children[k])
        return out
    for b in branches:
        g = ev[(ev.branch_child == b) & (ev.event == 'gain')]
        nodes = below(b); leaves = {T.names[k] for k in nodes if T.is_leaf[k]}
        clade = leaves | set(place.loc[place.panel_node.isin({T.names[k] for k in nodes}), 'species'])
        samp = rng.choice(g.family_id.unique(), size=50, replace=False)
        is_tip = bool(T.is_leaf[T.index[b]])
        if not is_tip:
            mem = res.members(b)
            child_of = {}
            for c in T.children[T.index[b]]:
                for k in below(T.names[c]):
                    if T.is_leaf[k]: child_of[T.names[k]] = T.names[c]
        q = {}
        for f in samp:
            row = g[g.family_id == f].iloc[0]
            ids = str(row.child_hog_ids).split(';')
            if is_tip:
                gl = ids
            else:
                gl = list(mem.loc[mem.hog_id.isin(ids), 'gene'])
            q[f] = gl
        allg = {x for v in q.values() for x in v}
        seqs = {}
        for sp in {x.split('_')[0] for x in allg}:
            seqs.update(read_fa(f'{P}/data/dense/{R}/proteomes/{sp}.fa', allg))
        reps = {}
        for f, gl in q.items():
            gl = sorted([x for x in gl if x in seqs], key=lambda x: (-len(seqs[x]), x))
            if is_tip: reps[f] = gl[:3]; continue
            chosen, used_sp, used_cl = [], set(), set()
            for x in gl:   # first pass: one per child clade
                s = sp_of.get(x, x.split('_')[0])
                if child_of.get(s) not in used_cl and s not in used_sp:
                    chosen.append(x); used_sp.add(s); used_cl.add(child_of.get(s))
                if len(chosen) == 3: break
            for x in gl:
                if len(chosen) == 3: break
                s = sp_of.get(x, x.split('_')[0])
                if s not in used_sp: chosen.append(x); used_sp.add(s)
            reps[f] = chosen
        qfa = f'{OUT}/{R}_{b}_queries.faa'
        with open(qfa, 'w') as fh:
            for f, rl in reps.items():
                for x in rl: fh.write(f'>{f}|{x}\n{seqs[x]}\n')
        hits = f'{OUT}/{R}_{b}_hits.tsv'
        if not os.path.exists(hits):
            subprocess.run(['nice', '-n', '10', DIAMOND, 'blastp', '-q', qfa, '-d', f'{P}/results/r2/panels/{R}/work/ladder/db.dmnd',
                            '-o', hits, '--very-sensitive', '-k', '0', '-e', '1e-5', '--query-cover', '30', '--max-hsps', '1',
                            '--threads', '8', '-f', '6', 'qseqid', 'sseqid', 'pident', 'evalue', '-b', '2', '-c', '1', '--quiet'], check=True)
        H = pd.read_csv(hits, sep='\t', header=None, names=['q', 's', 'pid', 'e'])
        H['fam'] = H.q.str.split('|').str[0]; H['ssp'] = H.s.map(sp_of)
        out_hits = H[~H.ssp.isin(clade)].groupby('fam').size()
        prod = gc[gc.branch_child == b].set_index('family_id')
        for f in samp:
            mine = int(out_hits.get(f, 0)) == 0
            famrows.append(dict(panel=R, branch=b, family_id=f, n_reps=len(reps[f]), n_outside_hits=int(out_hits.get(f, 0)),
                                my_seq_novel=mine, producer_seq_novel=bool(prod.loc[f, 'seq_novel']) if f in prod.index else None,
                                producer_split_repaired=bool(prod.loc[f, 'split_repaired']) if f in prod.index else None))
        fr = pd.DataFrame([r for r in famrows if r['branch'] == b])
        summary.append(dict(panel=R, branch=b, tip=is_tip, my_hog_gains=len(g), table_hog_gains=int(bt.loc[b, 'hog_gains']),
                            table_seq_novel=int(bt.loc[b, 'seq_novel']), table_frac_novel=round(bt.loc[b, 'seq_novel'] / bt.loc[b, 'hog_gains'], 3),
                            sample_frac_novel_mine=round(fr.my_seq_novel.mean(), 3), agree_with_producer=int((fr.my_seq_novel == fr.producer_seq_novel).sum()),
                            n_sample=len(fr), conservation_violations=viol, db_species_match_spec=db_ok,
                            n_db_species=len(dbsp), clade_round1_added=sorted(clade - leaves)))
        print(summary[-1], flush=True)
pd.DataFrame(famrows).to_csv(f'{OUT}/family_level.tsv', sep='\t', index=False)
json.dump(summary, open(f'{OUT}/summary.json', 'w'), indent=1, default=str)
