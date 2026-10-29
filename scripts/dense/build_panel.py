#!/usr/bin/env python3
"""Build one dense per-event panel (round 2) from its curated spec.

Inputs (curated, committed):
  scripts/dense/specs/<R>.tsv   one row per species: code, species, taxid, source, data_type, mc_state,
                                mc_conf (clear|grey), mc_ref, role (focal|control|outgroup), group, note
                                source = R1:<code>        reuse round-1 file data/proteomes/{euk,prok}/<code>.fa
                                         UniProt:UPxxxxxxxxx | NCBI:GCx_xxxxxxxxx.x | EukProt:EPxxxxx | LukProt:LPxxxxx/APxxxxx
  scripts/dense/panels.py       per panel: guide-tree edits (same grammar as scripts/data/build_guide_trees.py),
                                events (transitions), BUSCO lineage.
Steps (all re-runnable, cached):
  1. fetch    new raw files (scripts/data/fetch_proteomes.job: UniProt 2026_03 FTP, NCBI Datasets v2 zip, local EukProt/LukProt)
              -> manifest data/dense/download_manifest_dense.tsv
  2. process  scripts/data/preprocess_proteomes.process (one protein per gene, cleaning, >=50 aa, IDs <code>_<id>,
              CD-HIT 0.95 for transcriptomes) -> data/proteomes/dense/<code>.fa (+ idmaps/, stats/)
  3. link     data/dense/<R>/proteomes/<code>.fa -> relative symlink to the processed file
  4. busco    BUSCO 6.1.0 protein mode (scripts/data/run_busco.run) on new files; reused files take round-1 values
              -> data/dense/busco_summary_dense.tsv
  5. tree     ancestra.speciestree ncbi on <code>_<taxid> symlinks + panel edits -> data/dense/<R>/guide_tree.nwk,
              tree_edits.md, guide_tree_ncbi.nwk
  6. tables   data/dense/<R>/species.tsv, transitions.tsv
Usage: build_panel.py R1 [--steps fetch,process,link,busco,tree,tables] [--jobs 4] [--cpu 4]
"""
import argparse, csv, datetime, glob, hashlib, importlib.util, json, os, re, subprocess, sys, tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, f'{P}/scripts/data')
sys.path.insert(0, f'{P}/scripts/dense')
import fetch_proteomes as FP          # noqa: E402
import preprocess_proteomes as PP     # noqa: E402
import run_busco as RB                # noqa: E402
import build_guide_trees as GT        # noqa: E402
import panels as PANELS               # noqa: E402

PY = f'{P}/conda_envs/ancestra/bin/python'
MAN = f'{P}/data/dense/download_manifest_dense.tsv'
BUSCO_DENSE = f'{P}/data/dense/busco_summary_dense.tsv'
R1_TABLES = {'euk': f'{P}/data/species_final.tsv', 'prok': f'{P}/data/species_final_prok.tsv'}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def r1_index():
    idx = {}
    for k, f in R1_TABLES.items():
        d = pd.read_csv(f, sep='\t', dtype=str)
        for _, r in d.iterrows():
            idx[r['code']] = dict(kingdom=k, species=r['species'], taxid=r['ncbi_taxid'], source=r['source'],
                                  data_type=r['data_type_final'], busco=(r['busco_C'], r['busco_S'], r['busco_D'],
                                                                         r['busco_F'], r['busco_M']),
                                  lineage=r['busco_lineage'], md5=r['md5'], n=r['n_proteins_final'])
    return idx


def load_spec(R):
    d = pd.read_csv(f'{P}/scripts/dense/specs/{R}.tsv', sep='\t', dtype=str, comment='#').fillna('')
    d = d[d['code'] != '']
    assert d['code'].is_unique, d['code'][d['code'].duplicated()].tolist()
    r1 = r1_index()
    rows = []
    for _, r in d.iterrows():
        r = r.to_dict()
        if r['source'].startswith('R1:'):
            c = r['source'][3:]
            assert c == r['code'], f"R1 reuse must keep the round-1 code ({r['code']} vs {c})"
            m = r1[c]
            r['species'] = r['species'] or m['species']
            r['taxid'] = r['taxid'] or m['taxid']
            r['data_type'] = r['data_type'] or m['data_type']
            r['orig_source'] = m['source']
            r['file'] = f"{P}/data/proteomes/{m['kingdom']}/{c}.fa"
        else:
            assert r['code'] not in r1, f"{r['code']}: code already used in round 1 (choose another)"
            r['orig_source'] = r['source']
            r['file'] = f"{P}/data/proteomes/dense/{r['code']}.fa"
        assert r['mc_state'] in ('unicellular', 'simple', 'complex'), (r['code'], r['mc_state'])
        assert r['D'] in ('0', '1') and r['T'] in ('0', '1'), (r['code'], r['D'], r['T'])
        # PI rule (design review round 2): complex = D or T (D: germ-soma split or >=3 cell types; T: 3D tissue)
        assert (r['mc_state'] == 'complex') == (r['D'] == '1' or r['T'] == '1'), (r['code'], 'complex != D|T')
        assert r['role'] in ('focal', 'control', 'outgroup'), (r['code'], r['role'])
        assert r['taxid'].isdigit(), (r['code'], r['taxid'])
        rows.append(r)
    return rows, r1


# ----------------------------------------------------------------------------------------- fetch/process
def load_manifest():
    if os.path.exists(MAN):
        return pd.read_csv(MAN, sep='\t', dtype=str)
    return pd.DataFrame(columns=['code', 'source', 'url', 'local_path', 'md5', 'bytes', 'date'])


def step_fetch(rows, kingdom, threads):
    new = [r for r in rows if not r['source'].startswith(('R1:', 'AUGUSTUS:'))]
    man = load_manifest()
    have = set(man['code'])
    todo = [r for r in new if r['code'] not in have or not all(
        os.path.exists(f'{P}/{p}') for p in man.loc[man['code'] == r['code'], 'local_path'])]
    if not todo:
        print('fetch: nothing to do'); return
    FP.check_release()
    uptax = FP.up_taxid()
    res, errs = [], []
    with ThreadPoolExecutor(threads) as ex:
        futs = {ex.submit(FP.job, dict(r, lineage_group='Bacteria' if kingdom == 'prok' else 'Eukaryota'),
                          kingdom, uptax): r['code'] for r in todo}
        for fu in as_completed(futs):
            try:
                res += fu.result()
                print('fetched', futs[fu], flush=True)
            except Exception as e:
                errs.append((futs[fu], str(e))); print('ERROR fetch', futs[fu], e, flush=True)
    today = datetime.date.today().isoformat()
    add = pd.DataFrame(res, columns=['code', 'source', 'url', 'local_path'])
    add['md5'] = [md5(p) for p in add['local_path']]
    add['bytes'] = [str(os.path.getsize(p)) for p in add['local_path']]
    add['date'] = today
    add['local_path'] = [os.path.relpath(p, P) for p in add['local_path']]
    man = pd.concat([man[~man['local_path'].isin(add['local_path'])], add])
    man.sort_values(['code', 'local_path']).to_csv(MAN, sep='\t', index=False)
    if errs:
        print('FETCH ERRORS:', errs)


def process_augustus(r):
    """AUGUSTUS ab initio predictions (scripts/dense/augustus_predict.sh; one transcript per gene): same cleaning as
    round 1 (terminal stop trimmed, non-standard -> X, >= 50 aa, exact duplicates removed, IDs <code>_<id>)."""
    acc = r['source'].split(':', 1)[1]
    src = f'{P}/data/raw/augustus/{acc}/proteins.faa'
    out = r['file']
    if os.path.exists(out):
        return
    seen, n_raw, n_short, n_dup, kept = set(), 0, 0, 0, []
    for h, sq in PP.read_fasta(src):
        n_raw += 1
        sq = PP.clean(sq)
        if len(sq) < PP.MINLEN:
            n_short += 1; continue
        if sq in seen:
            n_dup += 1; continue
        seen.add(sq); kept.append((f"{r['code']}_{PP.sanitise(h.split()[0])}", sq))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out + '.tmp', 'w') as fh:
        for i, sq in kept:
            fh.write(f'>{i}\n' + '\n'.join(sq[k:k + 60] for k in range(0, len(sq), 60)) + '\n')
    os.replace(out + '.tmp', out)
    json.dump(dict(code=r['code'], source=r['source'], n_raw=n_raw, n_dropped_short=n_short, n_dropped_identical=n_dup,
                   n_final=len(kept), md5=md5(out), method='AUGUSTUS chlamy2011 ab initio; one transcript per gene'),
              open(f"{P}/data/proteomes/stats/{r['code']}.json", 'w'), indent=1)
    print(f"{r['code']}\traw={n_raw}\tshort={n_short}\tident={n_dup}\tfinal={len(kept)} (AUGUSTUS)")


def step_process(rows, jobs, cd_threads):
    for r in rows:
        if r['source'].startswith('AUGUSTUS:'):
            process_augustus(r)
    new = [r for r in rows if not r['source'].startswith(('R1:', 'AUGUSTUS:'))]
    man = load_manifest()
    manifest = defaultdict(list)
    for c, lp in zip(man['code'], man['local_path']):
        manifest[c].append(f'{P}/{lp}')
    errs = 0
    with ProcessPoolExecutor(jobs) as ex:
        futs = {ex.submit(PP.process, dict(code=r['code'], source=r['source'], data_type=r['data_type']),
                          'dense', manifest, cd_threads, False): r['code'] for r in new}
        for fu in as_completed(futs):
            c = futs[fu]
            try:
                st = fu.result()
                print(f"{c}\traw={st['n_raw']}\tgenes={st['n_after_gene_collapse']}\tshort={st['n_dropped_short']}"
                      f"\tident={st['n_dropped_identical']}\tcdhit={st['n_dropped_cdhit']}\tfinal={st['n_final']}", flush=True)
            except Exception as e:
                errs += 1; print(f'ERROR process {c}: {e}', flush=True)
    return errs


def step_link(R, rows):
    d = f'{P}/data/dense/{R}/proteomes'
    os.makedirs(d, exist_ok=True)
    want = {r['code'] for r in rows}
    for f in glob.glob(f'{d}/*.fa'):
        if os.path.basename(f)[:-3] not in want:
            os.remove(f)
    for r in rows:
        if not os.path.exists(r['file']):
            print('MISSING processed file', r['code']); continue
        dst = f"{d}/{r['code']}.fa"
        rel = os.path.relpath(r['file'], d)
        if os.path.islink(dst) and os.readlink(dst) == rel:
            continue
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(rel, dst)


def step_busco(rows, lineage, jobs, cpu):
    old = pd.read_csv(BUSCO_DENSE, sep='\t', dtype=str) if os.path.exists(BUSCO_DENSE) else pd.DataFrame()
    todo = [r['code'] for r in rows if not r['source'].startswith('R1:') and os.path.exists(r['file'])]
    out = []
    errs = 0
    with ThreadPoolExecutor(jobs) as ex:
        futs = {ex.submit(RB.run, c, 'dense', lineage, cpu): c for c in todo}
        for fu in as_completed(futs):
            try:
                code, res, how = fu.result()
                out.append(dict(code=code, **res))
                print(f"BUSCO {code} C={res['C']} D={res['D']} F={res['F']} ({how})", flush=True)
            except Exception as e:
                errs += 1; print('ERROR busco', futs[fu], e, flush=True)
    new = pd.DataFrame(out).astype(str)
    if len(old) and len(new):
        old = old[~((old['code'].isin(new['code'])) & (old['lineage'].isin(new['lineage'])))]
    allb = pd.concat([old, new]) if len(old) else new
    allb.sort_values(['code', 'lineage']).to_csv(BUSCO_DENSE, sep='\t', index=False)
    return errs


# --------------------------------------------------------------------------------------------- tree
def step_tree(R, rows, cfg):
    from ete3 import Tree
    d = f'{P}/data/dense/{R}'
    with tempfile.TemporaryDirectory(dir=f'{P}/work') as td:
        for r in rows:
            os.symlink(r['file'], f"{td}/{r['code']}_{int(r['taxid'])}.fa")
        out = f'{td}/raw.nwk'
        cmd = [PY, '-m', 'ancestra.speciestree', 'ncbi', '--proteomes', td, '--taxdump', f'{P}/data/raw/taxdump',
               '--leaf-naming', 'qfo', '--out', out]
        rr = subprocess.run(cmd, capture_output=True, text=True, cwd=P)
        if rr.returncode:
            sys.exit(rr.stderr)
        t = Tree(open(out).read(), format=1, quoted_node_names=True)
    for nd in t.traverse():
        nd.name = GT.sanitise(nd.name)
        if nd.is_leaf():
            nd.name = re.sub(r'_\d+$', '', nd.name)
    # internal nodes that coincide with a leaf name (NCBI species node above a strain leaf) are renamed
    leaves = {n.name for n in t.iter_leaves()}
    t = GT.collapse_unary(t, [])
    for nd in t.traverse():
        if not nd.is_leaf() and (nd.name in leaves or not nd.name):
            nd.name = (nd.name or 'node') + '_clade'
    open(f'{d}/guide_tree_ncbi.nwk', 'w').write(GT.newick(GT.canonical(t.copy())) + ';\n')
    log = []
    GT.apply(t, cfg['edits'], log)
    t = GT.collapse_unary(t, log)
    t = GT.canonical(t)
    codes = {r['code'] for r in rows}
    lv = [n.name for n in t.iter_leaves()]
    assert len(lv) == len(set(lv)) and set(lv) == codes, set(lv) ^ codes
    stems = {os.path.basename(f)[:-3] for f in glob.glob(f'{d}/proteomes/*.fa')}
    assert stems == codes, stems ^ codes
    internal = [n.name for n in t.traverse() if not n.is_leaf()]
    assert all(internal) and len(internal) == len(set(internal)), [x for x in internal if internal.count(x) > 1]
    assert not set(internal) & set(lv)
    assert all(re.fullmatch(r'[A-Za-z0-9_.\-]+', x) for x in internal + lv)
    for name in cfg.get('expect', []):
        GT.find(t, name)
    nwk = f'{d}/guide_tree.nwk'
    open(nwk, 'w').write(GT.newick(t) + ';\n')
    rr = subprocess.run([PY, '-c', f'from ancestra.newick import read_newick; t=read_newick("{nwk}"); print(len(t.leaves()))'],
                        capture_output=True, text=True, cwd=P)
    assert rr.stdout.strip() == str(len(lv)), rr.stderr
    n_poly = sum(1 for n in t.traverse() if len(n.children) > 2)
    md = [f'# {R} guide-tree edits ({cfg["title"]})', '',
          'Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy '
          '(`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on '
          '`<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; '
          'children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 '
          '(`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep '
          'the multifurcation.', '',
          f'Root: `{t.name}` with children ' + ', '.join(f'`{c.name}`' for c in t.children) +
          f'. {len(lv)} leaves, {len(internal)} named internal nodes, {n_poly} multifurcating nodes.', '',
          '| # | operation | nodes | rationale | reference |', '|---|---|---|---|---|']
    for i, op, args, why, ref in log:
        if op == 'collapse_unary':
            continue
        a = (f'{args[0]} := ({", ".join(args[1])})' if op == 'group' else ' -> '.join(map(str, args)))
        md.append(f'| {i} | {op} | {a} | {why} | {ref} |')
    if cfg.get('tree_notes'):
        md += ['', '## Open issues / notes', ''] + [f'- {x}' for x in cfg['tree_notes']]
    open(f'{d}/tree_edits.md', 'w').write('\n'.join(md) + '\n')
    print(f'{R}: tree leaves={len(lv)} internal={len(internal)} multifurcations={n_poly}')
    return t


# --------------------------------------------------------------------------------------------- tables
def step_tables(R, rows, cfg, r1):
    from ete3 import Tree
    d = f'{P}/data/dense/{R}'
    t = Tree(open(f'{d}/guide_tree.nwk').read(), format=1)
    bd = pd.read_csv(BUSCO_DENSE, sep='\t', dtype=str) if os.path.exists(BUSCO_DENSE) else pd.DataFrame(columns=['code', 'lineage'])
    lineage = cfg['busco_lineage']
    man = load_manifest()
    out = []
    state = {}
    for r in rows:
        c = r['code']
        if r['source'].startswith('R1:'):
            m = r1[c]
            C, S, D, F, M = m['busco']
            bl = m['lineage']
            src = f"{m['source']} (round-1 file data/proteomes/{m['kingdom']}/{c}.fa)"
            url = ''
        else:
            b = bd[(bd['code'] == c) & (bd['lineage'] == lineage)]
            if len(b):
                b = b.iloc[0]; C, S, D, F, M, bl = b['C'], b['S'], b['D'], b['F'], b['M'], b['lineage']
            else:
                C = S = D = F = M = ''; bl = ''
            src = r['source']
            u = man.loc[man['code'] == c, 'url'].tolist()
            url = u[0] if u else ''
        n = res = 0
        if os.path.exists(r['file']):
            for line in open(r['file']):
                if line.startswith('>'):
                    n += 1
                else:
                    res += len(line.strip())
        flags = []
        if r['data_type'] != 'genome':
            flags.append('TRANSCRIPTOME' if 'transcriptome' in r['data_type'] else r['data_type'].upper())
        try:
            if float(C) < 70:
                flags.append('BUSCO_C<70')
        except ValueError:
            flags.append('NO_BUSCO')
        ann = 'augustus_denovo' if r['source'].startswith('AUGUSTUS:') else 'curated'
        out.append(dict(code=c, species=r['species'], taxid=r['taxid'], source=src, url=url, data_type=r['data_type'],
                        annotation_source=ann,
                        busco_lineage=bl, busco_C=C, busco_S=S, busco_D=D, busco_F=F, busco_M=M,
                        mc_state=r['mc_state'], D=r['D'], T=r['T'], mc_conf=r['mc_conf'], mc_ref=r['mc_ref'], role=r['role'],
                        group=r['group'], n_proteins=n, n_residues=res, md5=md5(r['file']) if os.path.exists(r['file']) else '',
                        flags=';'.join(flags), note=r['note']))
        state[c] = r['mc_state']
    sp = pd.DataFrame(out)
    sp.to_csv(f'{d}/species.tsv', sep='\t', index=False)
    # transitions
    ev = []
    for e in cfg['events']:
        M = GT.find(t, e['M'])
        Pn = M.up
        sis = [c for c in Pn.children if c is not M]
        mleaves = sorted(M.get_leaf_names())
        sleaves = sorted(x for s in sis for x in s.get_leaf_names())
        ev.append(dict(event_id=e['id'], label=e['label'], kind=e['kind'], M=M.name, M_is_tip=int(M.is_leaf()),
                       P=Pn.name, sister_nodes=';'.join(s.name for s in sis), from_state=e['from'], to_state=e['to'],
                       confidence=e.get('conf', 'clear'), n_M=len(mleaves), n_sister=len(sleaves),
                       M_states=','.join(f'{k}:{v}' for k, v in sorted(pd.Series([state[x] for x in mleaves]).value_counts().items())),
                       sister_states=','.join(f'{k}:{v}' for k, v in sorted(pd.Series([state[x] for x in sleaves]).value_counts().items())) if sleaves else '',
                       M_species=','.join(mleaves), sister_species=','.join(sleaves), ref=e.get('ref', ''), note=e.get('note', '')))
    pd.DataFrame(ev).to_csv(f'{d}/transitions.tsv', sep='\t', index=False)
    print(f'{R}: species={len(sp)} proteins={sp["n_proteins"].sum():,} residues={sp["n_residues"].sum():,} states={sp["mc_state"].value_counts().to_dict()} '
          f'events={len(ev)} genome={sum(sp.data_type == "genome")} lowC={sum(sp["flags"].str.contains("BUSCO_C<70"))}')
    return sp, ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('panel')
    ap.add_argument('--steps', default='fetch,process,link,busco,tree,tables')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--cpu', type=int, default=4)
    a = ap.parse_args()
    assert a.jobs * a.cpu <= 16, 'CPU limit 16'
    cfg = PANELS.PANELS[a.panel]
    rows, r1 = load_spec(a.panel)
    kingdom = cfg.get('kingdom', 'euk')
    steps = a.steps.split(',')
    os.makedirs(f'{P}/work', exist_ok=True)
    if 'fetch' in steps:
        step_fetch(rows, kingdom, 6)
    if 'process' in steps:
        step_process(rows, a.jobs, a.cpu)
    if 'link' in steps:
        step_link(a.panel, rows)
    if 'busco' in steps:
        step_busco(rows, cfg['busco_lineage'], a.jobs, a.cpu)
    if 'tree' in steps:
        step_tree(a.panel, rows, cfg)
    if 'tables' in steps:
        step_tables(a.panel, rows, cfg, r1)


if __name__ == '__main__':
    main()
