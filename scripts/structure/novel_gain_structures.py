#!/usr/bin/env python3
"""Structural-novelty test for sequence-novel gains (and control-branch gains).

For each query UniProt accession: fetch its AFDB model (EBI per-protein file, User-Agent set, polite
rate limit, cached), keep models with mean pLDDT >= --min-plddt (default 70), run
`foldseek easy-search --exhaustive-search 1 -e 1e-3` against the reference DB of the *matching
transition's P-side species* (runs/structure/refdb/<transition_id>/db, built by build_refdb.py), and
call a query 'structurally old' if it has >=1 hit with E <= 1e-3 AND TM-score >= 0.5.

Input TSV (header required; tab-separated):
  accession      UniProt accession (aliases: uniprot, query, id)
  transition_id  e.g. T_metazoa (selects the reference DB)           [required unless --db is given]
  set            optional label, e.g. gain / control (default 'gain'); fractions are reported per (transition_id,set)
  sequence       optional; if given it is compared with the AFDB model's SEQRES (mismatch flagged)
  (any other columns are carried through to the output)

Outputs  <out>.queries.tsv  (one row per input row) and <out>.summary.tsv (per transition_id,set).

TM-score: foldseek alignment type 2 (3Di+AA) gives E-values; qtmscore/ttmscore/alntmscore are computed
from that alignment.  --tm-metric chooses which one is thresholded (default 'q' = normalised by the QUERY
length: "does the query's fold occur in the reference proteome"; 'max'/'min' of q,t; 'aln' = alntmscore).
All of q/t/aln of the best hit are always reported.  Best hit = highest TM metric among E<=cutoff hits
(ties: bits).  Offline use: pre-place models as <model-dir>/AF-<acc>-F1-model_v<N>.pdb.gz (they are not re-downloaded).
"""
import os, sys, csv, gzip, glob, argparse, subprocess, tempfile, shutil, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import afdb_fetch as A

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
       'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
FMT = 'query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,qcov,tcov,fident,alnlen'
ACC_COLS = ('accession', 'uniprot', 'query', 'id')


def model_seq(path):
    seq = []; seen = set()
    with gzip.open(path, 'rt') as f:
        for l in f:
            if l.startswith('ATOM') and l[12:16] == ' CA ':
                key = (l[21], l[22:27])
                if key not in seen:
                    seen.add(key); seq.append(AA3.get(l[17:20], 'X'))
    return ''.join(seq)


def tm_value(row, metric):
    q, t, aln = float(row['qtm']), float(row['ttm']), float(row['alntm'])
    return {'q': q, 't': t, 'max': max(q, t), 'min': min(q, t), 'aln': aln}[metric]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('queries_tsv'); ap.add_argument('out_prefix')
    ap.add_argument('--root', default=os.environ.get('P', '/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'))
    ap.add_argument('--db', help='override: a single Foldseek target DB for all queries (e.g. for tests)')
    ap.add_argument('--model-dir', help='where query models are cached (default <root>/runs/structure/query_models)')
    ap.add_argument('--min-plddt', type=float, default=70.0)
    ap.add_argument('--evalue', type=float, default=1e-3)
    ap.add_argument('--tm', type=float, default=0.5)
    ap.add_argument('--tm-metric', default='q', choices=['q', 't', 'max', 'min', 'aln'])
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--no-nice', action='store_true')
    ap.add_argument('--no-download', action='store_true', help='only use cached models')
    ap.add_argument('--keep-tmp', action='store_true')
    ap.add_argument('--foldseek', default=None)
    a = ap.parse_args()
    if not a.no_nice: os.nice(19)
    a.threads = min(a.threads, 4)
    FS = a.foldseek or os.environ.get('FOLDSEEK') or f'{a.root}/conda_envs/ancestra/bin/foldseek'
    mdir = a.model_dir or f'{a.root}/runs/structure/query_models'
    os.makedirs(mdir, exist_ok=True)

    rd = csv.DictReader(open(a.queries_tsv), delimiter='\t')
    accol = next((c for c in rd.fieldnames if c.lower() in ACC_COLS), None)
    if accol is None: sys.exit(f'no accession column (one of {ACC_COLS}) in {a.queries_tsv}')
    if not a.db and 'transition_id' not in rd.fieldnames: sys.exit('need transition_id column or --db')
    rows = list(rd)
    for r in rows:
        r['_acc'] = r[accol].strip()
        r['_tid'] = r.get('transition_id', 'DB') if not a.db else 'DB'
        r['_set'] = (r.get('set') or 'gain')
    # ---- fetch models + pLDDT filter
    accs = sorted({r['_acc'] for r in rows})
    def cached(acc):
        g = sorted(glob.glob(f'{mdir}/AF-{acc}-F1-model_v*.pdb.gz'))
        return g[-1] if g else None
    need = [x for x in accs if not cached(x)]
    status = {}
    if need and not a.no_download:
        print(f'[models] downloading {len(need)} of {len(accs)} models', file=sys.stderr)
        res = A.fetch_many(need, mdir, workers=8)
        for k, (st, p, v) in res.items(): status[k] = st
    info = {}
    for x in accs:
        p = cached(x)
        if p is None:
            info[x] = dict(path=None, plddt=None, st='no_model' if status.get(x) in (None, 'missing') else status[x])
        else:
            info[x] = dict(path=p, plddt=A.mean_plddt_from_pdb_gz(p), st='ok')
    # ---- search per transition
    hits = {}   # (tid, acc) -> best hit dict
    nh = {}
    groups = {}
    for r in rows:
        i = info[r['_acc']]
        if i['path'] and i['plddt'] >= a.min_plddt:
            groups.setdefault(r['_tid'], set()).add(r['_acc'])
    tmp = tempfile.mkdtemp(prefix='ngs_', dir=os.environ.get('TMPDIR'))
    for tid, qs in groups.items():
        db = a.db or f'{a.root}/runs/structure/refdb/{tid}/db'
        if not os.path.exists(db + '.dbtype'):
            print(f'[warn] no reference DB for {tid}: {db}', file=sys.stderr)
            for x in qs: hits[(tid, x)] = 'noDB'
            continue
        qd = f'{tmp}/q_{tid}'; os.makedirs(qd)
        for x in qs: os.symlink(os.path.abspath(info[x]['path']), f'{qd}/{x}.pdb.gz')
        out = f'{tmp}/hits_{tid}.tsv'
        cmd = [FS, 'easy-search', qd, db, out, f'{tmp}/ft_{tid}', '--exhaustive-search', '1', '-e', str(a.evalue),
               '--max-seqs', '1000', '--threads', str(a.threads), '-v', '1', '--format-output', FMT]
        if not a.no_nice: cmd = ['nice', '-n', '19'] + cmd
        print('[foldseek]', tid, len(qs), 'queries', file=sys.stderr, flush=True)
        subprocess.run(cmd, check=True)
        cols = FMT.split(',')
        for l in open(out):
            f = dict(zip(cols, l.rstrip('\n').split('\t')))
            q = re.sub(r'^AF-|-F1-model_v\d+$', '', f['query']); t = re.sub(r'^AF-|-F1-model_v\d+$', '', f['target'])
            h = dict(target=t, evalue=float(f['evalue']), bits=float(f['bits']), alntm=f['alntmscore'], qtm=f['qtmscore'],
                     ttm=f['ttmscore'], lddt=f['lddt'], qcov=f['qcov'], tcov=f['tcov'], fident=f['fident'])
            if h['evalue'] > a.evalue: continue
            nh[(tid, q)] = nh.get((tid, q), 0) + 1
            h['tm'] = tm_value(h, a.tm_metric)
            b = hits.get((tid, q))
            if not isinstance(b, dict) or (h['tm'], h['bits']) > (b['tm'], b['bits']): hits[(tid, q)] = h
        for x in qs: hits.setdefault((tid, x), None)
    if not a.keep_tmp: shutil.rmtree(tmp, ignore_errors=True)
    # ---- write outputs
    extra = [c for c in rd.fieldnames if c not in ('transition_id', 'set', 'sequence', accol)]
    hdr = ['accession', 'transition_id', 'set', *extra, 'model_status', 'mean_plddt', 'seq_matches_model', 'n_hits_E',
           'best_target', 'best_evalue', 'best_bits', 'best_qtm', 'best_ttm', 'best_alntm', 'best_lddt', 'best_qcov',
           'best_tcov', 'best_fident', 'structurally_old']
    summ = {}
    with open(a.out_prefix + '.queries.tsv', 'w') as w:
        w.write('\t'.join(hdr) + '\n')
        for r in rows:
            x, tid = r['_acc'], r['_tid']; i = info[x]
            seqok = ''
            if i['path'] and r.get('sequence'):
                seqok = int(model_seq(i['path']) == r['sequence'].strip())
            if i['path'] is None: ms, so = i['st'], ''
            elif i['plddt'] < a.min_plddt: ms, so = 'low_plddt', ''
            else:
                h = hits.get((tid, x))
                if h == 'noDB': ms, so = 'no_ref_db', ''
                else: ms, so = 'searched', ''
            h = hits.get((tid, x)) if ms == 'searched' else None
            if ms == 'searched':
                so = int(isinstance(h, dict) and h['tm'] >= a.tm)
            vals = ['', '', '', '', '', '', '', '', '', '']
            if isinstance(h, dict):
                vals = [h['target'], f"{h['evalue']:.2e}", h['bits'], h['qtm'], h['ttm'], h['alntm'], h['lddt'], h['qcov'], h['tcov'], h['fident']]
            w.write('\t'.join(str(v) for v in [x, tid, r['_set'], *[r.get(c, '') for c in extra], ms,
                    '' if i['plddt'] is None else f"{i['plddt']:.1f}", seqok, nh.get((tid, x), 0), *vals, so]) + '\n')
            s = summ.setdefault((tid, r['_set']), dict(n=0, model=0, plddt=0, searched=0, old=0))
            s['n'] += 1; s['model'] += i['path'] is not None; s['plddt'] += ms in ('searched', 'no_ref_db')
            if ms == 'searched': s['searched'] += 1; s['old'] += so
    with open(a.out_prefix + '.summary.tsv', 'w') as w:
        w.write('transition_id\tset\tn_queries\tn_with_model\tn_plddt_pass\tn_searched\tn_structurally_old\tfrac_old_of_searched\tfrac_old_of_all_queries\n')
        for (tid, st), s in sorted(summ.items()):
            fr = s['old'] / s['searched'] if s['searched'] else float('nan')
            w.write(f"{tid}\t{st}\t{s['n']}\t{s['model']}\t{s['plddt']}\t{s['searched']}\t{s['old']}\t{fr:.4f}\t{s['old']/s['n']:.4f}\n")
    print(open(a.out_prefix + '.summary.tsv').read())


if __name__ == '__main__':
    main()
