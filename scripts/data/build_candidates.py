#!/usr/bin/env python3
"""Resolve curated candidate species lists against proteome sources and write candidate tables.

Inputs : scripts/data/candidates_{euk,prok}_spec.txt (curated: species, transition, role, ...)
         data/raw/source_index.tsv (from build_source_index.py: UniProt ref proteomes, EukProt v3, LukProt v1.5.1)
         data/raw/uniprot/STATS_2026_03.txt (canonical FASTA entry counts per UniProt proteome)
         NCBI Datasets v2 API (live) for NCBI-only taxa
Outputs: data/species_candidates.tsv, data/species_candidates_prok.tsv

Source preference (auto): UniProt reference proteome (BUSCO C >= 60 or no better option)
  > EukProt/LukProt genome > EukProt/LukProt transcriptome; explicit overrides in the spec win.
Usage: build_candidates.py [euk|prok|both]
"""
import csv, os, re, sys, json, time, urllib.parse as up
import requests

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
RAW = os.path.join(ROOT, 'data', 'raw')
LIC = {'UniProt': 'CC BY 4.0 (UniProt)', 'EukProt': 'CC BY 4.0 (EukProt v3; original data-source terms apply)',
       'LukProt': 'CC BY 4.0 (LukProt v1.5.1; original data-source terms apply)',
       'NCBI': 'INSDC/NCBI: no use restrictions (cite submitter)', 'NONE': ''}
COLS = ['species', 'ncbi_taxid', 'lineage_group', 'transition_id', 'role', 'multicellularity_type', 'source',
        'data_type', 'completeness', 'n_proteins', 'license', 'notes',
        'tier', 'source_db', 'busco_detail', 'alt_sources', 'source_name', 'taxonomy']

def fnum(x):
    try: return float(str(x).rstrip('%'))
    except Exception: return -1.0

def load_index():
    rows = list(csv.DictReader(open(os.path.join(RAW, 'source_index.tsv')), delimiter='\t'))
    canon = {}
    for l in open(os.path.join(RAW, 'uniprot', 'STATS_2026_03.txt')):
        if l.startswith('UP'):
            f = l.split('\t'); m = re.match(r'\s*(\d+)', f[2]); canon[f[0]] = int(m.group(1)) if m else ''
    return rows, canon

def ncbi_report(acc=None, taxon=None):
    if acc:
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report'
        params = {}
    else:
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{up.quote(taxon)}/dataset_report'
        params = {'filters.has_annotation': 'true', 'page_size': 50}
    for _ in range(3):
        try:
            d = requests.get(url, params=params, timeout=60).json(); break
        except Exception:
            time.sleep(2)
    else:
        return None
    reps = [r for r in d.get('reports', []) if r.get('annotation_info')]
    if not reps: return None
    reps.sort(key=lambda r: (r['accession'].startswith('GCF'),
                             r['annotation_info'].get('stats', {}).get('gene_counts', {}).get('protein_coding', 0) or 0), reverse=True)
    r = reps[0]; a = r['annotation_info']
    return {'acc': r['accession'], 'name': r['organism']['organism_name'], 'taxid': str(r['organism']['tax_id']),
            'pc': a.get('stats', {}).get('gene_counts', {}).get('protein_coding', ''),
            'busco': a.get('busco', {}), 'provider': a.get('provider', ''),
            'level': r.get('assembly_info', {}).get('assembly_level', '')}

def resolve(q, override, rows, canon, prok):
    ql = q.lower()
    def hits(db):
        out = [r for r in rows if r['db'] == db and r['name'].lower().startswith(ql)]
        if db == 'UniProt':
            out = [r for r in out if (('Bacteria' in r['lineage'] or 'Archaea' in r['lineage']) if prok else ('Eukaryota' in r['lineage']))]
        return out
    cand = {db: hits(db) for db in ('UniProt', 'EukProt', 'LukProt')}
    def best(lst, genome_first=True):
        if not lst: return None
        return sorted(lst, key=lambda r: ((r['data_type'] == 'genome') if genome_first else 0, fnum(r['busco_C']),
                                          fnum(r['n_proteins'])), reverse=True)[0]
    def alt_str(chosen_id):
        alts = []
        for db, lst in cand.items():
            b = best(lst)
            if b and b['id'] != chosen_id:
                alts.append(f"{db}:{b['id']}({b['data_type'][:5]},B{b['busco_C']})")
        return ';'.join(alts)
    # explicit override
    if override and override not in ('AUTO',):
        if override == 'NONE':
            return {'source_db': 'NONE', 'source': 'NONE', 'data_type': 'none', 'completeness': '', 'n_proteins': '',
                    'busco_detail': '', 'taxid': '', 'name': q, 'alt': alt_str(None), 'tax': ''}
        db, acc = override.split(':', 1)
        if db == 'NCBI':
            rep = ncbi_report(taxon=q) if acc == 'AUTO' else ncbi_report(acc=acc)
            if not rep:
                return {'source_db': 'NONE', 'source': f'NCBI:{acc} (not found)', 'data_type': 'none', 'completeness': '',
                        'n_proteins': '', 'busco_detail': '', 'taxid': '', 'name': q, 'alt': '', 'tax': ''}
            b = rep['busco']
            return {'source_db': 'NCBI', 'source': rep['acc'], 'data_type': 'genome', 'completeness': b.get('complete', '') and round(100 * float(b['complete']), 1),
                    'n_proteins': rep['pc'], 'busco_detail': (f"NCBI annot BUSCO {b.get('lineage','')} C={b.get('complete','')}" if b else 'no BUSCO in NCBI report') + f"; provider={rep['provider']}; level={rep['level']}",
                    'taxid': rep['taxid'], 'name': rep['name'], 'alt': alt_str(None), 'tax': ''}
        dbn = {'UP': 'UniProt', 'EP': 'EukProt', 'LP': 'LukProt'}[db]
        r = [x for x in rows if x['db'] == dbn and x['id'] == acc]
        chosen = r[0] if r else None
    else:
        u = best(cand['UniProt'])
        g = best([x for x in cand['EukProt'] + cand['LukProt'] if x['data_type'] == 'genome'])
        t = best(cand['EukProt'] + cand['LukProt'])
        if u and (fnum(u['busco_C']) >= 60 or not g or fnum(g['busco_C']) <= fnum(u['busco_C']) + 10):
            chosen = u
        else:
            chosen = g or t or u
    if not chosen:
        if prok:
            rep = ncbi_report(taxon=q)
            if rep:
                return resolve(q, f"NCBI:{rep['acc']}", rows, canon, prok)
        return {'source_db': 'NONE', 'source': 'NOT FOUND', 'data_type': 'none', 'completeness': '', 'n_proteins': '',
                'busco_detail': '', 'taxid': '', 'name': q, 'alt': '', 'tax': ''}
    db = chosen['db']
    n = chosen['n_proteins']
    detail = ''
    if db == 'UniProt':
        n = canon.get(chosen['id'], n)
        m = re.search(r'BUSCO=([^;]*)', chosen['note']); detail = 'UniProt ' + (m.group(1) if m else '')
        m = re.search(r'CPD=([^;]*)', chosen['note']); detail += f"; CPD={m.group(1) if m else ''}"
    elif db == 'EukProt':
        detail = f"EukProt BUSCO v5 eukaryota_odb10 C={chosen['busco_C']}; " + chosen['note']
    else:
        detail = f"LukProt BUSCO C={chosen['busco_C']}; " + chosen['note']
    return {'source_db': db, 'source': chosen['id'], 'data_type': chosen['data_type'], 'completeness': chosen['busco_C'],
            'n_proteins': n, 'busco_detail': detail, 'taxid': chosen['taxid'], 'name': chosen['name'],
            'alt': alt_str(chosen['id']), 'tax': chosen['lineage']}

def build(kind):
    prok = kind == 'prok'
    rows, canon = load_index()
    spec = os.path.join(ROOT, 'scripts', 'data', f'candidates_{kind}_spec.txt')
    out = []
    for l in open(spec):
        if not l.strip() or l.startswith('#'): continue
        f = [x.strip() for x in l.rstrip('\n').split('|')]
        f += [''] * (8 - len(f))
        q, ov, lg, tr, role, mct, tier, notes = f[:8]
        r = resolve(q, ov, rows, canon, prok)
        sp = re.sub(r'\s*\(.*$', '', r['name']).strip() if r['source_db'] != 'NONE' else q
        dt = r['data_type']
        flag = []
        if 'transcriptome' in dt or dt in ('EST',): flag.append('TRANSCRIPTOME-ONLY')
        if r['source_db'] == 'NONE': flag.append('NO-PROTEOME-AVAILABLE')
        c = fnum(r['completeness'])
        if 0 <= c < 60: flag.append('LOW-COMPLETENESS')
        out.append({'species': sp, 'ncbi_taxid': r['taxid'], 'lineage_group': lg, 'transition_id': tr, 'role': role,
                    'multicellularity_type': mct, 'source': (r['source_db'] + ':' + r['source']) if r['source_db'] != 'NONE' else r['source'],
                    'data_type': dt, 'completeness': r['completeness'], 'n_proteins': r['n_proteins'],
                    'license': LIC.get(r['source_db'], ''), 'notes': ' '.join(f'[{x}]' for x in flag) + (' ' if flag else '') + notes,
                    'tier': tier, 'source_db': r['source_db'], 'busco_detail': r['busco_detail'], 'alt_sources': r['alt'],
                    'source_name': r['name'], 'taxonomy': r['tax']})
        print(f"{tier[:4]}\t{q[:35]:35s}\t{out[-1]['source']}\t{dt[:6]}\tB{r['completeness']}\tn{r['n_proteins']}\t{r['name'][:40]}", file=sys.stderr)
    fn = os.path.join(ROOT, 'data', 'species_candidates_prok.tsv' if prok else 'species_candidates.tsv')
    with open(fn, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(out)
    print(f'wrote {fn}: {len(out)} rows', file=sys.stderr)

if __name__ == '__main__':
    k = sys.argv[1] if len(sys.argv) > 1 else 'both'
    for kind in (['euk', 'prok'] if k == 'both' else [k]):
        build(kind)
