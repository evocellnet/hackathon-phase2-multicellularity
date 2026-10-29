#!/usr/bin/env python3
"""Candidate discovery for a dense panel: for each query taxon list the proteomes available in
 (1) UniProt reference proteomes 2026_03 (data/raw/uniprot/refprot_{euk,bac}.tsv, lineage match; canonical
     counts from STATS_2026_03.txt),
 (2) NCBI Datasets v2 annotated assemblies (best per species: RefSeq first, then most protein-coding genes),
 (3) LukProt v1.5.1 / EukProt v3 (metadata sheet; Taxonomy_LukProt/UniEuk match).
Usage: discover.py <panel> Taxon1 Taxon2 ...   -> data/dense/<panel>/candidates.tsv (appends/replaces rows of these taxa)
"""
import csv, os, sys, time, urllib.parse as up
import requests
P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
URL = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/dataset_report'  # v2 report endpoint returned HTTP 500 on 2026-10-29
COLS = ['query', 'db', 'acc', 'species', 'taxid', 'data_type', 'n_prot', 'busco', 'extra']


def ncbi(t):
    reps, tok = [], None
    while True:
        p = {'page_size': 1000, 'filters.has_annotation': 'true'}
        if tok:
            p['page_token'] = tok
        for i in range(4):
            try:
                d = requests.get(URL.format(up.quote(t)), params=p, timeout=180).json(); break
            except Exception:
                time.sleep(5)
        reps += d.get('reports', []); tok = d.get('next_page_token')
        if not tok or len(reps) > 20000:
            break
    best = {}
    for r in reps:
        sp = r['organism']['organism_name']; a = r.get('annotation_info', {})
        pc = a.get('stats', {}).get('gene_counts', {}).get('protein_coding', 0) or 0
        key = (r['accession'].startswith('GCF'), pc > 1000, a.get('busco', {}).get('complete', 0) or 0, pc)
        if sp not in best or key > best[sp][0]:
            best[sp] = (key, r)
    out = []
    for sp, (k, r) in sorted(best.items()):
        a = r.get('annotation_info', {}); b = a.get('busco', {})
        out.append(dict(query=t, db='NCBI', acc=r['accession'], species=sp, taxid=r['organism']['tax_id'],
                        data_type='genome', n_prot=k[3], busco=b.get('complete', ''),
                        extra=f"{a.get('provider','')}|{a.get('release_date','')}|{r.get('assembly_info',{}).get('assembly_level','')}|{b.get('busco_lineage','')}"))
    return out


def uniprot(t):
    stats = {}
    for line in open(f'{P}/data/raw/uniprot/STATS_2026_03.txt'):
        if line.startswith('UP'):
            f = line.rstrip('\n').split('\t'); stats[f[0]] = f
    out = []
    for fn in ('refprot_euk.tsv', 'refprot_bac.tsv'):
        for r in csv.DictReader(open(f'{P}/data/raw/uniprot/{fn}'), delimiter='\t'):
            lin = [x.strip() for x in r['Taxonomic lineage'].split(',')]
            if t in lin or r['Organism'].startswith(t):
                s = stats.get(r['Proteome Id'])
                out.append(dict(query=t, db='UniProt', acc=r['Proteome Id'], species=r['Organism'], taxid=r['Organism Id'],
                                data_type='genome', n_prot=(s[2].split()[0] if s else r['Protein count']),
                                busco=r['BUSCO'], extra=f"{r['Genome assembly ID']}|{r['Genome representation']}|{', '.join(lin[-4:-1])}"))
    return out


def lukprot(t):
    out = []
    for r in csv.DictReader(open(f'{P}/data/raw/lukprot/LukProt_v1.5.1_metadata.tsv'), delimiter='\t'):
        tx = (r.get('Taxonomy_LukProt') or '') + ';' + (r.get('Taxonomy_UniEuk_(only_for_EukProt)') or '')
        toks = [x.strip("' ") for x in tx.replace(',', ';').split(';')]
        if (t in toks or r['Genus_and_species'].startswith(t)) and r['LukProt_v1.5.1_status'] == 'included':
            lid = r['EukProt_ID']
            out.append(dict(query=t, db='LukProt' if lid.startswith('LP') else 'EukProt', acc=lid, species=r['Genus,_species_and_strain'],
                            taxid=r['NCBI_taxid'], data_type=r['Data_Source_Type'], n_prot=r['num_seqs'],
                            busco=r.get('BUSCO_Complete', ''), extra=f"{r['Data_Source_Name']}|{r.get('caution','')}"))
    return out


def main():
    panel, taxa = sys.argv[1], sys.argv[2:]
    d = f'{P}/data/dense/{panel}'; os.makedirs(d, exist_ok=True)
    fn = f'{d}/candidates.tsv'
    old = [r for r in csv.DictReader(open(fn), delimiter='\t')] if os.path.exists(fn) else []
    old = [r for r in old if r['query'] not in taxa]
    new = []
    for t in taxa:
        rows = uniprot(t) + ncbi(t) + lukprot(t)
        print(f'{t}: {len(rows)} candidates', file=sys.stderr)
        new += rows
    with open(fn, 'w') as fh:
        w = csv.DictWriter(fh, COLS, delimiter='\t'); w.writeheader()
        for r in old + new:
            w.writerow({k: r.get(k, '') for k in COLS})


if __name__ == '__main__':
    main()
