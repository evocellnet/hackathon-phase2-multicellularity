#!/usr/bin/env python3
"""Merge proteome metadata from UniProt reference proteomes (REST tsv), EukProt v3 (figshare + GitHub BUSCO)
and LukProt v1.5.1 (Zenodo metadata sheet) into one searchable index: data/raw/source_index.tsv
Columns: db, id, name, taxid, data_type, busco_C, n_proteins, lineage, note"""
import csv, re, os, sys
R = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'raw')
out = []
def busco_c(s):
    m = re.search(r'C:([\d.]+)%', s or '');  return m.group(1) if m else ''
# UniProt
for div in ('euk', 'bac', 'arc'):
    with open(f'{R}/uniprot/refprot_{div}.tsv') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            busco = r['BUSCO']
            lin = r['Taxonomic lineage'].replace('cellular organisms, ', '')
            out.append(['UniProt', r['Proteome Id'], r['Organism'], r['Organism Id'], 'genome', busco_c(busco),
                        r['Protein count'], lin, f"BUSCO={busco};CPD={r['CPD']};asm={r['Genome assembly ID']};repr={r['Genome representation']}"])
# EukProt
busco = {}
with open(f'{R}/eukprot/BUSCO_v5.2.2_eukaryota_odb10.EukProt_v03.2021_11_22.txt') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        busco[r['Input_file'].split('_')[0]] = r
tcs = set(l.split('_')[0] for l in open(f'{R}/eukprot/The_Comparative_Set_196_species.EukProt_v03.2021_11_22.txt') if l.strip())
lk = {}
with open(f'{R}/lukprot/LukProt_v1.5.1_metadata.tsv') as fh:
    rows = list(csv.reader(fh, delimiter='\t'))
hdr = rows[0]
for v in rows[1:]:
    v = v + [''] * (len(hdr) - len(v)); lk[v[0]] = dict(zip(hdr, v))
with open(f'{R}/eukprot/EukProt_included_data_sets.v03.2021_11_22.txt') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        e = r['EukProt_ID']; b = busco.get(e, {}); l = lk.get(e, {})
        out.append(['EukProt', e, r['Name_to_Use'].replace('_', ' ') + (f" ({r['Strain']})" if r['Strain'] else ''),
                    l.get('NCBI_taxid', ''), r['Data_Source_Type'], b.get('Complete', ''), l.get('num_seqs', ''),
                    r['Supergroup_UniEuk'] + '; ' + r['Taxogroup1_UniEuk'] + '; ' + r['Taxogroup2_UniEuk'],
                    f"TCS={'yes' if e in tcs else 'no'};src={r['Data_Source_Name']};action={r['Actions_Prior_to_Use']};"
                    f"OMArk_contam%={l.get('OMArk_Likely_Contamination_%','')}"])
# LukProt-only datasets (not in EukProt v3)
ep_ids = set(x[1] for x in out if x[0] == 'EukProt')
for e, l in lk.items():
    if e in ep_ids or l.get('LukProt_v1.5.1_status', '') != 'included':
        continue
    out.append(['LukProt', e, l['Genus,_species_and_strain'], l.get('NCBI_taxid', ''), l.get('Data_Source_Type', ''),
                l.get('BUSCO_Complete', '').rstrip('%'), l.get('num_seqs', ''),
                l.get('Supergroup_UniEuk', '') + '; ' + l.get('Taxogroup1_UniEuk', '') + '; ' + l.get('Taxonomy_LukProt', '')[-120:],
                f"src={l.get('Source','')};OMArk_contam%={l.get('OMArk_Likely_Contamination_%','')}"])
with open(f'{R}/source_index.tsv', 'w') as fh:
    fh.write('db\tid\tname\ttaxid\tdata_type\tbusco_C\tn_proteins\tlineage\tnote\n')
    for x in out:
        fh.write('\t'.join(str(y).replace('\t', ' ') for y in x) + '\n')
from collections import Counter
print(Counter(x[0] for x in out))
