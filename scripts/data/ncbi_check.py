#!/usr/bin/env python3
"""For each taxon name: (1) NCBI Datasets v2 annotated genomes (RefSeq/GenBank) with protein-coding gene counts,
(2) counts of TSA (transcriptome shotgun assembly) and SRA RNA-seq records via E-utilities.
Usage: ncbi_check.py "Genus species" ...   (prints TSV)"""
import sys, time, requests, urllib.parse as up
DS = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{}/dataset_report'
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
def genomes(name):
    try:
        r = requests.get(DS.format(up.quote(name)), params={'page_size': 20}, timeout=60).json()
    except Exception as e:
        return f'ERR {e}', ''
    reps = r.get('reports', [])
    ann = [x for x in reps if x.get('annotation_info')]
    best = sorted(ann, key=lambda x: (x['accession'].startswith('GCF'),
                  x['annotation_info'].get('stats', {}).get('gene_counts', {}).get('protein_coding', 0)), reverse=True)
    s = ';'.join(f"{x['accession']}:{x['organism']['organism_name']}:pc={x['annotation_info'].get('stats',{}).get('gene_counts',{}).get('protein_coding','?')}:tx={x['organism']['tax_id']}" for x in best[:3])
    return f"{len(reps)} asm/{len(ann)} annot", s
def esearch(db, term):
    time.sleep(0.4)
    try:
        return requests.get(EU, params={'db': db, 'term': term, 'retmode': 'json'}, timeout=60).json()['esearchresult']['count']
    except Exception as e:
        return 'ERR'
print('query\tgenomes\tbest_annotated\tTSA_nuccore\tSRA_RNAseq')
for n in sys.argv[1:]:
    g, s = genomes(n)
    tsa = esearch('nuccore', f'"{n}"[Organism] AND tsa-master[prop]')
    sra = esearch('sra', f'"{n}"[Organism] AND "transcriptomic"[Source]')
    print(f'{n}\t{g}\t{s}\t{tsa}\t{sra}', flush=True)
