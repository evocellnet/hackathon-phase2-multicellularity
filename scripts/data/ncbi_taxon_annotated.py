#!/usr/bin/env python3
"""List NCBI (Datasets v2) annotated genome assemblies (one per species, RefSeq preferred, then most protein-coding genes)
for each higher taxon given. Usage: ncbi_taxon_annotated.py Taxon1 Taxon2 ... > out.tsv"""
import sys, requests, urllib.parse as up
URL = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{}/dataset_report'
print('query_taxon\taccession\torganism\ttaxid\tprovider\tprotein_coding\tassembly_level\tbusco_annot')
for t in sys.argv[1:]:
    reps, tok = [], None
    while True:
        p = {'page_size': 500, 'filters.has_annotation': 'true'}
        if tok: p['page_token'] = tok
        d = requests.get(URL.format(up.quote(t)), params=p, timeout=120).json()
        reps += d.get('reports', []); tok = d.get('next_page_token')
        if not tok or len(reps) > 3000: break
    best = {}
    for r in reps:
        sp = r['organism']['organism_name']; a = r.get('annotation_info', {})
        pc = a.get('stats', {}).get('gene_counts', {}).get('protein_coding', 0) or 0
        key = (r['accession'].startswith('GCF'), pc)
        if sp not in best or key > best[sp][0]: best[sp] = (key, r)
    for sp, (k, r) in sorted(best.items()):
        a = r.get('annotation_info', {}); b = a.get('busco', {})
        print(f"{t}\t{r['accession']}\t{sp}\t{r['organism']['tax_id']}\t{a.get('provider','')}\t{k[1]}\t"
              f"{r.get('assembly_info',{}).get('assembly_level','')}\t{b.get('complete','')}", flush=True)
