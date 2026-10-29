#!/usr/bin/env python3
"""Write data/dense/R4/ANNOTATION_NOTES.txt: de novo (AUGUSTUS) proteomes in R4 vs annotated volvocine relatives.
Usage: annotation_notes.py [--runtime-note TEXT]"""
import argparse, glob, os, re
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
REL = ['Gonipec', 'Tetrsoc', 'Astrgub', 'Pleosta', 'Volvcar', 'Volvafr', 'Volvret', 'Chlarei']


def gff_stats(acc):
    g = glob.glob(f'{P}/data/raw/augustus/{acc}/chunks/chunk*.gff3')
    genes = cds = 0
    for f in g:
        for line in open(f):
            if '\tgene\t' in line:
                genes += 1
            elif '\tCDS\t' in line:
                cds += 1
    return genes, cds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runtime-note', default='')
    a = ap.parse_args()
    sp = pd.read_csv(f'{P}/data/dense/R4/species.tsv', sep='\t', dtype=str).fillna('')
    sp['n_proteins'] = sp['n_proteins'].astype(int)
    sp['mean_len'] = (sp['n_residues'].astype(int) / sp['n_proteins']).round(0)
    dn = sp[sp['annotation_source'] == 'augustus_denovo']
    L = ['R4 ANNOTATION NOTES - de novo gene prediction for unannotated volvocine assemblies', '',
         'Why: Yamagishiella unicocca and Eudorina sp. (undifferentiated Volvocaceae; needed as simple-state relatives of the '
         'differentiated Pleodorina/Volvox events) have public assemblies (Volvocales Genome Consortium) but no published gene '
         'annotation in NCBI, UniProt 2026_03, EukProt v3 or LukProt v1.5.1 (checked 2026-10-29).', '',
         'Method (scripts/dense/augustus_predict.sh): AUGUSTUS 3.5.0 (conda_envs/augustus, bioconda; libboost 1.85) ab initio, '
         'species model chlamy2011 (Chlamydomonas reinhardtii, shipped with AUGUSTUS), no hints/RNA-seq/protein evidence, '
         '--gff3=on --protein=on --uniqueGeneId=true --alternatives-from-sampling=false (one transcript per gene), both strands, '
         'partial genes allowed (default). Assemblies split into 7 length-balanced chunks run in parallel (nice 10). '
         'Proteins taken from the GFF3 protein comments, then the round-1 cleaning (terminal stop trimmed, non-standard -> X, '
         '>= 50 aa, exact duplicates removed, IDs <code>_<chunk>_<gene.t1>).', '']
    if a.runtime_note:
        L += ['Runtime: ' + a.runtime_note, '']
    L += ['Caveats: ab initio predictions without evidence miss/merge/split genes and miss lineage-specific genes poorly '
          'captured by the Chlamydomonas model; gene-content differences involving these two proteomes must be checked with a '
          'sensitivity run without them (species.tsv column annotation_source == augustus_denovo).', '',
          'Comparison (BUSCO 6.1.0 eukaryota_odb10, protein mode):', '',
          f"{'code':10s} {'species':34s} {'source':26s} {'annot':16s} {'n_prot':>7s} {'meanlen':>7s}  C     S     D    F    M"]
    for _, r in pd.concat([dn, sp[sp['code'].isin(REL)]]).iterrows():
        L.append(f"{r['code']:10s} {r['species'][:34]:34s} {r['source'][:26]:26s} {r['annotation_source']:16s} "
                 f"{r['n_proteins']:7d} {r['mean_len']:7.0f}  {r['busco_C']:5s} {r['busco_S']:5s} {r['busco_D']:4s} {r['busco_F']:4s} {r['busco_M']}")
    L.append('')
    for _, r in dn.iterrows():
        acc = r['source'].split(':')[1]
        g, c = gff_stats(acc)
        L.append(f"{r['code']}: assembly {acc}; AUGUSTUS genes={g}, CDS features={c}, proteins after cleaning={r['n_proteins']}")
    open(f'{P}/data/dense/R4/ANNOTATION_NOTES.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
