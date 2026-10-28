#!/usr/bin/env python3
"""Step 6: fill the per-species columns of data/species_final.tsv and data/species_final_prok.tsv.

Adds / overwrites: final_source_file (data/proteomes/<kingdom>/<code>.fa), raw_source_file, preprocessing,
md5 (of the final FASTA), n_proteins_raw, n_proteins_final, busco_lineage, busco_C/S/D/F/M, busco_version,
low_quality (BUSCO C < 50), flag_DF15 (transcriptome-like set with BUSCO D or F > 15 %), and the trait columns
from data/species_traits.tsv (mc_state, multicellular, reversal, transition_roles).
Also moves FASTA files of species no longer selected to data/proteomes/_dropped/ (so that the Ancestra input
directories contain exactly the selected species), and writes data/proteomes/summary.tsv.
"""
import hashlib, json, os, shutil
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def main():
    busco = pd.read_csv(f'{P}/data/qc/busco_summary.tsv', sep='\t').set_index('code')
    traits = (pd.read_csv(f'{P}/data/species_traits.tsv', sep='\t', dtype=str).set_index('code')
              if os.path.exists(f'{P}/data/species_traits.tsv') else None)
    summ = []
    for kingdom, f in (('euk', 'species_final.tsv'), ('prok', 'species_final_prok.tsv')):
        path = f'{P}/data/{f}'
        d = pd.read_csv(path, sep='\t', dtype=str)
        keep = set(d['code'])
        os.makedirs(f'{P}/data/proteomes/_dropped', exist_ok=True)
        for fa in os.listdir(f'{P}/data/proteomes/{kingdom}'):
            if fa.endswith('.fa') and fa[:-3] not in keep:
                shutil.move(f'{P}/data/proteomes/{kingdom}/{fa}', f'{P}/data/proteomes/_dropped/{fa}')
                print('moved unselected', fa, 'to data/proteomes/_dropped/')
        cols = {k: [] for k in ('final_source_file', 'raw_source_file', 'preprocessing', 'md5', 'n_proteins_raw',
                                'n_proteins_final', 'busco_lineage', 'busco_C', 'busco_S', 'busco_D', 'busco_F',
                                'busco_M', 'busco_version', 'low_quality', 'flag_DF15')}
        for _, r in d.iterrows():
            c = r['code']
            fa = f'data/proteomes/{kingdom}/{c}.fa'
            st = json.load(open(f'{P}/data/proteomes/stats/{c}.json'))
            h = hashlib.md5(open(f'{P}/{fa}', 'rb').read()).hexdigest()
            assert h == st['md5'], f'{c}: FASTA changed after preprocessing'
            b = busco.loc[c]
            assert b['input_md5'] == h, f'{c}: BUSCO was run on a different FASTA version'
            tl = r.get('data_type_final', r['data_type']) != 'genome' if kingdom == 'euk' else False
            vals = dict(final_source_file=fa, raw_source_file=st['source_file'], preprocessing=st['method'], md5=h,
                        n_proteins_raw=st['n_raw'], n_proteins_final=st['n_final'], busco_lineage=b['lineage'],
                        busco_C=b['C'], busco_S=b['S'], busco_D=b['D'], busco_F=b['F'], busco_M=b['M'],
                        busco_version=b['busco_version'], low_quality=int(b['C'] < 50),
                        flag_DF15=int(bool(tl) and (b['D'] > 15 or b['F'] > 15)))
            for k, v in vals.items():
                cols[k].append(v)
            summ.append(dict(code=c, kingdom=kingdom, species=r['species'], data_type=r.get('data_type_final', r['data_type']),
                             n_proteins_final=st['n_final'], busco_C=b['C'], low_quality=vals['low_quality']))
        for k, v in cols.items():
            d[k] = v
        if traits is not None:
            for k in ('mc_state', 'multicellular', 'multicellular_permissive', 'reversal', 'transition_roles'):
                d[k] = [traits.loc[c, k] if c in traits.index else '' for c in d['code']]
        d.to_csv(path, sep='\t', index=False)
        print(f"{kingdom}: {len(d)} species, {int(d['n_proteins_final'].astype(int).sum())} proteins, "
              f"low_quality={int(d['low_quality'].astype(int).sum())}, flag_DF15={int(d['flag_DF15'].astype(int).sum())}")
    pd.DataFrame(summ).to_csv(f'{P}/data/proteomes/summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
