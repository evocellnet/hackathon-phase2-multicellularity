#!/usr/bin/env python3
"""Step 1 of the dataset build: final species selection + species codes.

Input : data/species_candidates.tsv, data/species_candidates_prok.tsv (scout tables)
Output: data/species_final.tsv, data/species_final_prok.tsv
        (all candidate columns + selection_reason, code, final_source_file, md5,
         n_proteins_final; the last three are filled by finalize_tables.py)

Rules (DATA_SOURCES_SCOUTING.md section 3, DESIGN.md section 2):
  eukaryotes : all `core` rows that have a data source, plus the optional rows in
               EUK_OPTIONAL_ADD (each with a reason: fill = a transition had <2 taxa on
               one side; depth = more taxa for a transition already covered;
               thicken = under-sampled supergroup background). Target 190-210, cap 220.
  prokaryotes: all `core` rows plus optional rows only where a transition had <2
               taxa on a side (PROK_OPTIONAL_ADD).
  `exclude` rows and rows without data (source NONE) are never selected.
Re-runnable; deterministic.
"""
import os, re, sys
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# species name (as in the candidate table) -> reason
EUK_OPTIONAL_ADD = {
    # --- fill: transition side with <2 taxa among core rows
    'Pigoraptor chileana': 'fill: T_capsaspora_aggregation unicellular side (only Ministeria in core)',
    'Syssomonas multiformis': 'fill: T_ichthyosporea_coenocytic unicellular side (only Corallochytrium in core)',
    'Rhizamoeba saxonica': 'fill: T_copromyxa unicellular side (only Vermamoeba in core)',
    'Arcella intermedia': 'fill: T_copromyxa unicellular side (Tubulinea)',
    'Ceratiomyxa fruticulosa': 'fill: T_myxogastria_plasmodium unicellular sister (none in core)',
    'Bryopsis sp. KO-2023': 'fill: T_siphonous_coenocytic multicellular side (only Caulerpa in core)',
    'Ostreobium quekettii': 'fill: T_siphonous_coenocytic multicellular side',
    'Vaucheria litorea': 'fill: T_xanthophyceae_filaments multicellular side (only Tribonema in core)',
    # --- depth: more taxa for a covered transition
    'Amoebidium parasiticum': 'depth: T_ichthyosporea_coenocytic multicellular side',
    'Codosiga hollandica': 'depth: T_choano_colonies colonial craspedid',
    'Didymoeca costata': 'depth: T_choano_colonies / T_metazoa unicellular acanthoecid',
    'Cavostelium apophysatum': 'depth: T_dictyostelia unicellular protosteloid sister',
    'Stygamoeba regulata': 'depth: T_copromyxa outgroup (Discosea)',
    'Tetraselmis striata': 'depth: T_ulvophyceae chlorophyte backbone (Chlorodendrophyceae)',
    'Volvox africanus': 'depth: T_volvocales multicellular (second Volvox)',
    'Chlamydomonas schloesseri': 'depth: T_volvocales unicellular sister',
    'Rhodelphis marinus': 'depth: T_rhodophyta outgroup (Rhodelphidia, sister of red algae)',
    'Agarophyton chilense': 'depth: T_rhodophyta multicellular Florideophyceae',
    'Coleochaete scutata': 'depth: T_streptophyta grade (Coleochaetophyceae)',
    'Closterium sp. NIES-67': 'depth: T_streptophyta unicellular Zygnematophyceae',
    'Ceratopteris richardii': 'depth: T_streptophyta fern (fills land-plant gap)',
    'Piromyces finnis': 'depth: T_fungi_hyphae unicellular-rhizoidal sister (Neocallimastigomycota)',
    'Basidiobolus meristosporus CBS 931.73': 'depth: T_fungi_hyphae Zoopagomycota',
    'Coemansia umbellata': 'depth: T_fungi_hyphae Zoopagomycota (Kickxellomycotina)',
    'Mortierella hygrophila': 'depth: T_fungi_hyphae Mucoromycota (Mortierellomycotina)',
    'Botrytis cinerea': 'depth: T_fungi_pezizomycotina (Leotiomycetes)',
    'Zymoseptoria brevis': 'depth: T_fungi_pezizomycotina (Dothideomycetes)',
    'Taphrina deformans': 'depth: T_fungi_neolecta unicellular Taphrinomycotina sister',
    'Wallemia mellicola': 'depth: T_fungi_agaricomycetes unicellular sister side',
    'Tremella mesenterica': 'depth: T_fungi_agaricomycetes (Tremellomycetes, dimorphic jelly fungus)',
    'Oscarella pearsei': 'depth: T_metazoa homoscleromorph sponge',
    'Caenorhabditis elegans': 'depth: T_metazoa (Ecdysozoa model)',
    'Pythium oligandrum': 'depth: T_oomycete_hyphae',
    # --- thicken: under-sampled supergroups (control branches)
    'Anaeramoeba flamelloides': 'thicken: Metamonada (Anaeramoebidae)',
    'Monocercomonoides exilis PA203': 'thicken: Metamonada (Preaxostyla, LukProt genome)',
    'Trichomonas vaginalis': 'thicken: Metamonada (Parabasalia)',
    'Willaertia magna': 'thicken: Discoba (Heterolobosea; T_acrasis sister side)',
    'Mantamonas plastica': 'thicken: CRuMs',
    'Emiliania huxleyi': 'thicken: Haptista (Haptophyta)',
    'Palpitomonas bilix': 'thicken: Cryptista',
    'Symbiodinium microadriaticum': 'thicken: Alveolata (Dinoflagellata, none in core)',
    'Chromera velia': 'thicken: Alveolata (Chrompodellida)',
    'Reticulomyxa filosa': 'thicken: Rhizaria (Foraminifera)',
    'Lotharella globosa': 'thicken: Rhizaria (Cercozoa, Chlorarachnea)',
    'Brevimastigomonas motovehiculus': 'thicken: Rhizaria (Cercozoa)',
    'Hoilungia hongkongensis': 'depth: T_metazoa second placozoan (LukProt genome)',
    'Ephydatia muelleri': 'depth: T_metazoa second demosponge genome',
}

# rows not in the scout table (added after design review round 1, evo_biology Should 6;
# only taxa with a ready-made proteome were added - Hyphochytrium, Developayella, Yamagishiella,
# Eudorina, Interfilum: no annotated proteome in UniProt/EukProt/LukProt/NCBI on 2026-10-28)
EUK_EXTRA_ROWS = [
    dict(species='Paraphelidium tribonemae', ncbi_taxid='2548444', lineage_group='Holomycota-Aphelida',
         transition_id='T_fungi_hyphae', role='outgroup', multicellularity_type='',
         source='EukProt:EP00158', data_type='transcriptome', completeness='87.8', n_proteins='10439',
         license='CC BY 4.0 (EukProt v3; original data-source terms apply)',
         notes='[TRANSCRIPTOME-ONLY] aphelid, sister of Fungi (Torruella et al. 2018); added after review round 1',
         tier='optional', source_db='EukProt', busco_detail='EukProt BUSCO v5 eukaryota_odb10 C=87.8',
         alt_sources='', source_name='Paraphelidium tribonemae X-108 (Partr_v1)',
         taxonomy='Opisthokonta; Holomycota; Aphelida; Paraphelidium',
         selection_reason='review1: sister of Fungi; breaks the long fungal stem (Torruella et al. 2018)'),
]

# memory cap (PI rule, 2026-10-28, after design review round 1): eukaryote total <= 3.0 M proteins.
# Drop in order (1) tier=optional AND all roles == background, largest n_proteins_final first; (2) then other
# optional taxa, largest first, but only if every transition side the taxon is on (role multicellular or
# closest_unicellular_sister in the scout table) keeps >= 2 species; never drop the only representative of a side.
# Stop as soon as the total <= CAP. Needs the per-species stats of preprocess_proteomes.py (run select -> fetch ->
# preprocess -> select again); dropped taxa are listed in data/species_cap_drops.tsv.
EUK_PROTEIN_CAP = 3_000_000
EUK_DROP = {}   # no manual drops (the cap rule decides)

# data_type audit: EukProt 'genome' sets that are really transcript-derived
DATA_TYPE_FINAL = {
    'EukProt:EP00667': ('transcriptome', 'EukProt labels it genome, but the set is the translated transcriptome '
                        '(EG_transcript_N) of the genome paper; processed as a transcriptome'),
    'EukProt:EP00117': ('genome', 'Compagen SCIL_T-PEP: genome gene models (scgid gene IDs, 35.8k genes); '
                        'isoforms collapsed by scgid; kept as genome but check BUSCO D/F'),
}

PROK_OPTIONAL_ADD = {
    'Haloarcula marismortui': 'fill: T_haloferax_tissue unicellular side (only Halobacterium in core)',
    'Candidatus Electronema sp. V4': 'fill: T_cable_bacteria multicellular side (only Electrothrix in core)',
    'Desulfotalea psychrophila': 'fill: T_cable_bacteria unicellular side (only Desulfobulbus in core)',
}

USED = set()
SKIP_WORDS = {'candidatus', 'cand.'}


def make_code(species: str, taxid, used: set) -> str:
    toks = [t for t in re.split(r'[\s_]+', species.strip()) if t]
    toks = [t for t in toks if t.lower() not in SKIP_WORDS]
    genus = re.sub(r'[^A-Za-z]', '', toks[0])
    sp = re.sub(r'[^A-Za-z]', '', toks[1]) if len(toks) > 1 else 'sp'
    placeholder = sp.lower() in ('sp', 'spp', 'cf', 'aff', '')
    if placeholder:
        tail = ''.join(re.sub(r'[^A-Za-z0-9]', '', t) for t in toks[2:])[:6]
        code = genus[:4].capitalize() + 'sp' + (tail or str(taxid))
    else:
        code = genus[:4].capitalize() + sp[:3].lower()
        if code in used:  # collision: one more genus letter
            code = genus[:5].capitalize() + sp[:3].lower()
    base, k = code, 1
    while code in used:
        k += 1
        code = f'{base}{k}'
    used.add(code)
    return code


def select(path_in, path_out, add, kind, extra=(), drop=()):
    d = pd.read_csv(path_in, sep='\t', dtype=str)
    if extra:
        d = pd.concat([d, pd.DataFrame(list(extra))], ignore_index=True)
        add = dict(add, **{r['species']: r['selection_reason'] for r in extra})
    has_data = ~d['source'].fillna('NONE').str.upper().str.startswith('NONE')
    core = (d['tier'] == 'core') & has_data
    missing = [s for s in add if s not in set(d['species'])]
    if missing:
        sys.exit(f'[{kind}] optional species not found in candidate table: {missing}')
    bad = d[d['species'].isin(add) & ((d['tier'] != 'optional') | ~has_data)]
    if len(bad):
        sys.exit(f'[{kind}] selected optional rows that are not optional/with data: {bad.species.tolist()}')
    sel = d[(core | d['species'].isin(add)) & ~d['species'].isin(list(drop))].copy()
    sel['selection_reason'] = [('core' if t == 'core' else add[s]) for s, t in zip(sel['species'], sel['tier'])]
    sel['data_type_final'] = [DATA_TYPE_FINAL.get(s, (t, ''))[0] for s, t in zip(sel['source'], sel['data_type'])]
    sel['data_type_note'] = [DATA_TYPE_FINAL.get(s, (t, ''))[1] for s, t in zip(sel['source'], sel['data_type'])]
    used = USED
    sel['code'] = [make_code(s, t, used) for s, t in zip(sel['species'], sel['ncbi_taxid'])]
    for c in ('final_source_file', 'md5', 'n_proteins_final'):
        sel[c] = ''
    sel.to_csv(path_out, sep='\t', index=False)
    print(f'[{kind}] candidates={len(d)} core_with_data={int(core.sum())} '
          f'optional_added={len(add)} dropped={len(drop)} selected={len(sel)} -> {os.path.relpath(path_out, P)}')
    print(sel.groupby(['source_db', 'data_type']).size().to_string())
    return sel


def sides(row):
    """transition sides of a scout row: {(transition, 'mc'|'uni')}"""
    out = set()
    for t, r in zip(str(row['transition_id']).split(';'), str(row['role']).split(';')):
        if r == 'multicellular':
            out.add((t, 'mc'))
        elif r == 'closest_unicellular_sister':
            out.add((t, 'uni'))
    return out


def apply_cap(sel, cap, path_out):
    import json
    stats = {}
    for c in sel['code']:
        f = f'{P}/data/proteomes/stats/{c}.json'
        if not os.path.exists(f):
            print(f'[euk] cap NOT applied: no preprocessing stats for {c} (run preprocess, then select again)')
            return sel
        stats[c] = json.load(open(f))['n_final']
    sel = sel.copy()
    sel['n_tmp'] = sel['code'].map(stats)
    total = int(sel['n_tmp'].sum())
    side_count = {}
    for _, r in sel.iterrows():
        for sd in sides(r):
            side_count[sd] = side_count.get(sd, 0) + 1
    drops = []
    opt = sel[sel['tier'] == 'optional']
    is_bg = opt['role'].fillna('').apply(lambda x: all(y == 'background' for y in x.split(';')))
    step1 = opt[is_bg].sort_values('n_tmp', ascending=False)
    step2 = opt[~is_bg].sort_values('n_tmp', ascending=False)
    for step, cand in ((1, step1), (2, step2)):
        for _, r in cand.iterrows():
            if total <= cap:
                break
            sd = sides(r)
            if step == 2 and any(side_count[x] - 1 < 2 for x in sd):
                continue
            for x in sd:
                side_count[x] -= 1
            total -= int(r['n_tmp'])
            drops.append(dict(species=r['species'], code=r['code'], n_proteins_final=int(r['n_tmp']), rule_step=step,
                              role=r['role'], transition_id=r['transition_id'], total_after=total))
    dd = pd.DataFrame(drops)
    dd.to_csv(f'{P}/data/species_cap_drops.tsv', sep='\t', index=False)
    print(f'[euk] protein cap {cap:,}: dropped {len(dd)} taxa -> total {total:,}')
    if len(dd):
        print(dd.to_string(index=False))
    assert total <= cap, 'cap rule could not reach the cap'
    out = sel[~sel['code'].isin(set(dd.get('code', [])))].drop(columns='n_tmp')
    out.to_csv(path_out, sep='\t', index=False)
    return out


if __name__ == '__main__':
    e = select(f'{P}/data/species_candidates.tsv', f'{P}/data/species_final.tsv', EUK_OPTIONAL_ADD, 'euk',
               EUK_EXTRA_ROWS, EUK_DROP)
    e = apply_cap(e, EUK_PROTEIN_CAP, f'{P}/data/species_final.tsv')
    assert 180 <= len(e) <= 220, len(e)
    p = select(f'{P}/data/species_candidates_prok.tsv', f'{P}/data/species_final_prok.tsv', PROK_OPTIONAL_ADD, 'prok')
    allc = list(e['code']) + list(p['code'])
    assert len(set(allc)) == len(allc), 'duplicate codes across kingdoms'
    assert all(re.fullmatch(r'[A-Za-z0-9]+', c) for c in allc)
