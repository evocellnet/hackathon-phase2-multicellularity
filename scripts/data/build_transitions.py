#!/usr/bin/env python3
"""Step 5: transition table + per-species traits.

Each transition is a STATE CHANGE (design review round 1, evo_biology Must 1):
  from_state -> to_state, states: unicellular | simple | complex | grey (coenocytic, plasmodial, facultative
  colonies, culture-only aggregation) | reversal (= multicellular -> unicellular).
  M = guide-tree node (or tip) of the derived state; P = parent of M = last ancestor in the PRIOR state.
  analysis_tier: primary | secondary | exploratory | reversal | excluded | dropped_no_data;
  primary_eligible = 1 only for 'primary'.
Inputs : data/guide_tree_{euk,prok}.nwk, data/species_final*.tsv, data/qc/busco_summary.tsv (optional)
Outputs: data/transitions.tsv, data/species_traits.tsv
prior_state_ref_*: the nearest species outside M that are in the prior state (unicellular, non-reversal, for
U->S/X; simple for S->X; multicellular for reversals), found by walking up from M - the comparison set to use when
the tree-sister of M is itself multicellular (e.g. Tetrabaena vs Gonium+Volvocaceae).
Every M/P is derived from the tree; M given as leaf lists is resolved to their MRCA and flagged when the MRCA
contains species of a different state (non-clade).
"""
import os, sys
import pandas as pd
from ete3 import Tree

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# ---------------------------------------------------------------- per-species states (default: unicellular)
S, X, G = 'simple', 'complex', 'grey'
STATE = {}
for c in ('Allomac Rhizmic Mucolus Morthyg Basimer Coemumb Gonapro Ustimay Wallmel Taphdef '
          'Dictdis Dictpur Hetepal Tieglac Cavefas Acytsub Fontalb Acrakon Coprpro Guttvul '
          'Volvcar Volvafr Pleosta Astrgub Gonipec Tetrsoc Ulvamut Klebnit Colescu Pyroyez Porpumb '
          'Schiisc Phaecon Tribmin Ectosil Discmes '
          # prokaryotes
          'Planaga Lyngaes Phornig Micrvag Leptbor Pseubic Myxoxan Coracor Alloaur Cystfus Soracel Choncro Nannexe '
          'Halioch Strecoe Streven Kitaset Franaln Micrech Therfus Nocafar Saccery Magnmul MagnspHK1 Elecsca '
          'ElecspV4 Jeonsac Methace Methmaz Methbar').split():
    STATE[c] = S
for c in ('Amphque Ephymue Sycocil Oscapea Tricadh Hoilhon Mnemlei Nemavec Hydrvul Homosap Branbel Stropur '
          'Drosmel Caenele Capitel Pyroomp Neurcra Aspenid Botrcin Zymobre Coprcin Schicom Laccbic Dacrpri '
          'Tremmes Neolirr Marcpol Physpat Selamoe Ceraric Ambotri Arabtha Charbra Choncri Graccho Agarchi '
          'Sacclat Fucuser Dictdic '
          'Nostpun Nostsp103690 Tricvar Anabcyl Cylirac Noduspu Calopar Scythof Fiscthe Chlofri').split():
    STATE[c] = X
for c in ('Salpros Barrmon Choafle Codohol Capsowc Sphaarc Creofra Chroper Amoepar Coralim Syssmul Physpol '
          'Caullen BryospKO2023 Ostrque Vauclit Phytinf Saprpar Pytholi Plasbra LabyspSRHaC Chloatm Retifil '
          'Halovol Methsoe Methhun').split():
    STATE[c] = G
REVERSAL = set('Sacccer Yarrlip Liposta Schipom Saitcom Crypneo Rhodtor Spirmus Mesoend ClosspNIES67'.split())
STATE_NOTE = {
    'Ustimay': 'dimorphic (yeast + filamentous dikaryon)', 'Taphdef': 'dimorphic', 'Wallmel': 'filamentous (not a yeast)',
    'Yarrlip': 'Saccharomycotina reversal; dimorphic (forms hyphae)', 'Ectosil': 'filamentous brown alga',
    'Discmes': 'filamentous; sister to all other Phaeophyceae', 'Volvcar': 'differentiated (germ-soma)',
    'Volvafr': 'differentiated (germ-soma)', 'Pleosta': 'differentiated (partial soma)',
    'Chloatm': 'packets of cells (grey)', 'Retifil': 'multinucleate reticulose plasmodium (grey)',
    'Anaedeh': 'non-fruiting myxobacterium', 'Vulginc': 'non-fruiting myxobacterium',
    'Coralim': 'coenocytic stage (grey)', 'Syssmul': 'forms cell aggregates (grey)',
    'Mnemlei': 'ctenophore', 'Tremmes': 'jelly fruiting body; Tremellomycetes possibly a separate origin',
}

# ---------------------------------------------------------------- transitions
# id, kingdom, label, type, from, to, tier, M (node name | tip code | ('mrca', [...]) | ('tips', [...])), P override, notes
T = [
    ('T_metazoa', 'euk', 'Animals', 'clonal;complex;obligate', 'unicellular', 'complex', 'primary', 'Metazoa', None,
     'P state possibly facultatively colonial (choanoflagellate-like). Metazoa root = Porifera/Ctenophora/ParaHoxozoa '
     'polytomy: run sponge-first and ctenophore-first variants'),
    ('T_choano_colonies', 'euk', 'Choanoflagellate colonies', 'clonal(+aggregative in C. flexa);simple;facultative',
     'unicellular', 'grey', 'exploratory', ('tips', ['Salpros', 'Barrmon', 'Choafle', 'Codohol']), 'Choanoflagellata',
     'grey zone; coloniality scattered in Craspedida, no defensible internal origin node: tip-level only'),
    ('T_capsaspora_aggregation', 'euk', 'Capsaspora aggregation', 'aggregative;simple;facultative', 'unicellular',
     'grey', 'exploratory', 'Capsowc', None, 'single species; aggregative stage is a culture phenotype'),
    ('T_ichthyosporea_coenocytic', 'euk', 'Ichthyosporea/Teretosporea coenocytes', 'coenocytic;simple;obligate',
     'unicellular', 'grey', 'excluded', 'Teretosporea', None,
     'coenocytic development likely ancestral to Teretosporea (Corallochytrium coenocytic stage): P and M share the '
     'state; excluded from primary tests'),
    ('T_fungi_hyphae', 'euk', 'Fungal hyphae', 'hyphal(coenocytic/septate);simple;obligate', 'unicellular', 'simple',
     'primary', 'Blastocladiomycota_Amastigomycota', None,
     'M = Blastocladiomycota + terrestrial fungi; contains yeast reversals (coded separately); coenocytic hyphae are '
     'grey under a strict coding; Gonapodya = separate tip origin; Kiss et al. 2019: no kinase/receptor/adhesion expansion expected'),
    ('T_fungi_hyphae_gonapodya', 'euk', 'Gonapodya hyphae', 'hyphal;simple;obligate', 'unicellular', 'simple',
     'secondary', 'Gonapro', None, 'separate tip origin of hyphae within Chytridiomyceta'),
    ('T_fungi_pezizomycotina', 'euk', 'Pezizomycotina fruiting bodies', 'clonal;complex;facultative', 'simple', 'complex',
     'primary', 'Pezizomycotina', None,
     'P (saccharomyceta) is hyphal; Saccharomycotina = reversal (Nagy et al. 2014), not a unicellular sister'),
    ('T_fungi_agaricomycetes', 'euk', 'Agaricomycete fruiting bodies', 'clonal;complex;facultative', 'simple', 'complex',
     'primary', 'Agaricomycetes_Dacrymycetes', None,
     'M = Agaricomycetes + Dacrymycetes; Tremellomycetes possibly a separate origin (sensitivity: T_fungi_tremella); '
     'Cryptococcus = reversal'),
    ('T_fungi_tremella', 'euk', 'Tremella jelly fruiting body', 'clonal;complex;facultative', 'simple', 'complex',
     'exploratory', 'Tremmes', None, 'sensitivity coding (separate Tremellomycetes origin)'),
    ('T_fungi_neolecta', 'euk', 'Neolecta fruiting body', 'clonal;complex;facultative', 'simple', 'complex',
     'secondary', 'Neolirr', None, 'tip; BUSCO-limited; Taphrinomycotina internal order unresolved (polytomy)'),
    ('T_rev_saccharomycotina', 'euk', 'Saccharomycotina yeasts', 'reversal', 'simple', 'reversal', 'reversal',
     'Saccharomycotina', None, 'reversal to yeast growth (Yarrowia still dimorphic)'),
    ('T_rev_schizosaccharomyces', 'euk', 'Fission yeast', 'reversal', 'simple', 'reversal', 'reversal', 'Schipom', None,
     'tip reversal; P is a polytomy'),
    ('T_rev_saitoella', 'euk', 'Saitoella yeast', 'reversal', 'simple', 'reversal', 'reversal', 'Saitcom', None,
     'tip reversal; P is a polytomy'),
    ('T_rev_cryptococcus', 'euk', 'Cryptococcus yeast', 'reversal', 'simple', 'reversal', 'reversal', 'Crypneo', None,
     'tip reversal (dimorphic in sexual cycle)'),
    ('T_rev_rhodotorula', 'euk', 'Rhodotorula yeast', 'reversal', 'simple', 'reversal', 'reversal', 'Rhodtor', None,
     'tip reversal; ancestral state of Pucciniomycotina debated'),
    ('T_rev_zygnematophyceae', 'euk', 'Zygnematophyceae unicells', 'reversal', 'simple', 'reversal', 'reversal',
     'Zygnematophyceae', None,
     'all three sampled zygnematophytes are unicellular; filamentous Zygnematophyceae (Zygnema, Spirogyra) unsampled, '
     'so the node-level reversal is a sampling-dependent simplification'),
    ('T_dictyostelia', 'euk', 'Dictyostelid sorocarps', 'aggregative;simple;facultative', 'unicellular', 'simple',
     'primary', 'Dictyostelia', None, 'aggregative SIMPLE (Knoll 2011); Physarum (sister clade) is coenocytic/grey'),
    ('T_myxogastria_plasmodium', 'euk', 'Physarum plasmodium', 'coenocytic plasmodium;simple;obligate', 'unicellular',
     'grey', 'exploratory', 'Physpol', None, 'grey zone; transcriptome'),
    ('T_fonticula', 'euk', 'Fonticula sorocarp', 'aggregative;simple;facultative', 'unicellular', 'simple', 'secondary',
     'Fontalb', None, 'single tip; Nuclearia unavailable, Parvularia (transcriptome) is the sister'),
    ('T_acrasis', 'euk', 'Acrasis sorocarp', 'aggregative;simple;facultative', 'unicellular', 'simple', 'secondary',
     'Acrakon', None, 'P (Tetramitia) is a polytomy: Acrasis / Neovahlkampfia / Naegleriidae'),
    ('T_copromyxa', 'euk', 'Copromyxa sorocarp', 'aggregative;simple;facultative', 'unicellular', 'simple', 'secondary',
     'Coprpro', None, 'transcriptome; sister within Tubulinea contested (Brown et al. 2011); P is a polytomy'),
    ('T_guttulinopsis', 'euk', 'Guttulinopsis sorocarp', 'aggregative;simple;facultative', 'unicellular', 'simple',
     'excluded', 'Guttvul', None, 'presence-only (transcriptome BUSCO ~11%)'),
    ('T_volvocales_tetrabaena', 'euk', 'Volvocine origin 1 (Tetrabaenaceae)', 'clonal;simple;obligate', 'unicellular',
     'simple', 'primary', 'Tetrsoc', None,
     'Lindsey et al. 2021: Tetrabaenaceae arose independently (unsampled unicellular sister); P = TGV_clade inferred unicellular'),
    ('T_volvocales_gv', 'euk', 'Volvocine origin 2 (Gonium + Volvocaceae)', 'clonal;simple;obligate', 'unicellular',
     'simple', 'primary', 'Gonium_Volvocaceae', None,
     'Goniaceae paraphyletic (Lindsey et al. 2021); germ-soma differentiation arose >=4 times inside (Volvox, Pleodorina): '
     'coded "differentiated" at species level, not complex'),
    ('T_streptophyta', 'euk', 'Streptophyte multicellularity', 'clonal;simple;obligate', 'unicellular', 'simple',
     'primary', 'Klebsormidiophyceae_Phragmoplastophyta', None,
     'Chlorokybus packets make P grey; Zygnematophyceae inside M are reversals'),
    ('T_embryophyta', 'euk', 'Land plants', 'clonal;complex;obligate', 'simple', 'complex', 'primary', 'Embryophyta',
     None, 'nested S->X; P = Zygnematophyceae + Embryophyta'),
    ('T_chara', 'euk', 'Chara', 'clonal;complex;obligate', 'simple', 'complex', 'secondary', 'Charbra', None,
     'tip S->X'),
    ('T_rhodophyta', 'euk', 'Red algal multicellularity (Eurhodophytina)', 'clonal;simple;obligate', 'unicellular',
     'simple', 'primary', 'Eurhodophytina', None,
     'Stylonematophyceae/Compsopogonophyceae filaments = independent origins (not sampled as MC); Bangiales low BUSCO'),
    ('T_florideophyceae', 'euk', 'Florideophyceae', 'clonal;complex;obligate', 'simple', 'complex', 'primary',
     'Rhodymeniophycidae', None, 'nested S->X; all sampled florideophytes are Rhodymeniophycidae'),
    ('T_px_filaments', 'euk', 'PX-clade filaments', 'clonal;simple;obligate', 'unicellular', 'simple', 'secondary',
     'PX_clade', None,
     'filaments predate Phaeophyceae (Schizocladia, Phaeothamnion, Tribonema filamentous); P = Ochrophyta polytomy '
     '(incl. Raphidophyceae); Xanthophyceae also contain unicells'),
    ('T_phaeophyceae', 'euk', 'Brown algae', 'clonal;complex;obligate', 'simple', 'complex', 'primary', 'Phaeophyceae',
     None, 'nested S->X; Discosporangium and Ectocarpus are filamentous (coded simple); sensitivity M = '
     'Dictyotophycidae_Fucophycidae'),
    ('T_xanthophyceae_filaments', 'euk', 'Xanthophyceae', 'clonal;simple/siphonous;obligate', 'unicellular', 'grey',
     'exploratory', 'Xanthophyceae', None, 'subsumed by T_px_filaments; Vaucheria siphonous (grey)'),
    ('T_ulvophyceae', 'euk', 'Ulva', 'clonal;simple;obligate', 'unicellular', 'simple', 'secondary', 'Ulvamut', None,
     'no single Ulvophyceae M (Del Cortona et al. 2020); Ulva tip only'),
    ('T_siphonous_coenocytic', 'euk', 'Bryopsidales siphons', 'coenocytic;simple;obligate', 'unicellular', 'grey',
     'exploratory', 'Bryopsidales', None, 'grey zone (giant coenocytes)'),
    ('T_oomycete_hyphae', 'euk', 'Oomycete hyphae', 'coenocytic hyphae;simple;obligate', 'unicellular', 'grey',
     'exploratory', 'Oomycota', None, 'grey; hyphochytrid sisters unsampled (no proteome); P = Gyrista'),
    ('T_plasmodiophorid_plasmodium', 'euk', 'Plasmodiophora plasmodium', 'plasmodium;simple;obligate', 'unicellular',
     'grey', 'exploratory', 'Plasbra', None, 'grey; P = Rhizaria polytomy'),
    ('T_labyrinthula_nets', 'euk', 'Labyrinthula ectoplasmic nets', 'colonial net;simple;obligate', 'unicellular',
     'grey', 'exploratory', 'LabyspSRHaC', None, 'grey zone'),
    ('T_sorogena', 'euk', 'Sorogena sorocarp', 'aggregative;simple;facultative', 'unicellular', 'simple',
     'dropped_no_data', None, None, 'no sequence data for Sorogena; Platyophrya kept as background'),
    ('T_sorodiplophrys', 'euk', 'Sorodiplophrys sorocarp', 'aggregative;simple;facultative', 'unicellular', 'simple',
     'dropped_no_data', None, None, 'no sequence data; labyrinthulids kept as background'),
    # ---------------- prokaryotes
    ('T_cyano_filament', 'prok', 'Filamentous cyanobacteria', 'clonal;simple;obligate', 'unicellular', 'simple',
     'secondary', ('mrca', ['Planaga', 'Lyngaes', 'Phornig', 'Micrvag', 'Leptbor', 'Pseubic', 'Nostpun']), None,
     'MRCA of filamentous taxa = all non-Gloeobacter cyanobacteria: contains unicellular taxa (multicellularity early '
     'with reversals, Schirrmeister et al. 2011/2013, or multiple origins)'),
    ('T_cyano_heterocyst', 'prok', 'Heterocystous cyanobacteria', 'clonal;complex(>=3 cell types);obligate', 'simple',
     'complex', 'primary', 'Nostocales', None,
     'P = Nostocales + Chroococcidiopsis; P state uncertain (Chroococcidiopsis baeocytous/unicellular)'),
    ('T_cyano_branching', 'prok', 'True-branching cyanobacteria', 'clonal;complex;obligate', 'complex', 'complex',
     'exploratory', ('tips', ['Fiscthe', 'Chlofri']), 'Nostocales', 'tips; MRCA = Nostocales polytomy'),
    ('T_myxobacteria', 'prok', 'Myxobacterial fruiting bodies', 'aggregative;simple;facultative', 'unicellular',
     'simple', 'primary', ('mrca', ['Myxoxan', 'Coracor', 'Alloaur', 'Cystfus', 'Soracel', 'Choncro', 'Nannexe', 'Halioch']),
     None, 'M = Myxococcota contains non-fruiting Anaeromyxobacter and Vulgatibacter (losses or independent gains); '
           'P = Pseudomonadati polytomy'),
    ('T_streptomyces', 'prok', 'Streptomycete mycelium', 'clonal;simple;obligate', 'unicellular', 'simple', 'primary',
     'Streptomycetaceae', None,
     'P = Actinomycetes polytomy; other mycelial actinomycetes (Frankia, Micromonospora, ...) suggest mycelium may '
     'be ancestral to Actinomycetes with losses'),
    ('T_actino_mycelium_other', 'prok', 'Other mycelial actinomycetes', 'clonal;simple;obligate', 'unicellular',
     'simple', 'exploratory', ('tips', ['Franaln', 'Micrech', 'Therfus', 'Nocafar', 'Saccery']), 'Actinomycetes',
     'tips; no clade'),
    ('T_magnetotactic_mmp', 'prok', 'Multicellular magnetotactic bacteria', 'clonal;simple;obligate', 'unicellular',
     'simple', 'secondary', ('mrca', ['Magnmul', 'MagnspHK1']), None,
     'metagenomic drafts; MRCA = Desulfobacterales contains unicellular taxa'),
    ('T_cable_bacteria', 'prok', 'Cable bacteria', 'clonal;simple;obligate', 'unicellular', 'simple', 'primary',
     'Cable_bacteria', None, 'Electrothrix + Electronema'),
    ('T_jeongeupia', 'prok', 'Jeongeupia sacculi', 'clonal;simple;facultative', 'unicellular', 'simple', 'secondary',
     'Jeonsac', None, 'tip; P = Jeongeupia (polytomy of 3)'),
    ('T_methanosarcina', 'prok', 'Methanosarcina packets', 'clonal;simple;facultative', 'unicellular', 'simple',
     'primary', 'Methanosarcina', None, ''),
    ('T_haloferax_tissue', 'prok', 'Haloferax tissue-like growth', 'clonal;simple;facultative', 'unicellular', 'grey',
     'exploratory', 'Halovol', None, 'single species; condition-dependent'),
    ('T_archaea_filament', 'prok', 'Archaeal filaments', 'sheathed filaments;simple;obligate', 'unicellular', 'grey',
     'exploratory', ('tips', ['Methsoe', 'Methhun']), 'Methanomicrobia', 'tips'),
]


def leafset(nd):
    return [l.name for l in nd.iter_leaves()]


def main():
    trees = {k: Tree(open(f'{P}/data/guide_tree_{k}.nwk').read(), format=1) for k in ('euk', 'prok')}
    sp = pd.concat([pd.read_csv(f'{P}/data/species_final.tsv', sep='\t', dtype=str).assign(kingdom='euk'),
                    pd.read_csv(f'{P}/data/species_final_prok.tsv', sep='\t', dtype=str).assign(kingdom='prok')])
    codes = set(sp['code'])
    for c in [c for c in list(STATE) if c not in codes]:
        del STATE[c]   # species not (or no longer) in the final set
    state = {c: STATE.get(c, 'unicellular') for c in sp['code']}
    rows, member = [], {c: [] for c in codes}
    for (tid, k, label, typ, fr, to, tier, Mspec, Pover, notes) in T:
        t = trees[k]
        byname = {n.name: n for n in t.traverse()}
        flags = []
        if Mspec is None:
            rows.append(dict(transition_id=tid, kingdom=k, label=label, type=typ, from_state=fr, to_state=to,
                             analysis_tier=tier, primary_eligible=0, M='NA', M_kind='NA', P='NA',
                             sister_clade_nodes='NA', n_mc_species=0, n_sister_species=0, mc_species='',
                             sister_species='', prior_state_ref_node='NA', n_prior_state_ref=0,
                             prior_state_ref_species='', flags='no_data', notes=notes))
            continue
        if isinstance(Mspec, tuple) and Mspec[0] == 'tips':
            tips = Mspec[1]
            Mnode = byname[Pover]
            Pnode = Mnode
            M, Mkind, Pname = ','.join(tips), 'tips', Pover
            mc = tips
            sisters = [n for n in leafset(Pnode) if n not in tips]
            sister_nodes = 'tips_outside_M_within_P'
        else:
            if isinstance(Mspec, tuple) and Mspec[0] == 'mrca':
                Mnode = t.get_common_ancestor(*[byname[c] for c in Mspec[1]])
                Mkind = 'mrca_node'
            else:
                Mnode = byname[Mspec]
                Mkind = 'tip' if Mnode.is_leaf() else 'node'
            Pnode = byname[Pover] if Pover else Mnode.up
            M, Pname = Mnode.name, Pnode.name
            mc = leafset(Mnode)
            sisters = [n for n in leafset(Pnode) if n not in set(mc)]
            sib = [c for c in Pnode.children if c is not Mnode]
            sister_nodes = ';'.join(c.name for c in sib)
        # flags
        if len(Pnode.children) > 2:
            flags.append(f'P_polytomy({len(Pnode.children)})')
        if to != 'reversal':
            wrong = [c for c in mc if state[c] == 'unicellular']
            if wrong and Mkind != 'tips':
                flags.append('M_contains_unicellular(' + ','.join(sorted(wrong)) + ')')
            mcount = [c for c in mc if state[c] != 'unicellular']
        else:
            mcount = mc
            wrong = [c for c in mc if c not in REVERSAL]
            if wrong:
                flags.append('M_contains_non_reversal(' + ','.join(sorted(wrong)) + ')')
        if fr == 'unicellular':
            ms = [c for c in sisters if state[c] in (S, X)]
            if ms:
                flags.append('sister_contains_multicellular(' + ','.join(sorted(ms)) + ')')
        # nearest relatives in the PRIOR state (walk up from M until such leaves exist outside M)
        want = ({S, X} if to == 'reversal' else {'unicellular'} if fr == 'unicellular' else {S} if fr == 'simple'
                else {X} if fr == 'complex' else {G})
        mcs, ref_node, ref = set(mc), 'NA', []
        nd = Mnode if Mkind != 'tips' else Pnode
        while nd is not None:
            ref = [c for c in leafset(nd) if c not in mcs and state[c] in want
                   and not (fr == 'unicellular' and c in REVERSAL)]
            if ref:
                ref_node = nd.name
                break
            nd = nd.up
        tx = sp.set_index('code')
        tr = [c for c in mc + sisters if tx.loc[c, 'data_type_final'] != 'genome'] if 'data_type_final' in tx else []
        if tr:
            flags.append(f'transcriptome_taxa={len(tr)}')
        for c in mcount:
            member[c].append(f'{tid}:M')
        for c in sisters:
            member[c].append(f'{tid}:sister')
        rows.append(dict(transition_id=tid, kingdom=k, label=label, type=typ, from_state=fr, to_state=to,
                         analysis_tier=tier, primary_eligible=int(tier == 'primary'), M=M, M_kind=Mkind, P=Pname,
                         sister_clade_nodes=sister_nodes, n_mc_species=len(mcount), n_sister_species=len(sisters),
                         mc_species=','.join(sorted(mc)), sister_species=','.join(sorted(sisters)),
                         prior_state_ref_node=ref_node, n_prior_state_ref=len(ref),
                         prior_state_ref_species=','.join(sorted(ref)),
                         flags=';'.join(flags), notes=notes))
    tdf = pd.DataFrame(rows)
    tdf.to_csv(f'{P}/data/transitions.tsv', sep='\t', index=False)
    # traits
    busco = {}
    bp = f'{P}/data/qc/busco_summary.tsv'
    if os.path.exists(bp):
        b = pd.read_csv(bp, sep='\t')
        busco = dict(zip(b['code'], b['C']))
    mtype = {}
    deepest = {}
    for (tid, k, label, typ, fr, to, tier, Mspec, Pover, notes) in T:
        for c, roles in member.items():
            if f'{tid}:M' in roles and to != 'reversal':
                # most derived transition wins (later-listed nested ones are more derived when M is smaller)
                n = tdf.set_index('transition_id').loc[tid, 'n_mc_species']
                if c not in deepest or n < deepest[c][0]:
                    deepest[c] = (n, tid, typ)
    out = []
    for _, r in sp.iterrows():
        c = r['code']
        st = state[c]
        out.append(dict(code=c, kingdom=r['kingdom'], species=r['species'],
                        multicellular=int(st in (S, X)), multicellular_permissive=int(st in (S, X, G)),
                        mc_state=st, reversal=int(c in REVERSAL),
                        mc_type=(deepest[c][2] if c in deepest and st != 'unicellular' else ('reversal' if c in REVERSAL else 'unicellular')),
                        transition_id=deepest[c][1] if c in deepest and st != 'unicellular' else '',
                        transition_roles=';'.join(member[c]), state_note=STATE_NOTE.get(c, ''),
                        data_type=r.get('data_type_final', r['data_type']), busco_C=busco.get(c, '')))
    pd.DataFrame(out).to_csv(f'{P}/data/species_traits.tsv', sep='\t', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 120)
    print(tdf[['transition_id', 'analysis_tier', 'M', 'P', 'n_mc_species', 'n_sister_species',
               'prior_state_ref_node', 'n_prior_state_ref']].to_string())
    for _, r in tdf.iterrows():
        if r['flags']:
            print(r['transition_id'], '|', r['flags'])


if __name__ == '__main__':
    main()
