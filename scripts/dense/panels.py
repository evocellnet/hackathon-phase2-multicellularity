"""Per-panel configuration for the dense round-2 datasets (scripts/dense/build_panel.py).

edits : guide-tree edits, same grammar as scripts/data/build_guide_trees.py (op, args, rationale, reference);
        ops: group (new clade from existing siblings), move, dissolve (polytomy), rename. Member spec 'A|B' = first existing.
events: transitions inside the panel. M = node (or leaf) where the new state is first present; P (= parent of M) and
        sister nodes are derived from the tree. kind: origin (uni->simple), nested (simple->complex),
        origin_complex (uni->complex in one step at our resolution), reversal (loss), sensitivity (alternative coding).
"""
PANELS = {}

# ----------------------------------------------------------------------------------------------------------- R2 Fungi
PANELS['R2'] = dict(title='Fungi / Holomycota', busco_lineage='eukaryota_odb10', edits=[], events=[])
_SP16 = 'Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653'
_N18 = 'Nagy, Kovacs & Krizsan 2018 Biol Rev 93:1778'
PANELS['R2']['edits'] = [
    ('dissolve', ('Fungi_incertae_sedis',), 'NCBI bin is not a clade; replaced by the groups below', '-'),
    ('group', ('Mucoromycota', ['Mucoromycotina', 'Morthyg']), 'Mortierellomycotina + Mucoromycotina', _SP16),
    ('group', ('Zoopagomycota', ['Basimer', 'Coemumb']), 'Entomophthoromycotina + Kickxellomycotina', _SP16),
    ('group', ('Mucoromycota_Dikarya', ['Mucoromycota', 'Dikarya']), 'Mucoromycota sister to Dikarya', _SP16),
    ('group', ('Amastigomycota', ['Zoopagomycota', 'Mucoromycota_Dikarya']), 'terrestrial (non-flagellated) fungi', _SP16),
    ('group', ('Blastocladiomycota_Amastigomycota', ['Allomac', 'Amastigomycota']),
     'Blastocladiomycota sister to terrestrial fungi; M of the hyphal transition',
     'James et al. 2006 Nature 443:818; Chang et al. 2015 GBE 7:1590; Li et al. 2021 Curr Biol 31:1653'),
    ('group', ('Eumycota', ['chytridiomyceta', 'Blastocladiomycota_Amastigomycota']), 'Fungi sensu stricto', _SP16),
    ('group', ('Aphelida_Eumycota', ['Amoeocc|Amoepro', 'Eumycota']), 'Aphelida sister to Fungi s.s.; Rozellomycota branch first',
     'Torruella et al. 2018 Commun Biol 1:231'),
    ('group', ('Holomycota', ['Rotosphaerida', 'Fungi']), 'nucleariids + Fonticula sister to Fungi',
     'Brown et al. 2009 Mol Biol Evol 26:2699; Torruella et al. 2015 Curr Biol 25:2404'),
    ('group', ('Yarrowia_Saccharomyces', ['Yarrlip', 'Sacccer']), 'Lipomycetaceae branch first within Saccharomycotina',
     'Shen et al. 2018 Cell 175:1533'),
    ('dissolve', ('Agaricomycetes_incertae_sedis',), 'NCBI bin (Auriculariales, Sebacinales, Polyporales) is not a clade', '-'),
    ('group', ('Polyporales_Agaricales', ['Tramcoc', 'Agaricales']),
     'Polyporales + Agaricomycetidae form the late-diverging Agaricomycetes; Auriculariales and Sebacinales branch earlier (order unresolved)',
     'Hibbett et al. 2007 Mycol Res 111:509; Varga et al. 2019 Nat Ecol Evol 3:668'),
    ('group', ('Agaricomycetes_Dacrymycetes', ['Agaricomycetes', 'Dacrpri']),
     'Dacrymycetes sister to Agaricomycetes; Tremellomycetes sister to both', 'Nagy et al. 2016 MBE 33:959; He et al. 2019 Fungal Divers 99:105'),
    ('group', ('Wallemiomycotina_Agaricomycotina', ['Wallmel', 'Agaricomycotina']), 'Wallemiomycotina sister to Agaricomycotina',
     'Zhao et al. 2017 Fungal Divers 84:43; He et al. 2019'),
    ('group', ('Ustilaginomycotina_Agaricomycotina', ['Ustimay', 'Wallemiomycotina_Agaricomycotina']),
     'Pucciniomycotina sister to Ustilaginomycotina + Agaricomycotina', 'He et al. 2019; Li et al. 2021'),
]
PANELS['R2']['expect'] = ['Holomycota', 'Eumycota', 'Blastocladiomycota_Amastigomycota', 'Pezizomycotina',
                          'Agaricomycetes_Dacrymycetes', 'Saccharomycotina']
PANELS['R2']['tree_notes'] = [
    'Polytomies kept: Pezizomycotina (Orbiliomycetes / Pezizomycetes / leotiomyceta), leotiomyceta (Eurotiomycetes, '
    'Dothideomycetes, sordariomyceta), Taphrinomycotina (Neolecta, Saitoella, Schizosaccharomyces, Taphrina), '
    'Agaricomycetes (Auriculariales, Sebacinales, Polyporales+Agaricales).',
    'Nuclearia/Parvularia: no genome; Parvularia atlantis transcriptome (EukProt) is the only sister of Fonticula.',
]
PANELS['R2']['events'] = [
    dict(id='R2_hyphae', label='hyphal multicellularity (Blastocladiomycota + terrestrial fungi)', kind='origin', M='Blastocladiomycota_Amastigomycota',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Kiss et al. 2019 Nat Commun 10:4080; ' + _N18,
         note='stepwise hyphal origin; Kiss et al. place most hyphal genes earlier (pre-adaptation)'),
    dict(id='R2_gonapodya', label='Gonapodya hyphal-like thallus', kind='origin', M='Gonapro', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Chang et al. 2015 GBE 7:1590', note='tip event'),
    dict(id='R2_fonticula', label='Fonticula sorocarpic aggregation', kind='origin', M='Fontalb',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Brown et al. 2009 Mol Biol Evol 26:2699',
         note='aggregative; tip event; only sister is the Parvularia transcriptome'),
    dict(id='R2_pezizomycotina', label='Pezizomycotina fruiting bodies', kind='nested', M='Pezizomycotina', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref=_N18 + '; Nguyen et al. 2017 Nat Commun 8:14444',
         note='flag: the Pezizomycotina stem may hide several convergent origins (Nagy et al. 2018)'),
    dict(id='R2_neolecta', label='Neolecta fruiting body', kind='nested', M='Neolirr', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref='Nguyen et al. 2017 Nat Commun 8:14444',
         note='tip event; Taphrinomycotina ancestor yeast-like or dimorphic (P state uncertain)'),
    dict(id='R2_agaricomycetes', label='Agaricomycetes + Dacrymycetes fruiting bodies', kind='nested', M='Agaricomycetes_Dacrymycetes',
         **{'from': 'simple', 'to': 'complex'}, ref=_N18 + '; Krizsan et al. 2019 PNAS 116:7409'),
    dict(id='R2_tremella', label='Tremella basidiocarp', kind='nested', M='Tremmes', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref=_N18, note='tip event; sister Cryptococcus is a yeast'),
    dict(id='R2_endogone', label='Endogone sporocarp', kind='nested', M='Endosp', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref=_N18 + '; Chang et al. 2019 New Phytol 222:511', note='tip event'),
    dict(id='R2_rev_saccharomycotina', label='Saccharomycotina yeasts', kind='reversal', M='Saccharomycotina', conf='grey',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Nagy et al. 2014 Nat Commun 5:4471', note='Yarrowia keeps true hyphae'),
    dict(id='R2_rev_schizosaccharomyces', label='fission yeast', kind='reversal', M='Schipom',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Nagy et al. 2014 Nat Commun 5:4471'),
    dict(id='R2_rev_saitoella', label='Saitoella yeast', kind='reversal', M='Saitcom',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Nagy et al. 2014 Nat Commun 5:4471'),
    dict(id='R2_rev_cryptococcus', label='Cryptococcus yeast', kind='reversal', M='Crypneo', conf='grey',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Nagy et al. 2014 Nat Commun 5:4471'),
    dict(id='R2_rev_rhodotorula', label='Rhodotorula yeast', kind='reversal', M='Rhodtor',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Nagy et al. 2014 Nat Commun 5:4471'),
    dict(id='R2_rev_serendipita', label='Serendipita (no fruiting body)', kind='reversal', M='Sereind', conf='grey',
         **{'from': 'complex', 'to': 'simple'}, ref=_N18 + '; Zuccaro et al. 2011 PLoS Pathog 7:e1002290'),
]

# ------------------------------------------------------------------------------------------------------- R4 Volvocales
_L21 = 'Lindsey et al. 2021 BMC Biol 19:182'
PANELS['R4'] = dict(title='Volvocine algae (Chlamydomonadales)', busco_lineage='eukaryota_odb10', edits=[
    ('move', ('Chlaeus', 'Chlamydomonadales'), 'genus Chlamydomonas is polyphyletic; C. eustigma is outside Reinhardtinia',
     'Proschold et al. 2001 Protist 152:265'),
    ('move', ('Chlapri', 'Chlamydomonadales'), 'C. priscui (UWO241) is outside Reinhardtinia (Moewusinia)',
     'Proschold et al. 2001 Protist 152:265; Possmayer et al. 2016 J Phycol'),
    ('move', ('Chlamydomonas', 'Chlamydomonadales'), 'Chlamydomonas s.str. (C. reinhardtii group) belongs to Reinhardtinia, '
     'not with other NCBI Chlamydomonadaceae', 'Proschold et al. 2001 Protist 152:265; Craig et al. 2021 Plant Cell 33:1016'),
    ('dissolve', ('Chlamydomonadaceae',), 'NCBI family is polyphyletic', 'Proschold et al. 2001 Protist 152:265'),
    ('move', ('Gonipec', 'Chlamydomonadales'), 'NCBI places Gonium in Volvocaceae; Goniaceae is paraphyletic', _L21),
    ('dissolve', ('Volvox',), 'genus Volvox is polyphyletic; Volvocaceae internal order left unresolved',
     'Herron & Michod 2008 Evolution 62:436; ' + _L21),
    ('group', ('Volvox_africanus_reticuliferus', ['Volvafr', 'Volvret']), 'V. reticuliferus sister of V. africanus',
     'Nozaki et al. 2015 PLoS One 10:e0142632'),
    ('group', ('Astrephomene_Volvocaceae', ['Astrgub', 'Volvocaceae']), 'Astrephomene sister to Volvocaceae', _L21 + ' (MLBS 98)'),
    ('group', ('Gonium_Volvocaceae', ['Gonipec', 'Astrephomene_Volvocaceae']),
     'Gonium sister to Astrephomene + Volvocaceae (Goniaceae paraphyletic); M of the 2nd colonial origin', _L21),
    ('group', ('TGV_clade', ['Tetrsoc', 'Gonium_Volvocaceae']),
     'Tetrabaenaceae + Goniaceae + Volvocaceae among SAMPLED taxa; an unsampled unicellular lineage is sister of '
     'Tetrabaenaceae, so the ancestor here is inferred unicellular and colonial multicellularity arose twice', _L21),
    ('group', ('Chlamydomonas_TGV', ['Chlamydomonas', 'TGV_clade']),
     'core Reinhardtinia: C. reinhardtii/incerta/schloesseri sister to the colonial lineages', 'Craig et al. 2021 Plant Cell 33:1016; ' + _L21),
    ('group', ('Reinhardtinia', ['Edapdeb', 'Chlamydomonas_TGV']), 'Edaphochlamys outside core Reinhardtinia',
     'Craig et al. 2021; ' + _L21),
], events=[
    dict(id='R4_tetrabaena', label='Tetrabaena 4-celled colony (origin 1)', kind='origin', M='Tetrsoc',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_L21 + '; Featherston et al. 2018 Mol Biol Evol 35:855',
         note='tip event; Tetrabaena BUSCO C=32 (flag)'),
    dict(id='R4_gonium_volvocaceae', label='Goniaceae + Volvocaceae colonies (origin 2)', kind='origin', M='Gonium_Volvocaceae',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_L21 + '; Hanschen et al. 2016 Nat Commun 7:11370'),
    dict(id='R4_sens_single_origin', label='single colonial origin (alternative coding)', kind='sensitivity', M='TGV_clade', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Herron et al. 2009 PNAS 106:3254; Hanschen et al. 2016',
         note='older single-origin hypothesis; use as sensitivity coding only'),
    dict(id='R4_astrephomene', label='Astrephomene soma (independent germ-soma)', kind='nested', M='Astrgub',
         **{'from': 'simple', 'to': 'complex'}, ref='Yamashita et al. 2016 BMC Evol Biol 16:243; ' + _L21, note='complex = D only'),
    dict(id='R4_pleodorina', label='Pleodorina partial soma', kind='nested', M='Pleosta', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref='Herron & Michod 2008 Evolution 62:436; ' + _L21, note='complex = D only'),
    dict(id='R4_volvox_carteri', label='Volvox carteri germ-soma', kind='nested', M='Volvcar',
         **{'from': 'simple', 'to': 'complex'}, ref='Prochnik et al. 2010 Science 329:223; ' + _L21, note='complex = D only'),
    dict(id='R4_volvox_africanus', label='Volvox africanus/reticuliferus germ-soma', kind='nested', M='Volvox_africanus_reticuliferus',
         **{'from': 'simple', 'to': 'complex'}, ref='Herron & Michod 2008; ' + _L21, note='complex = D only'),
], expect=['Reinhardtinia', 'TGV_clade', 'Gonium_Volvocaceae', 'Astrephomene_Volvocaceae', 'Volvocaceae'],
    tree_notes=[
        'Volvocaceae internal order left as a polytomy (Volvox polyphyletic); Lindsey et al. 2021 infer 4-6 origins of '
        'cellular differentiation, coded here as separate tip/clade events.',
        'Undifferentiated Volvocaceae (Yamagishiella, Eudorina, Pandorina) have no published annotation; if present they are '
        'AUGUSTUS ab initio predictions (chlamy2011 model) on the unannotated Volvocales Genome Consortium assemblies (flagged).',
        'Chlamydomonadales outside Reinhardtinia (C. eustigma, C. priscui, Sanguina, Dunaliella) kept as a polytomy.'])

# ----------------------------------------------------------------------------------------------------- R3 Streptophyta
_W14 = 'Wickett et al. 2014 PNAS 111:E4859; One Thousand Plant Transcriptomes Initiative 2019 Nature 574:679'
PANELS['R3'] = dict(title='Streptophyta', busco_lineage='eukaryota_odb10', edits=[
    ('group', ('Chlorophyta_Streptophyta', ['Chlorophyta', 'Streptophyta']),
     'Prasinodermophyta (Prasinoderma) sister to Chlorophyta + Streptophyta', 'Li L. et al. 2020 Nat Ecol Evol 4:1220'),
    ('group', ('Mesostigmatophyceae_Chlorokybophyceae', ['Mesovir', 'Chloatm']),
     'Mesostigma + Chlorokybus form the earliest streptophyte clade', 'Wang S. et al. 2020 Nat Plants 6:95; Lemieux et al. 2007 BMC Biol 5:2'),
    ('group', ('Klebsormidiophyceae_Phragmoplastophyta', ['Klebnit', 'Streptophytina']),
     'grade: Klebsormidiophyceae, then Phragmoplastophyta (NCBI Streptophytina); M of streptophyte multicellularity', _W14),
    ('group', ('Anydrophyta', ['Zygnematophyceae', 'Embryophyta']), 'Zygnematophyceae sister to land plants',
     _W14 + '; Cheng et al. 2019 Cell 179:1057'),
    ('group', ('Coleochaetophyceae_Anydrophyta', ['Colescu', 'Anydrophyta']), 'Coleochaetophyceae sister to Anydrophyta', _W14),
    ('group', ('Setaphyta', ['Marcpol', 'Physpat']), 'liverworts + mosses (bryophyte monophyly; hornworts not sampled)',
     'Puttick et al. 2018 Curr Biol 28:733; Harris et al. 2020 Curr Biol 30:2001'),
], events=[
    dict(id='R3_streptophyte_filaments', label='streptophyte multicellularity (Klebsormidium + Phragmoplastophyta)', kind='origin',
         M='Klebsormidiophyceae_Phragmoplastophyta', **{'from': 'unicellular', 'to': 'simple'},
         ref='Umen 2014 Cold Spring Harb Perspect Biol 6:a016170; Bowles et al. 2020 Curr Biol 30:530',
         note='Chlorokybus sarcinoid packets coded unicellular (grey)'),
    dict(id='R3_embryophyta', label='land plants (embryophytes)', kind='nested', M='Embryophyta',
         **{'from': 'simple', 'to': 'complex'}, ref='Knoll 2011; Bowles et al. 2020 Curr Biol 30:530',
         note='tree sister (Zygnematophyceae) is secondarily unicellular in our sampling; simple-state references = Coleochaete, Chara, Klebsormidium'),
    dict(id='R3_chara', label='Chara secondary complexity', kind='nested', M='Charbra', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref='Nishiyama et al. 2018 Cell 174:448', note='tip event'),
    dict(id='R3_rev_zygnematophyceae', label='Zygnematophyceae unicellularity', kind='reversal', M='Zygnematophyceae', conf='grey',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Cheng et al. 2019 Cell 179:1057; Feng et al. 2024 Nat Genet 56:1018',
         note='filamentous Zygnematophyceae (Zygnema, Spirogyra) have no public annotated genome here; reversal timing uncertain'),
], expect=['Klebsormidiophyceae_Phragmoplastophyta', 'Anydrophyta', 'Embryophyta', 'Zygnematophyceae'],
    tree_notes=['No public annotated genome found (NCBI, UniProt 2026_03, EukProt, LukProt; 2026-10-29) for Zygnema, Penium, '
                'Spirogyra, Anthoceros/hornworts, Coleochaete (genomes of Feng et al. 2024 and Li et al. 2020 are hosted outside '
                'INSDC/UniProt or JGI-restricted). Coleochaete is a transcriptome (BUSCO 46).',
                'Zygnematophyceae internal order, and tracheophyte internal order kept from NCBI.'])

# ------------------------------------------------------------------------------------------ R6 Stramenopiles / brown algae
PANELS['R6'] = dict(title='Stramenopiles / Phaeophyceae', busco_lineage='eukaryota_odb10', edits=[], events=[])
_D24 = 'Denoeud et al. 2024 Cell 187:6943'
_SIL = 'Silberfeld et al. 2014 Cryptogam Algol 35:117'
PANELS['R6']['edits'] = [
    ('group', ('Discosporangiales', ['Chorten', 'Discmes']), 'Choristocarpus + Discosporangium = Discosporangiales', _SIL + '; ' + _D24),
    ('move', ('Discosporangiales', 'Phaeophyceae'), 'NCBI places Discosporangiales inside Sphacelariales; they are the earliest-diverging '
     'brown algal order', _SIL + '; ' + _D24),
    ('group', ('Dictyotophycidae', ['Dictdic', 'Sphacelariales|Halopan']), 'Dictyotales + Sphacelariales', _SIL + '; ' + _D24),
    ('group', ('Fucophycidae', ['Desmher', 'Ectocarpales', 'Fucaceae', 'Laminariales']),
     'brown algal crown radiation (Fucophycidae); internal order left unresolved', _SIL + '; ' + _D24),
    ('group', ('Dictyotophycidae_Fucophycidae', ['Dictyotophycidae', 'Fucophycidae']),
     'Discosporangiales sister to all other brown algae', _SIL + '; ' + _D24),
    ('group', ('Phaeophyceae_Schizocladia', ['Phaeophyceae', 'Schiisc']), 'Schizocladia sister to brown algae',
     'Kawai et al. 2003 J Phycol 39:556; ' + _D24),
    ('group', ('Gyrista', ['Ochrophyta', 'Oomycota']), 'Ochrophyta + Pseudofungi (oomycetes); Bigyra outside',
     'Derelle et al. 2016 MBE 33:2890; Burki et al. 2020 TREE 35:43'),
]
PANELS['R6']['expect'] = ['PX_clade', 'Phaeophyceae', 'Fucophycidae', 'Dictyotophycidae', 'Discosporangiales']
PANELS['R6']['tree_notes'] = [
    'Ochrophyta internal order (diatoms, Pelagophyceae, Eustigmatophyceae, Chrysophyceae, Raphidophyceae, PX clade) left as '
    'a polytomy; PX clade internal order (Chrysoparadoxa, Phaeothamnion, Tribonema, Schizocladia+Phaeophyceae) unresolved.',
    'Phaeoexplorer genomes vary strongly in completeness; Sphacelaria, Macrocystis, Saccorhiza x2, Himanthalia, Laminaria were '
    'rejected (BUSCO eukaryota C < 35); Discosporangiales (Discmes 46, Chorten 34) and Halopteris (48) kept as the only '
    'representatives of their orders (flagged). Ishigeales: no genome.']
PANELS['R6']['events'] = [
    dict(id='R6_px_filaments', label='PX clade filaments (Schizocladia, Phaeothamnion, Tribonema, brown algae)', kind='origin', M='PX_clade',
         conf='grey', **{'from': 'unicellular', 'to': 'simple'}, ref='Derelle et al. 2016 MBE 33:2890; ' + _D24,
         note='Chrysoparadoxa (coccoid) inside PX: coded as reversal; P = Ochrophyta polytomy'),
    dict(id='R6_fucophycidae', label='Fucophycidae parenchyma (PRIMARY complex event)', kind='nested', M='Fucophycidae',
         **{'from': 'simple', 'to': 'complex'}, ref='Bringloe et al. 2020 Crit Rev Plant Sci; ' + _SIL + '; ' + _D24,
         note='Ectocarpus/Pylaiella/Chordaria secondarily filamentous (reversals below)'),
    dict(id='R6_sens_phaeophyceae', label='Phaeophyceae stem (alternative coding)', kind='sensitivity', M='Phaeophyceae', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref='Cock et al. 2010 Nature 465:617; ' + _D24),
    dict(id='R6_dictyota', label='Dictyota parenchymatous thallus', kind='nested', M='Dictdic', conf='grey',
         **{'from': 'simple', 'to': 'complex'}, ref='Bringloe et al. 2020; ' + _SIL, note='tip event; sister Halopteris (polystichous)'),
    dict(id='R6_rev_ectocarpus', label='Ectocarpus filaments', kind='reversal', M='Ectosil', conf='grey',
         **{'from': 'complex', 'to': 'simple'}, ref='Cock et al. 2010 Nature 465:617; Bringloe et al. 2020'),
    dict(id='R6_rev_pylaiella', label='Pylaiella filaments', kind='reversal', M='Pylalit', conf='grey',
         **{'from': 'complex', 'to': 'simple'}, ref='Bringloe et al. 2020'),
    dict(id='R6_rev_chordaria', label='Chordaria pseudoparenchyma', kind='reversal', M='Chorlin', conf='grey',
         **{'from': 'complex', 'to': 'simple'}, ref='Bringloe et al. 2020'),
    dict(id='R6_rev_chrysoparadoxa', label='Chrysoparadoxa coccoid', kind='reversal', M='Chryaus', conf='grey',
         **{'from': 'simple', 'to': 'unicellular'}, ref='Wetherbee et al. 2019 J Phycol'),
]

# ------------------------------------------------------------------------------------------------------- R5 Rhodophyta
PANELS['R5'] = dict(title='Rhodophyta', busco_lineage='eukaryota_odb10', edits=[], events=[])
_M17 = 'Munoz-Gomez et al. 2017 Curr Biol 27:1677'
PANELS['R5']['edits'] = [
    ('dissolve', ('Bangiophyceae',), 'NCBI node lumps Bangiales with Porphyridiophyceae and Cyanidiophyceae; regrouped below', '-'),
    ('group', ('Cyanidiophyceae', ['Cyanmer', 'Galdieria']), 'Cyanidiophyceae sister to all other red algae',
     'Yoon et al. 2006 J Phycol 42:482; ' + _M17),
    ('group', ('Proteorhodophytina', ['Porppur', 'Compcae', 'Stylonemataceae']),
     'Porphyridiophyceae + Compsopogonophyceae + Stylonematophyceae (+ Rhodellophyceae, unsampled)', _M17),
    ('group', ('Eurhodophytina', ['Bangiaceae', 'Rhodymeniophycidae']), 'Bangiophyceae + Florideophyceae',
     'Yoon et al. 2006; ' + _M17),
    ('group', ('Rhodophytina', ['Proteorhodophytina', 'Eurhodophytina']), 'mesophilic red algae', _M17),
    ('group', ('Rhodophyta_Rhodelphidia', ['Rhodophyta', 'Rhodemar']), 'Rhodelphis sister to red algae', 'Gawryluk et al. 2019 Nature 572:240'),
    ('group', ('Archaeplastida', ['Cyanpar', 'Rhodophyta_Rhodelphidia', 'Chlorophyta']),
     'Archaeplastida; order of glaucophytes / red / green left unresolved', 'Burki et al. 2020 TREE 35:43'),
]
PANELS['R5']['expect'] = ['Cyanidiophyceae', 'Proteorhodophytina', 'Eurhodophytina', 'Rhodymeniophycidae']
PANELS['R5']['tree_notes'] = [
    'Florideophyceae are represented only by Rhodymeniophycidae (Chondrus, Gracilariopsis, Agarophyton, Delisea); no annotated '
    'genome of Hildenbrandiophycidae, Nemaliophycidae or Corallinophycidae passed QC (Calliarthron EukProt C 13, Hildenbrandia EST). '
    'The Florideophyceae stem and the Rhodymeniophycidae crown are therefore the same branch here.',
    'Proteorhodophytina internal order unresolved (polytomy); filamentous Bangiopsis and Compsopogon are MMETSP transcriptomes (flag).',
    'Bangiales genomes have low BUSCO completeness (Porphyra 44, Pyropia yezoensis 50, P. vietnamensis 61).']
PANELS['R5']['events'] = [
    dict(id='R5_eurhodophytina', label='Bangiophyceae + Florideophyceae multicellularity', kind='origin', M='Eurhodophytina', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_M17 + '; Yoon et al. 2006',
         note='contested: filamentous Proteorhodophytina (Bangiopsis, Compsopogon) imply either a single Rhodophytina origin or several'),
    dict(id='R5_florideophyceae', label='Florideophyceae (complex)', kind='nested', M='Rhodymeniophycidae',
         **{'from': 'simple', 'to': 'complex'}, ref='Knoll 2011; Collen et al. 2013 PNAS 110:5247',
         note='only Rhodymeniophycidae sampled: M is the crown of one subclass (stem/crown not separable)'),
    dict(id='R5_bangiopsis', label='Stylonematophyceae filaments (Bangiopsis)', kind='origin', M='Bangsp', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_M17, note='tip event; transcriptome'),
    dict(id='R5_compsopogon', label='Compsopogon filaments', kind='origin', M='Compcae', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_M17, note='tip event; transcriptome'),
    dict(id='R5_sens_rhodophytina', label='single origin at Rhodophytina (alternative coding)', kind='sensitivity', M='Rhodophytina',
         conf='grey', **{'from': 'unicellular', 'to': 'simple'}, ref=_M17,
         note='implies reversals in Porphyridium and Rhodosorus'),
]

# ---------------------------------------------------------------------------------------------------------- R1 Holozoa
PANELS['R1'] = dict(title='Holozoa', busco_lineage='eukaryota_odb10', edits=[], events=[])
_T15 = 'Torruella et al. 2015 Curr Biol 25:2404; Hehenberger et al. 2017 Curr Biol 27:2043'
PANELS['R1']['edits'] = [
    ('rename', ('Opisthokonta_incertae_sedis', 'Pluriformea'), 'NCBI bin holds exactly Corallochytrium + Syssomonas (Pluriformea)',
     'Hehenberger et al. 2017 Curr Biol 27:2043'),
    ('group', ('Teretosporea', ['Ichthyosporea', 'Pluriformea']), 'Pluriformea sister to Ichthyosporea', _T15),
    ('group', ('Choanozoa', ['Metazoa', 'Choanoflagellata']), 'choanoflagellates = sister of animals',
     'King et al. 2008 Nature 451:783; Brunet & King 2017 Dev Cell 43:124'),
    ('group', ('Filozoa', ['Filasterea', 'Choanozoa']), 'Filasterea sister to Choanozoa', _T15),
    ('group', ('Holozoa', ['Teretosporea', 'Filozoa']), 'Holozoa', _T15),
    ('group', ('Holomycota', ['Fontalb', 'Spizpun']), 'Fonticula (nucleariids) + Fungi', 'Brown et al. 2009 Mol Biol Evol 26:2699'),
    ('move', ('Tentaculata', 'Metazoa'),
     'NCBI Eumetazoa (ctenophores + ParaHoxozoa) is not supported; sponge-first vs ctenophore-first contested -> '
     'Porifera / Ctenophora / ParaHoxozoa polytomy (as round 1)', 'Simion et al. 2017 Curr Biol 27:958; Schultz et al. 2023 Nature 618:110'),
    ('rename', ('Eumetazoa', 'ParaHoxozoa'), 'remaining clade = Placozoa + Cnidaria + Bilateria (order unresolved)',
     'Ryan et al. 2010 EvoDevo 1:9; Laumer et al. 2019 eLife 8:e36278'),
    ('group', ('Calcarea_Homoscleromorpha', ['Sycocil', 'Oscalob']), 'Calcarea sister to Homoscleromorpha',
     'Simion et al. 2017; Whelan et al. 2015 PNAS 112:5773'),
    ('group', ('Silicea', ['Aphrvas', 'Heteroscleromorpha']), 'Hexactinellida + Demospongiae', 'Simion et al. 2017; Whelan et al. 2015'),
]
PANELS['R1']['expect'] = ['Holozoa', 'Choanozoa', 'Metazoa', 'Choanoflagellata', 'Teretosporea', 'ParaHoxozoa']
PANELS['R1']['tree_notes'] = [
    'Metazoa root kept as Porifera / Ctenophora / ParaHoxozoa polytomy (ctenophore-first: Schultz et al. 2023; run variants if needed).',
    'Choanoflagellata internal order (NCBI Salpingoecidae polytomy) and Filasterea internal order kept from NCBI.',
    'Teretosporea: Abeoforma and Pirum EukProt gene sets unusable (BUSCO C < 3); Ichthyophonus (C 54) and Tunicaraptor (C 56) dropped.']
PANELS['R1']['events'] = [
    dict(id='R1_metazoa', label='animals', kind='origin_complex', M='Metazoa',
         **{'from': 'unicellular', 'to': 'complex'}, ref='Knoll 2011; Brunet & King 2017 Dev Cell 43:124; Sebe-Pedros et al. 2017 Nat Rev Genet 18:498',
         note='P (Choanozoa) may have been facultatively colonial (simple): alternative coding simple->complex'),
    dict(id='R1_salpingoeca', label='Salpingoeca rosette/chain colonies', kind='origin', M='Salpingoeca', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Dayel et al. 2011 Dev Biol 357:73; Richter et al. 2018 eLife 7:e34226',
         note='facultative; exploratory tip-level coding'),
    dict(id='R1_choanoeca', label='Choanoeca cup colonies', kind='origin', M='Choanoeca', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Brunet et al. 2019 Science 366:326', note='facultative; transcriptomes'),
    dict(id='R1_barroeca', label='Barroeca rosettes', kind='origin', M='Barrmon', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Hake et al. 2021 mBio 12:e0182821'),
    dict(id='R1_codosiga', label='Codosiga colonial clusters', kind='origin', M='Codohol', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Richter et al. 2018 eLife 7:e34226', note='transcriptome'),
    dict(id='R1_capsaspora', label='Capsaspora aggregation', kind='origin', M='Capsowc', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Sebe-Pedros et al. 2013 eLife 2:e01287', note='aggregative, facultative'),
    dict(id='R1_teretosporea', label='coenocytic development (Ichthyosporea + Corallochytrium)', kind='origin', M='Teretosporea', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Dudin et al. 2019 eLife 8:e49801; Kozyczkowska et al. 2021 PNAS 118:e2025504118',
         note='grey zone (coenocyte -> cellularisation); Syssomonas coded unicellular inside'),
]

# ----------------------------------------------------------------------------------------------------- R8 Cyanobacteria
PANELS['R8'] = dict(title='Cyanobacteria', kingdom='prok', busco_lineage='bacteria_odb10', edits=[], events=[])
_S11 = 'Schirrmeister et al. 2011 BMC Evol Biol 11:45; Schirrmeister et al. 2013 PNAS 110:1791'
_SH13 = 'Shih et al. 2013 PNAS 110:1053'
PANELS['R8']['edits'] = [
    ('group', ('Cyanophyceae_excl_Gloeobacter', ['Colechth', 'Gloelit', 'Leptbor', 'Pseubic', 'Spirsub', 'Acaryochloridales',
                                                  'Chroococcidiopsidales', 'Nostocales', 'Oscillatoriophycideae', 'Pleurocapsales',
                                                  'Synechococcales']),
     'Gloeobacter sister to all other cyanobacteria', _SH13 + '; ' + _S11),
    ('group', ('Nostocales_Chroococcidiopsis', ['Nostocales', 'Chroococcidiopsidales']),
     'Chroococcidiopsidales sister to heterocystous Nostocales', _SH13 + '; ' + _S11),
    ('group', ('Branching_Nostocales', ['Fiscthe', 'Chlofri']), 'true-branching Fischerella + Chlorogloeopsis form a clade',
     _SH13 + '; Gugger & Hoffmann 2004 Int J Syst Evol Microbiol 54:349'),
    ('group', ('Melainabacteria_Cyanophyceae', ['Vampchl', 'Cyanophyceae']), 'Vampirovibrionia (Melainabacteria) sister to '
     'oxygenic cyanobacteria; Sericytochromatia outside', 'Soo et al. 2017 Science 355:1436'),
]
PANELS['R8']['expect'] = ['Nostocales', 'Nostocales_Chroococcidiopsis', 'Cyanophyceae_excl_Gloeobacter', 'Branching_Nostocales']
PANELS['R8']['tree_notes'] = [
    'Backbone below Cyanophyceae_excl_Gloeobacter left as a polytomy (orders per NCBI / Strunecky et al. 2023); internal '
    'Nostocales order (except the branching clade) unresolved.',
    'Ancestral cyanobacterial state is debated (Schirrmeister et al. 2011, 2013: early filamentous origin with reversals; '
    'Hammerschmidt et al. 2021: stepwise). Filament origins are therefore coded as tip/clade events with a sensitivity coding.']
PANELS['R8']['events'] = [
    dict(id='R8_heterocyst', label='heterocystous Nostocales (D)', kind='nested', M='Nostocales',
         **{'from': 'simple', 'to': 'complex'}, ref='Herrero et al. 2016 FEMS Microbiol Rev 40:831; ' + _S11,
         note='complex = D only; sister Chroococcidiopsidales unicellular (likely reversal, Schirrmeister 2011)'),
    dict(id='R8_branching', label='true branching / multiseriate (T, exploratory)', kind='T_event', M='Branching_Nostocales', conf='grey',
         **{'from': 'complex', 'to': 'complex'}, ref=_SH13 + '; Herrero et al. 2016', note='D=1 -> D=1,T=1'),
    dict(id='R8_oscillatoriales', label='Oscillatoriales filaments', kind='origin', M='Oscillatoriales', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_S11, note='sister Chroococcales unicellular'),
    dict(id='R8_leptolyngbya', label='Leptolyngbya filaments', kind='origin', M='Leptbor', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_S11, note='tip in backbone polytomy'),
    dict(id='R8_pseudanabaena', label='Pseudanabaena filaments', kind='origin', M='Pseubic', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_S11, note='tip in backbone polytomy'),
    dict(id='R8_coleofasciculus', label='Coleofasciculus filaments', kind='origin', M='Colechth', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_S11, note='tip in backbone polytomy'),
    dict(id='R8_spirulina', label='Spirulina filaments', kind='origin', M='Spirsub', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref=_S11, note='tip in backbone polytomy'),
    dict(id='R8_hyella', label='Hyella pseudofilaments (Pleurocapsales)', kind='origin', M='Hyelpat', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Rippka et al. 1979 J Gen Microbiol 111:1'),
    dict(id='R8_hydrococcus', label='Hydrococcus pseudoparenchymatous colonies (Pleurocapsales)', kind='origin', M='Hydrriv', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Rippka et al. 1979 J Gen Microbiol 111:1'),
    dict(id='R8_sens_early_origin', label='single early filamentous origin (alternative coding)', kind='sensitivity',
         M='Cyanophyceae_excl_Gloeobacter', conf='grey', **{'from': 'unicellular', 'to': 'simple'}, ref=_S11,
         note='implies reversals in unicellular orders'),
]

# ------------------------------------------------------------------------------------------------------- R7 Amoebozoa
PANELS['R7'] = dict(title='Amoebozoa', busco_lineage='eukaryota_odb10', edits=[], events=[])
PANELS['R7']['edits'] = [
    ('group', ('Holozoa', ['Capsowc', 'Monobre']), 'Holozoa (outgroup structure)', 'Torruella et al. 2015 Curr Biol 25:2404'),
    ('group', ('Myxogastria_Ceratiomyxa', ['Physpol', 'Cerafru']), 'Ceratiomyxa sister to Myxogastria',
     'Fiore-Donno et al. 2010 Protist 161:55; Kang et al. 2017 MBE 34:2258'),
]
PANELS['R7']['expect'] = ['Dictyostelia', 'Eumycetozoa', 'Myxogastria_Ceratiomyxa', 'Tubulinea']
PANELS['R7']['tree_notes'] = [
    'Physarum is a transcriptome: the Schaap et al. 2015 genome annotation is not distributed through NCBI/UniProt/EukProt/LukProt '
    '(checked 2026-10-29); most non-dictyostelid Amoebozoa (protosteliids, Tubulinea, Ceratiomyxa) also have transcriptome data only.',
    'Dictyostelid groups 1-4 sampled (Cavenderia; Acytostelium + Heterostelium; Tieghemostelium; Dictyostelium + Polysphondylium violaceum).']
PANELS['R7']['events'] = [
    dict(id='R7_dictyostelia', label='dictyostelid aggregative fruiting', kind='origin', M='Dictyostelia',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Kin & Schaap 2021 Genes 12:487; Gloeckner et al. 2016 Nat Commun 7:12085',
         note='aggregative, facultative; Acytostelium lacks stalk-cell death (contrast)'),
    dict(id='R7_myxogastria', label='plasmodial slime moulds (Physarum + Ceratiomyxa)', kind='origin', M='Myxogastria_Ceratiomyxa', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Schaap et al. 2015 Genome Biol Evol 8:109', note='coenocytic plasmodium; transcriptomes'),
    dict(id='R7_sens_macromycetozoa', label='single origin of fruiting in Macromycetozoa (alternative coding)', kind='sensitivity',
         M='Eumycetozoa', conf='grey', **{'from': 'unicellular', 'to': 'simple'}, ref='Fiore-Donno et al. 2010 Protist 161:55'),
    dict(id='R7_copromyxa', label='Copromyxa sorocarpic aggregation', kind='origin', M='Coprpro', conf='grey',
         **{'from': 'unicellular', 'to': 'simple'}, ref='Brown et al. 2011 Protist 162:277', note='tip; transcriptome'),
]
