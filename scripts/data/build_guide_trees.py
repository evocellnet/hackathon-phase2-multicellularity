#!/usr/bin/env python3
"""Step 4: guide species trees for Ancestra.

1. NCBI taxonomy tree with Ancestra's own builder:
      python -m ancestra.speciestree ncbi --proteomes <tmp>/<code>_<taxid>.fa ... --taxdump data/raw/taxdump
             --leaf-naming qfo --out <raw.nwk>
   ('qfo' naming reads the taxid from the file stem '<stem>_<taxid>', so symlinks <code>_<taxid>.fa -> final
   FASTA are made in a temporary directory; leaves are renamed back to <code> = the real file stem).
2. Node names sanitised to [A-Za-z0-9_.-] (spaces -> '_', quotes removed).
3. Curated edits (EUK_EDITS / PROK_EDITS below; each with rationale + reference) -> current consensus.
   Operations: group (new clade from existing siblings), move, dissolve (polytomy), rename.
4. Unary nodes collapsed; children sorted alphabetically by (sub)tree label for a canonical Newick
   (Ancestra's merge depends slightly on child order).
5. Checks: leaves == file stems of data/proteomes/<kingdom>/*.fa; internal names unique and disjoint
   from leaves; tree rooted (root has >= 2 children); a list of expected clades is verified; the result
   is parsed back with ancestra.newick.
Outputs: data/guide_tree_euk.nwk, data/guide_tree_prok.nwk, data/guide_tree_ncbi_{euk,prok}.nwk (unedited,
sanitised), data/guide_tree_edits.md (edit log with references).
"""
import glob, os, re, subprocess, sys, tempfile
import pandas as pd
from ete3 import Tree

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
PY = f'{P}/conda_envs/ancestra/bin/python'

# ----------------------------------------------------------------------------------------------- edits
# (op, args, rationale, reference). Node names refer to the sanitised NCBI names or to clades created earlier.
EUK_EDITS = [
    # ---- dissolve non-clades ('incertae sedis' bins of NCBI)
    ('rename', ('Opisthokonta_incertae_sedis', 'Pluriformea'),
     'NCBI bin holds exactly Corallochytrium + Syssomonas, which form the clade Pluriformea',
     'Hehenberger et al. 2017 Curr Biol 27:2043'),
    ('dissolve', ('Fungi_incertae_sedis',), 'NCBI bin is not a clade; replaced by the edits below', '-'),
    ('dissolve', ('CRuMs_incertae_sedis',), 'not a clade; internal CRuMs order left unresolved', 'Brown et al. 2018 GBE 10:427'),
    ('dissolve', ('Eukaryota_incertae_sedis',), 'not a clade (Telonema, Palpitomonas placed below)', '-'),
    # ---- Opisthokonta: Holozoa vs Holomycota
    ('group', ('Teretosporea', ['Ichthyosporea', 'Pluriformea']),
     'Pluriformea sister to Ichthyosporea', 'Torruella et al. 2015 Curr Biol 25:2404; Hehenberger et al. 2017'),
    ('group', ('Choanozoa', ['Metazoa', 'Choanoflagellata']), 'choanoflagellates = sister of animals',
     'King et al. 2008 Nature 451:783; Brunet & King 2017 Dev Cell 43:124'),
    ('group', ('Filozoa', ['Filasterea', 'Choanozoa']), 'Filasterea sister to Choanozoa',
     'Torruella et al. 2015; Hehenberger et al. 2017'),
    ('group', ('Holozoa', ['Teretosporea', 'Filozoa']), 'Holozoa', 'Torruella et al. 2015'),
    ('group', ('Aphelida_Eumycota', ['Paratri', 'chytridiomyceta', 'Allomac', 'Mucorineae', 'Morthyg', 'Basimer',
                                     'Coemumb', 'Dikarya']),
     'Aphelida (Paraphelidium) sister to Fungi s.s.; Rozellomycota branch before them',
     'Torruella et al. 2018 Commun Biol 1:231'),
    ('group', ('Eumycota', ['chytridiomyceta', 'Allomac', 'Mucorineae', 'Morthyg', 'Basimer', 'Coemumb', 'Dikarya']),
     'Fungi sensu stricto (zoosporic + terrestrial fungi)', 'Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653'),
    ('group', ('Mucoromycota', ['Mucorineae', 'Morthyg']), 'Mortierellomycotina + Mucoromycotina',
     'Spatafora et al. 2016'),
    ('group', ('Zoopagomycota', ['Basimer', 'Coemumb']), 'Entomophthoromycotina (Basidiobolus) + Kickxellomycotina',
     'Spatafora et al. 2016'),
    ('group', ('Mucoromycota_Dikarya', ['Mucoromycota', 'Dikarya']), 'Mucoromycota sister to Dikarya',
     'Spatafora et al. 2016; Li et al. 2021'),
    ('group', ('Amastigomycota', ['Zoopagomycota', 'Mucoromycota_Dikarya']),
     'terrestrial (non-flagellated) fungi; Zoopagomycota sister to Mucoromycota+Dikarya', 'Spatafora et al. 2016; Li et al. 2021'),
    ('group', ('Blastocladiomycota_Amastigomycota', ['Allomac', 'Amastigomycota']),
     'Blastocladiomycota sister to terrestrial fungi (Chytridiomyceta branch first); M of the hyphal transition',
     'James et al. 2006 Nature 443:818; Chang et al. 2015 GBE 7:1590; Li et al. 2021'),
    ('group', ('Holomycota', ['Rotosphaerida', 'Fungi']), 'Holomycota = Fonticula/nucleariids + Fungi',
     'Brown et al. 2009 PLoS One 4:e5656; Torruella et al. 2015'),
    # Basidiomycota
    ('group', ('Wallemiomycotina_Agaricomycotina', ['Wallmel', 'Agaricomycotina']),
     'Wallemiomycotina sister to Agaricomycotina', 'Zhao et al. 2017 Fungal Divers 84:43; He et al. 2019 Fungal Divers 99:105'),
    ('group', ('Ustilaginomycotina_Agaricomycotina', ['Ustimay', 'Wallemiomycotina_Agaricomycotina']),
     'Pucciniomycotina sister to Ustilaginomycotina + Agaricomycotina', 'He et al. 2019; Li et al. 2021'),
    ('group', ('Agaricomycetes_Dacrymycetes', ['Agaricales', 'Dacrpri']),
     'Dacrymycetes sister to Agaricomycetes; Tremellomycetes sister to both', 'Nagy et al. 2016 MBE 33:959; He et al. 2019'),
    # Ascomycota: Taphrinomycotina order left unresolved; Pezizomycotina structure = NCBI
    # ---- Metazoa
    ('move', ('Mnemlei', 'Metazoa'),
     'NCBI Eumetazoa (ctenophores + ParaHoxozoa, excluding sponges) is not supported; the sponge-first vs '
     'ctenophore-first question is contested, so Porifera / Ctenophora / ParaHoxozoa form a polytomy',
     'Simion et al. 2017 Curr Biol 27:958; Whelan et al. 2017 Nat Ecol Evol 1:1737; Schultz et al. 2023 Nature 618:110'),
    ('rename', ('Eumetazoa', 'ParaHoxozoa'), 'remaining clade = Placozoa + Cnidaria + Bilateria (order unresolved)',
     'Ryan et al. 2010 EvoDevo 1:9; Laumer et al. 2019 eLife 8:e36278'),
    ('group', ('Calcarea_Homoscleromorpha', ['Sycocil', 'Oscapea']), 'Calcarea sister to Homoscleromorpha',
     'Simion et al. 2017; Whelan et al. 2015 PNAS 112:5773'),
    # ---- Amoebozoa
    ('group', ('Myxogastria_Ceratiomyxa', ['Physpol', 'Cerafru']), 'Ceratiomyxa sister to Myxogastria',
     'Fiore-Donno et al. 2010 Protist 161:55; Kang et al. 2017 MBE 34:2258'),
    # ---- Amorphea / CRuMs
    ('group', ('Apusomonadida_Opisthokonta', ['Thectra', 'Opisthokonta']), 'apusomonads sister to opisthokonts',
     'Brown et al. 2013 Proc R Soc B 280:20131755; Burki et al. 2020 TREE 35:43'),
    ('group', ('Obazoa', ['Pygsbif', 'Apusomonadida_Opisthokonta']), 'Breviatea + Apusomonadida + Opisthokonta',
     'Brown et al. 2013'),
    ('group', ('Amorphea', ['Amoebozoa', 'Obazoa']), 'Amorphea', 'Adl et al. 2012 J Eukaryot Microbiol 59:429; Burki et al. 2020'),
    ('group', ('Podiata', ['CRuMs', 'Amorphea']),
     'CRuMs sister to Amorphea (not nested in Obazoa; the PI brief said "sister of Obazoa" - the published '
     'phylogenomic result is sister of Amorphea)', 'Brown et al. 2018 GBE 10:427; Burki et al. 2020'),
    # ---- Archaeplastida
    ('dissolve', ('Bangiophyceae',),
     'NCBI node lumps Bangiales with Porphyridiophyceae and Cyanidiophyceae; regrouped below', '-'),
    ('group', ('Cyanidiophyceae', ['Cyanmer', 'Galdsul']), 'Cyanidiophyceae sister to all other red algae',
     'Yoon et al. 2006 J Phycol 42:482; Munoz-Gomez et al. 2017 Curr Biol 27:1677'),
    ('group', ('Proteorhodophytina', ['Porppur', 'Rhodmar']), 'Porphyridiophyceae + Stylonematophyceae',
     'Munoz-Gomez et al. 2017 Curr Biol 27:1677'),
    ('group', ('Eurhodophytina', ['Bangiaceae', 'Rhodymeniophycidae']), 'Bangiophyceae + Florideophyceae',
     'Yoon et al. 2006; Munoz-Gomez et al. 2017'),
    ('group', ('Rhodophytina', ['Proteorhodophytina', 'Eurhodophytina']), 'mesophilic red algae',
     'Munoz-Gomez et al. 2017'),
    ('group', ('Rhodophyta_Rhodelphidia', ['Rhodophyta', 'Rhodemar']), 'Rhodelphis sister to red algae',
     'Gawryluk et al. 2019 Nature 572:240'),
    ('group', ('Chlorophyta_Streptophyta', ['Chlorophyta', 'Streptophyta']),
     'Prasinodermophyta (Prasinoderma) sister to Chlorophyta + Streptophyta', 'Li L. et al. 2020 Nat Ecol Evol 4:1220'),
    ('group', ('Archaeplastida', ['Cyanpar', 'Rhodophyta_Rhodelphidia', 'Viridiplantae']),
     'Archaeplastida; order of glaucophytes / red / green left unresolved', 'Burki et al. 2020; Strassert et al. 2021 Nat Commun 12:1879'),
    ('group', ('UTC_clade', ['CS_clade', 'Trebouxiophyceae', 'Ulvophyceae']),
     'Chlorodendrophyceae (Tetraselmis) sister to Ulvophyceae+Trebouxiophyceae+Chlorophyceae; UTC order unresolved',
     'Leliaert et al. 2012 Crit Rev Plant Sci 31:1'),
    # Volvocales (Lindsey et al. 2021: two origins of colonial multicellularity, Goniaceae paraphyletic)
    ('move', ('Gonipec', 'Chlamydomonadales'), 'NCBI places Gonium in Volvocaceae', 'Nozaki et al. 2000 MPE 17:256'),
    ('dissolve', ('Volvox',), 'genus Volvox is polyphyletic; Volvocaceae internal order left unresolved',
     'Herron & Michod 2008 Evolution 62:436; Lindsey et al. 2021 BMC Biol 19:182'),
    ('group', ('Astrephomene_Volvocaceae', ['Astrgub', 'Volvocaceae']), 'Astrephomene sister to Volvocaceae',
     'Lindsey et al. 2021 BMC Biol 19:182 (MLBS 98)'),
    ('group', ('Gonium_Volvocaceae', ['Gonipec', 'Astrephomene_Volvocaceae']),
     'Gonium sister to Astrephomene + Volvocaceae (Goniaceae paraphyletic); = M of the 2nd volvocine origin',
     'Lindsey et al. 2021'),
    ('group', ('TGV_clade', ['Tetrsoc', 'Gonium_Volvocaceae']),
     'Tetrabaenaceae + Goniaceae + Volvocaceae: a clade among SAMPLED taxa; Lindsey et al. 2021 place an '
     'unsampled unicellular lineage as sister of Tetrabaenaceae, so the ancestor at this node is inferred '
     'unicellular and multicellularity arose twice (Tetrabaena tip; Gonium_Volvocaceae)', 'Lindsey et al. 2021'),
    ('group', ('Chlamydomonas_TGV', ['Chlamydomonas', 'TGV_clade']),
     'core Reinhardtinia: Chlamydomonas (reinhardtii/incerta/schloesseri) sister to the colonial lineages',
     'Craig et al. 2021 Plant Cell 33:1016; Lindsey et al. 2021'),
    ('group', ('Reinhardtinia', ['Edapdeb', 'Chlamydomonas_TGV']), 'Edaphochlamys outside core Reinhardtinia',
     'Craig et al. 2021; Lindsey et al. 2021'),
    # Streptophyta
    ('group', ('Mesostigmatophyceae_Chlorokybophyceae', ['Mesovir', 'Chloatm']),
     'Mesostigma + Chlorokybus form the earliest streptophyte clade', 'Wang S. et al. 2020 Nat Plants 6:95; Lemieux et al. 2007 BMC Biol 5:2'),
    ('group', ('Klebsormidiophyceae_Phragmoplastophyta', ['Klebnit', 'Streptophytina']),
     'grade: Klebsormidiophyceae, then Phragmoplastophyta (NCBI Streptophytina = Charophyceae+Coleochaetophyceae+'
     'Zygnematophyceae+Embryophyta)', 'Wickett et al. 2014 PNAS 111:E4859; One Thousand Plant Transcriptomes 2019 Nature 574:679'),
    ('group', ('Anydrophyta', ['Zygnematophyceae', 'Embryophyta']), 'Zygnematophyceae sister to land plants',
     'Wickett et al. 2014; 1KP 2019; Cheng et al. 2019 Cell 179:1057'),
    ('group', ('Coleochaetophyceae_Anydrophyta', ['Colescu', 'Anydrophyta']), 'Coleochaetophyceae sister to Anydrophyta',
     'Wickett et al. 2014; 1KP 2019'),
    ('group', ('Setaphyta', ['Marcpol', 'Physpat']), 'liverworts + mosses (bryophyte monophyly)',
     'Puttick et al. 2018 Curr Biol 28:733; Harris et al. 2020 Curr Biol 30:2001'),
    # ---- TSAR
    ('group', ('Phaeophyceae_Schizocladia', ['Phaeophyceae', 'Schiisc']), 'Schizocladia sister to brown algae',
     'Kawai et al. 2003 J Phycol 39:556; Denoeud et al. 2024 Cell 187:6943 (Phaeoexplorer)'),
    ('group', ('Fucophycidae', ['Ectosil', 'Sacclat', 'Fucuser']), 'brown algal crown radiation (Fucophycidae)',
     'Silberfeld et al. 2014 Cryptogam Algol 35:117; Denoeud et al. 2024'),
    ('group', ('Dictyotophycidae_Fucophycidae', ['Dictdic', 'Fucophycidae']),
     'Discosporangiales sister to the rest of Phaeophyceae', 'Silberfeld et al. 2014; Denoeud et al. 2024'),
    ('group', ('Gyrista', ['Ochrophyta', 'Oomycota']), 'Ochrophyta + Pseudofungi (oomycetes)',
     'Derelle et al. 2016 MBE 33:2890; Burki et al. 2020'),
    ('rename', ('Colpodellida', 'Chrompodellida'), 'NCBI node = Vitrella + Chromera (chrompodellids)', 'Janouskovec et al. 2015 PNAS 112:10200'),
    ('group', ('Apicomplexa_Chrompodellida', ['Apicomplexa', 'Chrompodellida|Vitrbra|Chrovel']), 'chrompodellids sister to Apicomplexa',
     'Janouskovec et al. 2015'),
    ('group', ('Dinozoa', ['Perkmar', 'Symbmic']), 'Perkinsozoa + Dinoflagellata', 'Janouskovec et al. 2015; Burki et al. 2020'),
    ('group', ('Myzozoa', ['Apicomplexa_Chrompodellida', 'Dinozoa|Perkmar|Symbmic']), 'Myzozoa sister to ciliates', 'Janouskovec et al. 2015'),
    ('group', ('TSAR', ['Telosub', 'Sar']), 'Telonemia sister to SAR', 'Strassert et al. 2019 MBE 36:757; Burki et al. 2020'),
    ('group', ('Cryptista', ['Palpbil', 'Guilthe']), 'Palpitomonas sister to cryptophytes (Cryptista)',
     'Yabuki et al. 2014 Sci Rep 4:4641; Burki et al. 2020'),
    ('group', ('Diaphoretickes', ['Archaeplastida', 'TSAR', 'Haptista', 'Cryptista|Guilthe|Palpbil']),
     'Diaphoretickes; positions of Haptista and Cryptista unresolved (polytomy)', 'Burki et al. 2020'),
    # ---- Discoba
    ('dissolve', ('Acrasida',),
     'NCBI groups Acrasis with Naegleriidae to the exclusion of Neovahlkampfia; relationships among Tetramitia '
     'families are not resolved in published phylogenies -> polytomy', 'Brown et al. 2012 Eur J Protistol 48:103; Panek et al. 2017 (Handbook of the Protists: Heterolobosea)'),
    ('group', ('Discicristata', ['Heterolobosea', 'Euglenozoa']), 'Heterolobosea + Euglenozoa; Jakobida sister',
     'Hampl et al. 2009 PNAS 106:3859; Burki et al. 2020'),
    # ---- root
    ('group', ('Bikonta', ['Diaphoretickes', 'Discoba', 'Metamonada|Parapyr|Preaxostyla', 'Gefioke']),
     'ROOT between Podiata (Amorphea + CRuMs) and all other eukaryotes (unikont/bikont-type root, PI brief). '
     'Order inside is a polytomy. Alternatives (Opimoda/Diphoda root with Metamonada + Malawimonadida on the '
     'Amorphea side, Derelle et al. 2015; root near Discoba, Al Jewari & Baldauf 2023) are untested variants; '
     'they affect only backbone branches, not the transition nodes', 'Stechmann & Cavalier-Smith 2002 Science 297:89; '
     'Derelle et al. 2015 PNAS 112:E693; Burki et al. 2020; Al Jewari & Baldauf 2023 Sci Adv 9:eade4973'),
]

PROK_EDITS = [
    ('group', ('Cyanophyceae_excl_Gloeobacter', ['Nostocales', 'Oscillatoriophycideae', 'Leptbor', 'Pseubic', 'Stancya',
                                                  'Chrothe', 'Synechococcales', 'Acaryochloridales', 'Gloelit']),
     'Gloeobacter sister to all other cyanobacteria', 'Shih et al. 2013 PNAS 110:1053; Schirrmeister et al. 2013 PNAS 110:1791'),
    ('group', ('Nostocales_Chroococcidiopsis', ['Nostocales', 'Chrothe']),
     'Chroococcidiopsidales sister to heterocystous Nostocales', 'Shih et al. 2013; Schirrmeister et al. 2011 BMC Evol Biol 11:45'),
    ('group', ('Cable_bacteria', ['Elecsca', 'ElecspV4']), 'Electrothrix + Electronema form the cable-bacteria clade',
     'Kjeldsen et al. 2019 PNAS 116:19116'),
]

# expected clades (leaf sets) that must hold after editing
EUK_EXPECT = {
    'Holozoa': None, 'Holomycota': None, 'Amorphea': None, 'Diaphoretickes': None, 'Discoba': None,
    'Choanozoa': None, 'Dictyostelia': None, 'Phaeophyceae': None, 'Anydrophyta': None,
    'Eurhodophytina': None, 'TGV_clade': None, 'Gonium_Volvocaceae': None, 'Evosea': None,
}


def sanitise(n):
    n = n.replace("'", '').strip()
    return re.sub(r'[^A-Za-z0-9_.\-]', '_', n)


def find(t, name):
    hits = [n for n in t.traverse() if n.name == name]
    if len(hits) != 1:
        raise KeyError(f'node "{name}": {len(hits)} matches')
    return hits[0]


def exists(t, name):
    return any(n.name == name for n in t.traverse())


def resolve(t, spec):
    """'A|B|C' -> first name that exists in the tree (None if none)."""
    for x in spec.split('|'):
        if exists(t, x):
            return x
    return None


def apply(t, edits, log):
    """Apply edits. Robust to the taxon set: a member spec 'A|B' uses the first existing node; rename/dissolve of
    a node that is absent (collapsed because of sampling) and groups with < 2 resolvable members are skipped and
    logged as such."""
    for i, (op, args, why, ref) in enumerate(edits, 1):
        if op == 'rename':
            if not exists(t, args[0]):
                log.append((i, 'skipped:rename', args, 'node absent in this taxon set (unary, collapsed)', ref))
                continue
            find(t, args[0]).name = args[1]
        elif op == 'dissolve':
            if not exists(t, args[0]):
                log.append((i, 'skipped:dissolve', args, 'node absent in this taxon set (unary, collapsed)', ref))
                continue
            nd = find(t, args[0])
            assert not nd.is_leaf() and not nd.is_root()
            nd.delete(prevent_nondicotomic=False, preserve_branch_length=False)
        elif op == 'move':
            nd, dest = find(t, args[0]), find(t, args[1])
            nd.detach()
            dest.add_child(nd)
        elif op == 'group':
            name, members = args
            res = [resolve(t, m) for m in members]
            missing = [m for m, r in zip(members, res) if r is None]
            res = [r for r in res if r is not None]
            if len(res) < 2:
                log.append((i, 'skipped:group', (name, res), f'fewer than 2 members present ({res})', ref))
                continue
            nodes = [find(t, m) for m in res]
            parents = {id(n.up) for n in nodes}
            if len(parents) != 1:
                raise ValueError(f'edit {i} group {name}: members not siblings: '
                                 + ', '.join(f'{n.name}<{n.up.name}' for n in nodes))
            par = nodes[0].up
            if len(par.children) == len(nodes):
                # the group would equal its parent: rename the parent instead of creating a unary chain
                log.append((i, 'group=parent', (name, res), f'{why} [equals existing node {par.name}; renamed]', ref))
                par.name = name
                continue
            new = par.add_child(name=name)
            for n in nodes:
                n.detach()
                new.add_child(n)
            args = (name, res)
            if missing:
                why = why + f' [absent here: {", ".join(missing)}]'
        else:
            raise ValueError(op)
        log.append((i, op, args, why, ref))


def collapse_unary(t, log):
    for nd in list(t.traverse('postorder')):
        if not nd.is_leaf() and len(nd.children) == 1:
            log.append(('-', 'collapse_unary', (nd.name, nd.children[0].name), 'unary node removed', '-'))
            if nd.is_root():
                continue
            nd.delete(prevent_nondicotomic=False)
    while len(t.children) == 1:
        t = t.children[0]
        t.up = None
    return t


def canonical(t):
    key = {}
    for nd in t.traverse('postorder'):
        key[id(nd)] = nd.name if nd.is_leaf() else min(key[id(c)] for c in nd.children)
    for nd in t.traverse():
        if nd.children:
            nd.children.sort(key=lambda c: (c.name if not c.is_leaf() else '', key[id(c)]))
    return t


def newick(nd):
    if nd.is_leaf():
        return nd.name
    return '(' + ','.join(newick(c) for c in nd.children) + ')' + nd.name


def ncbi_raw(kingdom, table):
    d = pd.read_csv(f'{P}/data/{table}', sep='\t', dtype=str)
    with tempfile.TemporaryDirectory(dir=f'{P}/work' if os.path.isdir(f'{P}/work') else None) as td:
        for c, tx in zip(d['code'], d['ncbi_taxid']):
            src = f'{P}/data/proteomes/{kingdom}/{c}.fa'
            if not os.path.exists(src):
                sys.exit(f'missing proteome {src}')
            os.symlink(src, f'{td}/{c}_{int(tx)}.fa')
        out = f'{td}/raw.nwk'
        cmd = [PY, '-m', 'ancestra.speciestree', 'ncbi', '--proteomes', td, '--taxdump', f'{P}/data/raw/taxdump',
               '--leaf-naming', 'qfo', '--out', out]
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=P)
        if r.returncode:
            sys.exit(r.stderr)
        print(' '.join(cmd[1:]).replace(td, '<tmp>'), '->', r.stdout.strip())
        t = Tree(open(out).read(), format=1, quoted_node_names=True)
    for nd in t.traverse():
        nd.name = sanitise(nd.name)
        if nd.is_leaf():
            nd.name = re.sub(r'_\d+$', '', nd.name)
    return t, set(d['code'])


def check(t, codes, kingdom, expect=()):
    leaves = [n.name for n in t.iter_leaves()]
    assert len(leaves) == len(set(leaves)), 'duplicate leaves'
    stems = {os.path.basename(f)[:-3] for f in glob.glob(f'{P}/data/proteomes/{kingdom}/*.fa')}
    assert set(leaves) == codes == stems, (set(leaves) ^ stems)
    internal = [n.name for n in t.traverse() if not n.is_leaf()]
    assert all(internal), 'unnamed internal node'
    assert len(internal) == len(set(internal)), [x for x in internal if internal.count(x) > 1]
    assert not set(internal) & set(leaves)
    assert all(re.fullmatch(r'[A-Za-z0-9_.\-]+', x) for x in internal + leaves)
    assert len(t.children) >= 2
    for name in expect:
        find(t, name)
    return len(leaves), len(internal), sum(1 for n in t.traverse() if len(n.children) > 2)


def main():
    os.makedirs(f'{P}/work', exist_ok=True)
    md = ['# Guide-tree edits', '',
          'Built by `scripts/data/build_guide_trees.py`. Base: NCBI taxonomy (taxdump in `data/raw/taxdump`, '
          'downloaded 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on symlinks '
          '`<code>_<taxid>.fa`; leaves renamed back to the proteome file stems (`<code>`); node names sanitised '
          '(spaces to `_`, quotes removed). The unedited trees are `data/guide_tree_ncbi_{euk,prok}.nwk`.',
          '', 'Rules: resolve a polytomy or move a taxon only where recent phylogenomic work agrees; otherwise keep '
          'the multifurcation (Ancestra merges k-way). Every internal node has a unique name: NCBI names are kept '
          'where the leaf set matches; new clades get the literature name or an explicit compound name '
          '(e.g. `Gonium_Volvocaceae`). Children are sorted alphabetically (canonical Newick; Ancestra output '
          'depends slightly on child order). Unary nodes are collapsed (logged).', '']
    for kingdom, table, edits, expect, label in (
            ('euk', 'species_final.tsv', EUK_EDITS, EUK_EXPECT, 'Eukaryotes'),
            ('prok', 'species_final_prok.tsv', PROK_EDITS, (), 'Prokaryotes')):
        t, codes = ncbi_raw(kingdom, table)
        t = collapse_unary(t, [])
        open(f'{P}/data/guide_tree_ncbi_{kingdom}.nwk', 'w').write(newick(canonical(t.copy())) + ';\n')
        log = []
        apply(t, edits, log)
        t = collapse_unary(t, log)
        t = canonical(t)
        n_leaf, n_int, n_poly = check(t, codes, kingdom, expect)
        out = f'{P}/data/guide_tree_{kingdom}.nwk'
        open(out, 'w').write(newick(t) + ';\n')
        # parse back with Ancestra's reader
        r = subprocess.run([PY, '-c', f'from ancestra.newick import read_newick; t=read_newick("{out}"); '
                            f'print(len(t.leaves()))'], capture_output=True, text=True, cwd=P)
        assert r.stdout.strip() == str(n_leaf), r.stderr
        print(f'{kingdom}: leaves={n_leaf} internal={n_int} multifurcations={n_poly} edits={len(edits)} -> {out}')
        md += [f'## {label}', '',
               f'Root: `{t.name}` with children ' + ', '.join(f'`{c.name}`' for c in t.children) + '. '
               f'{n_leaf} leaves, {n_int} named internal nodes, {n_poly} multifurcating nodes.', '',
               '| # | operation | nodes | rationale | reference |', '|---|---|---|---|---|']
        for i, op, args, why, ref in log:
            a = (f'{args[0]} := ({", ".join(args[1])})' if op == 'group' else ' -> '.join(map(str, args)))
            md.append(f'| {i} | {op} | {a} | {why} | {ref} |')
        md.append('')
    md += ['## Changes after design review round 1 (reviews/round1/evo_biology.md, data_algorithms_compute.md)', '',
           '- **Transitions coded as state changes** (`data/transitions.tsv`, columns `from_state`, `to_state`, '
           '`analysis_tier`, `primary_eligible`; P = last ancestor in the prior state; `prior_state_ref_*` = nearest '
           'species in the prior state when the tree-sister of M is itself multicellular). States per species in '
           '`data/species_traits.tsv` (unicellular / simple / complex / grey; `reversal`).',
           '- **Reversals** (M->U): Saccharomycotina (Nagy et al. 2014 Nat Commun 5:4471), Schizosaccharomyces, '
           'Saitoella, Cryptococcus, Rhodotorula; Zygnematophyceae (Spirogloea, Mesotaenium, Closterium) nested in the '
           'streptophyte M1. Wallemia coded filamentous (simple), Ustilago/Taphrina dimorphic (simple).',
           '- **Volvocales: two origins** (Lindsey et al. 2021 BMC Biol 19:182, doi:10.1186/s12915-021-01087-0): '
           'Goniaceae paraphyletic, Astrephomene sister to Volvocaceae and Gonium sister to both (tree edits '
           '`Astrephomene_Volvocaceae`, `Gonium_Volvocaceae`); Tetrabaenaceae an independent origin whose unicellular '
           'sister is unsampled. M1 = Tetrsoc (tip), M2 = Gonium_Volvocaceae, both with P = TGV_clade (inferred '
           'unicellular); unicellular reference = Chlamydomonas (prior_state_ref).',
           '- **Streptophyta**: M1 = Klebsormidiophyceae_Phragmoplastophyta (U->S, P = Streptophyta); Embryophyta = '
           'nested S->X (P = Anydrophyta); Chara tip S->X (secondary).',
           '- **Rhodophyta**: M1 = Eurhodophytina (U->S); Florideophyceae (Rhodymeniophycidae node) nested S->X.',
           '- **Brown algae**: filaments predate Phaeophyceae (Schizocladia, Phaeothamnion, Tribonema filamentous): '
           'PX_clade = U->S (secondary, P = Ochrophyta polytomy); Phaeophyceae = S->X (primary).',
           '- **Dictyostelia** = aggregative simple (Knoll 2011 Annu Rev Earth Planet Sci 39:217); Physarum = grey.',
           '- **Grey / exploratory**: choanoflagellate colonies (tip-level, no internal M), Capsaspora aggregation, '
           'Physarum, Bryopsidales, Vaucheria/Xanthophyceae, oomycetes, Plasmodiophora, Labyrinthula. '
           '**Excluded from primary**: Ichthyosporea/Teretosporea (coenocytic state likely ancestral to Teretosporea), '
           'Guttulinopsis (BUSCO ~11 %, presence-only). **Dropped (no data)**: Sorogena, Sorodiplophrys.',
           '- **Fungal hyphae**: M = Blastocladiomycota + terrestrial fungi (edit `Blastocladiomycota_Amastigomycota`), '
           'P = Eumycota; Gonapodya = separate tip origin.',
           '- **Added taxa**: Paraphelidium tribonemae (EukProt EP00158, transcriptome; sister of Fungi, Torruella et al. '
           '2018). Not added (no annotated proteome in UniProt, EukProt v3, LukProt v1.5.1 or NCBI on 2026-10-28; only '
           'unannotated assemblies): Hyphochytrium catenoides (GCA_900088475.1), Developayella elegans (EukProt '
           '"not included"), Yamagishiella unicocca (GCA_003116995.1), Eudorina sp. (GCA_003117095.1), Interfilum (none).',
           '- **Memory cap (PI rule)**: eukaryote total <= 3.0 M proteins; optional taxa dropped in the written order '
           '(optional+background largest first, then other optional taxa that leave >= 2 species per transition side): '
           'Trichomonas, Tetraselmis, Symbiodinium, Emiliania, Chromera, Anaeramoeba, Mantamonas, Monocercomonoides, '
           'Palpitomonas, Reticulomyxa (`data/species_cap_drops.tsv`). Consequences for the tree: Metamonada = '
           'Paratrimastix only, Cryptista = Guillardia only, Dinozoa = Perkinsus only, Chrompodellida = Vitrella only '
           '(the corresponding edits are logged as skipped or resolved to the single remaining taxon).',
           '- **Transcriptomes**: `cd-hit -c 0.95 -G 0 -aS 0.9 -n 5` after longest-ORF-per-gene; sets with BUSCO D or '
           'F > 15 % flagged (`flag_DF15`). data_type audit: Euglena gracilis (EP00667) is a translated '
           'transcriptome (reprocessed as such); Sycon (EP00117, Compagen T-PEP) kept as genome gene models '
           '(isoforms collapsed by scgid).',
           '- **Canonical Newick**: children sorted alphabetically (internal-node name, then smallest leaf name).', '']
    md += ['## Known open issues', '',
           '- Eukaryote root: unikont/bikont-type root (Podiata vs Bikonta) as asked in the brief; alternative roots '
           '(Opimoda/Diphoda; Discoba) should be run as tree variants (review round 1, data/algorithms item 4).',
           '- Metazoa: Porifera / Ctenophora / ParaHoxozoa polytomy; run sponge-first and ctenophore-first variants.',
           '- Contested placements kept as in NCBI or as polytomies: Chromosphaera (within Ichthyosporea), '
           'Copromyxa (Tubulinea, Elardia), Telonemia/Haptista/Cryptista order, Rhizaria internal order, '
           'Choanoflagellata internal order (NCBI Salpingoecidae kept), Taphrinomycotina internal order, '
           'Bigyra monophyly, Ulvophyceae monophyly (Del Cortona et al. 2020 PNAS 117:2551).',
           '- Volvocales: Tetrabaena + (Gonium + (Astrephomene + Volvocaceae)) form a clade among the sampled taxa, '
           'but the unicellular sister of Tetrabaenaceae (Lindsey et al. 2021) is not sampled; transitions.tsv '
           'therefore codes two origins with the same parent node (TGV_clade, inferred unicellular).', '']
    open(f'{P}/data/guide_tree_edits.md', 'w').write('\n'.join(md))


if __name__ == '__main__':
    main()
