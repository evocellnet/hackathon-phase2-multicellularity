# R2 guide-tree edits (Fungi / Holomycota)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Opisthokonta` with children `Capsowc`, `Holomycota`. 40 leaves, 34 named internal nodes, 4 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | dissolve | Fungi_incertae_sedis | NCBI bin is not a clade; replaced by the groups below | - |
| 2 | group | Mucoromycota := (Mucoromycotina, Morthyg) | Mortierellomycotina + Mucoromycotina | Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653 |
| 3 | group | Zoopagomycota := (Basimer, Coemumb) | Entomophthoromycotina + Kickxellomycotina | Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653 |
| 4 | group | Mucoromycota_Dikarya := (Mucoromycota, Dikarya) | Mucoromycota sister to Dikarya | Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653 |
| 5 | group | Amastigomycota := (Zoopagomycota, Mucoromycota_Dikarya) | terrestrial (non-flagellated) fungi | Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653 |
| 6 | group | Blastocladiomycota_Amastigomycota := (Allomac, Amastigomycota) | Blastocladiomycota sister to terrestrial fungi; M of the hyphal transition | James et al. 2006 Nature 443:818; Chang et al. 2015 GBE 7:1590; Li et al. 2021 Curr Biol 31:1653 |
| 7 | group | Eumycota := (chytridiomyceta, Blastocladiomycota_Amastigomycota) | Fungi sensu stricto | Spatafora et al. 2016 Mycologia 108:1028; Li et al. 2021 Curr Biol 31:1653 |
| 8 | group | Aphelida_Eumycota := (Amoeocc, Eumycota) | Aphelida sister to Fungi s.s.; Rozellomycota branch first | Torruella et al. 2018 Commun Biol 1:231 |
| 9 | group | Holomycota := (Rotosphaerida, Fungi) | nucleariids + Fonticula sister to Fungi | Brown et al. 2009 Mol Biol Evol 26:2699; Torruella et al. 2015 Curr Biol 25:2404 |
| 10 | group | Yarrowia_Saccharomyces := (Yarrlip, Sacccer) | Lipomycetaceae branch first within Saccharomycotina | Shen et al. 2018 Cell 175:1533 |
| 11 | dissolve | Agaricomycetes_incertae_sedis | NCBI bin (Auriculariales, Sebacinales, Polyporales) is not a clade | - |
| 12 | group | Polyporales_Agaricales := (Tramcoc, Agaricales) | Polyporales + Agaricomycetidae form the late-diverging Agaricomycetes; Auriculariales and Sebacinales branch earlier (order unresolved) | Hibbett et al. 2007 Mycol Res 111:509; Varga et al. 2019 Nat Ecol Evol 3:668 |
| 13 | group | Agaricomycetes_Dacrymycetes := (Agaricomycetes, Dacrpri) | Dacrymycetes sister to Agaricomycetes; Tremellomycetes sister to both | Nagy et al. 2016 MBE 33:959; He et al. 2019 Fungal Divers 99:105 |
| 14 | group | Wallemiomycotina_Agaricomycotina := (Wallmel, Agaricomycotina) | Wallemiomycotina sister to Agaricomycotina | Zhao et al. 2017 Fungal Divers 84:43; He et al. 2019 |
| 15 | group | Ustilaginomycotina_Agaricomycotina := (Ustimay, Wallemiomycotina_Agaricomycotina) | Pucciniomycotina sister to Ustilaginomycotina + Agaricomycotina | He et al. 2019; Li et al. 2021 |

## Open issues / notes

- Polytomies kept: Pezizomycotina (Orbiliomycetes / Pezizomycetes / leotiomyceta), leotiomyceta (Eurotiomycetes, Dothideomycetes, sordariomyceta), Taphrinomycotina (Neolecta, Saitoella, Schizosaccharomyces, Taphrina), Agaricomycetes (Auriculariales, Sebacinales, Polyporales+Agaricales).
- Nuclearia/Parvularia: no genome; Parvularia atlantis transcriptome (EukProt) is the only sister of Fonticula.
