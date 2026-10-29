# R5 guide-tree edits (Rhodophyta)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Archaeplastida` with children `Cyanpar`, `Chlorophyta`, `Rhodophyta_Rhodelphidia`. 18 leaves, 14 named internal nodes, 3 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | dissolve | Bangiophyceae | NCBI node lumps Bangiales with Porphyridiophyceae and Cyanidiophyceae; regrouped below | - |
| 2 | group | Cyanidiophyceae := (Cyanmer, Galdieria) | Cyanidiophyceae sister to all other red algae | Yoon et al. 2006 J Phycol 42:482; Munoz-Gomez et al. 2017 Curr Biol 27:1677 |
| 3 | group | Proteorhodophytina := (Porppur, Compcae, Stylonemataceae) | Porphyridiophyceae + Compsopogonophyceae + Stylonematophyceae (+ Rhodellophyceae, unsampled) | Munoz-Gomez et al. 2017 Curr Biol 27:1677 |
| 4 | group | Eurhodophytina := (Bangiaceae, Rhodymeniophycidae) | Bangiophyceae + Florideophyceae | Yoon et al. 2006; Munoz-Gomez et al. 2017 Curr Biol 27:1677 |
| 5 | group | Rhodophytina := (Proteorhodophytina, Eurhodophytina) | mesophilic red algae | Munoz-Gomez et al. 2017 Curr Biol 27:1677 |
| 6 | group | Rhodophyta_Rhodelphidia := (Rhodophyta, Rhodemar) | Rhodelphis sister to red algae | Gawryluk et al. 2019 Nature 572:240 |
| 7 | group=parent | Archaeplastida -> ['Cyanpar', 'Rhodophyta_Rhodelphidia', 'Chlorophyta'] | Archaeplastida; order of glaucophytes / red / green left unresolved [equals existing node Eukaryota; renamed] | Burki et al. 2020 TREE 35:43 |

## Open issues / notes

- Florideophyceae are represented only by Rhodymeniophycidae (Chondrus, Gracilariopsis, Agarophyton, Delisea); no annotated genome of Hildenbrandiophycidae, Nemaliophycidae or Corallinophycidae passed QC (Calliarthron EukProt C 13, Hildenbrandia EST). The Florideophyceae stem and the Rhodymeniophycidae crown are therefore the same branch here.
- Proteorhodophytina internal order unresolved (polytomy); filamentous Bangiopsis and Compsopogon are MMETSP transcriptomes (flag).
- Bangiales genomes have low BUSCO completeness (Porphyra 44, Pyropia yezoensis 50, P. vietnamensis 61).
