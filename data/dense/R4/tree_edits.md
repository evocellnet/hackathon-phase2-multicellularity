# R4 guide-tree edits (Volvocine algae (Chlamydomonadales))

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Chlorophyta` with children `Ostrtau`, `core_chlorophytes`. 21 leaves, 13 named internal nodes, 3 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | move | Chlaeus -> Chlamydomonadales | genus Chlamydomonas is polyphyletic; C. eustigma is outside Reinhardtinia | Proschold et al. 2001 Protist 152:265 |
| 2 | move | Chlapri -> Chlamydomonadales | C. priscui (UWO241) is outside Reinhardtinia (Moewusinia) | Proschold et al. 2001 Protist 152:265; Possmayer et al. 2016 J Phycol |
| 3 | move | Chlamydomonas -> Chlamydomonadales | Chlamydomonas s.str. (C. reinhardtii group) belongs to Reinhardtinia, not with other NCBI Chlamydomonadaceae | Proschold et al. 2001 Protist 152:265; Craig et al. 2021 Plant Cell 33:1016 |
| 4 | dissolve | Chlamydomonadaceae | NCBI family is polyphyletic | Proschold et al. 2001 Protist 152:265 |
| 5 | move | Gonipec -> Chlamydomonadales | NCBI places Gonium in Volvocaceae; Goniaceae is paraphyletic | Lindsey et al. 2021 BMC Biol 19:182 |
| 6 | dissolve | Volvox | genus Volvox is polyphyletic; Volvocaceae internal order left unresolved | Herron & Michod 2008 Evolution 62:436; Lindsey et al. 2021 BMC Biol 19:182 |
| 7 | group | Volvox_africanus_reticuliferus := (Volvafr, Volvret) | V. reticuliferus sister of V. africanus | Nozaki et al. 2015 PLoS One 10:e0142632 |
| 8 | group | Astrephomene_Volvocaceae := (Astrgub, Volvocaceae) | Astrephomene sister to Volvocaceae | Lindsey et al. 2021 BMC Biol 19:182 (MLBS 98) |
| 9 | group | Gonium_Volvocaceae := (Gonipec, Astrephomene_Volvocaceae) | Gonium sister to Astrephomene + Volvocaceae (Goniaceae paraphyletic); M of the 2nd colonial origin | Lindsey et al. 2021 BMC Biol 19:182 |
| 10 | group | TGV_clade := (Tetrsoc, Gonium_Volvocaceae) | Tetrabaenaceae + Goniaceae + Volvocaceae among SAMPLED taxa; an unsampled unicellular lineage is sister of Tetrabaenaceae, so the ancestor here is inferred unicellular and colonial multicellularity arose twice | Lindsey et al. 2021 BMC Biol 19:182 |
| 11 | group | Chlamydomonas_TGV := (Chlamydomonas, TGV_clade) | core Reinhardtinia: C. reinhardtii/incerta/schloesseri sister to the colonial lineages | Craig et al. 2021 Plant Cell 33:1016; Lindsey et al. 2021 BMC Biol 19:182 |
| 12 | group | Reinhardtinia := (Edapdeb, Chlamydomonas_TGV) | Edaphochlamys outside core Reinhardtinia | Craig et al. 2021; Lindsey et al. 2021 BMC Biol 19:182 |

## Open issues / notes

- Volvocaceae internal order left as a polytomy (Volvox polyphyletic); Lindsey et al. 2021 infer 4-6 origins of cellular differentiation, coded here as separate tip/clade events.
- Undifferentiated Volvocaceae (Yamagishiella, Eudorina, Pandorina) have no published annotation; if present they are AUGUSTUS ab initio predictions (chlamy2011 model) on the unannotated Volvocales Genome Consortium assemblies (flagged).
- Chlamydomonadales outside Reinhardtinia (C. eustigma, C. priscui, Sanguina, Dunaliella) kept as a polytomy.
