# R3 guide-tree edits (Streptophyta)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Viridiplantae` with children `Prascol`, `Chlorophyta_Streptophyta`. 21 leaves, 20 named internal nodes, 0 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | group | Chlorophyta_Streptophyta := (Chlorophyta, Streptophyta) | Prasinodermophyta (Prasinoderma) sister to Chlorophyta + Streptophyta | Li L. et al. 2020 Nat Ecol Evol 4:1220 |
| 2 | group | Mesostigmatophyceae_Chlorokybophyceae := (Mesovir, Chloatm) | Mesostigma + Chlorokybus form the earliest streptophyte clade | Wang S. et al. 2020 Nat Plants 6:95; Lemieux et al. 2007 BMC Biol 5:2 |
| 3 | group | Klebsormidiophyceae_Phragmoplastophyta := (Klebnit, Streptophytina) | grade: Klebsormidiophyceae, then Phragmoplastophyta (NCBI Streptophytina); M of streptophyte multicellularity | Wickett et al. 2014 PNAS 111:E4859; One Thousand Plant Transcriptomes Initiative 2019 Nature 574:679 |
| 4 | group | Anydrophyta := (Zygnematophyceae, Embryophyta) | Zygnematophyceae sister to land plants | Wickett et al. 2014 PNAS 111:E4859; One Thousand Plant Transcriptomes Initiative 2019 Nature 574:679; Cheng et al. 2019 Cell 179:1057 |
| 5 | group | Coleochaetophyceae_Anydrophyta := (Colescu, Anydrophyta) | Coleochaetophyceae sister to Anydrophyta | Wickett et al. 2014 PNAS 111:E4859; One Thousand Plant Transcriptomes Initiative 2019 Nature 574:679 |
| 6 | group | Setaphyta := (Marcpol, Physpat) | liverworts + mosses (bryophyte monophyly; hornworts not sampled) | Puttick et al. 2018 Curr Biol 28:733; Harris et al. 2020 Curr Biol 30:2001 |

## Open issues / notes

- No public annotated genome found (NCBI, UniProt 2026_03, EukProt, LukProt; 2026-10-29) for Zygnema, Penium, Spirogyra, Anthoceros/hornworts, Coleochaete (genomes of Feng et al. 2024 and Li et al. 2020 are hosted outside INSDC/UniProt or JGI-restricted). Coleochaete is a transcriptome (BUSCO 46).
- Zygnematophyceae internal order, and tracheophyte internal order kept from NCBI.
