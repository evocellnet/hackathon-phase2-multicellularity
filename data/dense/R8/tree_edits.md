# R8 guide-tree edits (Cyanobacteria)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Cyanobacteriota_Melainabacteria_group` with children `Sericyt`, `Melainabacteria_Cyanophyceae`. 36 leaves, 18 named internal nodes, 6 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | group | Cyanophyceae_excl_Gloeobacter := (Colechth, Gloelit, Leptbor, Pseubic, Spirsub, Acaryochloridales, Chroococcidiopsidales, Nostocales, Oscillatoriophycideae, Pleurocapsales, Synechococcales) | Gloeobacter sister to all other cyanobacteria | Shih et al. 2013 PNAS 110:1053; Schirrmeister et al. 2011 BMC Evol Biol 11:45; Schirrmeister et al. 2013 PNAS 110:1791 |
| 2 | group | Nostocales_Chroococcidiopsis := (Nostocales, Chroococcidiopsidales) | Chroococcidiopsidales sister to heterocystous Nostocales | Shih et al. 2013 PNAS 110:1053; Schirrmeister et al. 2011 BMC Evol Biol 11:45; Schirrmeister et al. 2013 PNAS 110:1791 |
| 3 | group | Branching_Nostocales := (Fiscthe, Chlofri) | true-branching Fischerella + Chlorogloeopsis form a clade | Shih et al. 2013 PNAS 110:1053; Gugger & Hoffmann 2004 Int J Syst Evol Microbiol 54:349 |
| 4 | group | Melainabacteria_Cyanophyceae := (Vampchl, Cyanophyceae) | Vampirovibrionia (Melainabacteria) sister to oxygenic cyanobacteria; Sericytochromatia outside | Soo et al. 2017 Science 355:1436 |

## Open issues / notes

- Backbone below Cyanophyceae_excl_Gloeobacter left as a polytomy (orders per NCBI / Strunecky et al. 2023); internal Nostocales order (except the branching clade) unresolved.
- Ancestral cyanobacterial state is debated (Schirrmeister et al. 2011, 2013: early filamentous origin with reversals; Hammerschmidt et al. 2021: stepwise). Filament origins are therefore coded as tip/clade events with a sensitivity coding.
