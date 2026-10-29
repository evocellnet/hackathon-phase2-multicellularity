# R7 guide-tree edits (Amoebozoa)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Eukaryota` with children `Amoebozoa`, `Opisthokonta`. 25 leaves, 20 named internal nodes, 4 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | group | Holozoa := (Capsowc, Monobre) | Holozoa (outgroup structure) | Torruella et al. 2015 Curr Biol 25:2404 |
| 2 | group | Myxogastria_Ceratiomyxa := (Physpol, Cerafru) | Ceratiomyxa sister to Myxogastria | Fiore-Donno et al. 2010 Protist 161:55; Kang et al. 2017 MBE 34:2258 |

## Open issues / notes

- Physarum is a transcriptome: the Schaap et al. 2015 genome annotation is not distributed through NCBI/UniProt/EukProt/LukProt (checked 2026-10-29); most non-dictyostelid Amoebozoa (protosteliids, Tubulinea, Ceratiomyxa) also have transcriptome data only.
- Dictyostelid groups 1-4 sampled (Cavenderia; Acytostelium + Heterostelium; Tieghemostelium; Dictyostelium + Polysphondylium violaceum).
