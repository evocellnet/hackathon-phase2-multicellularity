# R6 guide-tree edits (Stramenopiles / Phaeophyceae)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Stramenopiles` with children `Caferoe`, `Gyrista`. 29 leaves, 17 named internal nodes, 5 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | group | Discosporangiales := (Chorten, Discmes) | Choristocarpus + Discosporangium = Discosporangiales | Silberfeld et al. 2014 Cryptogam Algol 35:117; Denoeud et al. 2024 Cell 187:6943 |
| 2 | move | Discosporangiales -> Phaeophyceae | NCBI places Discosporangiales inside Sphacelariales; they are the earliest-diverging brown algal order | Silberfeld et al. 2014 Cryptogam Algol 35:117; Denoeud et al. 2024 Cell 187:6943 |
| 3 | group | Dictyotophycidae := (Dictdic, Sphacelariales) | Dictyotales + Sphacelariales | Silberfeld et al. 2014 Cryptogam Algol 35:117; Denoeud et al. 2024 Cell 187:6943 |
| 4 | group | Fucophycidae := (Desmher, Ectocarpales, Fucaceae, Laminariales) | brown algal crown radiation (Fucophycidae); internal order left unresolved | Silberfeld et al. 2014 Cryptogam Algol 35:117; Denoeud et al. 2024 Cell 187:6943 |
| 5 | group | Dictyotophycidae_Fucophycidae := (Dictyotophycidae, Fucophycidae) | Discosporangiales sister to all other brown algae | Silberfeld et al. 2014 Cryptogam Algol 35:117; Denoeud et al. 2024 Cell 187:6943 |
| 6 | group | Phaeophyceae_Schizocladia := (Phaeophyceae, Schiisc) | Schizocladia sister to brown algae | Kawai et al. 2003 J Phycol 39:556; Denoeud et al. 2024 Cell 187:6943 |
| 7 | group | Gyrista := (Ochrophyta, Oomycota) | Ochrophyta + Pseudofungi (oomycetes); Bigyra outside | Derelle et al. 2016 MBE 33:2890; Burki et al. 2020 TREE 35:43 |

## Open issues / notes

- Ochrophyta internal order (diatoms, Pelagophyceae, Eustigmatophyceae, Chrysophyceae, Raphidophyceae, PX clade) left as a polytomy; PX clade internal order (Chrysoparadoxa, Phaeothamnion, Tribonema, Schizocladia+Phaeophyceae) unresolved.
- Phaeoexplorer genomes vary strongly in completeness; Sphacelaria, Macrocystis, Saccorhiza x2, Himanthalia, Laminaria were rejected (BUSCO eukaryota C < 35); Discosporangiales (Discmes 46, Chorten 34) and Halopteris (48) kept as the only representatives of their orders (flagged). Ishigeales: no genome.
