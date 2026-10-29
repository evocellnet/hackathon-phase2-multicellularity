# R1 guide-tree edits (Holozoa)

Built by `scripts/dense/build_panel.py` (edits in `scripts/dense/panels.py`). Base: NCBI taxonomy (`data/raw/taxdump`, 2026-10-28) via `python -m ancestra.speciestree ncbi --leaf-naming qfo` on `<code>_<taxid>.fa` symlinks; leaves renamed to file stems; node names sanitised; unary nodes collapsed; children sorted (canonical Newick). Unedited tree: `guide_tree_ncbi.nwk`. Same rules as round 1 (`data/guide_tree_edits.md`): resolve or move only where recent phylogenomic work agrees, otherwise keep the multifurcation.

Root: `Eukaryota` with children `Thectra`, `Opisthokonta`. 37 leaves, 28 named internal nodes, 5 multifurcating nodes.

| # | operation | nodes | rationale | reference |
|---|---|---|---|---|
| 1 | rename | Opisthokonta_incertae_sedis -> Pluriformea | NCBI bin holds exactly Corallochytrium + Syssomonas (Pluriformea) | Hehenberger et al. 2017 Curr Biol 27:2043 |
| 2 | group | Teretosporea := (Ichthyosporea, Pluriformea) | Pluriformea sister to Ichthyosporea | Torruella et al. 2015 Curr Biol 25:2404; Hehenberger et al. 2017 Curr Biol 27:2043 |
| 3 | group | Choanozoa := (Metazoa, Choanoflagellata) | choanoflagellates = sister of animals | King et al. 2008 Nature 451:783; Brunet & King 2017 Dev Cell 43:124 |
| 4 | group | Filozoa := (Filasterea, Choanozoa) | Filasterea sister to Choanozoa | Torruella et al. 2015 Curr Biol 25:2404; Hehenberger et al. 2017 Curr Biol 27:2043 |
| 5 | group | Holozoa := (Teretosporea, Filozoa) | Holozoa | Torruella et al. 2015 Curr Biol 25:2404; Hehenberger et al. 2017 Curr Biol 27:2043 |
| 6 | group | Holomycota := (Fontalb, Spizpun) | Fonticula (nucleariids) + Fungi | Brown et al. 2009 Mol Biol Evol 26:2699 |
| 7 | move | Tentaculata -> Metazoa | NCBI Eumetazoa (ctenophores + ParaHoxozoa) is not supported; sponge-first vs ctenophore-first contested -> Porifera / Ctenophora / ParaHoxozoa polytomy (as round 1) | Simion et al. 2017 Curr Biol 27:958; Schultz et al. 2023 Nature 618:110 |
| 8 | rename | Eumetazoa -> ParaHoxozoa | remaining clade = Placozoa + Cnidaria + Bilateria (order unresolved) | Ryan et al. 2010 EvoDevo 1:9; Laumer et al. 2019 eLife 8:e36278 |
| 9 | group | Calcarea_Homoscleromorpha := (Sycocil, Oscalob) | Calcarea sister to Homoscleromorpha | Simion et al. 2017; Whelan et al. 2015 PNAS 112:5773 |
| 10 | group | Silicea := (Aphrvas, Heteroscleromorpha) | Hexactinellida + Demospongiae | Simion et al. 2017; Whelan et al. 2015 |

## Open issues / notes

- Metazoa root kept as Porifera / Ctenophora / ParaHoxozoa polytomy (ctenophore-first: Schultz et al. 2023; run variants if needed).
- Choanoflagellata internal order (NCBI Salpingoecidae polytomy) and Filasterea internal order kept from NCBI.
- Teretosporea: Abeoforma and Pirum EukProt gene sets unusable (BUSCO C < 3); Ichthyophonus (C 54) and Tunicaraptor (C 56) dropped.
