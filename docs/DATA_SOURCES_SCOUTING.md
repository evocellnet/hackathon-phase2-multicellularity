# Data sources scouting (proteomes for the Ancestra runs)

Scouted 2026-10-28 (12:23 to 12:40 CET). Only metadata tables and 1–2 test files were downloaded, all under `data/raw/`.
Scripts are in `scripts/data/`. The outputs are `data/species_candidates.tsv` (eukaryotes) and `data/species_candidates_prok.tsv` (prokaryotes).

## 1. Sources

| Source | Version / URL | Licence | Size | Access (tested) | Notes / problems |
|---|---|---|---|---|---|
| **UniProt reference proteomes** | Release **2026_03** (2 Sep 2026). REST `https://rest.uniprot.org/proteomes/stream?query=reference:true AND taxonomy_id:<tx>&format=tsv&fields=upid,organism,organism_id,lineage,busco,cpd,genome_assembly,protein_count`. FTP `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/{Eukaryota,Bacteria,Archaea}/UPxxxx/UPxxxx_<taxid>.fasta.gz` | CC BY 4.0 | 3,521 eukaryotic, 17,972 bacterial and 641 archaeal reference proteomes. Per-proteome files are MB-sized; the full tarball is 348 GB, so do not download it | Tested: metadata saved to `data/raw/uniprot/refprot_{euk,bac,arc}.tsv`. Downloaded *Capsaspora* `UP000008743` (8,738 canonical proteins). `STATS` gives the canonical count per proteome | The `*.fasta` file is the canonical, gene-centric set (one protein per gene); isoforms are in `*_additional.fasta`, and `gene2acc` maps proteins to genes. **Problem 1:** in 2026_03 the human canonical file `UP000005640` holds **147,503** entries (20,652 in the previous release), so it must be reduced to one protein per gene with `gene2acc`, or Ensembl canonical used instead. A scan of `STATS` found only one other proteome that jumped (*Paspalum*, not used). **Problem 2:** the REST `protein_count` field includes isoforms, so the tables use the canonical counts from `STATS`. **Problem 3:** BUSCO is lineage-specific (n=255 eukaryota, 1519 chlorophyta, 758 fungi, …), so it cannot be compared across clades. CPD is "Unknown" for about 96% of proteomes and is useless as a filter. Two reference proteomes are fragmentary (*Oxytricha* UP000053232 with 810 proteins, *Andalucia* UP000799049 with 1,023), so EukProt is used for those. AlphaFold DB and InterPro are keyed on UniProt accessions, which is an advantage for the structure analyses |
| **EukProt v3** (Richter et al. 2022, PCJ) | figshare `10.6084/m9.figshare.12417881.v3` (2021_11_22). Proteins: `https://ndownloader.figshare.com/files/34434377` | CC BY 4.0. Original data-source terms still apply | `EukProt_proteins.v03.2021_11_22.tgz` is **5.49 GB** (single tgz, so the whole file is needed). TCS subset (196 spp.) is figshare 21586065 (1.17 GB). Pfam v34 annotations are figshare 24680811 (0.83 GB). Genome annotations tgz is 40 MB | Tested: HTTP range request works (`proteins/` tar header). Metadata in `data/raw/eukprot/`: included sets (993), not-included sets, **BUSCO v5.2.2 eukaryota_odb10 for all sets** (GitHub `beaplab/EukProt/annotations/v03`), TCS list | 375 genome, 497 transcriptome and 103 single-cell datasets. **The metadata table has no NCBI taxids and no protein counts**; LukProt's sheet fills both. Transcriptome sets are redundant (up to about 60k "proteins") and break Ancestra's one-protein-per-gene assumption. The data dates from 2021 |
| **LukProt v1.5.1.rev2** (Sobala 2024, GBE) | Zenodo `13829058` (24 Sep 2024). BLAST server: lukprot.hirszfeld.pl | CC BY 4.0 | `LukProt_v1.5.1_single_species_FASTA.7z` is **5.59 GB**. OMArk output 0.2 GB, BUSCO tables 13 MB, metadata sheet 0.9 MB | Tested: metadata sheet downloaded and converted to TSV (`scripts/data/ods_to_tsv.py`; output `data/raw/lukprot/LukProt_v1.5.1_metadata.tsv`) | 1,281 proteomes: the EukProt v3 sets plus 41 from AniProtDB and 250 new, mostly animals and holozoans. **The metadata gives NCBI taxid, sequence counts, BUSCO and OMArk contamination for all EukProt sets.** Needs `7z` to unpack (check the conda env). It holds improved genomes we use: *Ministeria vibrans* LP00007, *Sphaeroforma arctica* LP00010 (BUSCO 91, against 63 in UniProt), *Hoilungia*, *Monocercomonoides* |
| **"Eukaryoma"** | Not found (time-boxed to 15 min). Searched Zenodo, the web, and the GitHub and EukProt pages | – | – | – | No resource with this name could be located. Possibly an unreleased or internal dataset, or a misremembered name. **Ask the user.** Related resource found: *The Structural History of Eukarya* (SHE; bioRxiv 10.64898/2026.02.04.703859): Foldseek all-vs-all of AlphaFold v4 structures for 1,542 UniProt eukaryotic proteomes, at https://she-app.serve.scilifelab.se/. It could be useful later for the structure-remodelling analysis |
| **NCBI Datasets v2 (RefSeq/GenBank)** | `https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<GCx>/download?include_annotation_type=PROT_FASTA` (add `GENOME_GFF` for gene order). No `datasets` CLI is installed; plain HTTP works | INSDC: no use restrictions (cite the submitters) | Few MB per genome | Tested: *Schizocladia* GCA_964200225.2 `protein.faa` (21,272 proteins) | Main use: **Phaeoexplorer brown-algal genomes** (Denoeud et al. 2024; 42 Phaeophyceae plus *Schizocladia*, *Phaeothamnion*, *Heterosigma*), *Bryopsis*, *Ostreobium*, *Closterium*, extra Florideophyceae, *Willaertia*, and the MMP, cable-bacteria and *Jeongeupia* prokaryotes. RefSeq `.faa` files contain isoforms, so collapse to the longest protein per gene with the GFF. NCBI reports often carry **no BUSCO**, so completeness must be computed by us. `scripts/data/ncbi_check.py` and `ncbi_taxon_annotated.py` query genomes, TSA and SRA (output: `data/raw/ncbi_annotated_by_taxon.tsv`) |
| **Ensembl Genomes** | Release **63** (`https://ftp.ensemblgenomes.ebi.ac.uk/pub/{protists,fungi,plants,metazoa}/current/`). REST `rest.ensembl.org/info/genomes/division/EnsemblProtists` (236 genomes) | Open | Per genome | Tested: FTP and REST reachable | Proteomes overlap with UniProt, which imports from Ensembl. Most useful for **GFF3 gene order** (Ancestra `--gene-order` builder for Ensembl GFF3) and as a source of canonical-per-gene human data |
| **JGI (Mycocosm, Phycocosm, Phytozome)** | – | Login required. Many genomes are under a data-use moratorium | – | Not attempted | Avoid. **Ask the user** only if essential, e.g. *Zygnema circumcarinatum* or *Anthoceros* annotations, or unpublished fungal and algal genomes |

## 2. Candidate tables

`data/species_candidates.tsv` and `data/species_candidates_prok.tsv` are built with `scripts/data/build_source_index.py` followed by `scripts/data/build_candidates.py`. The curated inputs are `scripts/data/candidates_{euk,prok}_spec.txt`. To change the selection, edit the spec files and re-run.

- **Source rule.** UniProt reference proteome first (when BUSCO ≥ 60, or when no alternative is more than 10 points better). Then an EukProt/LukProt genome, then an EukProt/LukProt transcriptome, then NCBI. Explicit overrides are recorded.
- **`transition_id` and `role` are parallel `;` lists.** One species can inform several transitions. For example, *S. rosetta* is the sister for `T_metazoa` and a multicellular taxon for `T_choano_colonies`.
- **`tier`:** `core` = proposed for the main run, `optional` = add if the budget allows or for sensitivity analysis, `exclude` = documented only.
- **Flags in `notes`:** `[TRANSCRIPTOME-ONLY]`, `[LOW-COMPLETENESS]` (BUSCO < 60), `[NO-PROTEOME-AVAILABLE]`.
- **Extra columns:** `busco_detail` (which BUSCO lineage was used, CPD), `alt_sources`, `taxonomy`.

**Eukaryotes:** 260 rows, of which 151 are core (about 2.56 M proteins once human is reduced to 1 protein per gene). Core plus optional is 255 rows and about 4.6 M proteins. Core data types: 131 genome, 18 transcriptome, 2 missing. Core sources: UniProt 98, EukProt 41, LukProt 3, NCBI 7, none 2. Core plus optional sources: UniProt 149, EukProt 82, LukProt 5, NCBI 14.

**Prokaryotes:** 120 rows, of which 79 are core (about 0.40 M proteins). Core plus optional is 119 rows and about 0.61 M proteins. All are genomes: 111 from UniProt and 8 from NCBI.

Eukaryotic transitions covered (counts are core / available):

| Transition | Multicellular taxa | Unicellular sisters | Notes |
|---|---|---|---|
| T_metazoa | 11/18 | 6/9 | choanoflagellates. Filastereans and ichthyosporeans are outgroups |
| T_choano_colonies | 3/4 | 3/5 | |
| T_capsaspora_aggregation | 1 | 1/2 | |
| T_ichthyosporea_coenocytic | 3/5 | 1/2 | |
| T_fonticula | 1 | 1 | |
| T_fungi_hyphae | 4/10 | 2/3 | |
| T_fungi_pezizomycotina | 3/10 | 3/4 | |
| T_fungi_neolecta | 1 | 2/4 | |
| T_fungi_agaricomycetes | 4/7 | 2/4 | |
| T_dictyostelia | 6/8 | 2/4 | |
| T_copromyxa | 1 | 1/4 | |
| T_myxogastria_plasmodium | 1 | – | grey zone |
| T_volvocales | 5/6 | 3/4 | |
| T_ulvophyceae / T_siphonous_coenocytic | 2/4 | – | |
| T_streptophyta | 7/10 | 4/5 | plus 2 missing |
| T_rhodophyta | 4/7 | 2/2 | |
| T_phaeophyceae | 5/8 | 4/4 | |
| T_xanthophyceae_filaments | 1/2 | – | |
| T_oomycete_hyphae | 2/4 | – | |
| T_sorodiplophrys | **0 (no data)** | 3/5 | |
| T_sorogena | **0 (no data)** | 1/1 | *Platyophrya* only |
| T_guttulinopsis | 1 | 2/5 | transcriptome, BUSCO 11 |
| T_acrasis | 1 | 2/4 | |
| T_plasmodiophorid_plasmodium | 1 | – | |
| T_none (background) | – | – | 18 core / 33 available, plus all sisters and outgroups serve as control nodes |

Prokaryotic transitions covered (counts are core / available):

| Transition | Multicellular taxa | Unicellular sisters |
|---|---|---|
| T_cyano_filament | 14/19 | 8/9 |
| T_cyano_heterocyst | 10/13 | – |
| T_cyano_branching | 2/3 | – |
| T_myxobacteria | 8/12 | 2/4 (non-fruiting *Anaeromyxobacter*, *Vulgatibacter*) |
| T_streptomyces | 4/8 | 5/8 |
| T_actino_mycelium_other | 5/8 | – |
| T_magnetotactic_mmp | 2 (metagenomic drafts, completeness unknown) | 2/3 |
| T_cable_bacteria | 1/2 | 1/2 |
| T_jeongeupia | 1 | 2/3 |
| T_methanosarcina | 3/4 | 3/4 |
| T_haloferax_tissue | 1 | 1/2 |
| Other | Thiotrichales, Chloroflexi and archaeal filaments are optional | |

Transcriptome-only or problematic multicellular taxa:

| Taxon | Status | Source |
|---|---|---|
| *Sorogena stoianovitchae* | No data anywhere: no NCBI genome, TSA or SRA, and absent from EukProt and LukProt | – |
| *Sorodiplophrys stercorea* | No data anywhere (same checks) | – |
| *Guttulinopsis vulgaris* | Transcriptome, BUSCO 11 | EukProt EP00461 |
| *Copromyxa protea* | Transcriptome, BUSCO 83 | EukProt EP01127 |
| *Physarum polycephalum* | Transcriptome, BUSCO 81. Genome annotation not on NCBI | EukProt EP00026 |
| *Choanoeca flexa* | Transcriptome | EukProt EP00730 |
| *Fonticula alba* | **Genome** (not a transcriptome), BUSCO 71 | UniProt UP000030693 |
| *Acrasis kona* | **Genome**, BUSCO 86 | UniProt UP001431209 |

## 3. Recommended selection strategy

1. **Main eukaryotic run: about 190–210 proteomes.** This is the 151 core taxa plus about 40–60 optional ones. Choose optional taxa first where a transition has fewer than 2 taxa on either side (Volvocales sisters, Pezizomycotina and Agaricomycetes depth, Phaeophyceae from Phaeoexplorer, Tubulinea for Copromyxa, Cercozoa for Guttulinopsis), and second to thicken the background in under-sampled supergroups (Discoba, Metamonada, CRuMs, Haptista).
   - The total would be about 3.0–3.5 M proteins.
   - Runtime: extrapolating from the benchmarks (78 proteomes / 1 M proteins `accurate` in about 1.5 h at 12 threads, with roughly quadratic all-vs-all cost) gives about 10–15 h for `accurate` at 24 threads, or a few hours for `fast`.
   - **Benchmark a 40-proteome subset first.**
   - Going above about 250 proteomes adds mainly redundant fungi and animals and dilutes the rare transitions.
2. **Data source per taxon.**
   - **UniProt FTP:** canonical `*.fasta` from release 2026_03, pinned. Human must be reduced to one protein per gene with `gene2acc`.
   - **EukProt v3 (one 5.5 GB tgz):** 41 core taxa, 82 including optional.
   - **LukProt:** 3 core taxa, 5 including optional. Download the 5.6 GB 7z, or use the EukProt/other original sources for these 5 files.
   - **NCBI:** 7 core taxa, 14 including optional (Phaeoexplorer etc.). Collapse isoforms using the GFF.
3. **Uniform QC after download.** Run BUSCO eukaryota_odb10 (or OMArk) on every final FASTA, because the BUSCO values in the tables come from different lineage sets and are not comparable. Use OMArk to screen transcriptomes for contamination; the LukProt sheet already has OMArk results for EukProt sets.
   - Transcriptomes: collapse redundancy (CD-HIT about 0.95–1.0, or longest ORF per Trinity "gene") before Ancestra.
   - Prefix protein IDs with a species code (Ancestra needs unique IDs).
4. **Guide tree.** Start from NCBI taxonomy using the taxids in the tables, then fix it by hand to the current consensus (Holomycota/Holozoa, Pluriformea, CRuMs, Volvocales internal order, Phaeophyceae, Discoba). Record every edit.
5. **Transitions with one multicellular taxon** (*Capsaspora*, *Fonticula*, *Neolecta*, *Acrasis*, *Copromyxa*, *Guttulinopsis*) can only be analysed as "terminal branch" gains and losses. Interpret them with matched unicellular terminal branches as controls. Results for *Guttulinopsis* (BUSCO 11) should be reported for presence only.
6. **Prokaryotes:** a separate run of the 79 core proteomes (`accurate` should be cheap).

## 4. Not accessible / ask the user

| Item | Status | Question for the user |
|---|---|---|
| "Eukaryoma" | Not found | What is it? Is there a link? |
| *Sorogena stoianovitchae* | No sequence data in public databases | Does the user know of an unpublished or collaborator transcriptome? If not, this transition cannot be covered |
| *Sorodiplophrys stercorea* | Same as above | Same as above |
| JGI Mycocosm, Phycocosm, Phytozome | Login and moratoria | Should we use them at all? Candidates: *Zygnema circumcarinatum* (filamentous zygnematophyte), *Anthoceros* (hornwort), *Penium*, more Ulvophyceae and fungi |
| *Physarum polycephalum* genome annotation (Schaap et al. 2016) | Not on NCBI or UniProt | Is a genome-based set preferred over the EukProt transcriptome? |
| *Colpoda* (Colpodea, sister of *Sorogena*) | SRA RNA-seq only | Should we de novo assemble it? This costs hours and is optional |
| *Nuclearia* (sister of *Fonticula*) | No proteome; *Parvularia atlantis* is used instead | – |
| *Pocheina* | No proteome | – |
| Phaeoexplorer portal (Roscoff) | Proteomes are on ENA/NCBI already. The portal also hosts additional annotation (orthogroups) | Not needed unless required |
| `7z` binary | Needed for LukProt archives | Check it is in the conda env |
