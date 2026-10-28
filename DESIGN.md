# DESIGN — Genomic transitions to multicellularity across the tree of life

Version: v0 draft (pre-review). Owner: PI agent. Fixed by the task: use Ancestra to reconstruct ancestral gene content before and
after the documented transitions to multicellularity, run a broad set of analyses, compare with existing knowledge, and assess
whether shifts in evolutionary pressure can be detected in sequences or structures.

## 1. Questions

- **Q1 Gene-content change.** How many gene families are gained, duplicated and lost on branches that lead to multicellular clades,
  compared with branches that do not? Is a transition a burst of novelty, or mostly re-use (co-option) of genes that already existed?
- **Q2 Function.** Which functions (Pfam domains; pre-registered "toolkit" categories: adhesion, extracellular matrix, signalling
  kinases and receptors, transcription factors, programmed cell death, cell-cycle and cytoskeleton) are gained or expanded at transitions?
  Which are gained or expanded at more independent transitions than chance predicts (convergence)?
- **Q3 Pre-adaptation.** For families that characterise each multicellular ancestor, what fraction were already present in its last
  unicellular ancestor (co-option index)? How does this differ between clonal and aggregative, simple and complex transitions?
- **Q4 Domain architecture.** Are new domain combinations (architectures) created at transitions more than on control branches?
- **Q5 Evolutionary pressure.** Do multicellular lineages show signatures of changed selective constraint: (a) faster or slower
  relative evolutionary rates (RER) of conserved single-copy orthologs; (b) proteome-wide traits (protein length, predicted disorder,
  low complexity, domains per protein, paralog fraction); (c) individual families whose rates shift convergently at independent transitions?
- **Q6 Structure.** Are families "gained" at transitions truly new (no structural homolog in the pre-transition ancestor's proteome)
  or fast-diverged copies of older families (structural homolog present)? Do their structures change?
- **Q7 Literature concordance.** Which published claims (LITERATURE.md §7) does our data reproduce, contradict or not address?

## 2. Data

- **Main eukaryote run:** ~190–210 proteomes (`data/species_candidates.tsv`, core + selected optional), ~3.0–3.5 M proteins. Sources
  in order of preference: UniProt reference proteomes 2026_03 (genome-based; AlphaFold models exist), EukProt v3, LukProt v1.5.1, NCBI
  (Phaeoexplorer brown algae). Each transition: ≥2 multicellular species and ≥2 closest unicellular relatives where data exist,
  plus a background spanning all eukaryotic supergroups (provides the control branches).
- **Prokaryote run (separate):** ~79 core proteomes (filamentous/heterocystous cyanobacteria, Myxococcota, Streptomyces and other
  mycelial actinomycetes, cable bacteria, magnetotactic multicellular bacteria, Methanosarcina, haloarchaea). Separate because orthology
  across domains costs much and informs none of the transitions.
- **Pre-processing:** one protein per gene (UniProt gene2acc / NCBI GFF longest isoform; transcriptomes: CD-HIT 0.95 + longest ORF per
  Trinity gene where headers allow); species-prefixed IDs; length ≥ 50 aa; uniform BUSCO v5 (eukaryota_odb10 / bacteria / archaea)
  on the final files; recorded in DATA.md with checksums and licences.
- **Guide tree:** NCBI taxonomy (taxdump downloaded 2026-10-28) via `ancestra.speciestree`, then edited to current consensus
  (Holozoa/Holomycota, Pluriformea, Volvocales order, Phaeophyceae, Discoba, CRuMs). Every edit logged with its reference.
  Internal node names = levels. Transition nodes are labelled in a separate table (`data/transitions.tsv`): for each transition the
  multicellular crown node M, its parent P (last unicellular ancestor), the type (clonal/aggregative, simple/complex, facultative/obligate),
  and whether M is a single species (terminal-branch transition).
- **Questions to the user (not blocking):** what "Eukaryoma" is; data for *Sorogena*/*Sorodiplophrys*; JGI-restricted charophytes.

## 3. Methods

### 3.1 Ancestral gene content (Ancestra)
- Timing test on a 40-proteome subset (`accurate` and `fast`) before committing; full run with the preset that fits the budget
  (prefer `accurate`). No gene-order layer (not available for most taxa; it gave no measurable gain in Ancestra's benchmarks).
- Parse OrthoXML to per-branch events (module `mcell.hogs`, validated on Ancestra's simulator with known truth):
  gains (HOG originates on the branch), losses (HOG present at parent, absent in the child clade), duplications (one parent HOG →
  k child HOGs counts k−1). Conservation check per branch.
- Gains are reported in two flavours: all, and "supported" (present in ≥2 descendant species, or in ≥2 child lineages), because
  single-species gains are dominated by annotation and contamination artefacts.

### 3.2 Annotation
- Pfam-A 38.2 with pyhmmer hmmsearch (gathering cutoffs) on all proteins; HOG function = domains present in ≥50% of members.
- Toolkit categories defined *before* looking at results from Pfam accession lists drawn from the literature (file `data/toolkit_pfam.tsv`,
  each entry with its citation). Frozen at design time.
- Disorder (metapredict v3), low-complexity (SEG-like via pyhmmer/simple entropy), length, domain counts per protein.

### 3.3 Branch-level statistics (Q1, Q2, Q4)
- Unit = branch of the guide tree. Transition branches = P→M stems (and terminal branches for single-species transitions).
  Controls = all other branches.
- Branch lengths in substitutions/site from IQ-TREE (LG+G4 or LG+C20 if affordable) on a concatenated alignment of ~200 near-single-copy
  HOGs with the guide-tree topology fixed (gives the time-proxy covariate; it is NOT used to change the topology).
- Model: negative-binomial GLM `events ~ is_transition + log(branch_length) + log(n_descendant_species) + frac_transcriptome_desc +
  mean_BUSCO_desc + terminal`. Report the transition effect (rate ratio, 95% CI). Internal and terminal branches also analysed separately.
- Functional enrichment: for each transition, Pfam domains/toolkit categories among gained (and expanded) HOGs vs. all HOGs present
  at M (Fisher, BH-FDR). Convergence: number of transitions in which a domain/HOG is gained/expanded vs. a null that assigns transition
  labels to random control branches matched for depth class and descendant sampling (10,000 draws, seed list fixed).
- Expansion: per HOG, copy number per species; test multicellular vs sister unicellular species (per transition, Wilcoxon), then
  convergence across transitions as above.
- Domain architectures: ordered Pfam architecture per protein; "new architecture at M" = present in ≥2 descendant species of M,
  absent from all species outside M. Count per branch, same GLM.

### 3.4 Co-option index (Q3)
- For each transition: families significantly enriched/expanded in the multicellular clade or in toolkit categories present at M →
  fraction already present at P (and at the grandparent). Compare to the same fraction for matched control nodes.

### 3.5 Evolutionary pressure (Q5)
- **RER:** HOGs present as single copy in ≥70% of species (or of a clade-level subset for deeper coverage) → MAFFT, trimAl, IQ-TREE branch
  lengths on the pruned guide topology → relative rates per branch after RERconverge-style normalisation (residual vs. genome-wide expected
  rate). Tests: (i) genome-wide: are RERs on multicellular stem/crown branches higher than on matched unicellular branches?
  (ii) per gene: association of RER with the multicellular trait across transitions (Wilcoxon / correlation, BH-FDR), with trait-permutation null.
- **Proteome traits:** sister-clade contrasts per transition (multicellular clade mean vs sister unicellular clade mean); sign test /
  Wilcoxon across transitions; phylogenetic GLS on all species as a second view.
- **dN/dS (optional, focused):** only where divergence is low enough to avoid synonymous saturation and CDS are public (e.g. within
  Volvocales or within cyanobacterial genera) with HyPhy RELAX. Decided after the core results.

### 3.6 Structure (Q6)
- AlphaFold DB v4 models for UniProt proteomes (public per-proteome tar files) for one multicellular and one sister species per
  transition where available; Foldseek search of gained-HOG representatives against the pre-transition species' structures →
  fraction of "sequence-new" gains that have a structural homolog (E < 1e-3, TM ≥ 0.5). Control: the same for gains on control branches.
- Ancestral sequence reconstruction (IQ-TREE `-asr`) for a few convergent candidate families, with descriptive structural comparison.

### 3.7 Positive and negative controls (pre-registered)
- **Positive (literature recovery):** e.g., cadherin and receptor tyrosine kinase domains present at the Choanozoa node (pre-metazoan);
  integrin β present in Capsaspora/Filasterea; Wnt, TGF-β ligands, Hox-type homeobox and Notch/Delta complete set at the Metazoa node; volvocine
  pherophorin/VARL expansion; expansion of kinases in Dictyostelia. The full list (with expected direction) is written to
  `data/positive_controls.tsv` before results are inspected.
- **Negative:** (a) pseudo-transitions on random unicellular branches matched for depth/sampling; (b) trait-label permutation for RER and
  trait contrasts; (c) toy simulations with known events for the parser.
- **Robustness:** orthology-free check of Q1–Q2 with Pfam domain counts (no HOGs); `fast` vs `accurate` on the 40-proteome subset;
  genome-only species subset; excluding transcriptome taxa.

## 4. Reproducibility
- Seeds: [11, 23, 37, 41, 53] (first is default). All scripts in `scripts/`, library in `mcell/`, tests in `tests/`.
- DATA.md: source, version, checksum, licence per dataset. Raw data not committed. Results (tables, figures) committed.
- Every key result re-derived by an agent that did not produce it before being cited in the report.

## 5. Compute plan (24 CPU, 150 GB, 48 h allocation; target ≤36 h)
| t (h) | milestone |
|---|---|
| 0–3 | literature, data scouting, design + reviews |
| 3–7 | dataset build, QC, guide tree, positive-control list; Ancestra timing test |
| 7–15 | Ancestra main run (12 threads) ‖ Pfam + disorder annotation (12 threads); prokaryote run after |
| 15–28 | analyses Q1–Q6 in parallel agents |
| 28–34 | verification agents, figures |
| 34–40 | report |
Each stage writes checkpoints (Ancestra stages are resumable).
