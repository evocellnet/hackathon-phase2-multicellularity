# DESIGN — Round 2: resolving single events and the step from simple to complex multicellularity

Version: v0 (2026-10-29 15:45, round-2 t ≈ 0.6 h; pre-review). Round 1 is frozen at git tag `phase1-final`.
Inputs: PROMPT.md (rules unchanged), PROMPT_ROUND2.md (new scientific points), the department-head review
(reviews/round2_department_head_review.txt; public page https://claude.ai/artifact/NfrexbqRBFmS3A97VvTqcR).
Resources: 48 CPU, 500 GB RAM, 1× A100-40GB, ~30 h soft deadline (START_ROUND2.txt: 2026-10-29 15:07).

## 0. What changes and why

Round 1 asked whether seven very different origins share one genomic signature. With 7 stems against 30 controls it had
13–35% power at a 1.5-fold effect, so its gene-content nulls are absence of evidence. Its most useful products were a map of where
Ancestra's ancestral gene content cannot be trusted (deep-node splits, duplications placed one node too deep, domains gained inside older
families) and two leads: more sequence-novel genes at the three nested simple→complex transitions (p = 0.016, outside the Holm family)
and many new domain architectures on the animal stem. Round 2 therefore (i) fixes round-1 errors, (ii) validates round-1 events with an
independent orthology resource, (iii) resolves single events with dense, focused Ancestra runs, (iv) asks how complex multicellularity
arises from simple multicellularity, and (v) replaces class-level aggregates with specific, pre-defined protein families.

## 1. Workstreams

| id | question | data / method | priority |
|---|---|---|---|
| W0 | Errata: depth SMD, power table on real branch lengths, positive-control count, H9 orientation dependence | round-1 files; exact enumeration over 1,024 orientations | done first |
| W1 | Do round-1 stem gains/losses/duplications agree with an independent HOG set? | EukaryOMA (1,328 genomes, OMA/FastOMA HOGs) projected on our tree; event-by-event agreement, stems vs controls; C03/C06/C07 | high |
| W2 | What happens at each event when it is densely sampled? | 8 focused Ancestra runs (`accurate`), R1 Holozoa, R2 Fungi, R3 Streptophyta, R4 Volvocales, R5 Rhodophyta, R6 Phaeophyceae/stramenopiles, R7 Amoebozoa, R8 Cyanobacteria; 25–60 genome-based proteomes each; in-run controls | high |
| W3 | How does complex multicellularity arise from simple multicellularity? | nested simple→complex events in R2 (Pezizomycotina, Agaricomycetes), R3 (Embryophyta), R4 (Volvox), R5 (Florideophyceae), R6 (complex brown algae), R8 (heterocystous); per-event and pooled | high |
| W4 | Which specific protein families expand, recurrently? | pre-defined sets: ELM motif-binding domains, cell-death machinery, cell-type-specification TF families, adhesion/ECM/signalling; plus Pfam-wide per-family recurrence with transposon masking; architectures and GO functions reported | medium–high |
| W5 | Does structure help Ancestra at deep nodes, and are stem novelties structurally new? | ProstT5 3Di on GPU for all proteins; Ancestra `--family-rescue 3di` vs default on R1 and R2; Foldseek (3Di+AA) search of stem-gained families against all proteins outside the clade | medium |
| W6 | Is there a structure-level signature of changed constraint? | per-branch 3Di vs amino-acid divergence of single-copy orthologs (structural substitutions per sequence substitution), multicellular vs sister branches | medium (exploratory) |
| W7 | Do stem novelties and retained duplicates join ancient protein complexes? | Cox et al. 2024/2026 LECA complex map (CF-MS; Zenodo 10.5281/zenodo.17716736) mapped via UniProt accessions | low |

## 2. Definitions

- **Multicellularity state** per species: unicellular / simple (colonies, filaments, mycelia, aggregation; little differentiation, no 3D
  tissue) / complex (3D tissues with differentiated cell types or a germ–soma split), following Knoll 2011 and Nagy et al. 2018; each assignment
  carries a reference; ambiguous cases flagged and tested in sensitivity analyses.
- **Event**: a state change on a branch: uni→simple, uni→complex (animals, where no simple intermediate is sampled), simple→complex.
  Coded per panel in `data/dense/<R>/transitions.tsv` before the run.
- **Gains**: as in round 1, the primary gain measure is *sequence-novel* families (no DIAMOND hit outside the clade in all 187 round-1
  proteomes plus the panel), because raw Ancestra gains are inflated 2.1–2.6× on simulations. Duplications and losses use Ancestra HOGs and are
  reported with the W1 agreement rate for that kind of event.
- **In-run controls**: all non-event internal branches of the same panel tree with ≥2 descendant species; each event is compared with
  controls by the negative-binomial residual model of round 1 (covariates: log branch length, log descendants, transcriptome fraction, BUSCO),
  fitted per panel. Depth imbalance is reported as a standardised mean difference for every comparison.

## 3. Tests

Round 2 is mostly exploratory. A small confirmatory family is fixed here, before any round-2 output exists, and Holm-corrected together:

- **R2-P1 (W3)** Simple→complex events carry more sequence-novel families than in-run controls. Statistic: random-effects meta-analysis of per-event
  log rate ratios (residual model) across the nested events in R2–R6, R8. Caveat: round 1 saw p = 0.016 on three of these events; the dense
  panels add new species but are not independent data, so this is a partial replication.
- **R2-P2 (W4)** Motif-binding domains (ELM interaction-domain Pfam list, 200 families, downloaded 2026-10-29) are expanded in multicellular
  species: phylogenetic regression of per-species copy number (log, offset = log proteome size) on state (uni/simple/complex) over all species,
  with a phylogenetic covariance from the round-1 tree (PGLS); and the sign of multicellular-minus-sister contrasts across independent events.
- **R2-P3 (W4)** Cell-death machinery (Pfam list fixed in `data/r2_sets/cell_death.tsv` with citations before results) — same test.
- **R2-P4 (W4)** Cell-type-specification TF families (`data/r2_sets/tf_celltype.tsv`, with citations) — same test.
- **R2-P5 (W6)** Ratio of structural (3Di) to sequence (amino-acid) divergence differs between multicellular and sister branches — sister-clade
  contrasts across events, sign/Wilcoxon test.

Everything else (per-family recurrence scans, architectures, W1, W5, W7, single-event tables) is exploratory and reported as effect sizes with
confidence intervals. Pre-specified decision rule for W5: the 3Di-rescued Ancestra output replaces the default for R1/R2 only if it improves the
deep-node positive controls (C03/C06/C07) and does not lower EukaryOMA agreement; otherwise structure is used only in post-hoc analyses.

## 4. Controls and checks

- Positive controls: round-1 controls re-scored on each relevant panel; now separated into Ancestra-event controls and per-species-count controls.
- Negative controls: per-panel label permutation of events among matched internal branches; transposon-domain mask in all Pfam scans
  (list in `data/r2_sets/transposon_mask.tsv`).
- Method checks: ProstT5 GPU 3Di agrees 97.1% per residue with Foldseek's own ProstT5 implementation (5,000 proteins); Ancestra panel runs
  pass `verify_orthoxml`; a pilot on each panel measures the largest family-graph component before the full run (LESSONS rule 2).
- Every headline number re-derived by an agent that did not produce it.

## 5. Compute plan (deadline ≈ 2026-10-30 21:00)

| round-2 t (h) | milestone |
|---|---|
| 0–2 | checkpoint, review, errata, design + one review round; GPU ProstT5 started (all 269 round-1 proteomes, ~17 h at 48 proteins/s) |
| 1–6 | dense panel data (agent), W1 EukaryOMA cross-validation, W4 set definitions, Pfam/disorder for new species |
| 4–14 | Ancestra panel runs (`accurate`, 16 threads, 2–3 concurrent), 3Di-rescue variants for R1/R2 when their 3Di is complete |
| 10–22 | W2–W7 analyses in parallel agents |
| 22–27 | verification, figures |
| 27–30 | report (round-2 sections integrated into the technical report; errata section) |
