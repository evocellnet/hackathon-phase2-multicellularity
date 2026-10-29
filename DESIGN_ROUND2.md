# DESIGN — Round 2: resolving single events and the step from simple to complex multicellularity

Version: v1 (2026-10-29 16:15, round-2 t ≈ 1.1 h), revised after one review round (statistics, evolutionary biology,
compute/AI; reviews/round2_design_review_*.txt). Change log in §7. v0 text below is superseded where §6 says so. Round 1 is frozen at git tag `phase1-final`.
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

## 6. Revisions after review (v1; these override §1–§5 where they conflict)

**Confirmatory family (Holm over 3, one-sided α = 0.05; frozen at tag `prereg-r2` before any round-2 Ancestra output):**
- **R2-P1** Simple→complex events carry more sequence-novel families than matched in-run controls. Unit = event; *confirmatory events are
  only those without any round-1 novelty value* (round 1 already saw Pezizomycotina, Agaricomycetes, Embryophyta in its nested set and
  Florideophyceae, Phaeophyceae, heterocystous cyanobacteria in exploratory tables). Expected confirmatory events: Volvox germ–soma
  (Volvocaceae), Fucophycidae (parenchymatous brown algae), Astrephomene (independent soma) — the exact list is fixed in
  `data/r2_prereg/p1_events.tsv` when the panels are built and before their Ancestra output exists. Statistic: per-event one-sided
  permutation p-value against the in-run controls (event label permuted among internal non-crown branches matched for depth class),
  combined by Stouffer's method. Effect sizes: per-event rate ratio (95% CI) and a pooled ratio with modified Hartung–Knapp CI
  (q* = max(1, q); Röver et al. 2015), τ² and prediction interval; nested events handled by one effect per panel in a sensitivity analysis.
  All simple→complex events (seen and unseen) are pooled as a *secondary* estimate. Rough power with 3 events: 18% at 2-fold, 41% at 3-fold —
  we expect an estimate, not a decision.
- **R2-P2** Motif-binding domains (ELM LIG/DOC/DEG/TRG interaction domains mapped to Pfam 38.2: 115 primary families; WD40/TPR/ARM in a
  sensitivity set; 41 ELM entries given only as InterPro/SMART IDs not used) are more abundant in multicellular than unicellular species.
- **R2-P3** Cell-death machinery (22 Pfam families incl. non-animal forms: metacaspases, VPE, HET/HeLo, NB-ARC, gasdermins, MazF) — same test.
- Statistic for P2/P3: PGLS (phylolm, Pagel's λ) of the equal-weight mean of standardised per-family log(count + 0.5) on a one-sided
  multicellular-vs-unicellular contrast, with log proteome size and proteome disorder fraction as covariates, on genome-based proteomes with
  BUSCO ≥ 50% of the round-1 187-species set on the round-1 phylogenomic tree. Sign must survive leave-one-origin-out; sensitivity:
  all proteomes, phylogenetic NB GLMM, complex-vs-simple contrast, sets including repeats. Lineage-restricted families are analysed only within
  panels where they occur.
- Moved to exploratory: **P4** (cell-type TF families; 26 Pfam incl. KNOX, RWP-RK, velvet, HetR: TF expansion in animals and plants is already
  known and round 1 used these families) and **P5** (3Di vs amino-acid divergence: ProstT5 3Di is a function of sequence, so the ratio partly
  measures the model; to be calibrated against AlphaFold-derived 3Di, restricted to confident Pfam-domain regions and compared within matched
  amino-acid distance bins; one-sided expectation = higher on multicellular branches, read together with round-1 H9).

**States.** Two binary axes coded per species with references: D (germ–soma split or ≥3 cell types) and T (3D tissue with cell–cell
connections); complex = D or T (primary), T-only as sensitivity. Volvox and Nostocales are D-only; Fucophycidae, Florideophyceae, land plants,
fruiting bodies are T. Fruiting bodies arose 8–11 times in fungi (Nagy et al. 2018): the Pezizomycotina stem may hide several origins.

**Panels.** Ancestra preset `fast` (as round 1; the 2.1–2.6× gain-inflation calibration and the 2.2× rule were made on `fast`), with the verified
patched Leiden wrapper (scripts/ancestra_patch/), a 30-min family-graph pilot per panel, ≤40 proteomes / ≤700k proteins per panel. Priority
R2, R4, R3, R6, R5, R1, R8, R7. Go/no-go: a panel not launched by round-2 t = 8 h is dropped and recorded.

**Controls/models.** One NB GLM fitted on the control branches of all panels (panel intercepts, shared slopes; covariates log branch length,
log descendants, transcriptome fraction, BUSCO, distance to nearest outgroup, which proxies homology-detection failure; Weisman et al. 2020);
branches inside multicellular crowns excluded from controls; the sister-clade stem reported as the depth-matched comparison; SMD reported.

**Scans.** Per-family NB GLM, multicellular clade vs sister, log-proteome offset, BH q < 0.1 within event after transposon masking (224 Pfam,
`data/r2_sets/transposon_mask.tsv`); recurrence (significant in ≥2 events) calibrated by permuting event labels; fold-changes shrunk (ashr)
before ranking; architectures and GO terms (pfam2go) reported for top families.

**Structure (W5).** ProstT5 3Di (own GPU implementation; windowing 1,000/250) validated against AlphaFold-derived 3Di of round-1 AFDB models
by length/disorder bin and checked for integrity (one record per ID, 3Di length = sequence length). `--family-rescue 3di` tested on R1 only, and only
if R1's 3Di is complete by round-2 t = 14 h; decision rule as §3. Foldseek on CPU (its GPU build fell back to CPU here); no MMseqs2-GPU (would
change the novelty definition).

**Independence.** EukaryOMA shares 122 input proteomes with us, so W1 is a cross-method check, not independent data; W1 is reported separately
for overlapping and non-overlapping species.

**Process.** Round-2 outputs only under results/r2/ and runs/r2/<panel>/, with refuse-if-exists guards and DONE.json (md5); every new analysis
gets a D-entry with its time cost; report skeleton by t = 10 h; analysis freeze t = 24 h; numbers in the report come from script output files and
are re-derived by a second agent.

## 7. Change log
- v1: confirmatory family reduced from 5 to 3 tests (P4, P5 exploratory); P1 restricted to unseen events; Stouffer + modified HK; pooled NB
  controls; PGLS specification; D/T state axes; ELM restricted to LIG/DOC/DEG/TRG; curated cell-death and TF lists; `fast` preset; panel caps,
  pilots, priority and go/no-go; 3Di validation against AFDB; W1 independence caveat; output guards.
