# DESIGN — Genomic transitions to multicellularity across the tree of life

Version: **v3.1** (after review rounds 1–3 and pre-results calibration of the statistics engine: reviews/round1/*.md, reviews/round2/*.md; earlier versions in docs/design_history/). §4.6 is the binding statistical specification and overrides §4.1–4.2 where they differ. Owner: PI agent. Change log in §9.
Fixed by the task: use Ancestra to reconstruct ancestral gene content before and after documented transitions to multicellularity,
run a broad set of analyses, compare with existing knowledge, and assess whether shifts in evolutionary pressure are detectable
in sequences or structures.

## 0. The main risk, and how the design answers it
The most likely misleading headline is *"transitions are bursts of new gene families, with convergent adhesion/signalling/TF gains"*
(external PI review). Four mechanisms push towards it: (i) Ancestra splits deep families, so each split looks like a gain on a long
stem branch (our simulation: gains inflated 2.1–2.6×, losses/duplications accurate; results/validation_hogparser/summary.md); its
fixed family-size cap (7,500 genes) will split large kinase/TF/adhesion superfamilies; (ii) unicellular sisters are more often
transcriptomes; (iii) toolkit lists come from animal literature; (iv) many transitions rest on one species.
Answers: a gain-evidence ladder with sequence-level novelty as the primary gain statistic (§3.2); per-transition comparison against
matched control branches of the same class (§4.1); rarefaction and data-type sensitivity (§4.4); a Pfam-wide analysis before any
toolkit analysis and a toolkit that includes non-animal families (§3.3); single-species transitions never reported as raw gains.

## 1. Questions
- **Q1 Gene-content change.** How much gain (new families), duplication and loss happens on branches where multicellularity arises,
  relative to matched branches where it does not? Burst of novelty or re-use (co-option)?
- **Q2 Function.** Which functions are gained/expanded/lost at transitions, and which recur across independent transitions beyond chance?
- **Q3 Pre-adaptation (co-option).** What fraction of the families that characterise each multicellular clade already existed in its
  last ancestor in the prior state?
- **Q4 Domain architecture.** Do transitions create new domain combinations more often than control branches?
- **Q5 Evolutionary pressure.** Do multicellular lineages evolve their proteins faster or slower (relaxed or intensified constraint),
  and do protein traits (length, disorder) shift — measured on shared single-copy families to avoid annotation artefacts?
- **Q6 Structure.** Are sequence-novel gains structurally novel, or diverged copies of older folds?
- **Q7 Literature concordance.** A numeric scorecard, written before results, of which published claims our data reproduce.

## 2. Transitions (state changes, not "origins")
Each transition is coded in `data/transitions.tsv` as a **state change** `from_state → to_state` with states
{unicellular, simple, complex} and type labels (clonal/aggregative; obligate/facultative), plus `reversal` events
(e.g. Saccharomycotina yeasts, Zygnematophyceae). **M** = crown node (MRCA of sampled species in the new state; a leaf if only one
species); **P** = parent of M = last sampled ancestor in the prior state. `primary_eligible` = 1 only if M is an internal node with
≥2 sampled species in the new state AND ≥2 sampled species in the sister lineage(s), the state assignment is not grey-zone, and data
are not low quality. Grey-zone (choanoflagellate rosettes, *Capsaspora* aggregation, *Physarum* plasmodium, ichthyosporean coenocytes),
single-species (*Fonticula*, *Acrasis*, *Copromyxa*, *Neolecta*) and very-low-quality (*Guttulinopsis*) transitions are analysed as
**exploratory** only. *Sorogena* and *Sorodiplophrys* have no public proteomes and are not analysed.
Corrections adopted from review (evo_biology.md): two Volvocales origins (Lindsey et al. 2021; Goniaceae not monophyletic); streptophyte
M1 = Klebsormidiophyceae + Phragmoplastophyta, Embryophyta = nested simple→complex; filamentous Schizocladia/Phaeothamnion/Tribonema
→ Phaeophyceae = simple→complex; Dictyostelia = aggregative simple; fungal complex fruiting bodies = nested simple→complex.
The table and the guide tree are frozen in PREREG (§7) before any Ancestra output is parsed.

## 3. Data and core computation
### 3.1 Proteomes and guide tree
- **Eukaryotes:** ~200 proteomes, ≤ 3.0 M proteins (150 GB memory limit), from UniProt 2026_03 reference proteomes, EukProt v3,
  LukProt v1.5.1, NCBI (Phaeoexplorer). **Prokaryotes:** ~80 UniProt proteomes, separate run (also the pipeline dry run).
- One protein per gene; transcriptomes collapsed with `cd-hit -c 0.95 -G 0 -aS 0.9`; ≥50 aa; species-prefixed IDs; uniform BUSCO 6.1
  (odb10) on final files; transcriptomes with BUSCO D or F > 15% flagged.
- Guide tree: NCBI taxonomy edited to current consensus, every edit referenced (data/guide_tree_edits.md); children sorted
  deterministically. Two to three variants (root position; contested placements) are re-merged from the Ancestra checkpoints if time allows.
- Resolved species tree with branch lengths: IQ-TREE (LG+G4, topology constrained to the guide tree) on concatenated BUSCO
  eukaryota_odb10 single-copy markers. Used for branch-length covariates and relative-rate tests; never to relabel transitions.

### 3.2 Ancestral gene content (Ancestra) and the gain-evidence ladder
- **`fast` is the pre-registered primary preset.** Proteins capped at ≤ 3.0 M by the written drop rule in DATA.md; the run is killed
  at 145 GB RSS. `accurate` re-uses the `fast` S0/S2 stages (the proteome directory must not change) and is a sensitivity analysis.
- Events per branch from `mcell.hogs` (validated: exact on true HOGs; on Ancestra HOGs losses/dups r ≥ 0.98, gains 2.1–2.6× inflated).
- **Gain-evidence ladder** for every HOG gained on an internal branch:
  1. *HOG-level gain* (raw Ancestra; upper bound, never a headline);
  2. *split-repaired gain*: gained HOG re-merged with an older HOG "lost" on the same branch when the novelty search (3) finds
     reciprocal best hits between them (Ancestra's within-family scores cannot be used: split HOGs fall in different families);
     accepted only if on the simulations the gain ratio falls to ≤ 1.3× with loss/dup correlations unchanged;
  3. **sequence-novel gain (primary)**: one DIAMOND search (`--ultra-sensitive -k 0 -e 1e-5 --query-cover 30`; `-k 0` is required,
     otherwise in-clade paralogues fill the hit list) of up to 3 representatives per gained HOG against all eukaryote + prokaryote
     proteomes, hits then filtered by clade membership; a gain is sequence-novel if no representative hits any proteome outside the
     clade. Run on stems AND on every control branch in the pool U (§4.6), so the filter acts identically on both;
  4. *architecture-novel gain*: Pfam architecture absent outside the clade (secondary).
- Duplications and losses use Ancestra events directly. Presence supported only by transcriptomes requires ≥2 species; losses ignore
  transcriptome-only taxa (treated as unknown).
- **Calibration on our tree:** Ancestra's simulator on the guide topology pruned to ~60 taxa, to measure gain/loss/dup bias per branch
  class (stem vs other, root children), reported with results.

### 3.3 Annotation
- Pfam-A 38.2 (pyhmmer, gathering cutoffs) on all proteins; metapredict v3 disorder; protein length; domains per protein; ordered
  Pfam architectures.
- **Pfam-wide analysis first** (all domains, no prior list). **Toolkit categories** (adhesion/ECM, signalling, kinases/pTyr, TF DBDs,
  programmed cell death, cell cycle/cytokinesis) from review Appendix A2, including plant/algal TF and non-animal adhesion families;
  immune LRR proteins excluded; C2H2 reported separately. Frozen in `data/toolkit_pfam.tsv`.
- Transmembrane/signal peptide (H7, exploratory): a documented open predictor if installable within 30 min, else a hydropathy heuristic.

## 4. Statistics
### 4.1 Unit, classes and the primary comparison
- Unit = branch of the guide tree. Three classes: **transition stems** (P→M), **within-multicellular**, **within-unicellular**
  branches. Nested transitions (e.g. Embryophyta, Agaricomycetes complex) form a separate set. Internal branches are the primary
  analysis; terminal branches are analysed separately ("supported" gains are undefined on them).
- **Per-transition test:** for each primary-eligible transition, the observed statistic on its stem is compared with ≥200 matched
  within-unicellular control branches (matched on depth class, log branch length tertile, descendant species tertile, fraction
  transcriptome descendants; matching relaxed stepwise and reported if <200 exist); report observed/expected with prediction interval
  and the empirical percentile.
- **Across transitions:** random-effects summary of log(observed/expected) (forest plot), sign test and Stouffer's Z, with
  leave-one-transition-out and leave-one-supergroup-out estimates. Model-based view: NB-GLM
  `events ~ class + log(branch_length) + log(n_desc) + frac_transcriptome_desc + mean_BUSCO_desc + root_child`, with offset
  log(families at parent) for losses and duplications; sensitivity without branch length.
- **Primary p-values** from stratified relabelling: transition labels re-assigned to matched control branches (10,000 draws, seed 11,
  without replacement), recomputing the full statistic.
- **Smallest effect of interest:** 1.5-fold. Minimum detectable effects are reported; non-significant results are worded
  "not detected at this power", never "absent".

### 4.2 Primary hypotheses (Holm correction across the 5 tests; one statistic each)
| id | hypothesis (LITERATURE.md §7) | primary statistic |
|---|---|---|
| H1/H2 | co-option vs burst of novelty | sequence-novel gains per stem, observed/expected vs matched controls (§4.1); co-option index (§4.3) as companion estimate |
| H4 | TF families expand at transitions | compositional: share of stem duplications in TF-DBD HOGs vs matched controls; complex×transition interaction as estimate only |
| H5 | adhesion/ECM/communication expand convergently | per class (Appendix A2 adhesion/ECM + signalling), stem duplication share percentile vs controls; statistic = sum over transitions of percentiles, vs relabelling null |
| H9 | drift / relaxed constraint in eukaryotes | per-transition relative-rate test (§4.5); across transitions Wilcoxon + Hodges–Lehmann interval |
| H12 | loss of unicellular functions (re-scoped) | flagellar-module losses on stems with **pre-registered direction per transition** (loss expected: fungal hyphae, Dictyostelia, Zygnematophyceae+Embryophyta lineage; not expected: Metazoa, Volvocales, Phaeophyceae, Ulva, aggregative amoebae) — concordance count vs relabelling |

H8 (clonal > aggregative) is **demoted to exploratory** (aggregative transitions are mostly single-species; confounded).
Secondary/exploratory (labelled; BH within the logged exploratory family): H3, H6, H7, H8, H10, H11, H13, the prokaryote gene-count
side of H9, and reversals (yeasts, Zygnematophyceae: genes gained at the earlier transition lost preferentially).
Cut for time: H14 (HGT), H15 (site convergence).

### 4.3 Co-option index (Q3)
Families defined independently of gain timing: Pfam domains and toolkit HOGs with significantly higher copy number in the
multicellular clade than in its sister (per transition, family-aware test) → fraction present at P and at P's parent, compared with
the same fraction for matched control nodes.

### 4.4 Robustness (all reported)
Rarefaction without re-running Ancestra (multicellular clades subsampled to 2–3 species; sisters and outgroups dropped one at a time,
following Richter et al. 2018); genome-only subset; `fast` vs `accurate`; guide-tree variants; orthology-free Pfam-count version of
Q1–Q2; positive controls kept out of the test families.

### 4.5 Evolutionary pressure (Q5)
- **Genome-wide rate shift (H9):** per transition, outgroup-anchored relative-rate test: from gene-wise trees of shared single-copy
  BUSCO/HOG markers, distance from the outgroup-anchored P node to multicellular vs sister-lineage tips; per-transition median
  log-ratio with bootstrap CI over genes. (RER normalisation removes genome-wide shifts, so RER is not used for this question.)
- **Gene-level shifts:** RERconverge-style relative rates (python re-implementation, validated on simulated trees with planted
  shifts) for up to ~3,000 near-single-copy HOGs on the fixed topology; foreground = stems + multicellular crowns; significance by
  permulations; BH-FDR; enrichment of top genes by Pfam/toolkit (exploratory).
- **Protein traits:** length and disorder on shared single-copy families (same HOGs on both sides), sister-clade contrasts per
  transition, Wilcoxon + Hodges–Lehmann; whole-proteome versions secondary (confounded by data type).
- Optional: Kr/Kc radical/conservative substitution ratio; burial-stratified rates with AlphaFold models; Volvocales-only RELAX.

### 4.6 Binding statistical specification (v3; from reviews/round2/statistics.md)
**Control pool.** U = internal branches b=(u→v) with u, v unicellular; |tips(v)| ≥ 2; |tips(sibling(v))| ≥ 2; v not a root child; b not
adjacent to any stem (not the edge into P, not P→sister); not inside multicellular or grey-zone clades. Expected |U| ≈ 40–60.
Covariates x = z-scored [depth, log branch length, log n_desc, log n_sister_tips, frac_transcriptome_desc, frac_transcriptome_sister]
(+ log d for H4/H5). For each primary stem i, N_i = the k = min(20, |U|) Mahalanobis-nearest branches in U; standardised mean
differences reported, |SMD| > 0.5 flagged; optional regression-adjusted expectation for poorly matched deep stems.
**One relabelling statistic for every primary test** (equal weight per transition):
```
u(b,N) = (#{c∈N\b: s_c<s_b} + 0.5·#{c∈N\b: s_c=s_b} + 0.5) / (|N\b| + 1)
T_obs  = mean_i u(stem_i, N_i)
10,000 draws (seed 11): draw c_i ∈ N_i without replacement across i (reject duplicates); T_d = mean_i u(c_i, N_i \ c_i)
p (one-sided where a direction is stated) = (1 + #{T_d ≥ T_obs}) / 10,001 ; two-sided uses |T−½|
```
Effect size per transition: log O/E_i = log((y_i+0.5)/(mean_{N_i} y + 0.5)); pooled by REML random effects (within-transition variance
from the controls); "prediction interval" = 5–95% range of the controls. Stouffer and sign tests are descriptive only.
**Statistics per primary test (s_b):**
- *H1/H2*: s_b = number of sequence-novel gains on b. Two-sided (burst vs. deficit).
- *H4*: D_b = set of HOGs with ≥1 duplication on b; d = |D_b|; a = |D_b ∩ C| with C = HOGs whose members carry a frozen TF-DBD domain
  (excluding C2H2) in ≥50% of members; s_b = logit((a+0.5)/(d+1)) − logit(π_u), π_u = share of HOGs present at u in C. d ≥ 10 required
  (stems with d < 10 are listed and excluded); branches flagged with whole-genome duplication excluded (sensitivity with them). One-sided.
- *H5*: as H4 with C = adhesion/ECM ∪ signalling (one union statistic; per-class results secondary). One-sided.
- *H9*: per marker gene g, IQ-TREE `-te` (LG+G4) on the guide tree pruned to the taxa present; r_ig = log(L̄(P→M tips)/L̄(P→S tips)),
  L̄ = equal-split weighted mean path length; gene included if ≥1 outgroup, ≥2 M tips with MRCA M, ≥1 S tip, ≥50% site coverage;
  s = r_i = median_g r_ig (2,000 gene bootstraps for its CI); controls = nodes in U with two child clades of ≥2 tips, random orientation
  (seed 11). Two-sided. Reported as a **rate shift**; "relaxed constraint" only if a within-protein ratio (e.g. Kr/Kc or surface/core)
  supports it. Rerun on genome-only taxa. Caveat: r covers the whole post-P history; longer generation times bias r below 0.
- *H12*: directions re-derived on the frozen tree as "flagellar module present at P, absent at M" (data/h12_expected_direction_frozen.tsv;
  Chara/bryophytes/ferns lie inside streptophyte M1, so a stem loss is impossible there; with the fungal path stem of D7 all 26
  Amastigomycota are non-flagellate, so loss is expected on the fungal and Dictyostelia stems: 2 loss-expected stems → H12 primary). f_b = mean over the 4 modules of
  the fraction of module HOGs present at u that are lost on b (modules with < 50% of HOGs present at u skipped); s_b = w_i·f_b with
  w = +1 where loss is expected, −1 otherwise. **H12 is primary only if ≥ 2 loss-expected primary stems remain after re-derivation;
  otherwise exploratory** (decided when PREREG is frozen, before results).

**v3.1 amendments (made before any eukaryote output was parsed; decisions D6–D7 in RESEARCH_TREE.md):**
- *Relabelling scheme (D6).* The draw-without-replacement null of v3 is anti-conservative because the matched sets N_i overlap
  (|U| = 30): 14.5% / 4.3% rejections at α = 0.05 / 0.01 under the null on the real tree (results/validation_branchstats/null_calibration.tsv).
  The primary scheme is therefore `permute`: statistic values are permuted over stems ∪ all controls and u, T recomputed with the
  original matched sets (4.9% / 0.9%). The v3 scheme is reported as a sensitivity analysis.
- *Residual count statistic for H1/H2 (D6).* Matching is poor (stems are deep, controls shallow; depth SMD 1.77), so s_b for H1/H2 is
  log((y_b+0.5)/(ŷ_b+0.5)) with ŷ from an NB-GLM fitted on U only: y ~ log(branch length) + log(n_desc) + log(n_sister_tips) +
  frac_transcriptome_desc (leave-one-out for controls; Poisson/quasi fallback logged). Compositional (H4/H5), rate (H9) and
  module-loss (H12) statistics are unchanged.
- *Path stems (D7, audit F1/F2).* A stem may be a chain of consecutive edges from P to M; events and branch lengths are summed over the
  edges, and the edges are excluded from U. Volvocine: Chlamydomonas_TGV>TGV_clade>Gonium_Volvocaceae (Lindsey et al. 2021: TGV_clade
  inferred unicellular; it exists only through the 32%-complete Tetrabaena). Fungal hyphae: Eumycota>Blastocladiomycota_Amastigomycota>
  Amastigomycota (Kiss et al. 2019 place hyphae at Zoopagomycota+Dikarya; Blastocladiomycota mostly unicellular).
- *Prokaryotes.* Their tree yields |U| = 2 control branches, so prokaryotic transitions are analysed descriptively (no relabelling test).
- *Primary family.* 7 eukaryotic origins (data/transitions.tsv analysis_set = primary); Holm across the 5 primary tests (H1/H2, H4, H5,
  H9, H12); the attainable-p and MDE table is results/prereg/power_table.tsv. `analysis_set` supersedes the older `analysis_tier` column.

**Multiplicity.** Holm across the primary tests that remain (≤ 5). PREREG lists per test n, k, the minimum attainable p and a simulated
minimum detectable effect (planted multiplicative effect; 80% power at α = 0.01 and 0.05).
**Reporting.** Per test: n, k, estimate with 95% CI, p with Monte-Carlo SE, Holm p, minimum attainable p. Wording: "no ≥ 1.5-fold effect"
only if the CI lies within (1/1.5, 1.5); otherwise "inconclusive; compatible with [a, b]".
**Co-option index (descriptive).** F_i = domains with log2((x̄_M+1)/(x̄_S+1)) ≥ 1 present in ≥ 50% of M species; C_i = fraction of F_i
present at P (Pfam: ≥1 genome or ≥2 transcriptomes in S ∪ outgroups; HOGs: split-repaired presence); report C_i − median(C_{N_i}) and u.
Also a toolkit-based version (toolkit families present in the M clade). Descriptive only (presence in S implies presence at P by construction).
**Expansions (H4/H5) must agree with the orthology-free Pfam-count version** to be called robust.

## 5. Structure (Q6, exploratory, bounded)
Sequence-novel gains with an AlphaFold DB model (EBI per-protein files; mean pLDDT ≥ 70) in ~10 transitions with UniProt species on
both sides; Foldseek against AFDB models of the P-side species' proteins (TM ≥ 0.5, E ≤ 1e-3) → fraction "structurally old"; same for
matched control gains. ASR for ≤ 2 families (IQ-TREE -asr), descriptive.

## 6. Positive and negative controls (frozen in PREREG)
- **Positive:** review Appendix A1 (16 eukaryotic + 2 prokaryotic controls: gain, co-option, loss and expected-null kinds) with expected
  node, direction, confidence and source, in `data/positive_controls.tsv`. Pass criterion fixed in advance: ≥ 11/16 eukaryotic controls
  recovered and every high-confidence failure explained; every failure reported.
- **Negative:** relabelling of transitions onto matched unicellular branches; trait permutation for rates; planted-effect and
  null-simulation tests for every analysis module; toy data for parsers.
- **Literature scorecard (Q7):** expected outcome per hypothesis × transition, written before results (`data/literature_scorecard.tsv`).

## 7. Process safeguards (from ai_failure_modes.md)
- **PREREG:** DESIGN.md, transitions.tsv, guide trees, positive_controls.tsv, toolkit_pfam.tsv, h12_sets.tsv, literature_scorecard.tsv,
  primary-test list and seeds are hashed in PREREG.md and tagged `prereg-v1` before any eukaryote Ancestra output is parsed;
  `scripts/check_prereg.py` verifies hashes. Later changes go to a dated *Deviations* section, never silently.
- **Numbers in the report** come only from committed `results/**/*.tsv` listed in `results/MANIFEST.tsv`; a lint flags literal numbers
  in report prose not traceable to the manifest.
- **Tests:** every analysis module gets a planted-effect test, a null test and a hand-computable case, reviewed by an agent that did
  not write it. `scripts/check_ids.py` asserts species codes agree across FASTA, trees and tables. Jobs write `DONE.json` only on success.
- **Citations:** every citation resolves to LITERATURE.md; a fresh agent audits ≥ 20 claim–citation pairs in the report.
- **Analysis registry:** every test run is appended to `results/analysis_registry.tsv` (hypothesis, statistic, n, p, primary/exploratory).
- **Re-derivation:** 5–8 headline claims re-derived by an agent that did not produce them.
- Seeds: [11, 23, 37, 41, 53]; 11 default.

## 8. Compute plan (24 CPU, 150 GB; START 12:21; target finish ≤ 36 h)
| t (h) | work |
|---|---|
| 0–2 | design rounds 1–2, dataset build, PREREG freeze |
| 2–3 | prokaryote Ancestra run (dry run); 40-proteome eukaryote timing/memory test |
| 3–10 | eukaryote `fast` run (12–16 threads) ‖ Pfam + metapredict (8–12 threads) ‖ BUSCO-marker concatenation; analysis code written and unit-tested on prokaryote outputs |
| 10–15 | `accurate` merge from checkpoints (if budget) ‖ gain-ladder DIAMOND ‖ RER gene trees |
| 15–28 | analyses Q1–Q6 |
| 28–34 | verification, re-derivation, figures |
| 34–40 | report |

**Drop order if behind schedule at t = 20 h:** (1) structure/ASR (§5); (2) `accurate` sensitivity run; (3) guide-tree variants;
(4) gene-level RER; (5) exploratory hypotheses beyond H3/H7/H8/H11. Primary tests, positive controls and verification are never dropped.

## 9. Change log
- v0 → v1: §1b mapping to literature hypotheses.
- v3 → v3.1: D6 (permute null, NB-residual H1/H2 statistic), D7 (path stems; fungal and volvocine stems), prokaryotes descriptive, H12 primary with 2 loss-expected stems (round-3 audit reviews/round3/prereg_audit.md).
- v2 → v3 (review round 2): binding statistical spec §4.6 (control pool U with sister covariates, k=20 nearest controls, one
  mean-percentile relabelling statistic, exact H4/H5/H9/H12 statistics, attainable-p and MDE table in PREREG, reporting wording);
  `fast` pre-registered primary; ≤3.0 M protein cut rule; DIAMOND `-k 0` single search with clade post-filter, also on control branches;
  re-merge via novelty-search hits; H12 conditional primary with re-derived directions; H9 called "rate shift"; AFDB downloads need a
  User-Agent header; drop order at t = 20 h (§8).
- v1 → v2 (review round 1, all five reviews): gain-evidence ladder with sequence-novel gains primary; split re-merge; transitions as
  state changes with primary eligibility; per-transition matched-control comparison + relabelling null; three branch classes,
  internal primary; compositional H4/H5; H8 demoted; H12 re-scoped with directions; H9 via relative-rate test, RER only gene-level;
  traits on shared single-copy families; Holm across 5 primary tests; cut H14/H15/LG+C20; positive controls replaced (Appendix A1);
  toolkit A2; memory cap 3.0 M proteins; `fast` first; PREREG hashing and report-number manifest.
