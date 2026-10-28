# DESIGN — Genomic transitions to multicellularity across the tree of life

Version: **v2** (after review round 1: reviews/round1/*.md; v1 archived in docs/design_history/). Owner: PI agent. Change log in §9.
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
- Run `fast` first (timing/memory gated by a 40-proteome test: cut to ≤ 2.6 M proteins if predicted peak > 110 GB), then `accurate`
  from the same S0/S2 checkpoints if time allows; compare per-branch counts between presets.
- Events per branch from `mcell.hogs` (validated: exact on true HOGs; on Ancestra HOGs losses/dups r ≥ 0.98, gains 2.1–2.6× inflated).
- **Gain-evidence ladder** for every HOG gained on an internal branch:
  1. *HOG-level gain* (raw Ancestra; upper bound, never a headline);
  2. *split-repaired gain*: gained HOG re-merged with a reciprocal-best older HOG "lost" on the same branch (`mcell.remerge`; accepted
     only if on the simulations the gain ratio falls to ≤ 1.3× with loss/dup correlations unchanged);
  3. **sequence-novel gain (primary)**: no DIAMOND hit (e ≤ 1e-5, query cover ≥ 30%; ultra-sensitive) from up to 3 representatives
     against all proteomes outside the clade;
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

## 9. Change log
- v0 → v1: §1b mapping to literature hypotheses.
- v1 → v2 (review round 1, all five reviews): gain-evidence ladder with sequence-novel gains primary; split re-merge; transitions as
  state changes with primary eligibility; per-transition matched-control comparison + relabelling null; three branch classes,
  internal primary; compositional H4/H5; H8 demoted; H12 re-scoped with directions; H9 via relative-rate test, RER only gene-level;
  traits on shared single-copy families; Holm across 5 primary tests; cut H14/H15/LG+C20; positive controls replaced (Appendix A1);
  toolkit A2; memory cap 3.0 M proteins; `fast` first; PREREG hashing and report-number manifest.
