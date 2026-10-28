# PREREG — pre-registration of the primary analyses

Frozen before any eukaryote Ancestra output was parsed (the eukaryote run was still in its family stage at freeze time;
the prokaryote run, which is descriptive only, was also unparsed). Git tag: `prereg-v1`. Hashes of every frozen input:
`PREREG_HASHES.tsv`; verify with `python scripts/prereg/check_prereg.py`. Design: DESIGN.md v3.1 (§4.6 binding).

## What is frozen
| item | file |
|---|---|
| design and binding statistical specification | DESIGN.md (v3.1) |
| transitions, states, M/P nodes, path stems, analysis sets | data/transitions.tsv (produced by scripts/prereg/apply_eligibility.py) |
| guide trees | data/guide_tree_euk.nwk, data/guide_tree_prok.nwk |
| species sets | data/species_final.tsv, data/species_final_prok.tsv |
| positive controls, node mapping, definitions | data/positive_controls.tsv, data/positive_controls_nodes.tsv, data/prereg_definitions.tsv |
| toolkit categories | data/toolkit_pfam.tsv |
| H12 modules and directions | data/h12_sets.tsv, data/h12_expected_direction_frozen.tsv |
| literature scorecard template | data/literature_scorecard_template.tsv |
| power / attainable p / minimum detectable effects | results/prereg/power_table.tsv |

## Primary family (Holm across these 5 tests; one statistic each; DESIGN §4.6)
Primary transitions (7 eukaryotic origins): Metazoa; fungal hyphae (path stem Eumycota → Amastigomycota); Dictyostelia;
Streptophyta (Klebsormidiophyceae + Phragmoplastophyta); Rhodophyta (Eurhodophytina); PX clade filaments (brown-algal lineage);
volvocine colonial lineage (path stem Chlamydomonas_TGV → Gonium + Volvocaceae).

| test | statistic | direction | null |
|---|---|---|---|
| H1/H2 | sequence-novel gains on the stem, NB-GLM residual vs control pool U | two-sided | permute (10,000, seed 11) |
| H4 | logit share of stem-duplicated HOGs that are TF-DBD HOGs minus baseline at u (d ≥ 10) | one-sided (expansion) | permute |
| H5 | as H4 for adhesion/ECM ∪ signalling | one-sided | permute |
| H9 | median per-gene log rate ratio, P→M vs P→S tips (BUSCO markers, fixed topology) | two-sided | permute (controls with two ≥2-tip child clades) |
| H12 | signed flagellar-module loss on stems; loss expected for fungal hyphae and Dictyostelia | one-sided | permute |

Secondary (same statistics, not in the Holm family): nested simple → complex transitions (Pezizomycotina, Agaricomycetes,
Embryophyta). Prokaryotic transitions (Myxococcota, Streptomyces, magnetotactic Desulfobacterales, Methanosarcina): descriptive only
(|U| = 2). Everything else is exploratory and labelled as such.

## Power (from results/prereg/power_table.tsv; placeholder branch lengths)
With 7 stems and |U| = 30 controls the design detects only large effects: H1/H2 minimum detectable effect at 80% power and α = 0.05 is
≈ 2.4-fold (NB dispersion 0.3) to ≈ 5.9-fold (dispersion 1.0); H4 ≈ 2.8-fold; H5 ≈ 2.2-fold; H9 ≈ 1.26–1.69-fold rate ratio. Power to detect
the 1.5-fold smallest effect of interest is 13–35% for H1/H2, H4 and H5. Non-significant results will therefore be reported as
"not detected at this power" with the compatible interval, never as absence.
Known calibration issue: in the power simulations the H1/H2 residual statistic rejected 8.5% of null datasets at α = 0.05 (200 simulations;
≈ 2 SE above nominal). H1/H2 p-values will be reported with this caveat and with the literal-scheme sensitivity.

## Positive-control pass criterion
≥ 11 of 16 eukaryotic controls recovered and every high-confidence failure explained; every failure reported.

## Deviations
Any change after the freeze is added here with date, reason and whether it was made before or after seeing results.

- **Dev-1 (2026-10-28 ~14:45, BEFORE any eukaryote output existed).** H1/H2 null distribution: parametric bootstrap from the NB-GLM
  fitted on the control pool (simulate counts for all branches from the fitted model, refit, recompute T) replaces the `permute` scheme
  for H1/H2 only. Reason: final calibration of the statistics engine showed that permuting residuals is anti-conservative when stems
  extrapolate beyond the controls' covariate range (false-positive rate 0.083 at a 1-SD covariate shift and 0.283 at 2.5 SD, α = 0.05),
  and our stems are deeper than the controls (depth SMD 1.77); the parametric bootstrap stays at 0.040–0.053
  (results/validation_branchstats/null_calibration.tsv). The `permute` H1/H2 p-value is reported as a sensitivity analysis.
  Also corrected: the null-calibration numbers quoted in DESIGN §4.6 v3.1 (14.5%/4.3%, 400 replicates) are superseded by
  13.8%/4.5% for the literal scheme and 5.1%/1.3% for `permute` (1,000 replicates); DESIGN.md is left unchanged because it is hashed.
- **Dev-2 (2026-10-28 ~18:15 (commit 18:19), BEFORE any eukaryote output existed).** Gain-ladder DIAMOND sensitivity `--very-sensitive` (two passes)
  instead of `--ultra-sensitive`: on 500 real proteins it recovers 99.65% of ultra-sensitive hits at about one third of the cost, and an
  ultra-sensitive recheck of 200 prokaryote 'novel' gains found 0 outside-clade hits (results/gainladder/diamond_timing.json). Reason:
  the eukaryote search would otherwise cost an estimated 150–700 CPU-h. The eukaryote ladder is run on the branches the tests use
  (all stems of primary, secondary and exploratory transitions, reversals, and every branch of the control pool U) rather than on all
  internal branches; the filter is identical on stems and controls.
- **Dev-3 (same time; outcome of a pre-registered acceptance test, not a choice).** The split-repair level (DESIGN §3.2 level 2)
  failed its acceptance criterion on the simulations (gain ratio stays 2.12–2.54× vs the ≤1.3× threshold; loss/dup correlations
  unchanged), so it is NOT used; the ladder is HOG gain → sequence-novel gain → architecture-novel gain. Simulation precision of
  sequence-novel gains is 0.42–0.56 (recall 0.91–0.97) because simulated deep families diverge beyond DIAMOND detection; novel-gain
  counts are therefore reported as upper bounds, and H1/H2 conclusions rest on the stem-vs-control comparison under the same filter.
- **Dev-4 (2026-10-28 19:47, BEFORE any eukaryote output existed; interpretation rule only, tests unchanged).** The pre-registered
  calibration on our tree (DESIGN §3.2; results/calibration_tree/summary.md; one draw per setting, 7 stems vs 9 controls, simulator
  without true gene birth) found spurious HOG gains inflated 1.7–2.2× more on stems than on controls after adjusting for branch length
  and clade size (losses/duplications unaffected; shift on the test's percentile scale +0.06 to +0.11, within noise). Interpretation
  rule fixed now: a "burst of novelty" (H2) is claimed only if H1/H2 is significant after Holm AND the pooled observed/expected
  effect for sequence-novel gains exceeds 2.2 (the upper calibrated bias); a significant effect below 2.2 is reported as "not
  distinguishable from method bias". A deficit or a null result is not explained by this bias and is interpreted as usual.
- **Clarification (2026-10-28 20:09, corrective entry; nothing changed).** "Before any eukaryote output existed" in Dev-1–Dev-4 means:
  no eukaryote HOGs, branch events or test statistics existed or had been inspected. Precisely: Dev-1 (~14:45) preceded all eukaryote
  Ancestra stages; Dev-2/Dev-3 (18:19) and Dev-4 (19:47) were made after the eukaryote *family* stage had finished (18:41 is when its
  completion was noticed; families were not parsed or inspected) and while scoring was running. The primary H9 rate test (BUSCO markers,
  independent of Ancestra) and the exploratory gene-level RER/Kr/Kc had already been run (15:21–15:40); Dev-1–Dev-4 concern only the
  gene-content tests (H1/H2) and gain ladder, not H9. Earlier timestamp edits (Dev-2, Dev-4) corrected values written without checking
  the clock; they are visible in git history.
