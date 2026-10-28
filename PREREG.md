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
(none yet — any change after the freeze is added here with date, reason and whether it was made before or after seeing results)
