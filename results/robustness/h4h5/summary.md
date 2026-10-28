# H4/H5 robustness: orthology-free version (pre-registered) and duplication re-placement (post hoc)

Written by the PI from the analysing agent's report; all numbers are in the TSVs of this directory
(`pfam_counts_tests.tsv`, `pfam_counts_per_transition.tsv`, `replacement_validation_metrics.tsv`, `replacement_tests.tsv`,
`replacement_per_transition.tsv`, `replacement_branch_values.tsv`) and in results/analysis_registry.tsv.

## Main point
Neither check changes the primary conclusion: transcription-factor (H4) and adhesion/ECM + signalling (H5) expansions are not
detected as a convergent feature of the seven primary transitions. The animal stem (Metazoa) is the one transition that is extreme
for both classes in the orthology-free version.

## A. Orthology-free version (pre-registered robustness, DESIGN §4.6)
Statistic per branch: s_b = log((mean class-protein count over species below the branch + 0.5) / (mean over the sister lineage + 0.5))
minus the same log ratio for all proteins (removes proteome-size, whole-genome-duplication and transcriptome effects). Child and sister
sets equal the frozen sets of data/transitions.tsv. Same matched control sets and mean-percentile permutation test (10,000 draws,
seed 11, one-sided) as the primary analysis.

| test | orthology-free T / p / pooled effect [95% CI] | HOG-based (primary) p / pooled | agreement (same conclusion at α=0.05 and same sign) |
|---|---|---|---|
| H4 TF-DBD | 0.609 / 0.184 / +0.27 [−0.35, 0.88] | 0.361 / +0.15 | yes |
| H5 adhesion/ECM ∪ signalling | — / 0.421 / +0.07 [−0.37, 0.52] | 0.772 / −0.18 | conclusion yes, sign no → H5 direction not robust |

Sensitivities: equal-weight clade means H4 p = 0.343, H5 p = 0.883; genome-only species H4 p = 0.272, H5 p = 0.323.
Secondary (nested complex) set: H4 p = 0.485 (sign differs from HOG-based), H5 p = 0.581 (agrees).
Metazoa is the only stem with an extreme percentile for both classes (0.976 each).

## B. Duplication re-placement (POST HOC sensitivity; motivated by positive-control C07 diagnosis)
Each duplication is moved to the most recent common ancestor of the species carrying its daughter copies (removing the
"one node too deep plus compensating losses" artefact). Validation on Ancestra simulations: correlation of per-branch duplications
with truth increased in 6 of 8 simulations (e.g. 0.53 → 0.58; 0.926 → 0.947) and decreased by ≤ 0.001 in 2; error decreased in
8 of 8. Caveat: only 2–6% of duplications move in simulations versus 38% of control-branch duplications in the real data, so the
simulations support the direction of the correction, not its magnitude.

Results: H4 p = 0.527 (pooled +0.08); H5 p = 0.482 (pooled −0.11) — conclusions unchanged. Notable shifts: Metazoa-stem TF-DBD
duplications 21 → 30 (percentile 0.64 → 0.98); Rhodophyta H5 percentile 0.071 → 0.690. The artefact mostly affects controls
(control branches keep a median 72% of their duplicated HOGs after re-placement, stems 83–107%). Other variants (matched sets fixed;
secondary set re-matched) also non-significant (H4 p = 0.483 / 0.843; H5 p = 0.372 / 0.253).
