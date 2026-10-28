# Round 1 — Statistics review (phylogenetic comparative methods)

**Main point.** All inference rests on about 20 eukaryotic transition branches. Of these:
- about 9 are single-species terminal branches (Acrasis, Fonticula, Copromyxa, Guttulinopsis, Capsaspora, Neolecta, Myxogastria, plasmodiophorid, choanoflagellate colonies);
- some are nested (Agaricomycetes and Pezizomycotina inside hyphal fungi);
- 4 have no unicellular sister in `species_final.tsv` (oomycete, siphonous, Xanthophyceae, plasmodiophorid).

The GLM as specified treats about 390 branches as independent and compares against a mixed control pool. Three tests are confounded by design:
- **H8:** aggregative status is aliased with being a terminal branch.
- **H9:** transcriptomes make up 33% of unicellular sisters but only 7% of multicellular species.
- **Genome-wide RER:** the test cannot work after RERconverge normalisation.

**Power.** Simulated NB-GLM, 370 controls, α = 0.05, independence assumed (optimistic); NB2, log BL ~ N(−2.5, 1), log μ = 4 + 0.7·log BL, 300 replicates, seed 11.

| transition branches | NB size θ | power at RR = 1.5 / 2 / 3 |
|---|---|---|
| 10 | 1 | 0.20 / 0.51 / 0.88 |
| 15 | 1 | 0.29 / 0.68 / 0.96 |
| 20 | 2 | 0.64 / 0.98 / 1.00 |

Only effects of about twofold or more are detectable. The sign test needs 9/10, 10/12 or 12/15 concordant contrasts. If each contrast has a 0.7 chance of pointing the expected way, its power is only 0.15–0.30. The Wilcoxon signed-rank test reaches 80% power only for standardised effects d ≥ 0.8–1.0.

---

## MUST

**M1. Three branch classes and a frozen eligibility table.**
- Replace `is_transition` with `state ∈ {transition_stem, within_multicellular, within_unicellular}`.
- Primary contrast: stem vs within_unicellular. Branches inside Metazoa or Streptophyta are not unicellular controls. Report within_multicellular vs within_unicellular separately: this is the "state" effect H9 is about.
- Nested transitions (simple→complex fungi) go in a separately labelled set, never in the same model as their parent transition.
- Before any results, freeze `transitions.tsv` with columns `eligible_internal`, `eligible_terminal`, `has_unicellular_sister`, `nested_in` and `type`, and print the n per hypothesis into DESIGN.md.

**M2. Fit terminal and internal branches separately.**
- "Supported gains" (≥2 descendant species) is undefined on terminal branches.
- **Primary (H1/H2):** internal branches only:
  `glm.nb(supported_gains ~ state + log(BL) + log(n_desc) + log(n_sister_desc) + min_BUSCO_sister + frac_tx_desc)`
- **Secondary:** terminal branches, with unicellular tips as the only controls:
  `gains ~ state + log(BL) + BUSCO_self + is_tx_self + log(n_proteins_self)`
- False stem gains come from sparse or incomplete **sister** sampling (Richter et al. 2018), hence the sister covariates.
- Keep ≤5 covariates. The effect is estimated from about 10–12 branches.

**M3. Exposure.**
- **Gains:** keep log(BL) as an estimated covariate, not an offset. Concatenated-tree BL measures substitutions/site (rate × time). H9 predicts that multicellularity raises the rate, so BL is partly a mediator. Floor BL at 1e‑4.
- **Losses and duplications:** add `offset(log(n_HOGs_at_parent))`, because only existing genes can be lost or duplicated.
- **Sensitivity:** refit (a) without BL and (b) with LSD2 relative-time lengths (`iqtree --date`, 3–5 logged calibrations). Report whether the sign of the effect is stable.
- Assert that the branch-length tree equals the Ancestra guide topology.

**M4. Inference that respects phylogeny.**
- Placement errors are anticorrelated between parent and child branches, and data quality clusters by clade. The Wald p-value therefore cannot be the primary p-value.
- Primary statistic: the fitted rate ratio (RR).
- Primary p-value: **stratified relabelling**.
  - Assign the n_T transition labels to random within_unicellular branches from the same strata as the real transitions.
  - Strata: terminal/internal × BL tertile × n_desc {1, 2–4, ≥5}.
  - Refit for each draw: 10,000 draws, seed 11, p = (1 + #{RR* ≥ RR}) / (N + 1).
- Also report:
  - the Wald 95% CI;
  - leave-one-transition-out and leave-one-supergroup-out RRs;
  - θ;
  - simulated quantile residuals.
- If θ < 0.5 or residuals cluster by supergroup, add a supergroup random intercept (glmmTMB `nbinom2`).

**M5. Main display: per-transition observed/expected (O/E).**
- Fit the GLM to controls only.
- For each transition, report O/E with a 95% **prediction** interval from NB simulation (seed 11).
- Pool with a random-effects summary of log(O/E), I², and the count above the 97.5th percentile.

**M6. Exactly one statistic per primary hypothesis.** Correct the six primary p-values with Holm (a confirmatory family) and also report BH.

- **H1/H2:** the M4 RR, two-sided because the two hypotheses predict opposite directions. The co-option index is linked algebraically (1 − index = toolkit gains / toolkit at M), so it is descriptive, not a second test.
- **H4, H5 and H12** use compositional outcomes, which separates them from the general burst in H2:
  - **H4:** `TF_dups ~ state + … + offset(log(total_dups))`.
  - **H5:** H5-class (gains + dups) with offset log(total gains + dups).
  - **H12:** binomial `cbind(lost, present − lost) ~ state + log(BL)`, using flagellar/phagotrophy HOGs present at P. Include only transitions whose P has such HOGs.
- **H4 complexity interaction** (about 6 vs 10 transitions; minimum detectable ratio of RRs ≈ 3): estimate with CI only.
- **H8:** aggregative ≈ terminal, so compare log O/E of aggregative vs clonal terminal transitions, and Dictyostelia vs internal clonal transitions, by Mann–Whitney. State the n per arm (about 5 vs 3–4) and the minimum detectable effect. Otherwise demote H8 to exploratory.

**M7. H9 contrasts.**
- Measure length and disorder **on HOGs that are single-copy in both clades of a contrast**. Take the median across HOGs of log(trait_M / trait_S).
- Do not use proteome-wide means. Transcriptome fragments shorten sister proteins, which favours H9 artefactually.
- Repeat with genome-only species.
- Contrasts must be phylogenetically independent: each species is used once, and nested pairs are dropped or merged.
- Report:
  - the Wilcoxon signed-rank p-value;
  - the Hodges–Lehmann median with its exact CI;
  - n;
  - the number of concordant signs.
- Treat PGLS with a binary predictor as secondary, because of phylogenetic pseudoreplication (Maddison & FitzJohn 2015). Give the number of origins as the effective n.

**M8. RER.**
- After normalisation, mean RER per branch across genes is ≈ 0 by construction, so test (i) cannot detect a genome-wide rate shift.
- Replace (i) with a **per-transition relative-rate test** on the concatenated tree:
  `r = log(mean path P→tips(M) / mean path P→tips(S))`
  - P is the common ancestor of M and its sister S.
  - Both paths span equal time, so r is a log rate ratio.
  - Summarise across transitions as in M7.
- Per-gene RER (ii):
  - R is in no project env: install `r-base` + RERconverge now.
  - Settings: `transform="sqrt"`, `weighted=TRUE`, `scaleForPproj=TRUE`.
  - Minimum 15 species per gene, with ≥2 transitions represented.
  - Foreground for H9 (a state hypothesis): `foreground2Paths(clade="weighted")` over multicellular clades, so Metazoa does not dominate. Stem-only foreground is a sensitivity analysis.
  - p-values from **permulations** (Saputra et al. 2021; 1,000 per gene), because plain correlations are anti-conservative.
  - BH correction; label as exploratory.

**M9. Convergence null (H5, per-domain).**
- Report the eligible control-pool size per transition. Require ≥5 distinct controls per stratum, or else widen the caliper and log the change. Check the large-clade stems (Metazoa, Fungi, Streptophyta) first: their pools will be small.
- Sample **without replacement**. Controls must not be nested in one another, inside multicellular clades, or on P's own lineage.
- Recompute the statistic on each pseudo-set with the identical pipeline, thresholds included.
- Use a continuous statistic (Σ per-branch z of the compositional proportion), not a count of significant transitions.

**M10. Report everything, including null results.**
- For every primary test, report: estimate, 95% CI, permutation p, Holm p, n, and the minimum detectable effect.
- Smallest effect of interest: RR = 1.5.
- Wording:
  - Write "no burst ≥1.5-fold" only if the upper CI is < 1.5.
  - Otherwise write "inconclusive; compatible with RR ∈ [a, b]".

## SHOULD

- **S1. Outgroup rarefaction.** On internal stems, drop 0, 25 and 50% of sister species and record how stem gains move.
- **S2. Duplicated BUSCOs.** Add BUSCO D% as a covariate in the duplication models, because transcriptome isoforms that survive filtering inflate duplications.
- **S3. Fisher enrichment.**
  - HOGs of one domain family are not independent.
  - Count each domain once per HOG, and report the number of distinct families behind each enrichment.
  - Apply BH within each transition; label as exploratory.
- **S4. Co-option set.** Define the set as frozen-toolkit HOGs present at M. Defining it by enrichment of gains makes the index "new" by construction; defining it by expansion makes it "old". Compare with a beta-binomial `present_at_P ~ state + log(BL_stem)`.
- **S5.** Positive controls are calibration: report the recovery rate, outside every test family.
- **S6. Seeds.** Use seed 11 for all primary results. The other four seeds only show Monte Carlo stability of the permutation p-values (report the range). Write the seed into each output's header.

## COULD

- **C1.** `scripts/power_sim.py`: rerun the power table using the real control BL/n_desc distribution and the fitted θ once Ancestra finishes. Store the project-specific minimum detectable effects.
- **C2.** Phylogenetic logistic regression (Ives & Garland 2010) of presence vs state for top candidate HOGs.
