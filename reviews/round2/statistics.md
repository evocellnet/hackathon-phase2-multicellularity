# Round 2 — Statistics review of DESIGN v2 (§2, §4, §6)

**Verdict.** The structure is sound but not yet implementable: the control pool is too small, the relabelling statistic is undefined, two tests cannot reach Holm significance, and the H12 directions contradict the sampled taxa.

## MUST

**M1. "≥200 matched controls" is impossible.** There are about 110 unicellular tips, which gives at most about 90 within-unicellular internal branches; about 40–60 remain after eligibility. Define the pool and the matched sets as follows:
```
U = internal b=(u→v): u,v unicellular; |tips(v)|≥2; |tips(sib(v))|≥2; v not root child;
    b not adjacent to any stem (edge into P, edge P→S); not inside multicellular/grey clades
x = z[depth, log BL, log n_desc, log n_sister, frac_tx_desc, frac_tx_sister]
N_i = k=min(20,|U|) Mahalanobis-nearest b∈U to stem i; report |SMD|, flag >0.5
```
The sister covariates are new. False stem gains come from sparse or transcriptomic sister sampling, and v2 matches only on descendants.

**M2. Relabelling statistic.** One rule for every relabelling test, giving each transition equal weight:
```
u(b,N) = (#{c∈N\b: s_c<s_b} + 0.5·#{c∈N\b: s_c=s_b} + 0.5)/(|N\b|+1)
T_obs  = mean_i u(stem_i, N_i)
repeat 10,000 (seed 11): draw c_i∈N_i without replacement across i, reject nested draws;
                         T_d = mean_i u(c_i, N_i)          # leave-self-out
p = (1+#{|T_d−½|≥|T_obs−½|})/10001                         # one-sided where stated
```
- **Why it is coherent.** It is a stratified permutation of the per-transition percentiles (Edgington). The reported per-transition percentile is the same u.
- **No double counting**, provided that no stem, nested stem or within-multicellular branch enters any N_i.
- **Effect size.** log O/E_i = log((y_i+0.5)/(mean_{N_i}y+0.5)). Pool with REML random effects, using the within-transition variance taken from the controls. The "prediction interval" is the 5–95% range of the controls.
- **Demoted.** Stouffer and the sign test become descriptive.

**M3. Attainable p.** About 6–9 transitions are primary-eligible. Holm's first threshold is 0.01. With n ≤ 7, exact Wilcoxon and sign tests have p ≥ 0.0156.
- Test H9 by relabelling (M6).
- Print n, the minimum attainable p and the MDE per test in PREREG.

**M4. H4/H5.** For each branch b=(u→v):
- **Counts.** D_b = HOGs with ≥1 duplication on b (HOGs, not events). d = |D_b|; a = |D_b ∩ C|.
- **Class membership.** A HOG is in class C if ≥50% of its members carry a frozen class domain.
- **Statistic.** s_b = logit((a+0.5)/(d+1)) − logit(π_u), where π_u is the share of HOGs present at u that are in C.
- **Classes.** H4: TF-DBD, excluding C2H2. H5: adhesion/ECM ∪ signalling as one union statistic, so Holm stays at 5 tests; per-class results are secondary.
- **Few duplications.** Require d ≥ 10 for stems and controls. List excluded stems, and add log d to the matching variables.
- **WGD.** Exclude WGD-flagged branches, with a sensitivity analysis that includes them. TF retention after WGD would fake H4.
- **Direction.** One-sided (expansion).

**M5. H12 directions are inconsistent.** Flagellate taxa sit inside both "loss expected" crowns:
- the hyphal crown M contains *Allomyces* and *Gonapodya*;
- streptophyte M1 contains *Chara*, the bryophytes, *Selaginella* and *Ceratopteris*.

Stem losses are therefore impossible there. Re-derive the directions as "present at P, absent at M" on the frozen tree; likely only Dictyostelia remains loss-expected. If fewer than 2 loss-expected stems remain, make H12 exploratory.
- **Statistic.** f_b = mean over the 4 modules of the fraction of module HOGs present at u that are lost on b. Skip modules with fewer than 50% of their HOGs present at u.
- **Signed score.** s_b = w_i·f_b, with w = +1 where loss is expected and −1 otherwise. Test one-sided via M2.

**M6. H9 is well defined only once the following are fixed.**
- **Branch lengths.** Per marker: IQ-TREE `-te` (LG+G4) on the species tree pruned to the taxa present. Fit once, so controls reuse the same fits.
- **Per-gene statistic.** r_ig = log(L̄(P→M tips)/L̄(P→S tips)), where L̄ is the equal-split weighted mean path length.
- **Gene inclusion.** P is identifiable only if the gene has ≥1 outgroup, ≥2 M tips whose MRCA is M, ≥1 S tip, and ≥50% site coverage.
- **Per transition.** r_i = median_g r_ig, with a 2,000× gene-bootstrap CI.
- **Null.** M2 with s = r, using control nodes with two child clades of ≥2 tips in random orientation. Two-sided.
- **Caveats.** r covers the whole post-P history, not only the stem. Longer generation times bias r below 0. Rerun on genome-only taxa.

## SHOULD

- **S1. Co-option index.**
  - **Family set.** d ∈ F_i if log2((x̄_M+1)/(x̄_S+1)) ≥ 1 and d is present in ≥50% of M.
  - **Index.** C_i = fraction of F_i present at P.
  - **Presence at P.** For Pfam, ≥1 genome or ≥2 transcriptomes in S ∪ outgroups. For HOGs, split-repaired presence.
  - **Report** C_i − median(C_{N_i}) and u. Descriptive only: anything present in S is present at P by construction.
- **S2. Report per test:**
  - n, k, the estimate with 95% CI, p with its Monte Carlo SE, the Holm p and the minimum attainable p;
  - the simulated MDE (planted multiplicative effect; 80% power at α = 0.01 and 0.05).
  - **Wording:** "no ≥1.5-fold effect" only if the CI lies within (1/1.5, 1.5); otherwise "inconclusive; compatible with [a, b]".
- **S3.** Run the DIAMOND gain ladder on every branch in U, not only stems.

## COULD

- **C1.** Regression-adjusted matching for poorly matched deep stems: E_i = mean_c y_c·exp(β̂ᵀ(x_i−x_c)).

## None needed

- **Holm across 5.** Valid under any dependence.
- **10,000 draws.** Sufficient.
