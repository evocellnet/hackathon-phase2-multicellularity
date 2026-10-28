# mcell.branchstats — branch statistics (DESIGN §4.6 + v3.1 amendments D6/D7)

Code: `mcell/branchstats.py`. Tests: `tests/test_branchstats.py` (38 tests, about 90 s).
Scripts:
- `scripts/validate_branchstats_null.py` → `results/validation_branchstats/null_calibration.tsv`
- `scripts/prereg_power_table.py` → `results/prereg/power_table.tsv` (+ `power_curves.tsv`)

Every function works on per-branch values. A branch id is the name of its
child node. Event tables come from `mcell.hogs` (`events`, `hogs`).

## API

| step | function | spec element |
|---|---|---|
| tree | `load_tree(newick)` | rooted Newick with named, unique internal nodes |
| branch table | `branch_table(tree, traits, branch_lengths=None)` | depth, n_desc, n_sister, is_terminal, root_child, branch_length/log_bl, frac_tx_desc/sister, mean_busco_desc. `branch_lengths` = None / `"tree"` / dict / a tree with lengths (path, Newick or ete3) |
| lengths | `branch_lengths_from_tree(guide, bl_tree)` | maps lengths by tip set (or its complement). Before `data/species_tree_bl_euk.nwk` exists, the placeholder 1.0 is used (with a warning) |
| classes + U | `build_branch_table(tree, traits, transitions, branch_lengths)` | branch_class, transition_id, role, primary / primary_transitions, `stem_path`, adjacent_to_stem, `in_U`, `U_exclusion` |
| stems | `stem_branches(bt)`, `path_edges(bt, b)` | transition_id → stem branch (= M). Primary = `analysis_set == 'primary'`, falling back to `primary_eligible == 1` |
| matching | `match_controls(bt, stems, pool, covariates, extra, k=20)` → `MatchResult` | z-scored Mahalanobis, k = min(20, \|U\|); `.smd`, `.flagged()` |
| relabelling | `relabel_test(values, match, alternative, n_draws=10000, seed=11, scheme='permute', signs)` | u, T_obs, p, p_mcse, `min_attainable_p` |
| floors | `min_attainable_p`, `count_relabellings` | 1/(n_draws+1) and the exact enumeration floor |
| effects | `effect_table`, `random_effects` (REML, DL fallback, z or HKSJ CI), `holm`, `leave_one_out` | §4.6 effect size, pooling, multiplicity, LOTO/LOSO |
| H1/H2 | `count_statistic(events, 'gain', hogs=seq_novel, branches=bt)` → `count_residual_statistic(y, bt)`; `fit_count_glm`; `residual_bootstrap_test` | D6 residual statistic |
| H4/H5 | `compositional_statistic(events, hogs, C, bt, min_d=10, exclude)` | s_b = logit((a+.5)/(d+1)) − logit(π_u) |
| H12 | `module_loss_statistic(events, module_hogs, ref_size, bt, hogs=hogs)`, `h12_signs` | f_b, with w_i passed through `signs=` |
| H9 | `h9_values`, `h9_control_branches` | plug-in r_i / r_c with random orientation |
| power | `simulate_power`, `mde` | NB planted multiplicative effect |
| registry | `append_registry(path, row)` | `results/analysis_registry.tsv` |

H1/H2 in v3.1:
```python
bl = bs.branch_lengths_from_tree(tree, "data/species_tree_bl_euk.nwk")   # placeholder 1.0 until it exists
bt = bs.build_branch_table(tree, traits, transitions, branch_lengths=bl)
m  = bs.match_controls(bt)
y  = bs.count_statistic(hr.events, "gain", hogs=seq_novel, branches=bt)  # path stems summed
s  = bs.count_residual_statistic(y, bt)["s"]
res = bs.relabel_test(s, m, alternative="two-sided")                     # scheme='permute'
```

## v3.1 amendments as implemented

**D6 relabelling scheme.** `scheme='permute'` is the default and the primary
scheme. It permutes the statistic values over the stems and every control in
any N_i, then recomputes u and T with the original matched sets. The literal v3
scheme (`'spec'`: c_i ∈ N_i, distinct across transitions, T_d = mean
u(c_i, N_i \ c_i)) remains as a sensitivity option. 'spec' is anti-conservative
because the N_i overlap. Draws without replacement make T_d less variable than
T_obs: the draws are negatively correlated with each other, while the stems'
u values are positively correlated through shared controls. See the
calibration section below.

**D6 H1/H2 residual statistic.** s_b = log((y_b+0.5)/(ŷ_b+0.5)), with ŷ from
an NB2 GLM fitted on U only. The model is y ~ log_bl + log_n_desc +
log_n_sister + frac_tx_desc; covariates are standardised on U.
- Each control gets a leave-one-out fit. Stems and all other branches use the
  full-U fit.
- The NB fitter is our own: IRLS alternating with an ML estimate of α. It agrees
  with `statsmodels` NegativeBinomial to within 1e-4 (tested).
- If the fit fails (no convergence, non-finite values, α at its bound), it falls
  back to Poisson with quasi-dispersion, raises a warning and records the
  fallback in `.attrs['fits']`.
- Covariates that are constant or missing in U are dropped. This includes the
  1.0 branch-length placeholder.

**Caveat.** All stems share the prediction error of the GLM. When stems lie
outside U's covariate range (extrapolation), their residuals move together and
the permutation test becomes anti-conservative. Rejection rates rise from 0.045
with no shift to 0.07 at a 1 SD shift and 0.245 at a 2.5 SD shift (α = 0.05;
see the table below). `residual_bootstrap_test` is the calibrated alternative:
it simulates from the fit on U, refits with LOO in every replicate and
recomputes T. It needs about 0.2 s per replicate. Recommendation: report it next
to the permutation p for H1/H2 once the branch lengths are known, and use it as
primary if stems extrapolate. The PI should decide this.

**D7 path stems.** `stem_path` "P>X>M" must be a parent→child chain ending at M.
An invalid path raises a warning and falls back to parent(M)→M.
- The stem row (M) gets `stem_path` "X;M", and its `branch_length` is the sum
  over the edges (`edge_branch_length` keeps the single edge). Its other
  covariates are those of the last edge.
- The intermediate edges get the stem's class with role `path_edge`.
- Excluded from U: the edge into P, every child edge of P, and every child edge
  of each intermediate node.
- Counts are summed over the path edges (`count_statistic(..., branches=bt)`, `path_sum`).
- H4/H5: D is the union of the duplicated HOGs over the path edges, each mapped
  to its ancestor at P through `hogs.parent_hog_id` (HOGs born on the path keep
  their own id). π_u is taken at P.
- H12: module HOGs present at P that are lost on any path edge (needs `hogs=`).
- On the current tables, T_fungi_hyphae and T_volvocine_tgv are 2-edge paths.
  `T_volvocales_gv` shares its M with `T_volvocine_tgv` (warning; only the
  primary one is used).

**Prokaryotes.** |U| = 2 on the prokaryote tree, so prokaryotic transitions are
descriptive only (no relabelling test).

## Null calibration (`results/validation_branchstats/null_calibration.tsv`)

The real eukaryote tree has 7 primary stems, |U| = 30, k = 20, and 28 distinct
controls. Values are iid under H0. Rejection rate at α = 0.01 / 0.05 / 0.10
(1,000 replicates, 2,000 draws, seed 11):

| scheme | two-sided | one-sided (greater) |
|---|---|---|
| spec (v3) | 0.045 / 0.138 / 0.220 | 0.036 / 0.113 / 0.166 |
| permute (v3.1) | 0.013 / 0.051 / 0.099 | 0.012 / 0.054 / 0.114 |

The earlier figures of 14.5% / 4.3% came from 400 replicates. The table has
1,000 replicates and gives the same picture.

For the H1/H2 residual statistic, two-sided, rejection rates at α = 0.01 /
0.05 / 0.10. Counts are NB with α = 0.3. Each value is followed by
(replicates; MC SE at α = 0.05).

| setting | permute on raw counts | permute on residual | parametric bootstrap on residual |
|---|---|---|---|
| real tree, constant mean 20 (bl placeholder) | — | 0.013 / 0.073 / 0.125 (400; 0.011) | 0.010 / 0.043 / 0.097 (300; 0.013) |
| synthetic, stems shifted 0 SD | 0.003 / 0.030 / 0.080 (400) | 0.008 / 0.048 / 0.090 (400) | 0.020 / 0.053 / 0.087 (300) |
| synthetic, 1 SD | 0.115 / 0.295 / 0.415 (400) | 0.025 / 0.083 / 0.158 (400) | 0.00 / 0.040 / 0.080 (50; 0.031) |
| synthetic, 2.5 SD | 0.815 / 0.945 / 0.970 (400) | 0.128 / 0.283 / 0.355 (400) | 0.00 / 0.040 / 0.100 (50; 0.031) |

The synthetic setting has |U| = 30, n = 7 and k = 20. Stems are shifted on
log bl and log n_desc, and counts depend on both.
- The residual statistic removes most of the bias that raw counts suffer.
- Permuting residuals is still anti-conservative once stems extrapolate.
- The parametric bootstrap stays calibrated.

Two bootstrap rows had to be cut down because the node was oversubscribed:
they use 50 replicates and 49 resamples (MC SE 0.031 at α = 0.05). The other
bootstrap rows use 99 resamples, so their α = 0.01 entries carry no
information.

Other schemes tried and rejected: a stem↔control swap, which was
anti-conservative one-sided (0.107 at α = 0.05).

## Power table (`results/prereg/power_table.tsv`)

Setup: the real eukaryote tree, the 7 primary stems, the permute scheme,
200 simulations per effect and 1,000 draws. Minimum attainable p =
1/10001 for every test. Assumptions, stated per row:

- **H1/H2**: NB counts with mean 20 on every branch; dispersion 0.3 or 1.0.
- **H4/H5**: d = 10 + NB(30, α 0.5); a ~ Bin(d, p_b) with logit p_b = logit π_u
  + N(0, 0.5²); π_u = 0.04 (H4) and 0.12 (H5). The effect is an odds ratio.
  Matching includes log d.
- **H9**: r ~ N(0, σ²) with σ = 0.15 or 0.30. The pool is U's nodes with two
  child clades of ≥2 tips: 10 nodes, so k = 10.
- **H12**: 40 module HOGs; loss probability ~ Beta(0.5, 12). The effect applies
  only to the 2 loss-expected stems (frozen directions); T_rhodophyta is
  not_applicable, so n = 6.

| test | assumption | MDE α=0.01 | MDE α=0.05 | power at 1.5× (α=0.05) | size at α=0.05 |
|---|---|---|---|---|---|
| H1/H2 | NB α = 0.3 | 3.06× | 2.42× | 0.35 | 0.085 |
| H1/H2 | NB α = 1.0 | > 8× | 5.89× | 0.13 | 0.085 |
| H4 | logit SD 0.5 | OR 3.59 | OR 2.77 | 0.28 | 0.05 |
| H5 | logit SD 0.5 | OR 2.98 | OR 2.20 | 0.35 | 0.05 |
| H9 | σ = 0.15 | 1.39× | 1.26× | 1.00 | 0.07 |
| H9 | σ = 0.30 | 1.89× | 1.69× | 0.66 | 0.055 |
| H12 | Beta(0.5, 12) | not reached | not reached (power ≤ 0.33 even at complete loss) | — | 0.055 |

Reading the table:
- No primary test has 80% power at the 1.5-fold smallest effect of interest,
  except H9 when σ ≤ 0.15. Non-significant results therefore have to be worded
  "not detected at this power".
- H12 with only 2 loss-expected stems out of 6 cannot reach 80% power under
  T = mean u. This argues for keeping H12 exploratory.
- The H1/H2 size of 0.085 at α = 0.05 is the mild extrapolation inflation
  described above. The stems differ from U on log n_sister (SMD −1.1).

## Choices where the specification is silent

- **Node states**: the deepest crown wins; nodes below no crown are unicellular.
  An uncoded multicellular tip or a grey tip is flagged in `state_note`. If M
  lists several tips, each tip is its own crown.
- **Adjacency**: covers all coded transitions of every tier, including grey and
  reversal. It covers the P named in the table as well as parent(M).
- **n_sister**: tips summed over all other children (polytomies).
- **SMD**: (x_stem − mean x_N) / SD_U.
- **Within-transition variance**: var_c(e_c)·(1+1/k). A zero variance is
  floored at the smallest positive one.
- **REML**: bounded optimiser with a check at the τ² = 0 boundary. It matches
  metafor on the BCG data. The CI is Wald z by default; `ci_method='hksj'` is
  available and is advisable with n = 7.
- **π_u of 0 or 1**: smoothed to (n_C+0.5)/(n_u+1) and flagged.
- **H12 module size**: the denominator of the 50% rule is supplied by the
  caller (`module_ref_size`).
- **min_attainable_p**: the larger of 1/(n_draws+1) and 1/(M+1), where M is
  the number of relabellings. For a single transition with k = 20 it is 1/21,
  so always report it next to p.
