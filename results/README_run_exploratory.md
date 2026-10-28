# Running the exploratory analyses (for the PI)

Everything here is **EXPLORATORY or descriptive**. Nothing is in the primary Holm family (PREREG.md). Every test
row goes to `results/analysis_registry.tsv` with `tier = exploratory` and a note that starts with `label=exploratory`.
Prokaryote rows also carry `DRY RUN prok`.
Code: `mcell/explore.py` (shared layer) and `scripts/analysis/{pfam_wide,cooption,architectures,protein_traits,reversals}.py`.
Tests: `tests/test_exploratory.py` (planted toy data plus hand-computable cases; about 1 s).

## Prerequisites
1. `runs/analysis/<run>/DONE.json` (`scripts/analysis/build_branch_events.py --run <run>`). `architectures.py` and
   `pfam_wide.py --no-hog` do not need it.
2. Full Pfam in `data/annotations/pfam/<run>/<code>.pfam.tsv.gz`. Species without a finished file are listed and
   left out; they are never counted as zero.
   - `architectures.py` refuses to run below 95 % species coverage unless `--allow-partial` is given.
   - `pfam_wide.py` and `cooption.py` drop non-evaluable contrasts. A contrast is evaluable when at least
     min(2, |M|) M species and at least 1 sister are annotated, and at least 80 % of the species outside M are annotated.
3. Disorder in `data/annotations/disorder/<run>/` (only finished species are used). Without it,
   `protein_traits.py` runs the length tests only.
4. Optional: `results/gainladder/<run>/gains_classified.tsv.gz` enables split-repaired presence in the HOG co-option version.
5. `scripts/prereg/check_prereg.py` must pass. Every script checks this and aborts on a mismatch.

## Commands (single line each; ≤ 1 thread, niced; eukaryotes after the Ancestra run and Pfam/disorder are complete)

```bash
cd /scratch/spascare/agentic-ai/workspaces/phase2-multicellularity && export PATH=$PWD/conda_envs/ancestra/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 python scripts/analysis/pfam_wide.py      --run euk
nice -n 19 python scripts/analysis/cooption.py       --run euk
nice -n 19 python scripts/analysis/architectures.py  --run euk
nice -n 19 python scripts/analysis/protein_traits.py --run euk
nice -n 19 python scripts/analysis/reversals.py      --run euk
nice -n 19 python scripts/analysis/exploratory_summary.py --run euk      # -> results/exploratory/euk/summary.md
```
All five in one go:
`for s in pfam_wide cooption architectures protein_traits reversals; do nice -n 19 python scripts/analysis/$s.py --run euk || break; done && python scripts/analysis/exploratory_summary.py --run euk`

Prokaryote dry run (already done; the results are in `results/exploratory/prok/`):
`... --run prok` for each script, with `reversals.py --run prok --pairs Synechococcales:T_cyano_filament,Chroococcales:T_cyano_filament,Acaryochloridales:T_cyano_filament`.
The prokaryote tree has no reversal rows, so this runs pseudo-reversals: unicellular cyanobacterial clades nested in the filamentous M clade, labelled `pseudo:`.

Common options:
- `--transitions default|primary|secondary|all|T_a,T_b`. Default: the 7 primary origins (euk) or the 4 PREREG
  descriptive transitions (prok).
- `--n-perm` / `--n-draws` (10,000), `--seed` (11), `--no-registry`.
- Each script rewrites its own output folder. Re-running appends new registry rows, so delete stale rows by hand if a run is repeated.

Runtime on prok: each script takes under 1 min once the Pfam tables are cached in `runs/analysis/<run>/pfam_full/`.
The first euk run builds those tables in about 5-15 min (overlap resolution for architectures). Everything else should take minutes.

## Definitions

**Contrasts.**
- Transition: M = tips below M; S = `sister_species` of data/transitions.tsv; P = the first node of the
  transition's own `stem_path` (else P, else parent(M)).
- Outside = every species not below M. Outside P = every species not below P.
- Pseudo-transitions (controls): for each branch b = (u → v) of the control pool, M = tips(v), S = tips(u) − tips(v), P = u.
- The pool is U (DESIGN §4.6; |U| = 30 for euk). If |U| < 10, as for prok (|U| = 2), a relaxed pool is used: internal
  within-unicellular branches, not root children, with ≥ 2 tips and ≥ 2 sister tips, adjacency to stems allowed (13 branches).
  It is labelled `relaxed(...)` everywhere and is not interpretable as a test.
- Matched sets N_i: `branchstats.match_controls` (Mahalanobis, k = min(20, |pool|)).

**1. pfam_wide.py.**
- Copy number = the number of proteins of a species carrying the domain (gathering cutoffs).
- Per transition and domain:
  - log2fc = log2((mean_M + 1)/(mean_S + 1)).
  - One-sided Mann-Whitney (M > S) on per-species counts when each side has ≥ 2 species. It is exact by enumerating
    labellings when C(n_M + n_S, n_M) ≤ 2000, otherwise the normal approximation. BH is applied within the transition.
  - **expanded**: q ≤ 0.05 with log2fc > 0, or log2fc ≥ 1 with presence in ≥ 50 % of M.
  - **gained**: present in ≥ 50 % of M and absent from every annotated species outside M.
  - **lost**: absent from all M and present in ≥ 50 % of S.
  - **contracted**: log2fc ≤ −1 with presence in ≥ 50 % of S.
  - With ≤ 3 species per side the minimum one-sided Mann-Whitney p is 0.05-0.17 (e.g. 1/C(5,2) = 0.1), so for small clades the fold rule decides.
- HOG view: per domain, the HOGs gained, duplicated (non-origin) or lost on the stem edges that carry the domain
  in ≥ 50 % of their members.
  - `hog_expanded`: gain + dup HOGs ≥ 2 and percentile u ≥ 0.95 among the matched controls.
  - Controls are compared with their own matched sets inside the pool.
- Convergence:
  - k = the number of transitions in which a domain is flagged; G_k = the number of domains with k ≥ 2 (or ≥ 3).
  - Null `permute` (primary, the analogue of D6): transitions and pseudo-transitions are exchangeable; each draw takes n
    distinct columns and recomputes k per domain and G_k. 10,000 draws, seed 11. BH over domains with k ≥ 1.
  - The `matched` null (one matched control per transition, distinct) is reported as a sensitivity only. Per domain it
    is anti-conservative for domains rarely flagged in controls (toy test: 17 % of null domains at p ≤ 0.01).
- `toolkit_by_transition.tsv` is the orthology-free Pfam-count view of the toolkit categories. DESIGN §4.6 requires it to
  agree with the H4/H5 expansions before those are called robust.

**2. cooption.py.** DESIGN §4.6.
- F_i = families with log2((mean_M + 1)/(mean_S + 1)) ≥ 1 that are present in ≥ 50 % of M.
- C_i = the fraction of F_i present at P, i.e. carried by ≥ 1 genome or ≥ 2 transcriptome species outside M.
- C_i^parent = the same fraction among species outside P.
- Reported: C_i − median C over N_i, the percentile u, and T = mean u. There is no p-value: this is descriptive by design.
- Versions:
  - `pfam_expanded`: all domains.
  - `toolkit_expanded`: toolkit domains only.
  - `toolkit_present`: toolkit domains present in ≥ 50 % of M.
  - `hog_toolkit`: Ancestra families with ≥ 50 % toolkit genes. With the gain ladder, a split-repaired gain on a branch
    below P counts as present at P (an upper bound).

**3. architectures.py (H11).**
- Architecture = Pfam accessions ordered by alignment start, after greedy overlap resolution by score, with tandem
  repeats collapsed (the same rule as `mcell/gainladder.py`).
- An architecture is new on branch p → c if c is the MRCA of the annotated species that carry it and ≥ 2 species carry it.
- Statistics per branch, summed over path-stem edges:
  - `new_arch_multi` (≥ 2 distinct domains): the H11 statistic.
  - `new_arch_recombined`: multi-domain, with every component domain also occurring outside the clade.
  - `new_arch_all`.
- Test: NB-GLM residual (`count_residual_statistic`) and raw counts, relabelling scheme `permute`, one-sided greater,
  10,000 draws, seed 11; effect = pooled log O/E. Exploratory.

**4. protein_traits.py.**
- HOGs are taken at P. A HOG is used when it is present in ≥ 50 % of M and of S species and is single-copy in every
  M ∪ S species that has it.
- Per HOG, d = log2(median length M / median length S); for disorder, d = the difference in median frac_disordered.
- Test: Wilcoxon signed-rank (two-sided) with paired Hodges-Lehmann and a 95 % CI.
- TF-DBD HOGs (H10) are tested single-copy (`tf_disorder_sc`) and at any copy number (`tf_disorder_all`, `tf_length_all`).
- Secondary tests compare per-species medians (whole proteome, all species and genome-only) with Mann-Whitney and the
  two-sample HL. `tf_disorder_excess` = TF disorder minus proteome disorder per species.
- For the paired length and disorder tests, u = the percentile of the transition's HL among its matched pseudo-transitions.
- BH is applied across all rows of the script.

**5. reversals.py.**
- Reversal branch = rev_parent → M of every `to_state = reversal` row.
- Preceding transition = the nearest ancestral origin (unicellular → simple/complex). Mapping: all five yeast reversals →
  T_fungi_hyphae (path Eumycota > Blastocladiomycota_Amastigomycota > Amastigomycota); Zygnematophyceae → T_streptophyta.
- HOGs at rev_parent are classified by whether an ancestor (or the HOG itself) was gained, or was the parent HOG of a
  non-origin duplication, on a stem edge. Lost = no descendant at rev_child.
- Test: Fisher exact (one-sided greater and two-sided), with Haldane/Woolf OR CIs.
- Sister-branch controls: ratio of ORs relative to the sisters that kept the multicellular state (majority of tips with
  `mc_state` simple/complex).
- Reversals of a single species whose data are transcriptome-only are flagged: their losses are uncertain.

## Choices where the specification is silent
1. Copy number counts proteins per domain, not domain hits (`n_hits` is kept in the cache).
2. Expansion uses an OR of a test and a fold rule (as in the brief), because with ≤ 3 species per side the Mann-Whitney
   cannot reach q ≤ 0.05. `expanded_by_test` and `expanded_by_fold` are reported separately.
3. Convergence uses the `permute` null, not draws from matched controls only (see the calibration note under 1).
4. **Path stems.** `mcell.branchstats.build_branch_table` keeps the first `stem_path` among transitions that share an M.
   `T_volvocales_gv` (exploratory) is listed before the primary `T_volvocine_tgv`, so the branch table gives
   Gonium_Volvocaceae a 1-edge stem, and the TGV_clade edge is neither summed nor in U. The exploratory scripts take each
   transition's own `stem_path` and patch a local copy of the branch table (`explore.patch_path_stems`): branch length
   0.0070 → 0.0178.
   **This affects the primary pipeline (H1/H2, H4, H5 inputs for T_volvocine_tgv). Reported to the PI; not changed here.**
5. HOG co-option presence follows the same species rule as Pfam (≥ 1 genome or ≥ 2 transcriptomes outside M), plus split repair when available.
6. Reversal sister reference: sisters with a majority multicellular tip state, or all sisters if there are none.
