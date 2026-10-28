# Running the eukaryote primary tests (for the PI)

## Prerequisites
1. `runs/euk_fast/DONE.json` exists (Ancestra `fast`, eukaryotes).
2. `results/gainladder/euk/branch_gain_ladder.tsv` exists (`scripts/gain_ladder.py ... --run euk`). H1/H2 needs it; `primary_tests.py`
   aborts without it (`--allow-missing-ladder` runs H4/H5/H12 only; H1/H2 must then be run later with `--tests H1H2`).
3. `python scripts/prereg/check_prereg.py` prints `PREREG OK` (every analysis script checks this itself and aborts on mismatch).

## The command (single line; ≤ 1 thread, niced; about 45-60 min, mostly the H1/H2 bootstrap)

```bash
cd /scratch/spascare/agentic-ai/workspaces/phase2-multicellularity && export PATH=$PWD/conda_envs/ancestra/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 && nice -n 19 python scripts/analysis/build_branch_events.py --run euk && nice -n 19 python scripts/analysis/primary_tests.py --run euk && nice -n 19 python scripts/analysis/score_positive_controls.py --run euk
```

- Step 1 parses `runs/euk_fast/hogs.orthoxml` once into `runs/analysis/euk/` (skipped if the cache exists; `--force` rebuilds).
- Step 2 writes `results/primary/euk/` and appends every test row to `results/analysis_registry.tsv`. `results/h9/h9_tests.tsv`
  already exists, so Holm over the primary family is computed in the same run (`holm.tsv`; `holm_p` in `primary_tests.tsv`
  and in the registry). If H9 is ever re-run: `python scripts/analysis/holm_primary.py`.
- Step 3 writes `results/positive_controls/euk/{controls.tsv,summary.md}`. It can be re-run at any time. C10 is `undetermined` until
  `runs/gainladder/euk/per_query.parquet` (or `results/gainladder/euk/gains_classified.tsv.gz`) exists.
- Prokaryote dry run (already done; re-run when `results/gainladder/prok/` exists so that the seq-novel column is filled):
  `nice -n 19 python scripts/analysis/primary_tests.py --run prok` (after `score_positive_controls.py --run prok`).

Runtime: the H1/H2 parametric bootstrap costs about 0.17 s per replicate on the eukaryote tree, so about 28 min for the
primary set (9,999 replicates; `--n-boot`) and about 6 min for the secondary set (1,999 replicates; `--n-boot-secondary`).
Everything else takes a few minutes.

## What is computed (DESIGN §4.6 + v3.1, PREREG incl. Dev-1)
| test | s_b | null / p | effect (pooled REML) |
|---|---|---|---|
| H1/H2 | log((y+.5)/(ŷ+.5)), y = sequence-novel gains (ladder column `seq_novel`, summed over path-stem edges), ŷ from the NB-GLM on U (LOO for controls) | **primary: parametric bootstrap** (Dev-1, two-sided); sensitivities: permute, v3 `spec` | log O/E of raw counts vs the matched controls |
| H4 | logit((a+.5)/(d+1)) − logit(π_u), C = TF-DBD HOGs (no C2H2), d ≥ 10, WGD-flagged branches excluded; matching adds log d | permute, one-sided; sensitivities: `spec`, WGD branches included | difference in adjusted logit (log odds ratio) |
| H5 | as H4, C = adhesion/ECM ∪ signalling | as H4 | as H4 |
| H12 | w_i · f_b, 4 flagellar modules, w from `data/h12_expected_direction_frozen.tsv`, not_applicable (Rhodophyta) excluded | permute, one-sided; primary only if ≥ 2 loss-expected stems keep a value (else the row is labelled exploratory and leaves the Holm family) | signed difference |

All tests: k = min(20, |valid U|) Mahalanobis-nearest controls; 10,000 draws; seed 11; per-transition table with the estimate,
5-95 % range of the controls and the percentile u; REML pooled effect (Wald and HKSJ CIs); leave-one-transition-out and
leave-one-supergroup-out (`loo.tsv`; LOO p-values use `permute` for every test, including H1/H2). Secondary nested transitions
(Pezizomycotina, Agaricomycetes, Embryophyta) get the same statistics, labelled `secondary`, outside Holm; H12 has no frozen
direction for them, so only f_b is reported.

## Specification gaps and how they were resolved (fixed before any eukaryote output existed)
1. **Class membership** (H4/H5 and the descriptive compositions): a HOG at any level is in class C if ≥ 50 % of its member genes
   carry a class domain. For H5 the union is taken per protein (a protein carries an adhesion/ECM or a signalling domain), then the
   50 % rule. The LRR rule follows `toolkit_pfam.tsv`: LRR_1/LRR_8 count as adhesion only on proteins without NB-ARC, NACHT or TIR;
   proteins with another adhesion domain stay adhesion. TF-DBD excludes C2H2 (its own category). Pfam source: the pre-registered HMM
   subset (`data/annotations/pfam_prereg/`), which contains every Pfam the tests and controls use; full Pfam is not needed.
2. **WGD exclusion** (H4/H5): an internal branch is flagged if its duplication rate (HOGs at the parent with ≥ 1 non-origin
   duplication on the branch / HOGs at the parent) is ≥ 3 × the median rate over all internal, non-root-child branches, or if
   `results/qc/flags.tsv` names it with a duplication/WGD/polyploidy issue. The current QC flags name only terminal branches (Spirmus,
   Selamoe, Allomac), which are never stems or controls. A flagged edge of a path stem excludes the whole stem. Flags are listed in
   `results/primary/euk/wgd_flags.tsv`; the "WGD included" sensitivity is reported.
3. **H12 "modules with < 50 % of HOGs present at u skipped"**: a module is present at u if ≥ 50 % of its components (h12_sets.tsv
   members; the three BBS2 domain models count as one component) are carried by ≥ 50 % of the members of some module HOG at u.
   Only module HOGs at nodes where the module is present enter f_b. Coverage per node: `h12_module_coverage.tsv`.
4. **Pool per test**: U restricted to branches with a defined statistic (d ≥ 10 and not flagged, for H4/H5; ≥ 1 module present, for
   H12) before matching, so k stays min(20, |valid U|).
5. **LOO-supergroup**: supergroup = majority first token of `lineage_group` of the multicellular species (Holozoa/Holomycota →
   Opisthokonta; Stramenopiles/Alveolata/Rhizaria → SAR; Archaeplastida, Amoebozoa as is).
6. **Bootstrap size**: Dev-1 gives no replicate count; 9,999 (primary) and 1,999 (secondary) are used. The minimum attainable p is
   1/(n_boot + 1).
7. **Positive controls** (`mcell/controls.py`): `origin(H)` is the level of the top of H's ancestor chain; `present_at` follows the
   lineage up or down; "loss on stem(C)" means present at parent(C) and no lineage copy at C. C13 identifies Hex1 by its UniProt
   accession in the proteome ids (`Neurcra_P87252`) instead of a reciprocal best hit, and eIF5A as the other N. crassa PF01287
   protein(s). H_hex is the top ancestor HOG of Hex1. "≥ 90 % Pezizomycotina species" is read with the leak tolerance
   (≤ 2 outside species). C09: a P node is not_evaluable if one side (M clade or sister lineage) has only low-quality species. C15:
   O/E = exp(adjusted-logit difference vs matched controls); controls follow the H5 WGD exclusion. C10(c): any representative of
   the gained family whose best outside DIAMOND hit is prokaryotic (gain-ladder search, e ≤ 1e-5). `undetermined` means an input
   is missing; the scorer can be re-run.

## Observation for the PI (no change made)
The frozen H4/H5 statistic gives s > 0 whenever a = 0 and π_u < 0.5/(d+1), because of the continuity correction. In the prokaryote
run, where TF-DBD domains are almost absent, every stem has s(H4) ≈ 2-5 with a = 0. For eukaryotes π_u(TF-DBD) should be a few
percent, so this only matters for stems with small d.

Tests: `tests/test_analysis.py` (synthetic OrthoXML with planted events; an end-to-end run of the eukaryote driver; a smoke test of all
16 eukaryote control rules on the real guide tree).
