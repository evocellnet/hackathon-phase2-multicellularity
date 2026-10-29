# Independent verification report

```text
INDEPENDENT VERIFICATION REPORT - primary eukaryote results
=============================================================
Verifier: independent agent (did not produce any of the results). Date: 2026-10-29.
Repo HEAD at verification: ff430fb (2026-10-29 00:57). Env: conda_envs/ancestra, nice 10, <= 6 threads.
All re-derivations use the verifier's own code: results/verification/scripts/ (outputs in results/verification/outputs/).
mcell.branchstats / mcell.analysis / mcell.hogs were NOT used for re-derivation. lxml was used for OrthoXML, ete3 for
Newick parsing only. The one exception is check 4a, which runs the project's OWN pre-fix and post-fix branchstats code.
That run diagnoses the bug; it does not re-derive anything.

OVERALL: 7 checks. 6 PASS, 1 PASS-with-flags (check 7: the registry has duplicates the report must handle).
No numerical discrepancy was found in any verified quantity.

---------------------------------------------------------------------------------------------------------------
CHECK 1  HOG event counts from OrthoXML                                                              PASS
---------------------------------------------------------------------------------------------------------------
Script: scripts/events_from_orthoxml.py, a one-pass lxml iterparse of runs/euk_fast/hogs.orthoxml.
  - Parsed 74,920 top-level groups and 2,969,261 genes. The taxonomy block was used as the species tree.
  - Semantics (OrthoXML/pyHAM):
      orthologGroup = speciation at its taxonId.
      paralogGroup  = duplication on the branch directly above its taxonId.
      ancestral gene at X = a gene-tree edge u->v with v at or below X and u strictly above X.
  - Counted per branch: gain = family whose top position lies on the branch; loss = ancestral gene at the parent with
    no descendant at the child; dup = k-1 for k descendants; dup_hogs = ancestral genes with k >= 2.

branch                          gains  losses  dups  dups@origin  dup_hogs  n_parent  n_child  | parquet + branch_values
Choanozoa>Metazoa (stem)          760   11593  5677      42         1338     35878    30764   | identical
Eumycetozoa>Dictyostelia (stem)   598    6726  1510      18          694     15542    10942   | identical
Alveolata>Myzozoa (control U)     186   10152  3618      44         1248     21129    14825   | identical
Amoebozoa>Evosea (control U)      261    7840  2913     126         1051     25705    21165   | identical
Volvocine path edges (bonus):
Chlamydomonas_TGV>TGV_clade        78    2886  2050      67         1017     17731    17040   | identical (events.parquet)
TGV_clade>Gonium_Volvocaceae       87    2173  2809      15         1583     17040    17778   | identical (events.parquet)
  - Sets of gained family IDs are identical to events.parquet for all 6 branches.
  - Conservation law holds, e.g. Metazoa: 35878 - 11593 + 5677 + 42 + 760 = 30764.
  - Path-stem summation in branch_values.tsv (row Gonium_Volvocaceae) is correct: gains 165 = 78+87, losses 5059,
    dups 4859, dups@origin 82, dup_hogs 2600 = 1017+1583.
  - H4_d = 2558 (< 2600) is consistent with taking the union of HOGs over both edges. The overlap was not
    independently checked.

---------------------------------------------------------------------------------------------------------------
CHECK 2  Sequence-novel gains (clade rule applied to raw DIAMOND hits)                               PASS
---------------------------------------------------------------------------------------------------------------
Script: scripts/seqnovel.py.
  - Inputs:
      gained families = the verifier's own OrthoXML parse (check 1);
      representatives = runs/gainladder/euk/reps.parquet;
      hits = raw runs/gainladder/euk/hits_pass1.tsv.gz + hits_pass2.tsv.gz.
  - Rule: a gain is sequence-novel if no hit of any searched representative lies in a species outside tips(child),
    including prokaryotes. Species = ID prefix; all 269 DB prefixes are valid species codes.
  - Every gained family had at least 1 searched representative. No representative lay outside its clade.

stem / branch                    HOG gains  seq-novel (mine)  reported (per_transition_H1H2 / ladder)
Metazoa                          760        594               594
Dictyostelia                     598        507               507
Myzozoa (control)                186        152               152 (branch_values)
Volvocine path (2 edges)         78+87      32+56 = 88        88

---------------------------------------------------------------------------------------------------------------
CHECK 3  Holm, T_obs, H12 direction weights                                                          PASS
---------------------------------------------------------------------------------------------------------------
Holm recomputed from the five p-values (m = 5):
  H12  p=0.020498 -> 0.10249   H1H2 p=0.113 -> 0.452   H9 p=0.26937 -> 0.80811
  H4   p=0.361064 -> 0.80811   H5   p=0.772223 -> 0.80811     (all identical to holm.tsv)
  No test is significant after Holm. This holds even at the most favourable Monte-Carlo p-value.

T_obs recomputed with the verifier's own u-score from branch_values.tsv + matched_sets.tsv:
  H12 0.738095 (reported 0.738095)
    u: Dicty .952, Fungi .929, Metazoa .905, PX .595, Strepto .548, Volvox .500.
  H4  0.547619 (=)    H5 0.411565 (=)
  H1H2 0.710884 (= mean of per-transition u)    H9 0.655844 (= from r_i and r_oriented)

Permutation p from the verifier's own implementation (independent RNG, 10,000 draws, scheme = permute):
  H12 0.0227 +- 0.0015 (reported 0.0205)   H4 0.357 (0.361)   H5 0.774 (0.772)   H9 0.277 (0.269)
  All agree within about 1.5 Monte-Carlo SE.

H12 weights: identical to data/h12_expected_direction_frozen.tsv.
  w = +1: T_dictyostelia, T_fungi_hyphae (loss expected).
  w = -1: T_metazoa, T_streptophyta, T_px_filaments, T_volvocine_tgv.
  T_rhodophyta = not_applicable and correctly excluded (n = 6).

---------------------------------------------------------------------------------------------------------------
CHECK 4  H9 and the stem-path bug                                                                    PASS (H9 unaffected)
---------------------------------------------------------------------------------------------------------------
Timeline:
  - H9 ran at 15:19:33 (h9_run.json; registry commit 9ccbf1c; results committed 507b4eb 15:21).
  - The fix (5a646f2) came at 18:31. So H9 was computed BEFORE the fix.

(a) The project's own build_branch_table, run at the pre-fix (507b4eb) and post-fix (HEAD) versions:
  - Exactly one row differs: Gonium_Volvocaceae.
      pre-fix:  stem_path "Gonium_Volvocaceae", branch_length 0.007028 (single edge);
      post-fix: stem_path "TGV_clade;Gonium_Volvocaceae", branch_length 0.017848.
  - The transition_id order of that row also changes.
  - in_U, U_exclusion, adjacent_to_stem and branch_class are identical, so |U| = 30 in both versions.
  - The H9 control pool is identical (the same 10 nodes: Diaphoretickes ... Holomycota).
  - k = |pool| = 10, so every stem's matched set equals the whole pool, and the branch-length covariate cannot matter.
  => The control pool does NOT change after the fix.

(b) H9 does not use branchstats' stem_path for P/M/S:
  - It takes P, M and S from data/transitions.tsv via rates.stem_triples.
  - The gene table records P=Chlamydomonas_TGV, M=Gonium_Volvocaceae, S=Chlamydomonas for T_volvocine_tgv.

(c) Independent recomputation (scripts/h9_volvox.py; own equal-split path code on runs/phylo/euk/genetrees/*.bl.nwk):
  - L_M includes both path edges. Tetrabaena belongs to neither side. Rule: >= 2 M tips whose MRCA is M, >= 1 S tip,
    >= 1 outgroup tip.
  - T_volvocine_tgv: 177 genes, r_i = 0.16297, identical to the reported values. The gene set is identical and the
    per-gene max |diff| is 5e-6. Excluded: 53 genes where MRCA != M (Gonium absent), 20 with < 2 M tips.
  - Control check, T_metazoa: 250 genes, r_i = -0.12254 (identical).
  - Difference attributable to the bug: 0. The bug never reached H9.
  - For contrast, the buggy single-edge coding (P=TGV_clade, S=Tetrabaena) would give r = -0.108 (64 genes). H9 did
    not use it.
  - Bootstrap CI depends on the RNG: the verifier's seed gives [0.055, 0.233]; the report gives [0.084, 0.233].
  - Sensitivity, not binding: with a looser rule that only requires the M tips to be monophyletic (230 genes),
    r = 0.152.

---------------------------------------------------------------------------------------------------------------
CHECK 5  Positive controls C01, C05, C10 (Pfam tables + OrthoXML + raw DIAMOND)                       PASS
---------------------------------------------------------------------------------------------------------------
Script: scripts/controls_mine.py.
  - Pfam data: data/annotations/pfam/euk/*.pfam.tsv.gz (187 species, GA cutoff).
  - PF-HOG = >= 2 members, >= 50% with a hit. Counted per (level, ancestral gene) row.

C01 Wnt (PF00110), gain on the Metazoa stem:
  - 69 PF-HOG rows in 4 families (HOG0003847/3900/3915/4103) with origin on Choanozoa>Metazoa.
  - 0 species outside Metazoa have a hit. All 4 families have 0 outside-clade DIAMOND hits (4/4 sequence-novel).
  - Identical to the report (69; 4 families; 0; 4/4). RECOVERED.

C10 GRAS (PF03514), gain on the Anydrophyta stem:
  - 57 PF-HOG rows in 1 family (HOG0038019). 0 outside euk species have a hit.
  - Best outside-clade hit of the representatives is Soracel (Sorangium cellulosum, a bacterium), bitscore 167,
    e = 4e-43, i.e. prokaryotic.
  - Identical to the report (57; 0; best hits 'prok'). RECOVERED.

C05 SH2+kinase:
  - (a) SH2 (PF00017) before PK_Tyr_Ser-Thr (PF07714) by envelope start: Metazoa 212, Choanoflagellata 97,
    Filasterea 17 (identical).
  - (b) PF-HOG(PF00102) present at Choanozoa with Metazoa + non-metazoan Filozoa members: 140 rows in 10 families
    (identical). This counts HOG rows at Choanozoa or any ancestral level whose lineage is present at Choanozoa, as
    in mcell.controls.
  - Restricted to rows AT the Choanozoa level only: 41 rows in 8 families.
  - Either reading satisfies '>= 1'. RECOVERED.

Minor flag: results/positive_controls/euk/failure_diagnosis.md section 'Summary' still says "pass count stays at
10/16 (C10 is still undetermined)". summary.md and controls.tsv give the final 11/16. Update the text before
reporting.

---------------------------------------------------------------------------------------------------------------
CHECK 6  PREREG integrity and chronology                                                             PASS
---------------------------------------------------------------------------------------------------------------
Hashes:
  - scripts/prereg/check_prereg.py -> "PREREG OK" (exit 0).
  - Recomputing sha256 independently: all 15 hashed files match PREREG_HASHES.tsv.
  - PREREG_HASHES.tsv is byte-identical to the version at the tag.
  - git diff prereg-v1..HEAD on the 15 frozen files is empty. PREREG.md itself changed only by appending
    Dev-1..Dev-4 and the Clarification.

Tag:
  - prereg-v1 = e775f54, 2026-10-28 14:20:48 +0100. It is an annotated tag and an ancestor of HEAD.

Eukaryote outputs (runs/ is not tracked; filesystem mtimes):
  - Ancestra families stage first file 14:50:13; scores 21:26 onward; hogs.orthoxml 22:30:51; DONE 22:36:41;
    runs/analysis/euk/events.parquet 22:39:39.
  - So no eukaryote family/HOG/event output existed before the tag. stages/input (13:20) holds inputs only.

First commit containing results/primary/euk/: 1d106ed 2026-10-28 22:43:16. This is 8 h 22 min after the tag.

Primary-test-related commits (+0100):
  e775f54 10-28 14:20:48  PREREG freeze (tag prereg-v1)
  9ccbf1c 10-28 14:29:41  statistics engine final; PREREG Dev-1 (H1/H2 parametric bootstrap)
  507b4eb 10-28 15:21:23  phylogenomics + primary H9 (p = 0.27)
  4177757 10-28 15:40:00  RER / Kr-Kc (exploratory); primary_tests.py, rates.py committed
  f8ba008 10-28 15:40:07  primary-analysis driver, positive-control scorer (prok dry run)
  c6e93c2 10-28 18:19:19  gain ladder (prok); PREREG Dev-2, Dev-3
  5a646f2 10-28 18:31:38  FIX stem_path priority (volvocine path stem)
  58fa05f 10-28 19:47:27  tree calibration; PREREG Dev-4 (15a4eb9 19:47:31 timestamp fix)
  50e5a7e 10-28 20:08:58  PREREG clarification (5ac3946 20:09:02 timestamp fix)
  eee4d1a 10-28 22:37:59  eukaryote Ancestra run complete
  1d106ed 10-28 22:43:16  primary H4/H5/H12 (euk) + controls 10/16 provisional
  214ded2 10-28 23:03:10  positive-control failure diagnosis
  cc2e03c 10-28 23:28:15  H4/H5 robustness (orthology-free; post-hoc re-placement)
  a43f118 10-29 00:55:22  primary family complete, Holm, controls 11/16, gain ladder euk
  0ae59b4 / 1413621 / ff430fb 10-29 00:55-00:57  research tree, structure (exploratory)

Minor flag: the Dev-1 entry says "~14:45", but it was committed at 14:29:41 (9ccbf1c). The stated time is later than
the commit. This is harmless, since no eukaryote output existed before 14:50. Correct it as was done for Dev-2/Dev-4.

---------------------------------------------------------------------------------------------------------------
CHECK 7  Analysis registry multiplicity                                                              PASS with FLAGS
---------------------------------------------------------------------------------------------------------------
results/analysis_registry.tsv has 227 rows.

Raw counts by registry tier:
  exploratory 178, secondary 27, primary 14, exploratory_posthoc 8.
  (There is no 'sensitivity' tier; sensitivity analyses carry tier exploratory/secondary with 'sensitivity' in notes.)

Duplicates:
  - 22 exact duplicate rows: same hypothesis, statistic, n, k, scheme, tier, p and notes, written by re-runs.
      * H4/H5/H12 rows from 21:41 UTC were rewritten at 23:54 UTC.
      * H1H2 rows from 23:11 UTC were rewritten at 23:54 UTC.
      * H4/H5 Pfam rows at 14:38 and 14:39.
      * positive_controls at 21:42 and 23:54.
  - 5 additional 'Holm across 5 primary tests' annotation rows repeat the 5 primary p-values.
  - Of the 14 primary-tier rows, therefore, only 5 are distinct tests (H1H2, H4, H5, H9, H12). The 9 extra rows are 4
    exact duplicates and 5 Holm annotations.

Distinct analyses after removing the 27 duplicate/annotation rows (200 rows):
  category      distinct  with p-value  descriptive
  primary            5          5            0
  secondary          8          5            3
  sensitivity       40         40            0     (notes contain 'sensitivity'; includes secondary-set sensitivities)
  posthoc            8          8            0
  exploratory      139        121           18     (includes 8 H9 seed-stability rows)
  TOTAL            200        179           21
=> Multiplicity to disclose:
  - 5 primary tests (the Holm family);
  - 174 further non-primary p-values in total (5 secondary, 40 sensitivity, 8 post hoc, 121 exploratory).
Recommendation: de-duplicate the registry, or add a 'superseded' flag, before it is cited as the multiplicity count.

---------------------------------------------------------------------------------------------------------------
Observations (not failures)
---------------------------------------------------------------------------------------------------------------
- H1/H2, T_volvocine_tgv: the raw O/E is < 1 (est -0.61), yet u = 0.929. The reason is that the NB-GLM gives a very low
  yhat (14.1) for this short 2-edge path stem (summed branch length 0.0178). The u for this stem therefore depends
  on the path-stem branch length, which is the quantity the bug had affected. The primary run used the post-fix value
  (branch_values shows stem_path TGV_clade;Gonium_Volvocaceae).
- H12, Monte-Carlo p: the verifier got 0.0227 and the report 0.0205. The Holm-adjusted value is about 0.10-0.11 in
  either case, so the conclusion is unchanged.
```
