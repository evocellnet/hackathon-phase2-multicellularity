# External PI review — round 2 (DESIGN.md v2)

## Verdict
v2 handles the round-1 headline risk well (sequence-novel gains primary, three branch classes, state coding, PREREG, manifest). No round-1 must was ignored outright. Four problems remain:
- the control design cannot be carried out;
- the primary tests are under-powered and one has wrong directions;
- the sister covariates were dropped;
- the compute plan is stale.

## 1. Round-1 musts dropped or mis-implemented
- **Sister covariates dropped** (stats M2, PI M1). Matching and the GLM use only descendant-side covariates. Spurious stem gains come from sparse *sister* sampling (Richter 2018).
- **No n per hypothesis and no power statement** (stats M1, PI M1).
- **Contamination screen dropped** (data M4), with no reason given.
- **Pfam-count expansion is only a robustness check** (data S6), although H4/H5 test exactly the superfamilies that the 7,500-gene cap splits.

## 2. Consistency and feasibility
- **"≥200 matched within-unicellular controls" do not exist.**
  - About 105 of the 196 tips are not multicellular, and they are scattered. That gives roughly 60–90 internal within-unicellular branches.
  - Few match the Metazoa, Fungi or Streptophyta stems in depth or length, and four-way tertile strata will hold 0–3 branches each.
- **Holm across 5 tests with ~6–7 eligible U→S transitions.** The smallest exact two-sided Wilcoxon p is 0.031 at n = 6 and 0.016 at n = 7. Both exceed the first Holm threshold (0.01), so H9 can never be significant.
- **H12 "loss expected" is wrong in 2 of 3 cases** (inherited from evo-review A3):
  - *Allomyces*, which is inside the hyphal M, has flagellated zoospores.
  - Bryophytes and *Chara* have flagellated sperm. Flagella were lost in Zygnematophyceae (a reversal) and in seed plants, not on the streptophyte stems.
- **Proteome size and memory.** The frozen set has 3.30 M proteins, above the 3.0 M cap. Predicted `fast` peak is ~125 GB (exponent 1.57, anchored at 221 GB for 4.75 M). §3.1 and §3.2 also give different caps (3.0 vs 2.6 M).
- **Pfam takes ~11 h at 12 threads**, not the 7 h slot. It gates H4/H5, Q3 and Q4.
- **Code volume** (~12 new modules, incl. a Python RERconverge) is the main schedule risk; otherwise feasible in ~44 h.

## 3. Ranked changes
**Must**
1. **Replace "≥200 matched controls".**
   - Primary statistic: per-transition observed/expected (O/E) from an NB-GLM fitted to controls only, with a prediction interval.
   - Relabel within coarse strata (depth × branch-length tertile), and report the eligible-control count per transition.
   - Use nearest-neighbour matching as a sensitivity analysis.
2. **Add sister covariates** to matching and the GLM: number of sister species, minimum sister BUSCO, sister transcriptome fraction.
3. **Before PREREG, tabulate n and the smallest attainable p for each hypothesis.** A test that cannot reach Holm becomes estimate-only (H9: Hodges–Lehmann estimate with CI).
4. **Correct H12.** For fungi, expect loss on the Zoopago+Mucoro+Dikarya branch; on streptophyte stems, expect no loss. If fewer than 3 loss-expected transitions remain, make H12 exploratory.
5. **Cut to ≤3.0 M proteins now, by a written rule.** Drop only over-sampled multicellular or background taxa, never a sister, P-side or eligibility-determining species. Start Pfam immediately.
6. **Report H9 as a rate shift.** Generation time and mutation rate also drive it. Claim "relaxed constraint" only from a within-protein ratio, by promoting Kr/Kc or buried/exposed rates.

**Should**
7. **H4/H5:** a headline requires sign agreement from the Pfam copy-number (Wagner) version. Report the share of toolkit genes in Leiden-split families.
8. **Homology-detection failure** inflates sequence-novel gains on fast lineages (Weisman et al. 2020), which confounds H1/H2 with H9. Add a ≥100 aa or Pfam-bearing subset and a gene-rate covariate.
9. **Contamination screen:** reuse the ladder's DIAMOND hits (near-zero extra cost) and report novelty without the flagged proteins.
10. **Pre-register a drop order at t = 20 h:** ASR → gene-level RER → tree variants → `accurate` → re-merge → structure. Use R RERconverge if it installs in 30 min, else drop gene-level RER.
11. **Co-option index:** also compute it on frozen toolkit HOGs present at M (stats S4). The expansion-defined set relies on an underpowered per-family test.
12. **Pre-state what a positive-control failure implies:** gain claims become descriptive.

**Could**
13. Add size-matched random Pfam sets as negative controls for H4/H5.
