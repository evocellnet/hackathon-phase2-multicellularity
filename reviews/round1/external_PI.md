# External PI review — round 1 (DESIGN.md v1)

## Overall verdict

Questions and the §1b literature mapping are strong. The design is **not yet the most informative one**: (1) the central statistic pools heterogeneous branches so that one
or two well-sampled transitions (Metazoa, fungi) will drive it, and the "control" class mixes unicellular branches with branches inside
multicellular clades; (2) it has no defence against the documented Ancestra failure mode that most directly fakes the headline result
(HOG splits at deep levels, README "Limitations"); (3) it tries to do too much: 7 questions, 16 hypotheses, two Ancestra runs,
structure, ASR and dN/dS in ≈40 h. Fix (1)–(2), cut (3).

## Most likely misleading headline

**"Transitions to multicellularity are marked by a burst of new gene families, with convergent gains of adhesion/signalling/TF
families."** This is the result the pipeline is biased to produce:

1. **Deep-merge splits.** Ancestra's own README reports merge splits at deep levels (Bilateria, Eukaryota) and S0 family
   completeness ≈0.72. A split HOG is a spurious *gain*. Deep stems (Metazoa, Fungi,
   Streptophyta) are inflated more than shallow controls.
2. **Sampling asymmetry.** Multicellular sides have 71 core taxa (4 transcriptomes); unicellular sisters have 45 (11 transcriptomes;
   3 of 6 Metazoa sisters). A family missing from an incomplete sister
   gets placed as a gain on the multicellular stem. More descendant species also means more gains (union). A GLM covariate does not remove this (cf. Richter 2018: ≈372 families moved).
3. **Ascertainment.** The toolkit categories (cadherin, integrin, RTK, Wnt, Notch, homeobox) are animal-defined, and Pfam models
   are built mostly from animal, plant and fungal sequences. "Enrichment" will be largest for animals by construction.
4. **Pooling.** About half of the 24 eukaryotic "transitions" are single-taxon terminal branches (including all five aggregative
   ones except Dictyostelia), some transcriptome-only (*Guttulinopsis* has BUSCO 11). Their gains are mostly artefacts.

## Ranked changes

### MUST

**M1. Make the transition, not the branch, the unit of inference; define controls properly.**
- Classify every branch as (a) *transition stem* P→M, (b) *unicellular control*, where both endpoints are unicellular, or
  (c) *within-multicellular*. The current "all other branches" mixes (b) and (c).
- Primary statistic per transition: the percentile of the stem among ≥200 matched (b)-branches. Match on stem length (subs/site),
  number of descendant species, mean descendant BUSCO, transcriptome fraction, and depth class.
- Combine across transitions with a sign test and Stouffer's Z, plus leave-one-transition-out. Forest plot per transition.
- Keep the NB-GLM as a secondary view, with transition as a random effect.
- Pre-register a **strict set** of transitions with ≥2 genomes on each side: Metazoa, Dictyostelia, Volvocales, Streptophyta,
  Rhodophyta, Phaeophyceae, Pezizomycotina, Agaricomycetes, and choanoflagellate colonies if both sides qualify. All other
  transitions form the "broad" sensitivity set.
- State the power up front. With 8 strict transitions, only 8/8 concordant gives a two-sided p < 0.05 (0.008); 7/8 gives p = 0.07.

**M2. Audit gains before reporting them (cost: hours).** For every HOG gained on transition and matched-control branches, run DIAMOND
(sensitive) of its members against all proteomes outside the clade. Add Foldseek where AlphaFold models exist.
- Classify each gain as **novel** (no hit), **diverged** (only a structural or weak sequence hit), or **split** (a confident hit to a
  HOG already present at P, which means a merge failure).
- Headline gains = novel set; also report S0-family-level gains.
- Apply the same audit to control branches, so that the split rate becomes a measured covariate.

**M3. Rarefy post hoc, without re-running Ancestra.** Recompute the placement of each gain as the LCA of the species carrying the
deepest HOG while:
- (a) subsampling each multicellular clade to k = 2 or 3 species, matched to its sister;
- (b) dropping each sister or outgroup taxon in turn (the Richter 2018 test);
- (c) excluding transcriptome taxa.

H1/H2 is only claimed if it survives (a) and (b). For terminal single-species transitions, report losses and domain-bearing gains
only, never raw gains.

**M4. Fix the RER test.** RERconverge-style residuals remove the genome-wide branch effect by construction, so test 3.5(i) cannot
detect a genome-wide shift.
- For genome-wide constraint (H9), use **outgroup-anchored relative-rate contrasts**: d(outgroup, multicellular) − d(outgroup,
  sister) on concatenated near-single-copy alignments, one contrast per transition, combined with a sign test.
- Run gene-level RER within clade subsets: Holozoa, Fungi, Chlorophyta+Streptophyta, Stramenopiles, Amoebozoa. A eukaryote-wide
  ≥70% single-copy filter leaves too few genes and saturated deep branches.

**M5. Cut to fit the time budget.**
- **Drop:** H15 (site-level convergence); H14 (the 79 multicellular-biased prokaryote genomes cannot support HGT calls); LG+C20
  (use LG+G4); the per-HOG Wilcoxon expansion test (2–5 species per side has no power; use Ancestra duplications with the M1 null);
  ASR + structure beyond ≤2 families.
- **Demote H8** (clonal > aggregative) to exploratory. Only Dictyostelia is well sampled; other aggregative lineages are
  confounded with terminal-branch/transcriptome status.
- **Prokaryote run:** run it first, because it is cheap (0.4 M proteins). Use it as the pipeline dry run and for H9/H12 in
  cyanobacteria only.
- **Size cap:** cap the eukaryote run at about 2.6–2.9 M proteins (core plus ≤25 optional taxa that fill sister gaps). Ancestra
  `fast` used 221 GB at 4.75 M proteins and the node has 150 GB, so measure RSS in the 40-proteome timing test before choosing
  `accurate`.
- Audit the `data_type` labels. *Sycon* EP00117 ("genome", 50,731 proteins, source `SCIL_T-PEP`) looks transcriptome-derived.
  *Telonema* and *Diphylleia* (41–59 k proteins) need redundancy collapse.

### SHOULD

**S1. Unbiased first, toolkit second.**
- For H5/Q2 the primary analysis is Pfam-wide: every domain, convergence across strict transitions against the M1 null, with BH
  correction.
- Extend the pre-registered toolkit with non-animal effectors: pherophorins/VARL, alginate and fucoidan enzymes, sulfotransferases,
  hydrophobins, cellulose synthases, dictyostelid polyketide synthases.
- Report how many hits fall in animal-defined versus lineage-specific domains.

**S2. Use reversals as converse tests.** Label reversal branches explicitly: *S. cerevisiae*, *Yarrowia*, *S. pombe*
(Saccharomycotina and Taphrinomycotina losses of hyphal multicellularity), *Cryptococcus* (yeast within Agaricomycotina). Test
whether HOGs gained at the corresponding transition are lost there more often than matched HOGs of the same age (H12b).

**S3. Literature concordance scorecard, written before results.**
Numeric expectations with tolerance (Richter 2018/Paps & Holland 2018 Metazoa vs Choanozoa novelty; Kiss 2019 no kinase/adhesion
expansion at hyphae; Denoeud 2024 brown-algal gains; Hanschen 2016 few volvocine gains), scored reproduced/contradicted/not addressed.

**S4. Analysis code ready before Ancestra finishes.** Build and validate every analysis script on the 40-proteome timing run and on
the prokaryote run while the main run executes.

**S5. Benchmark Pfam annotation time** on 5 proteomes before committing 3 M proteins. If it is too slow, annotate the full proteomes
for architectures and HOG representatives (≤10 per HOG) for HOG function.

### COULD

**C1. Structure-aware selection.** For ~300 clade-level single-copy HOGs, estimate per-site rates (IQ-TREE `-wsr`, or branch-partitioned rates) for the
  multicellular and the sister subtree. Stratify by AlphaFold burial and pLDDT. Relaxed constraint predicts a disproportionate increase at buried, high-pLDDT
  sites.

**C2. Radical-to-conservative substitution ratio (Kr/Kc)** as a protein-only proxy for the strength of selection, using the same
sister contrasts as M4.

**C3. Volvocales-only dN/dS** with HyPhy RELAX (Gonium/Tetrabaena/Pleodorina/Volvox vs Chlamydomonas), only for genes with
dS < 2.

## Feasibility summary

Feasible with M5 (Ancestra ≈6–10 h on ≤2.9 M proteins; M2 ≈2 h; M3 minutes). Without M5, analyses will still be running at report
time and the over-sampled animal and fungal transitions will dominate by default.
