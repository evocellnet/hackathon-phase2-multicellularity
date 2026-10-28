# Structure of sequence-novel gains (DESIGN §5 / Q6) -- EXPLORATORY

Everything below is **exploratory** (tier `exploratory` in `results/analysis_registry.tsv`, hypothesis labels `Q6_structure_*`); not in the Holm family.
Code: `scripts/structure/{sample_gains,novel_gain_structures,analyze_structure}.py`; data/versions: `docs/STRUCTURE_DATA.md`.

## Design
- Gains: `results/gainladder/euk/gains_classified.tsv.gz`, `ladder == seq_novel`. Stem = path-stem edges of each of the 7 primary transitions
  (`stem_path`, all edges after the first node). Query = the first (lowest rep_rank) UniProt-sourced representative of the HOG (accession = id after `<code>_`).
- Sample: <=60 stem gains per transition (seed 11) from gains with a UniProt representative, and an equal-sized sample pooled over the 20 H4 matched control
  branches of that stem (`results/primary/euk/matched_sets.tsv`; gain-weighted, so large control clades dominate). 840 queries (420 stem, 420 control; `results/structure/query_set.tsv`, `sampling_counts.tsv`).
- AFDB v6 models (EBI per-protein files), mean pLDDT >= 70. Foldseek exhaustive 3Di+AA search, E <= 1e-3, TM (query-normalised) >= 0.5 = "structurally old".
- Reference A (**Swiss-Prot**, 590k AFDB models): hit must be to a taxon OUTSIDE the gaining clade. Stem gains: transition's M (`data/structure_ref_exclusion.tsv`).
  Control gains: all taxa under the NCBI LCA of the control clade's species (`exclude_taxids`; 154/420 controls are on deep branches whose LCA is Eukaryota, so for them only prokaryote hits count).
  (A first version that excluded only M for controls inflated the control fraction; corrected, the earlier registry rows were removed and replaced.)
- Reference B (**P-side proteome DB**): only transitions whose reference species are fully downloaded: T_metazoa, T_fungi_hyphae, T_dictyostelia, T_streptophyta (4/7).
  Coverage of AFDB models: Monobre 8406/9155, Salpros 10341/11614, Gonapro 13341/13795, Planfun 15481/16631, Ostrtau 7172/7723, Micrcom 9634/10064 (89-97%).
  T_rhodophyta / T_px_filaments / T_volvocine_tgv await the remaining downloads (Porppur, Rhodmar, Nannsal, Phaetri, Chlarei). Controls whose gaining clade contains a reference species are dropped (65 of 240).
- Stats: Wilson CIs; per-stem Fisher; Mantel-Haenszel OR / CMH (continuity-corrected) stratified by stem; Breslow-Day for heterogeneity.

## Usable models (both references share this)
| set | sampled | no AFDB model | pLDDT < 70 | usable (searched) | median mean-pLDDT |
|---|---|---|---|---|---|
| stem | 420 | 94 | 263 | 63 (15.0%) | 57.4 |
| control | 420 | 58 | 273 | 89 (21.2%) | 60.7 |

No usable model: stem 357/420 = 85.0% (Wilson 81.3-88.1) vs control 331/420 = 78.8% (74.6-82.4); MH OR 1.51 (1.06-2.15), CMH p = 0.025 (Fisher pooled p = 0.025; Breslow-Day p = 0.001: heterogeneous across stems, strongest in T_rhodophyta 93% vs 75%, T_px_filaments 93% vs 78%).
"No model" = new UniProt accessions absent from AFDB (stem 22%, control 14%); "low pLDDT" dominates (63% / 65%). Sequence-novel gains are therefore mostly disordered/low-confidence or unmodelled; stem gains slightly more so.

## Fraction structurally old (of usable-model queries)
Swiss-Prot reference (outside own clade):
| stem | stem old | control old |
|---|---|---|
| T_metazoa | 0/10 (0%; 0-28) | 3/14 (21%; 8-48) |
| T_fungi_hyphae | 3/14 (21%; 8-48) | 3/9 (33%; 12-65) |
| T_dictyostelia | 4/17 (24%; 10-47) | 5/8 (63%; 31-86) |
| T_streptophyta | 2/8 (25%; 7-59) | 2/17 (12%; 3-34) |
| T_rhodophyta | 3/4 (75%; 30-95) | 1/15 (7%; 1-30) |
| T_px_filaments | 2/4 (50%; 15-85) | 3/13 (23%; 8-50) |
| T_volvocine_tgv | 1/6 (17%; 3-56) | 2/13 (15%; 4-42) |
| **pooled** | **15/63 = 23.8% (15.0-35.6)** | **19/89 = 21.3% (14.1-31.0)** |

MH OR (stem/control) = 0.96 (0.45-2.07), CMH p = 0.92; pooled Fisher p = 0.84; per-stem nominal p < 0.05 only for T_rhodophyta (n = 4) and none after any correction; Breslow-Day p = 0.009.
Sensitivity (controls with LCA = Eukaryota dropped; 36 controls left): control 9/36 = 25.0% (13.8-41.1); MH OR 0.92 (0.37-2.30), CMH p = 0.95.

P-side proteome DB (4 transitions): pooled stem 5/49 = 10.2% (4.4-21.8) vs control 6/22 = 27.3% (13.2-48.2) after dropping controls whose clade contains a reference species; MH OR 0.36 (0.09-1.40), CMH p = 0.26. Per stem stem/control: metazoa 0/10 vs 3/7, fungi 3/14 vs 2/8, dicty 1/17 vs 0/3, streptophyta 1/8 vs 1/4.
(Full tables: `per_stem_*.tsv`, `cmh_*.tsv`, `model_availability_*.tsv`.)

## Interpretation (descriptive)
- About 3/4 of sequence-novel gains have no confident structure to test, so any "fold is old" statement is about a ~15-20% subset biased to ordered, well-modelled proteins.
- Among testable gains roughly a fifth to a quarter are structurally old by either reference; stems do not differ from matched controls (CI for OR spans 0.45-2.1). Power is very low (n_stem testable = 4-17 per stem).
- Stem gains are somewhat more often without a usable model (OR 1.5), a weak hint of more disorder/novel short proteins on stems.
- Caveats: query = one UniProt representative per gain; AFDB lacks models for recent UniProt accessions; Swiss-Prot is eukaryote-sparse for protists (few outside hits possible for protist clades); LCA-based control exclusions are very conservative for deep clades; proteome reference runs for 3 transitions pending.
