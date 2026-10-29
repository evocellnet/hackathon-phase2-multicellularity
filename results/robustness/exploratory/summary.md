# Robustness of the three positive-looking exploratory results (POST HOC checks)

Written by the PI from the analysing agent's report (the agent could not write .md files). Numbers: TSVs in this directory, code
scripts/robustness/exploratory_posthoc.py, and registry rows whose note starts `posthoc` (86 in results/analysis_registry.tsv, 66 in the de-duplicated view; corrected 2026-10-29 from an earlier "102" that did not match the file).
All checks were designed after seeing the exploratory results; they can only weaken, not strengthen, a claim.

## Main point
Of three exploratory signals, one disappears (convergent domain loss), one survives most but not all checks (more new domain
architectures on transition stems), and one survives subsampling but rests on Pfam only (less co-option in animals and streptophytes).

## 1. Convergent Pfam-domain loss — DISAPPEARS
Original: 421 domains lost in ≥ 2 primary transitions vs 192 expected (p = 0.043 permute; 0.005 matched); reproduced exactly.
- Leave out Dictyostelia: 87 vs 117 expected (p = 0.60); leave out Dictyostelia and Rhodophyta: 43 vs 76 (p = 0.68).
- Genome-based species only: 88 vs 193 (p = 0.77); genome + BUSCO ≥ 80%: 157 vs 170 (p = 0.42).
- Loss counts fall with the number of multicellular species (Spearman ρ = −0.58, p = 0.0008), not with sister species (ρ = −0.17) or BUSCO.
- The whole excess is one pair: fungal hyphae × Dictyostelia share 177 lost domains vs 61 expected; the two reduced proteomes
  (Dictyostelia × Rhodophyta) share fewer than expected (69 vs 96).
- Function: cilium/flagellum domains strongly enriched among convergently lost domains (15% vs 1.9%, OR 17); metabolism depleted
  (OR 0.47); DNA repair not enriched. Excluding 358 cilium-linked domains: 258 vs 170 (p = 0.17); with a null preserving each
  transition's number of losses: 258 vs 264 (p = 0.70).
- Interpretation: the only convergent loss is the known flagellum loss in fungi and dictyostelids, i.e. the same signal as H12.

## 2. New Pfam architectures on stems (H11) — SURVIVES leave-one-out and transcriptome exclusion; WEAKENS under completeness filter
Original: residual T = 0.80, p = 0.005.
- Leave-one-transition-out: all p ≤ 0.017 (without Metazoa p = 0.017); 5 of 7 stems at percentiles 0.83–0.98.
- With BUSCO completeness added to the model: p = 0.017. Genome-based species only: p = 0.0004–0.0065 across versions.
- Restricting to species with BUSCO ≥ 80%: p = 0.16 (0.34 without Metazoa).
- Interpretation: suggestive that transition stems create new domain combinations more often than comparable branches, but sensitive
  to proteome completeness; not a robust finding.

## 3. Co-option index lower in Metazoa and Streptophyta — SURVIVES subsampling (slightly weaker)
- Subsampling the multicellular clade to 3–5 species (100 subsamples, seed 11): Metazoa C = 0.76–0.79, Streptophyta 0.84–0.87, vs
  0.97–0.99 for similarly subsampled controls; all subsamples below the control median; percentile 0.024 → 0.071.
- Caveats: Pfam over-represents animal- and plant-specific families (annotation bias toward well-studied lineages); the HOG-based
  version of the index shows no effect (results/exploratory/euk/summary.md §2).
