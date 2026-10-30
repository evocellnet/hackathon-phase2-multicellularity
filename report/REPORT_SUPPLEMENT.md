# Supplement to "From many origins to single events": detail moved from the main report

This supplement belongs to report/REPORT_v2_final.txt and is AI-written, like the report itself (see the human-written disclaimer at the top of the main report). It holds the detail a reader of the main text can skip: errata, data and reconstruction detail, verification detail, full family lists, the engineering of the structure work, extended methods and the full figure captions. Section numbers in brackets refer to the main report.

Contents: S1 Round-1 errata · S2 Panel construction and reconstruction · S3 Cross-method comparison: detail and verification · S4 Simple-to-complex novelty: notes, one-effect-per-panel sensitivity, further sensitivity analyses, positive controls · S5 Family sets: support, architectures, recurrence and per-event families · S6 Structure: implementation, validation, structure mode, structural novelty · S7 Protein complexes · S8 Elapsed time and compute · S9 Extended methods · S10 Full figure captions

---

## S1. Round-1 errata (main §2)

A separate agent re-checked the round-1 report's specific claims with new scripts (results/round1_errata/ERRATA.txt). None changes a round-1 conclusion.

Table S1. Corrections to the round-1 report.

| item | round-1 statement | corrected statement | effect |
|---|---|---|---|
| 1. Depth imbalance | stems 1.77 control SDs deeper than matched controls | 2.06 (mean depth 8.7 vs 5.0 edges from the root); the 1.77 came from an earlier stem definition with placeholder branch lengths and was never recomputed | strengthens the report's own warning |
| 2. Power table | effects "planted on the real tree" | the frozen script used branch length 1.0 everywhere (the branch-length tree did not exist yet); rerun with real branch lengths | no power value changed by more than 0.025 |
| 3. Positive controls | "Ancestra recovered 11 of 16 expected events", read as a validation | of 16 controls, only 10 need Ancestra to place an event on a specific branch; 7 passed, but one of these (C15, "no expansion of stem duplications") is a null control a method that misses duplications would also pass: 6 of 9 excluding it (5 of 8 high-confidence). All recent-stem gains (Wnt, TGF-beta, GRAS, Hex1) passed; all three deep-node placements failed; two of the passes did not use Ancestra at all | Ancestra's event placement is validated on recent stems and fails at deep nodes |
| 4. Rate test | p reported from one random seed | each control needed a random orientation; p ranged 0.19–0.61 across seeds; exact test over all 1,024 orientations: orientation-averaged p = 0.31, no orientation gives p ≤ 0.05; orientation-free pooled rate ratio 1.06 (1.10 reported) | conclusion (no 1.5-fold shift) unchanged; 9 of 10 rate-test controls contain multicellular species on one side, a limitation that remains |

Further review points: 36 proteomes were transcriptomes; five proteomes below 50% BUSCO completeness were kept; three complex transitions (florideophyte red algae, complex brown algae, heterocystous cyanobacteria) had been excluded for lack of sister species; and the round-1 record declared "Eukaryoma" not found, although the database exists (eukaryoma.unil.ch: 1,328 genomes, 33.8 million proteins) and the request had been misread.

Other round-1 results summarised in main §2: the planned null distribution (relabelling controls as stems) called 13.8% of no-effect datasets significant at nominal 5% and was replaced; the smallest Holm-adjusted p was 0.10, for flagellar loss in fungi and dictyostelids; positive controls were recovered in 11 of 16 cases (failures: the cadherin cytoplasmic domain, gained inside an older gene lineage; integrins with EGF-domain proteins, over-merged; T-box duplications placed one branch too deep); the nested simple-to-complex lead was p = 0.016 (1.9–3.3-fold; Pezizomycotina, Agaricomycetes, land plants), and the animal stem carried 1,131 new multi-domain architectures against a control median of 34.

---

## S2. Panel construction and reconstruction (main §4)

*Sources.* Fifty-three proteomes were new relative to round 1. Sources were the round-1 sources plus additional genomes from UniProt, NCBI, EukProt and LukProt. EukaryOMA was not used as a sequence source for the panels: its downloadable family file carries only internal identifiers, and its sequences were used in round 2 only to map proteins for the cross-method comparison (S3), where 0.50 million EukaryOMA sequences were retrieved in bulk through its bulk_retrieve interface. Known gaps: no annotated genomes were available for *Physarum*, *Zygnema*, *Penium*, hornworts, *Coleochaete* or *Nuclearia*, and the Holozoa and Amoebozoa panels still contain 10 and 9 transcriptomes.

*De novo annotation caveat.* Two volvocine genomes, *Yamagishiella unicocca* and *Eudorina* sp., had assemblies but no published protein annotation. The data agent predicted their genes *de novo* with AUGUSTUS (Stanke et al. 2006), without RNA or protein evidence, using the *Chlamydomonas* gene model. Their completeness (BUSCO 71% and 69%) is similar to the annotated *Gonium* (71%), but de novo gene models differ systematically from curated ones. A sensitivity panel without these two genomes (R4_curated) was agreed before either reconstruction was run.

*Reconstruction.* All panels used Ancestra's fast preset. A patched, output-identical version of Ancestra's community-detection step (verified byte-identical on test data) let large components be split in parallel. Measured run times were 15–45 minutes per panel (cyanobacteria 15, red algae 16, fungi 19, volvocines 22, stramenopiles 25, amoebozoa 36, streptophytes 37, holozoa 45), with peak memory of 2.6–16.9 GB. Every panel passed the format validator and the parser's conservation check (on every branch, genes at the parent plus gains minus losses plus duplications equal genes at the child). Because panels were this cheap, all eight ran, and the planned go/no-go cut-off at 8 hours never applied. The panels were built in about 2.4 hours.

*Sequence novelty.* Sequence-novel gains were called as in round 1, searching each panel's proteomes plus all round-1 proteomes outside the panel (eukaryotic, or prokaryotic for cyanobacteria). Across panels, 41–80% of Ancestra's HOG gains were sequence-novel (lowest in Amoebozoa, where transcriptomes are common; highest in Fungi). An independent verifier recounted raw gains from the OrthoXML for five branches (all matched) and re-called novelty for 200 sampled families with its own representatives (196 of 200 calls agreed; results/r2/verification/panels/VERIFICATION.txt).

---

## S3. Cross-method comparison: detail and verification (main §5)

*Protein mapping.* EukaryOMA's OrthoXML carries only internal OMA identifiers. Proteins were mapped by identical sequence (MD5 of our sequence against the sequence checksum that EukaryOMA's interface reports for all 2.31 million proteins of the 125 shared genomes), which mapped 79.8% of our proteins. For the 28 species with < 90% checksum mapping (different assemblies or gene-model versions), 0.50 million EukaryOMA sequences were retrieved in bulk (bulk_retrieve interface) and reciprocal best DIAMOND hits were accepted at ≥ 95% identity and ≥ 50% coverage of both proteins, for a total of 87.1% mapped (1,757,405 of 2,016,931). Fifteen species remained < 50% mapped (results/eukaryoma_xval/summary.txt).

*Duplications on the single-species stems.* The main text gives duplication agreement for the six internal stems (4–8%, as on controls). On the two single-species (tip) stems that also map exactly, agreement is higher: 37% (Ancestra to EukaryOMA, 440 of 1,203) and 67% (EukaryOMA to Ancestra, 449 of 670), computed as the difference between the pooled eight exact stems (11% and 22%) and the six internal stems. A duplication on a tip is two copies in one genome, which both methods detect more easily.

*Expected agreement: two estimators.* Both are binomial regressions (GLMs) of same-branch agreement on node depth and log number of descendant tips over non-stem branches of the pruned tree. The producer's fit uses 86 non-stem branches and predicts 22.0% for the six stems (11.9% for EukaryOMA-to-Ancestra gains); the independent verifier's fit uses 97 non-stem internal branches and predicts 24.7% for the stems and 13.1% for the controls. The difference is the choice of baseline branches. Controls of the same depth as the stems (depth 6–9, n = 8) agree at 21.5%.

*Independent verification* (results/r2/verification/eukaryoma/VERIFICATION.txt). A separate agent re-implemented the comparison without using the producer's code. It confirmed the protein mapping (150 of 150 checksum matches identical; 50 of 50 reciprocal-best-hit matches from the correct genome at ≥ 95.5% identity; 87.13% mapped) and reproduced all 20 per-branch agreement values and all pooled values within 2 percentage points. Gains matched exactly and losses within 3.4%. Three of the 10 duplication counts differed by 12–15% (for example, Ancestra's streptophyte stem, 2,651 against 2,343), because the original filtering step kept groups left with a single child after other species were removed, which made the parser merge some genuine duplications away. Agreement barely changed when restricted to events whose members were at least 80% mapped (stem gains 28.9% to 32.5%), so the disagreement lies in the methods' family structures, not in unmapped proteins.

*Overlapping versus non-overlapping species.* DESIGN_ROUND2 §6 promised the comparison split by species shared between the two resources and species present in only one. This split is not possible: both HOG sets were restricted to the 125 shared species before placement, by construction, so the comparison has no non-overlapping species.

*Round-1 direction and controls.* Averaged over five primary stems, the stems' percentile among controls was 0.63 with EukaryOMA and 0.72 with Ancestra on the same pruned tree (0.5 = typical); the sign agreed on four of five stems. With branch labels stripped, Ancestra places 9 T-box duplications on the animal stem, as expected.

---

## S4. Simple-to-complex novelty (main §6)

*The branch-length floor (deviation D-PI-1).* The analysis code applied log(max(branch length, 10⁻⁴)) to branch length and to distance to the outgroup; the written specification (results/r2/panels/METHODS_DEFINITIONS.txt) says only "log(branch length)". The floor changes exactly one control, an amoebozoan branch (Cavosteliaceae) of length 2.8 × 10⁻⁵; through the shared slopes (log branch length 0.37 → 0.42) it changes the ranking of the volvocine control residuals, moving the *V. africanus* p from 0.429 to 0.571 and the combined p from 0.0504 to 0.072. Both versions are reported, the as-specified one as confirmatory.

*One-effect-per-panel sensitivity* (scripts/r2/panels_sensitivity_onepanel.py; results/r2/panels_sensitivity_onepanel/summary.txt). Inputs: the independent verifier's per-branch leave-one-out residuals, as specified and as implemented.

Table S2. Dependence-aware sensitivity for the confirmatory test.

| quantity | as specified | as implemented |
|---|---|---|
| volvocine events, log RR (RR) | 2.61 (13.6), 2.70 (14.9), 0.24 (1.27) | (13.2), (14.2), (1.19) |
| six volvocine control residuals | 1.23, 0.29, 0.23, 0.06, −1.14, −3.01 | 1.19, 0.29, 0.28, 0.03, −1.08, −2.90 |
| volvocine panel effect: mean log RR (RR) | 1.85 (6.4) | 1.80 (6.1) |
| control triples (of 20) with mean ≥ observed; largest triple mean | 0; 0.59 | 0; 0.59 |
| volvocine panel p | 1/21 = 0.048 | 0.048 |
| Fucophycidae RR, p (5 controls) | 2.01, 0.167 | 1.84, 0.167 |
| two-panel Stouffer p (Holm, ×2) | 0.031 (0.062) | 0.031 (0.062) |
| pooled RR, one effect per panel: volvocine variance assuming perfectly correlated events | 3.86 (95% CI 0.00–10,040) | 3.63 |
| same, assuming independent events | 4.76 (0.01–2,742) | 4.40 |
| same, inverse-variance-weighted panel effect (convention of the secondary analysis) | 4.77 (0.01–2,893) | 4.41 |

Reading: both panel p-values sit at their floors, so 0.031 is the smallest two-panel p the design can produce, above the 0.025 that Holm's second step requires. The four-event version could in principle reach 0.019. The volvocine panel effect is carried by the two tip events; ordinary unicellular tips in the volvocine panel already score about 6.6-fold against the same internal-control model, about the size of the panel effect. With two panels, the pooled effect is uninformative.

*Tip artefact detail.* Ordinary non-event unicellular tips score about 8-fold (volvocine panel 6.6, brown-algal panel 9.9) under the internal-control model; origin tips about 19-fold, simple-to-complex tips about 41-fold, reversal tips 20.7 (10.7–40.0; nine events). With the implemented floor, the tip-control comparison gives combined p = 0.57 and RR 0.80 (as specified 0.49 and 0.81, 0.32–2.03). Simple-to-complex tips against tip controls give 0.99 (0.47–2.1); tip controls exist only for the volvocine and brown-algal panels. The problem was found by the analysis agent while inspecting interim output from five panels; the design had frozen tip events against internal controls before the data were seen.

*Further sensitivity analyses.* Removing the two de novo annotated volvocine genomes lowered the combined p to 0.033 (pooled RR 4.84), because removing close relatives raises the apparent novelty of *V. africanus*; with tip controls it was 0.32. All 15 nested simple-to-complex events (seen and unseen, tip and internal) give RR 11.2 (one effect per panel, k = 6: 10.3), dominated by tips. Internal origins pool to 0.83 (0.41–1.71) as implemented, internal reversals to 0.63 (k = 2). The pre-registered primary result remains the eight-panel run.

*Positive controls on panels.* Fourteen of 18 control scorings were recovered. The T-box control (gain at Holozoa, loss on the choanoflagellate stem, duplication on the animal stem), which failed in round 1, passed in the dense Holozoa panel. The cadherin and integrin controls still failed, as did a Notch-related control and the pherophorin threshold control. The GRAS transcription-factor gain on the stem of land plants plus Zygnematophyceae, which passed in round 1, failed on the streptophyte panel: no gain is placed on that stem. Dense sampling therefore fixes some deep-node placements and breaks at least one recent one.

---

## S5. Family sets (main §7)

*Cell-death test, further support.* 21 of 22 families have positive individual estimates. Without the NB-ARC and NACHT pair, p = 0.018. Within the lineages where they occur, lineage-restricted families also point upward: Bcl-2, Death, CARD, DED, BIR and gasdermin are expanded in animals relative to choanoflagellates; HeLo, HET and NACHT in Pezizomycotina fungi; NB-ARC in red algae; NB-ARC, TIR and MACPF in land plants. Of 20 significant within-origin changes, 19 are expansions; the exception is metacaspase, which is reduced in animals, where caspases take its role. The family list is data/r2_sets/cell_death.tsv; the inclusion rule results/r2/w4/P3_rule.txt; the independent re-derivation results/r2/verification/w4/.

*Architectures, detail.* Land plants carry on average 84 more NB-ARC proteins per species than the Zygnematophyceae, mostly NB-ARC–LRR with or without N-terminal TIR or coiled-coil domains. Red algae: NB-ARC–WD40 in all five multicellular species. Filamentous fungi: HET +37 per species in Pezizomycotina against budding yeasts. Animals: a sponge, a lancelet and a sea urchin carry hundreds of Death or NACHT proteins, the fly and nematode only a handful. In the dense panels, NB-ARC–WD40 families were duplicated at the origin of florideophyte red algae (8 duplicated families against 1.6 expected), and many apoptosis and immune families on the animal stem. NLR families evolve fast and are hard to annotate (fragments, split gene models), so lineage means are driven by a few genomes. Sources: results/r2/w4_celldeath/.

*Motif-binding sensitivity.* A Bayesian phylogenetic Poisson-lognormal model (MCMCglmm) gave a small positive estimate (posterior probability of no excess 0.03). We consider that model anti-conservative for many-family sets, because it treats lineage-specific expansions of each family as independent evidence. The only positive leave-one-origin-out estimate (without volvocines) is +0.005; without animals, −0.063. Within origins: animals 29 families up, 7 down; hyphal fungi 0 up, 12 down. Transcription factors: LOB also expanded in streptophytes.

*Recurrence across origins.* For each of the seven round-1 origins and three nested events, every Pfam family (transposon domains masked) was tested for a copy-number difference between the multicellular clade and its sister clade (4,226–9,523 families per event, q < 0.1 within event). Animals had the most significant families (1,661 expanded, 694 contracted); volvocines almost none (0 expanded, 10 contracted). The pre-registered null (permuting species labels within events) was too permissive: mixing multicellular and sister clades destroys the divergence that makes many families significant in any clade comparison, so every observed count appeared significant. A conditional null that keeps each origin's number of significant families but shuffles which families they are gives: expansions at two or more origins 42 against 54 expected (p = 0.98); contractions 157 against 90 (p = 0.001; 20 at three or more against 3.1). Most recurrent expansions: Cu/Zn superoxide dismutase, CAP (PR-1-like secreted proteins) and a ZBED1 RNase-like domain (three origins each). Most recurrent contractions: cAMP phosphodiesterases and Hint/intein-like domains (four origins each). Several strongly "contracted" families are bacterial-like glycoside hydrolases, probably bacterial prey or contaminant sequences in the transcriptomes of unicellular sisters.

*Domains at simple-to-complex events.* For the seven internal simple-to-complex events, Pfam domains of families gained or duplicated on the event branch were compared with those of control branches. Domains appear at two or more events more often than expected (1,556 vs 1,057), but this is a size effect: events carry more duplicated families, and after rarefying each branch to at most 50 annotated families the excess disappears (p = 0.77). The top "recurrent" domains are large generic families (peptidyl-prolyl isomerases, subtilases, P450s, ABC transporters). The protein kinase domain is the only domain in sequence-novel families at five of seven events (p = 0.023, q = 0.06).

*Per-event families* (results/r2/panels/<panel>/w3_event_families.tsv; transposon families excluded).
- Volvocine germ–soma events (*Astrephomene*, *Pleodorina*, *V. carteri*, *V. africanus/reticuliferus*): all four duplicate pherophorin-like DUF3707 families (6, 6, 21 and 10 duplicated families), matrix metalloprotease Peptidase_M11 families (2, 5, 7 and 11) and guanylate-cyclase families (5, 4, 5 and 15). All are duplications of existing families, not new families. (An earlier summary said that only the *Volvox* events carry M11 and three of four carry guanylate cyclases; that came from a top-domain shortlist and is corrected here.)
- Heterocystous Nostocales: two-component signalling (response regulators, histidine kinases, GAF), pentapeptide repeats, Uma2.
- Land plants: leucine-rich-repeat receptor kinases and F-box families (Shiu & Bleecker 2001).
- Florideophyceae: C2H2 zinc-finger and bZIP transcription factors.
- Fucophycidae: kinases, WD40, ankyrin and TPR families, and the largest duplication count of the seven internal events (5,924; the *Chara* tip event has 9,889).
- Complex fungi (Pezizomycotina, Agaricomycetes): transporters, Zn2Cys6 transcription factors, short-chain dehydrogenases, cellulose-binding domains.

Table S3. Share of sequence-novel families with any Pfam domain (results/r2/panels/pooled/final_8panels/w3/w3_event_summary.tsv).

| event | novel families | with Pfam | % |
|---|---|---|---|
| *V. carteri* (tip) | 1,300 | 18 | 1.4 |
| *Astrephomene* (tip) | 892 | 29 | 3.3 |
| *V. africanus/reticuliferus* | 104 | 4 | 3.8 |
| *Pleodorina* (tip) | 1,621 | 167 | 10.3 |
| Fucophycidae | 493 | 63 | 12.8 |
| Florideophyceae | 709 | 120 | 16.9 |
| heterocystous Nostocales | 308 | 53 | 17.2 |
| Agaricomycetes | 139 | 28 | 20.1 |
| Pezizomycotina | 417 | 160 | 38.4 |
| land plants | 259 | 132 | 51.0 |

---

## S6. Structure (main §8)

*Motivation.* In round 1, only about a fifth of sampled novel gains had a usable AlphaFold model, and of those, roughly a tenth to a quarter had a known fold outside the clade. The reviewer warned that Ancestra's structure-assisted family rescue had lowered precision in its own benchmark and that its structure-based scoring backend was not implemented, and advised against whole reconstructions in structure mode; round 2 therefore planned two targeted uses only.

*GPU implementation.* A pilot tested three ways of running ProstT5 on the A100. Foldseek's own GPU build silently ran the model on the CPU, at about 15 proteins per second; the job looked successful. The language-model library's standard text-generation mode was far too slow, because it predicts letters one at a time. An own implementation running only ProstT5's encoder followed by its small convolutional prediction head ran at 48 proteins per second and agreed with Foldseek's ProstT5 output at 97.1% of residues (4,996 proteins). Long proteins were cut into windows of 1,000 residues overlapping by 250, which raised throughput from about 21,000 to 34,000 residues per second. Windowing is exact for proteins up to 1,000 residues (about 92–93% of all proteins) and agrees at 90.3% of residues with full-length prediction for proteins of 1,000–3,000 residues. The queue was ordered by panel, animals and their relatives first, because the full set of 322 proteomes might not finish (the first 265 completed proteomes alone hold 3.80 million proteins and 1.71 billion residues; an early estimate of 1.43 billion residues for all 322 was too low). By the end of the round, 263 of 322 proteomes were complete; each structural analysis used the proteomes complete when it started, which is why the divergence analysis covers only three panels.

*Validation.* For 24 completed species, every protein had exactly one 3Di record of the right length, with no invalid characters (108 proteins over 3,000 residues in one sponge were skipped by a length cap). For 4,981 proteins with an AlphaFold model, agreement was 65.5% of 2.55 million residues. Accuracy was similar across model-confidence and disorder classes, slightly lower for windowed long proteins (64.9%) and for proteins from a protist without Swiss-Prot entries (62.6%, versus 66.5% for mostly animal Swiss-Prot proteins, which are probably in ProstT5's training data).

*Structure mode of Ancestra.* Small families are searched by Foldseek (3Di plus amino-acid letters) against representatives of large families and attached if they match; the guard re-checks each attachment. Holozoa panel (37 proteomes, 628,756 proteins; 99.98% with predicted 3Di letters). The rescue attached 25,035 small families (41,914 genes) to 2,400 large families; the guard confirmed 17%. Deep-node controls: T-box passed in all three runs; cadherin and integrin failed in all three (a domain gained inside an older family; an integrin that is truly present in one choanoflagellate). Agreement with EukaryOMA: 23.3% (default), 20.7% (structure mode; −2.6 percentage points, 95% CI −5.6 to +0.4), 21.5% (with guard; −1.8, −4.8 to +1.2); every comparison points down, including the reverse direction and other internal branches. Side effects: losses +7.5% on internal branches and duplications +24% on the animal-plus-choanoflagellate stem without the guard, the same loss of precision Ancestra's own benchmark had shown. Structure mode cost about 10–14% more CPU time on top of several GPU hours for the letters. Verdict: results/r2/w5/r1/verdict.txt.

*Structural novelty, detail.* Up to 200 sequence-novel families per branch (seed 11), searched through the same up to three representatives used for the sequence search. Families that are not sequence-novel were found again in 77–96% of cases; shuffled queries gave 0 hits in 9,747 families. The structural search covers only the 18–40 proteomes of each panel, whereas sequence novelty was defined against the panel plus all 187 round-1 proteomes. The first analysis (43.7% vs 33.3% unmatched; Mantel–Haenszel odds ratio 1.44; branch permutation p = 0.004) was concentrated at simple-to-complex events (1.87) rather than origins (1.15), and its reported adjustment covered only search-space size. The verifier found: unmatched families median 211 vs 519 amino acids; 49% vs 32% at least half disordered; 92% of families under 100 amino acids unmatched; events carry 8.1% vs 4.2% families under 100 amino acids. Exhaustive search rescued 21% of unmatched families (23.0% on events, 19.5% on controls); raw odds ratio then 1.23 (p = 0.07). At simple-to-complex events: adjusted 1.73 (1.20–2.51), exhaustive 1.79, among longer, mostly ordered proteins 2.25, branch-level 1.23 (0.95–1.58); in red algae the effect goes the other way. Origins adjusted 1.04. Sources: results/r2/w5/, results/r2/verification/w5/.

*Divergence analysis, detail.* For 300 single-copy genes per panel, the 3Di substitution model could not be estimated with IQ-TREE (GTR20 failed), so a panel-specific exchange matrix was estimated from pairs of closely related sequences (documented deviation; Poisson sensitivity). Gene-level one-sided p-values for the crown effect: 6 × 10⁻¹² (animals and relatives), 1.6 × 10⁻⁶ (fungi), 0.020 (volvocines); with branches as units only fungi reached p < 0.05 (0.020), and not in a ratio-of-sums version. The 21 event stems were not shifted (12 of 21 positive; sign test p = 0.33). Split-half control: −0.031 (animals), −0.002 (fungi), −0.276 (volvocines; five crown branches, regression to the mean on a short-branch panel). ProstT5 letters are themselves predicted from sequence, so this ratio partly measures the model. Sources: results/r2/w6/.

---

## S7. Protein complexes (main §9)

Only the published LECA complex map of Cox et al. (2026) was used; the raw co-fractionation elution profiles were not reanalysed. No unicellular sister of animals was fractionated. Proteins were mapped through UniProt accessions for 13 of our species. Odds ratios for complex-subunit involvement among stem duplications: animals 1.14, fungal hyphae 1.48, brown algae 0.82, volvocines 0.79; dictyostelids and streptophytes at control level. The one complex enriched on the animal stem traced to a single expanded family. Across 125 complexes, co-retention of subunits (keeping all or none) relative to a prevalence-matched null was slightly lower in multicellular species (paired p = 0.008, without phylogenetic correction).

---

## S8. Elapsed time and compute

Round 1: started 12:21 CET 28 October 2026 (START.txt), about 13.6 hours on 24 CPUs. Round 2: started 15:07:47 CET 29 October 2026 (START_ROUND2.txt) on 48 CPUs, 500 GB and one A100-40GB; the pre-registration tag prereg-r2 is stamped 15:48:45 (about 40 minutes in). Post-freeze deviations of the panel analysis: D-P1 to D-P11 plus D-PI-1 (12 entries; results/r2/panels/DEVIATIONS.txt).

---

## S9. Extended methods

**Data.** Round-1 data: 187 eukaryotic and 82 prokaryotic proteomes (UniProt reference proteomes 2026_03; EukProt v3, Richter et al. 2022; LukProt v1.5.1, Sobala 2024; NCBI annotations), one protein per gene, proteins ≥ 50 amino acids, transcriptomes reduced to the longest ORF and clustered at 95% identity with cd-hit (Fu et al. 2012). Panel data: data/dense/<panel>/ (proteomes, species.tsv, guide_tree.nwk, tree_edits.md, transitions.tsv); two de novo proteomes in the volvocine panel predicted with AUGUSTUS 3.5.0 ab initio, *Chlamydomonas* model, no hints (Stanke et al. 2006; data/dense/R4/ANNOTATION_NOTES.txt). Completeness: BUSCO 6.1.0, odb10 sets (Manni et al. 2021). EukaryOMA release All.Oct2026, HOG OrthoXML md5 296586b7eaae71fb1ce5676c5cc3736e. Cox et al. LECA complex map: Zenodo 10.5281/zenodo.17716736. Sources, versions, checksums and licences: DATA.md.

**Orthology.** Ancestra (github.com/evocellnet/ancestra-mvp-ai, commit ffb0bbe, tag freeze-v1), preset `fast`, with the output-identical parallel Leiden wrapper (Traag et al. 2019; scripts/ancestra_patch/); DIAMOND 2.2.8 (Buchfink et al. 2021) inside Ancestra. Events parsed with mcell.hogs (conservation law checked on every branch).

**Branch lengths.** Panel markers: root HOGs single-copy in ≥ 80% of species (≤ 150 markers), MAFFT --auto (Katoh & Standley 2013), trimAl -automated1 (Capella-Gutiérrez et al. 2009), IQ-TREE 3 (Wong et al. 2026) LG+G4 on the fixed topology.

**Sequence novelty.** DIAMOND blastp --very-sensitive -k 0 -e 1e-5 --query-cover 30 --max-hsps 1, two passes, ≤ 3 representatives per family; target database = panel proteomes plus round-1 proteomes outside the panel (eukaryotic, or prokaryotic for cyanobacteria).

**Panel statistics.** Pooled negative-binomial GLM on 55 internal non-crown controls (panel intercepts; shared slopes for log branch length, log descendants, transcriptome fraction, mean BUSCO, log distance to nearest outgroup); leave-one-out control residuals; depth classes (≤ 3, 4–6, ≥ 7 edges) with widening to ≥ 5 controls; one-sided permutation p; Stouffer combination; REML random effects with modified Hartung–Knapp CI (q* = max(1, q)). Specification: results/r2/panels/METHODS_DEFINITIONS.txt; deviations: results/r2/panels/DEVIATIONS.txt.

**Family tests.** Pfam-A 38.2 (Paysan-Lafosse et al. 2025) with pyhmmer (Larralde & Zeller 2023; gathering thresholds); disorder with metapredict 3.1 (Emenecker et al. 2021). Motif-binding domain list from ELM (Kumar et al. 2024). PGLS with phylolm (Ho & Ané 2014; Pagel's λ, Pagel 1999) on the round-1 IQ-TREE species tree, one-sided t-test, parametric bootstrap CI (1,000 replicates). Sensitivity: phylogenetic Poisson-lognormal GLMM (MCMCglmm, Hadfield 2010; two chains). Scans: per-family negative-binomial GLM with log proteome-size offset, likelihood-ratio test, BH q < 0.1 within event, fold-changes shrunk with ashr (Stephens 2017); recurrence nulls by label permutation and by conditional family shuffling (B = 1,000). Pfam-to-GO: pfam2go mapping of 2026-07-06 (Mitchell et al. 2015; Ashburner et al. 2000; Gene Ontology Consortium 2026). Rules fixed before testing: results/r2/w4/P3_rule.txt.

**Cross-method comparison.** Sequence-checksum and reciprocal-best-hit protein mapping (≥ 95% identity, ≥ 50% coverage), topological placement of both HOG sets on the 125-species pruned tree, event matching by member proteins (results/eukaryoma_xval/summary.txt; tests in tests/test_eukaryoma_xval.py). Independent re-implementation: results/r2/verification/eukaryoma/VERIFICATION.txt.

**Structure.** ProstT5 (Heinzinger et al. 2024) encoder plus CNN head, own implementation (PyTorch 2.6.0+cu124, transformers 5.22.0, sentencepiece 0.2.2) on one A100-40GB; windows of 1,000 residues with 250 overlap; validation against Foldseek 3Di (Foldseek 10.941cd33; van Kempen et al. 2024) of AlphaFold DB models (Varadi et al. 2024), 5,000 randomly sampled proteins from 24 completed species (seed 11; results/r2/w5_3di_validation/).

**Structure mode of Ancestra.** Holozoa panel, preset fast, --family-rescue 3di with and without --family-rescue-guard; family and scoring stages reused from the default run (identical input checksums), rescue, merge and output recomputed. Decision rule fixed before any output (results/r2/w5/decision_operationalisation.txt): deep-node positive controls must improve and stem-gain agreement with EukaryOMA (same code as the cross-method comparison, 29 shared species) must not fall. Verdict: results/r2/w5/r1/verdict.txt.

**Structural novelty.** Up to 200 sequence-novel families per event or control branch (seed 11), ≤ 3 representatives each; Foldseek databases built from our amino-acid and ProstT5 3Di files; search against panel proteins of species outside the clade, --alignment-type 2 (3Di+AA), E < 10⁻³, ≥ 50% query coverage, max-seqs 1,000. Positive control: sequence-old gains; negative control: shuffled queries. Statistics: Mantel–Haenszel, within-panel branch permutation, logistic regression with branch-clustered errors (results/r2/w5/). Independent verification (results/r2/verification/w5/): own re-search of 104 families, metapredict disorder and low complexity, logistic mixed models (lme4) with random intercepts for panel and branch, branch-level quasibinomial regression, covariate-adjusted permutation, and an exhaustive re-search (--exhaustive-search 1) of all 3,280 unmatched event and control families with a shuffled negative control.

**Structure versus sequence divergence.** Three panels (Holozoa; Fungi; Volvocales without the two de novo proteomes); root-level families single-copy in ≥ 80% of species, 300 drawn per panel; Pfam envelopes masked per taxon; MAFFT --auto alignments with 3Di projected onto the same columns; IQ-TREE 3.1.4 branch lengths on the fixed topology (amino acids LG+G4; 3Di with a panel-specific exchange matrix estimated by pair counting, because GTR20 estimation failed; Poisson sensitivity). Statistic: log ratio of 3Di to amino-acid branch length, residual against amino-acid-length deciles of unicellular reference branches; crown versus unicellular contrasts with genes and with branches as units; controls: constant-ratio simulation (AliSim), column-shuffled 3Di, amino-acid split halves (results/r2/w6/).

**Pre-registration, seeds and verification.** Round 1: tag prereg-v1. Round 2: tag prereg-r2 (design and frozen family lists, sha256 in data/r2_prereg/PREREG_R2_HASHES.tsv) and tag prereg-r2-p1 (confirmatory events, data/r2_prereg/p1_events.tsv). Seeds 11 (default), 23, 37, 41, 53. Every confirmatory number was re-derived by an agent that did not produce it (results/r2/verification/). Software environments: envs/ (Python 3.11.18; numpy, scipy and other versions in the pip-freeze and lock files).

**One-effect-per-panel sensitivity.** scripts/r2/panels_sensitivity_onepanel.py reads the verifier's residual tables (results/r2/verification/panels/ and its clip1e-4/ variant), enumerates all C(6,3) = 20 triples of volvocine controls, uses p = (1 + #triples with mean ≥ observed)/(20 + 1), combines the two panels with Stouffer's method, and pools one effect per panel with the same REML/modified Hartung–Knapp code as the verifier. Output: results/r2/panels_sensitivity_onepanel/{summary.txt,onepanel_sensitivity.json}.

---

## S10. Full figure captions

These are the long captions (also in report/figures_v2/CAPTIONS.txt); the main text carries shortened versions.

Colours are the same in every figure: grey = unicellular, green = simple multicellular (and origins of simple
multicellularity), blue = complex multicellular (and simple-to-complex events), orange = reversals,
pink = the four pre-registered (confirmatory) simple-to-complex events.

-----------------------------------------------------------------------------------------------------------
Figure 1. Eight focused reconstructions, each built around one or a few origins of multicellularity.
(a) Species (proteomes) per panel by multicellular state; numbers inside bars are species counts, numbers at
the right are panel totals (227 proteomes in all). (b) Coded events per panel by kind (54 events; the five
alternative codings used only as sensitivity analyses are not shown, giving 59 coded events in total).
Hatched: the origin of animals, a direct change from unicellular to complex because no simple intermediate is
sampled. "Other" is the exploratory true-branching event inside heterocystous cyanobacteria. (c) Species tree
of each panel drawn as a cladogram (tips aligned; branch lengths not to scale). Dots give the state of each
species; coloured branches are event branches, coloured by the state they lead to. Thick pink branches are the
four confirmatory simple-to-complex events: Astrephomene and Volvox carteri (single-species branches),
V. africanus/reticuliferus (Volvocales) and the Fucophycidae stem (Stramenopiles). Most of each tree lies inside
multicellular clades, which leaves 5-9 control branches per panel. Tips are labelled with abbreviated species
names and internal event branches with clade names (text coloured like the event). Sources: data/dense/<panel>/species.tsv,
transitions.tsv, guide_tree.nwk; results/r2/panels/<panel>/event_table.tsv.

-----------------------------------------------------------------------------------------------------------
Figure 2. Two orthology methods, Ancestra and EukaryOMA, often disagree on when gene families arose and
almost never on duplications. Both sets of hierarchical orthologous groups were restricted to the 125 shared
species and placed on the same pruned tree by the same parser. (a) New gene families and (b) duplications
on the six transition stems that map exactly onto the pruned tree (log scale). On the animal stem Ancestra
reports 804 new families and 6,617 duplications, EukaryOMA 2,502 and 196 (3.1-fold and 34-fold differences in
opposite directions). (c) Share of events that the other method places on the same branch. Dark bars: the six
stems pooled; grey bars: 26 control branches inside unicellular lineages; hatched bars: the share expected for
branches of the stems' depth and size (producer's binomial GLM fitted on 86 non-stem branches). Points are single
branches. Ancestra's stem gains are confirmed in 29% of cases (540 of 1,859 mapped) against 11% on controls and
22% expected (producer's binomial GLM on depth and log number of tips over 86 non-stem branches; the
independent verifier's equivalent fit on 97 branches gives 24.7%); EukaryOMA's stem gains in 9% (542 of 5,750). Duplications agree in 4-8% of cases on stems and
controls alike. Agreement is consistency between methods, not accuracy. Sources: results/eukaryoma_xval/
counts_per_branch.tsv, agreement_per_branch.tsv, agreement_stems_depth_adjusted.tsv.

-----------------------------------------------------------------------------------------------------------
Figure 3. Simple-to-complex events do not carry significantly more sequence-novel gene families; their
large per-event ratios come from comparing single-species branches with internal branches. The rate ratio
(RR) is the observed number of sequence-novel families on a branch divided by the number predicted by a
negative-binomial regression fitted on 55 internal control branches of all eight panels (panel intercepts;
shared effects of branch length, number of descendant species, transcriptome fraction, BUSCO completeness and
distance to the nearest outgroup). Bars are 95% confidence intervals (CI); diamonds are random-effects pooled
estimates (95% CI); open circles are tip (single-species) branches; the dashed line marks RR = 1.
(a) The pre-registered confirmatory test, computed exactly as specified (no branch-length floor; deviation
D-PI-1): RR 13.6, 14.9, 1.27 and 2.01 for Astrephomene, V. carteri, V. africanus/reticuliferus and Fucophycidae,
with one-sided permutation p against the 6, 6, 6 and 5 controls of their panels. Combined (Stouffer) p = 0.0504,
pooled RR 4.90 (95% CI 0.63-38.4); not significant after Holm correction (p < 0.025 required). The implemented
variant with a branch-length floor gives combined p = 0.072, pooled RR 4.63. Source: results/r2/verification/
panels/verify_stats.json (scripts/r2/verification/panels/verify_stats.py).
(b) Every tip scores high against internal controls, whatever its biology: 23 unicellular tips that are not
events (geometric mean RR 8.1), 15 origin tips (pooled 19.0), 8 simple-to-complex tips (40.8; pink: the two
confirmatory tips) and 9 reversal tips (20.7, 95% CI 10.7-40.0).
(c) Against the 23 unicellular tips of the volvocine and brown-algal panels (tip controls), Astrephomene scores
RR 0.56 (p = 0.85) and V. carteri 0.80 (p = 0.54). The confirmatory test recomputed with tip controls for its two
tip events gives pooled RR 0.81 (0.32-2.03), combined p = 0.49. Simple-to-complex tips (n = 4) pool to 0.99 and
reversal tips (n = 4) to 1.16.
(d) Internal branches only (exploratory, not pre-registered): the seven internal simple-to-complex events pool
to RR 2.70 (95% CI 1.29-5.63), against 0.83 (0.41-1.71) for 12 internal origins of simple multicellularity and
0.63 for 2 internal reversals (CI 0.00-634). Three of the seven events were seen in round 1, and the lower bound
does not clear the 2.2-fold stem inflation measured in round-1 simulations (dashed line). Values in (b)-(d) use
the implemented model (the floor changes one amoebozoan control only); computed exactly as specified, the
internal simple-to-complex estimate is 2.81 (1.34-5.89) and the internal origins 0.82 (0.38-1.77). Sources:
results/r2/panels/pooled/final_8panels/{event_residuals,branch_residuals,event_class_comparators}.tsv.

-----------------------------------------------------------------------------------------------------------
Figure 4. Among three pre-defined family sets, only cell-death families are more abundant in multicellular
genomes, mostly in complex ones, and the effect depends partly on proteome size. Estimates are coefficients
from phylogenetic generalised least squares (PGLS; Pagel's lambda) on 148 genome-based proteomes (68
multicellular, 80 unicellular or grey-zone), with log proteome size and proteome disorder fraction as
covariates. The response is the equal-weight mean over families of the standardised log copy number, so the
unit is one standard deviation (SD) of a family's log copy number; +0.167 SD for cell death corresponds to a
geometric-mean 1.21-fold excess (95% CI 1.04-1.42). Bars: 95% CI; p: one-sided t-test.
(a) Multicellular vs unicellular (black) and, in a three-state model, simple vs unicellular (green) and complex
vs unicellular (blue), for 12 widespread cell-death families (confirmatory; Holm-adjusted p = 0.014), 115
motif-binding domain families (confirmatory; p = 0.72) and 11 widespread cell-type transcription-factor families
(exploratory; p = 0.30). (b) Leave one of the seven primary origins out: the cell-death estimate stays positive
in all seven (one-sided p 0.0002-0.078; weakest without streptophytes); motif-binding estimates are negative in
six of seven. (c) Cell death under proteome-size sensitivities (independent re-derivation): removing the ten
largest proteomes (six of them multicellular) lowers the estimate by about 28% (p = 0.025, n = 138), and
modelling copy number as a share of the proteome gives +0.083 (p = 0.14). Sources: results/r2/w4/pgls/
pgls_results.tsv; results/r2/verification/w4/v2_pgls_output.txt.

-----------------------------------------------------------------------------------------------------------
Figure 5. The extra cell-death proteins of multicellular lineages are mostly immune-type receptors, and the
receptor type differs between lineages. (a) Mean number of proteins per species carrying each cell-death domain,
in the multicellular clade (MC) of each lineage and in its closest unicellular or simpler sister clade (sis.);
Pfam 38.2 gathering-threshold hits; dots = none. "Caspase" excludes proteins that also carry the metacaspase
domain (outside animals most Peptidase_C14 hits are metacaspases). Left labels group domains by function. Land
plants, streptophytes and red algae gain mostly NB-ARC and TIR proteins; filamentous fungi gain HET (heterokaryon
incompatibility) and NACHT proteins; animals gain death-fold adaptors (Death, CARD, DED), NACHT, caspases, Bcl-2
and BIR proteins; dictyostelids and volvocines show no excess. Several sister clades are transcriptomes (animals,
dictyostelids) or secondarily unicellular (yeasts for Pezizomycotina, Zygnematophyceae for land plants).
(b) NB-ARC and NACHT proteins (STAND ATPases) by architecture: C-terminal LRR, WD40, TPR or ankyrin repeats
(NLR-like immune / non-self receptors), an N-terminal effector domain without repeats (e.g. TIR, coiled coil,
CARD, HET), or neither (often fragments). The expansions use different repeat types: LRR in land plants (plant
NLRs), WD40 in red algae, TPR in brown algae, WD40/ankyrin with NACHT in fungi; animals are mixed. (c) Gene
families carrying a cell-death domain that arose (hatched, upper bar) or were duplicated (plain, lower bar) on
simple-to-complex branches of the eight focused reconstructions (and the unicellular-to-complex origin of
animals); dark = NB-ARC/NACHT families. Ticks: number expected from the share of cell-death families among
Pfam-annotated families on the panel's control branches. q < 0.1: one-sided Fisher test, Benjamini-Hochberg over
32 branch x kind rows (exploratory). * single-species branch (tip). Only the animal origin, Florideophyceae and
Pleodorina (doubtful ARM-repeat proteins with weak cell-death domain matches) and two tip branches exceed
expectation. Sources: results/r2/w4_celldeath/ (fig4b_matrix_families.tsv, fig4b_matrix_stand.tsv,
panels_event_celldeath.tsv); scripts/r2/w4_celldeath/04_figure.py (output fig5_celldeath_architectures.*; the
matrix files keep their historical fig4b_ prefix).

-----------------------------------------------------------------------------------------------------------
Figure 6. The families gained or duplicated at simple-to-complex events fit known biology, but only a small
annotatable minority can be described. (a) Number of gene families gained (sequence-novel) or duplicated on
each event branch that carry the given Pfam domain (transposon-domain families excluded; square-root colour
scale; dots = 0). Columns: four volvocine germ-soma events (three single-species branches) and six internal
simple-to-complex events; pink labels mark confirmatory events. All four volvocine events duplicate families of
pherophorin-like extracellular-matrix glycoproteins (DUF3707), volvocine matrix metalloproteases
(Peptidase_M11) and guanylate cyclases; in all four events these are duplications, not new families.
Heterocystous Nostocales carry two-component signalling families (response regulators, histidine kinases, GAF)
and pentapeptide repeats; land plants leucine-rich-repeat and F-box families; Florideophyceae C2H2 and bZIP
transcription factors; complex fungi (Pezizomycotina, Agaricomycetes) Zn2Cys6 transcription factors,
transporters, short-chain dehydrogenases and cellulose-binding domains. (b) Families with a protein-kinase
domain at the seven internal events: one sequence-novel kinase family at each of 5 of 7 events (darker bars;
recurrence permutation p = 0.023, BH q = 0.06), against 7-41 duplicated kinase families per event (lighter bars).
(c) Percentage of sequence-novel families with any Pfam domain, from 1.4% (V. carteri, 18 of 1,300) to 51%
(land plants, 132 of 259); below 21% at eight of the ten events (e.g. 63 of 493 in Fucophycidae). Sources: results/r2/panels/<panel>/w3_event_families.tsv (primary panels; R4, not the curated
sensitivity panel); results/r2/panels/pooled/final_8panels/w3/{w3_event_summary,w3_recurrent_domains}.tsv.

-----------------------------------------------------------------------------------------------------------
Figure 7. Structure predicted from sequence is accurate enough to search with, but it did not improve Ancestra
and gave no robust evolutionary signal. (a) Share of residues at which the 3Di structural letter predicted from
sequence by our GPU implementation of ProstT5 equals the letter computed by Foldseek from the AlphaFold
database model, by protein length (4,981 proteins, 2.55 million residues; overall 0.655). Blue: Swiss-Prot
proteins, mostly from animals (0.65–0.68 per length bin); orange: other proteins, mostly from the protist
Salpingoeca (0.62–0.63). Proteins of 1,000–3,000 residues were predicted in overlapping windows. Dotted line:
agreement expected by chance from letter frequencies (0.17); dashed line: always predicting the commonest
letter (0.31); grey band: 41–65%, the range reported for an equivalent sequence-based predictor (Heinzinger et
al. 2024). (b) Ancestra's optional structure mode on the animals-and-relatives panel (37 proteomes). Bars:
share of Ancestra's new families on the animal stem, on the stem of animals plus choanoflagellates
("choanoflag."), on both stems pooled, and on all other internal branches that EukaryOMA places on the same
branch (95% Wilson intervals). Grey: default, sequence only; light blue: structure mode; dark blue: structure
mode with its guard. The deep-node positive controls (text, top) were recovered 1 of 3 in every run. The pre-
registered rule required better deep-node controls and no loss of agreement; neither was met, so the default
was kept. (c) Odds that a sequence-novel family has no structural match outside its clade, on event branches
versus control branches (19 event and 55 control branches in eight panels; 8,967 families), from the
independent verifier's mixed models with random effects for panel and branch: unadjusted (grey circles);
adjusted for protein length, disorder, low complexity, family size, branch depth and search-space size (blue
squares); adjusted, with an exhaustive search that finds homologues the default search misses (orange
diamonds); adjusted, with branches rather than families as the units (pink triangles; not available for
origins). Lines: 95% CI. The general effect (top) disappears after adjustment; a lead remains at the seven
simple-to-complex events in family-level models but not with branches as units. (d) Structural change per unit
of sequence change on branches inside multicellular clades minus unicellular branches (log scale; median over
300 genes per panel, 95% bootstrap CI). Blue: real predicted letters; orange: negative control with
letters shuffled among species within each alignment column; grey: control comparing two random halves of the
sequence alignment (near zero in animals and fungi, but -0.28 in volvocines, whose five crown branches cannot
support the test). The shuffled control gives an effect comparable to or larger than the real data, so the
small positive values cannot be distinguished from prediction noise combined with tree shape. Sources:
results/r2/w5_3di_validation/per_protein_accuracy.tsv; results/r2/w5/r1/{xval_agreement_summary,controls_<run>}.tsv;
results/r2/verification/w5/work/glmm_results{,_exh}.txt; results/r2/w6/panel_tests.tsv; script
scripts/r2/figures/fig7_structure.py.
