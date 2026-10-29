# Department-head review of round 1 and assessment of the round-2 proposals

Independent senior reviewer (Fable agent), 2026-10-29; numbers checked against the result files.

## Summary judgement

This is an honest, pre-registered study whose headline numbers are supported by the files, and whose main product is methodological: it shows where Ancestra's reconstruction of ancestral gene content cannot be trusted. Its question, one genomic signature shared by seven very different origins, was one the design could not answer: with 7 transition branches and 30 controls, power at a 1.5-fold effect was 13-35% (results/prereg/power_table.tsv), so the null results are absence of evidence, as the report states; only the rate test genuinely excludes a 1.5-fold average shift. Round 2 should move from "all origins" to well-resolved single events and the simple-to-complex step.

## What was wrong or weak

1. The framing averages over heterogeneous events: facultative, aggregative dictyostelids beside obligate clonal lineages, animals coded as one unicellular-to-complex step, 26 fungal species against 4 chytrids. A mean percentile over seven such values (REPORT §4.3) is nearly uninformative.

2. Control matching failed, and the report understates it. Stems are deeper than controls by 2.06 control standard deviations (results/primary/euk/smd.tsv, ALL row), not 1.77 (§4.3). The volvocine "6-fold excess" is extrapolation (14 predicted gains at depth 13; controls near depth 5). The power table used placeholder branch lengths of 1.0 (power_table.tsv), yet §4.3 says "effects planted on the real tree".

3. The rate test has 10 controls, 9 containing multicellular species, and its p-value moves from 0.19 to 0.61 with the random orientation seed (results/h9/summary.md); average over all 1,024 orientations instead. The largest value (dictyostelids) rests on two transcriptome sisters although Physarum has a genome.

4. "11 of 16 positive controls" overstates the validation: five passes (C05, C09, C12, C14, C16) are per-species domain counts that never use Ancestra's events, while all three controls of deep-node HOG placement (C03, C06, C07) failed (results/positive_controls/euk/controls.tsv).

5. The flagellum test is a confirmation: four of its six stems are scored on "no loss expected" (data/h12_expected_direction_frozen.tsv), and it spends Holm correction the other tests needed.

6. Data: 36 of 187 proteomes are transcriptomes; five below 50% BUSCO completeness were kept; the complex transitions of red algae, brown algae and heterocystous cyanobacteria were excluded only for lacking two sister species.

7. "Eukaryoma" was declared not found (RESEARCH_TREE D2); EukaryOMA exists (eukaryoma.unil.ch: 1,328 genomes, 33.8 million proteins, 1.57 million deepest HOGs).

## What was good and deserves more effort

Hashed pre-registration, a calibrated null (the planned null rejected 13.8% of no-effect data), simulated gain inflation (2.1-2.6-fold), a bias rule fixed before results, independent re-derivation and a citation audit set a high standard. The most valuable science is the positive-control failure diagnosis (results/positive_controls/euk/failure_diagnosis.md): domains gained inside older gene lineages are invisible to HOG gains; large multi-domain families are over-merged; duplications are placed one node too deep. Two exploratory signals merit follow-up: more novel genes at the three nested simple-to-complex transitions (p = 0.016, 1.9-3.3-fold; Merényi et al. 2020, Bowles et al. 2020) and 1,131 new domain architectures on the animal stem against a control median of 34 (Ekman et al. 2007).

## What was not evaluated; the lead's proposals

- EukaryOMA (HIGH). An independent precomputed HOG set with 7-fold more genomes: first cross-validate Ancestra's stem gains and duplications event by event, then use its genomes. Pitfalls: different species tree, HOG levels and species mapping.
- Per-transition dense Ancestra runs (HIGH). Small runs avoid the 73% giant component and deep-node artefacts and allow the accurate preset and dense sampling (LukProt, EukaryOMA; genome sisters). They sharpen each event, not the cross-event sample size; controls must come from within each run.
- Complex from simple (HIGH). Three events already show a signal; add Florideophyceae, Phaeophyceae and heterocystous cyanobacteria with proper sisters, and date domain-level innovations. Pitfall: nested events are not independent of their simple ancestors (Pereira Lobo et al. 2025; Nagy et al. 2018).
- Specific Pfam families, architectures, functions, recurrence (MEDIUM-HIGH, cheap). Round 1 tested class aggregates and a domain-wide scan topped by transposon domains; per-transition tables of expanded families with annotation are missing. Mask transposons; use per-species counts in a phylogenetic regression over all 187 species.
- Motif-binding domains, cell death, cell-type transcription factors (MEDIUM). Sound; define Pfam sets before looking, avoid animal-centric lists; same power limits.
- 3D structures on the GPU (MEDIUM). Ancestra's 3Di rescue lowered benchmark precision (F -5.0) and its 3Di scoring backend is unimplemented, so do not run the whole pipeline in structure mode. Use ProstT5 3Di plus Foldseek-GPU on targeted sets to finish the structure check for all 7 transitions; AlphaFold 3 at proteome scale is infeasible.
- Marcotte co-fractionation data (LOW-MEDIUM). Cox et al. 2024 (bioRxiv 10.1101/2024.05.26.595818; 31 species, LECA-centred) has too few species per transition to test lineage-specific complex conservation; use it only to ask whether stem gains and retained duplicates join ancient complexes.

## Recommended priorities for the next 30 h

1. Checkpoint; correct the SMD, power-table and positive-control statements; fix the rate-test null.
2. EukaryOMA cross-validation of every stem event.
3. Dense per-event Ancestra runs (accurate preset), animals and fungi first.
4. Simple-to-complex analysis with the extended event set and domain-level dating.
5. Per-transition family, architecture and function tables; recurrent expansions; motif, cell-death and TF sets.
6. GPU structure checks for all 7 transitions.
7. More calibration replicates with a simulator that includes gene birth.
8. Complex-membership annotation from Cox et al.
