# Research tree — genomic transitions to multicellularity (Ancestra)

Legend: `[ ]` planned · `[~]` in progress · `[x]` done · `[-]` abandoned/negative · `[!]` key finding

START: 2026-10-28T12:21:22+01:00 (see START.txt). Budget: 48 h allocation; target finish ≤ 36 h after start, leaving margin.

## Status summary (latest first)

### M2 — 2026-10-28 14:25 (t ≈ 2.1 h): pre-registration frozen
Design went through three review rounds (5 + 3 reviewers + a final pre-registration audit) and is frozen as DESIGN v3.1 / PREREG.md
(tag `prereg-v1`, 15 hashed inputs) before any eukaryote Ancestra output existed. The reviews changed the project substantially: raw
Ancestra gains are inflated 2.1–2.6× on simulations, so the primary gain statistic is sequence-level novelty (DIAMOND, all
proteomes outside the clade); transitions are coded as state changes with an explicit eligibility rule (7 primary eukaryotic origins,
3 nested complex transitions secondary, 4 prokaryotic descriptive); two stems became path stems after the audit; the relabelling null
was replaced by a calibrated permutation scheme after it proved anti-conservative on the real tree.

Honest expectation, written before results: with 7 primary stems and 30 control branches the design only detects ≈2–6-fold effects for
gene-content tests (power at 1.5-fold 13–35%), but is sensitive for the rate test (H9). Running: eukaryote Ancestra `fast` (family
stage, ~1 h in), prokaryote run (scores stage), Pfam (full + pre-registered subset), gain-ladder and phylogenomics modules being
built. Next: wait for Ancestra; parse; positive controls; primary tests; then exploratory analyses.

### M1 — 2026-10-28 12:45 (t ≈ 0.4 h)
Literature review done (LITERATURE.md, ~9k words, 249 tool-verified references; 16 testable hypotheses H1–H16). Data scouting done:
151 core + 104 optional eukaryotic candidates, 79 core prokaryotic; sources UniProt 2026_03, EukProt v3, LukProt v1.5.1, NCBI
(Phaeoexplorer). "Eukaryoma" could not be found (closest: EUKARYOME, an rRNA database) and *Sorogena*/*Sorodiplophrys* have no public
proteomes — flagged for the user, not blocking. EukProt and LukProt archives downloaded; Pfam 38.2, BUSCO lineages and three conda envs ready.

DESIGN.md v1 integrates the hypotheses into a pre-registered primary family (H1/H2, H4, H5, H8, H9, H12) and is under review round 1
by five independent reviewers (external PI, evolutionary biology, statistics, data/algorithms/compute, AI failure modes). In parallel
the dataset build agent is assembling ~200 eukaryotic + ~80 prokaryotic proteomes and guide trees, and the HOG-parser agent is
validating gain/loss/duplication extraction on simulated data. Decision: start dataset build before design review ends (cheap to adjust taxa).

### M0 — 2026-10-28 12:30 (t = 0.2 h)
Project initialised: PROMPT.md saved verbatim, git repo created, Ancestra cloned (commit ffb0bbe, frozen method `freeze-v1`),
conda environment building, NCBI taxonomy and Pfam-A downloading. Two scouting agents are running in parallel: a literature
reviewer (→ LITERATURE.md) and a data scout (→ species candidate tables and source inventory).

Next: draft DESIGN.md from the two scouting reports, run up to three rounds of external review, then build the dataset and
launch the Ancestra run(s). Decision so far: separate eukaryote and prokaryote Ancestra runs, because orthology across domains
adds cost without informing the transitions (all prokaryotic transitions are compared within their own phyla).

## Tree

- [~] 0. Setup: environment, record files, generic resources (taxdump, Pfam)
- [x] 1. Literature review → LITERATURE.md (16 hypotheses; 249 refs)
- [x] 2. Data scouting → data/species_candidates*.tsv, docs/DATA_SOURCES_SCOUTING.md
- [x] 2b. HOG parser (mcell.hogs): exact on true HOGs; Ancestra gains inflated 2.1–2.6× on sims
- [x] 2c. Statistics engine (mcell.branchstats), calibrated permutation null
- [~] 2d. Gain-evidence ladder (mcell.gainladder)
- [~] 2e. Phylogenomics: BUSCO-marker tree, gene trees, H9 and RER (mcell.rates)
- [x] 4. Dataset: 187 euk (2.97 M proteins) + 82 prok (0.41 M); BUSCO; guide trees (62 edits); transitions.tsv
- [x] 3. Design → DESIGN.md v3.1 after 3 review rounds; PREREG frozen (tag prereg-v1)
- [~] 5. Ancestra runs: euk `fast` and prok `fast` running (no 40-proteome timing test: skipped since dataset was under the 3.0 M cap and the watchdog guards memory)
- [ ] 6. Analyses
- [ ] 7. Verification by independent agents
- [ ] 8. Report

## Decision log
- D0 (t=0): Eukaryote and prokaryote runs kept separate (see M0).
- D1 (t=0.3h): Dataset build launched before design review completes; taxon changes from review are cheap before the Ancestra run.
- D2 (t=0.4h): "Eukaryoma" not found → proceed with UniProt + EukProt + LukProt + NCBI; question logged for the user.
- D3 (t=1.2h): Eukaryote Ancestra `fast` launched (12 threads, 140 GB watchdog) with the guide tree as built; if the round-3 audit changes the tree, only the merge stage is re-run (S0/S2 cached). Prokaryote run (4 threads) and Pfam (8 threads) in parallel.
- D4 (t=1.3h): Volvocine primary transition coded at the TGV stem (all sampled colonial species descend from it); single- vs two-origin interpretation is a reported caveat; Tetrabaena-tip and Gonium+Volvocaceae rows kept exploratory.
- D5 (t=1.3h): Eligibility rule refined before results: low-quality species (BUSCO C<50%) do not count toward the ≥2 per side, rather than disqualifying a transition; M must be a single clade. Result: 7 primary eukaryotic origins, 7 secondary (3 nested complex, 4 prokaryotic). Script: scripts/prereg/apply_eligibility.py.
- D6 (t=1.6h): Statistics engine calibration on the real tree (before results): literal relabelling null anti-conservative (14.5% at α=0.05) because |U|=30 and matched sets overlap → primary scheme `permute` (4.9%); H1/H2 statistic = NB-GLM residual fitted on controls (matching poor: depth SMD 1.77). Prokaryotes |U|=2 → descriptive only.
- D7 (t=1.7h): Round-3 audit F1/F2 adopted on phylogenetic grounds (Lindsey 2021; Kiss 2019): volvocine and fungal-hyphae stems become 2-edge path stems ending at Gonium_Volvocaceae and Amastigomycota. Consequence (not the reason): fungal stem becomes flagellum-loss-expected, so H12 has 2 loss-expected stems and stays primary.
