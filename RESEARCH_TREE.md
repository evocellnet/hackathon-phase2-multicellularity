# Research tree — genomic transitions to multicellularity (Ancestra)

Legend: `[ ]` planned · `[~]` in progress · `[x]` done · `[-]` abandoned/negative · `[!]` key finding

START: 2026-10-28T12:21:22+01:00 (see START.txt). Budget: 48 h allocation; target finish ≤ 36 h after start, leaving margin.

## Status summary (latest first)

### M1 — 2026-10-28 ~13:50 (t ≈ 1.5 h)
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
- [~] 2b. HOG parser (mcell.hogs) + simulation validation
- [~] 4. Dataset build (started early, in parallel with design review)
- [~] 3. Design → DESIGN.md v1, review round 1 running (5 reviewers)
- [ ] 5. Ancestra runs (toy test → timing test → full eukaryote run; prokaryote run)
- [ ] 6. Analyses
- [ ] 7. Verification by independent agents
- [ ] 8. Report

## Decision log
- D0 (t=0): Eukaryote and prokaryote runs kept separate (see M0).
- D1 (t=1.3h): Dataset build launched before design review completes; taxon changes from review are cheap before the Ancestra run.
- D2 (t=1.5h): "Eukaryoma" not found → proceed with UniProt + EukProt + LukProt + NCBI; question logged for the user.
