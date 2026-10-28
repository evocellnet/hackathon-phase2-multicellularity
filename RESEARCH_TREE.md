# Research tree — genomic transitions to multicellularity (Ancestra)

Legend: `[ ]` planned · `[~]` in progress · `[x]` done · `[-]` abandoned/negative · `[!]` key finding

START: 2026-10-28T12:21:22+01:00 (see START.txt). Budget: 48 h allocation; target finish ≤ 36 h after start, leaving margin.

## Status summary (latest first)

### M0 — 2026-10-28 12:30 (t = 0.2 h)
Project initialised: PROMPT.md saved verbatim, git repo created, Ancestra cloned (commit ffb0bbe, frozen method `freeze-v1`),
conda environment building, NCBI taxonomy and Pfam-A downloading. Two scouting agents are running in parallel: a literature
reviewer (→ LITERATURE.md) and a data scout (→ species candidate tables and source inventory).

Next: draft DESIGN.md from the two scouting reports, run up to three rounds of external review, then build the dataset and
launch the Ancestra run(s). Decision so far: separate eukaryote and prokaryote Ancestra runs, because orthology across domains
adds cost without informing the transitions (all prokaryotic transitions are compared within their own phyla).

## Tree

- [~] 0. Setup: environment, record files, generic resources (taxdump, Pfam)
- [~] 1. Literature review → LITERATURE.md
- [~] 2. Data scouting → data/species_candidates*.tsv, docs/DATA_SOURCES_SCOUTING.md
- [ ] 3. Design → DESIGN.md, ≤3 review rounds
- [ ] 4. Dataset build (download, QC, one-protein-per-gene, guide tree) → DATA.md
- [ ] 5. Ancestra runs (toy test → timing test → full eukaryote run; prokaryote run)
- [ ] 6. Analyses
- [ ] 7. Verification by independent agents
- [ ] 8. Report

## Decision log
- D0 (t=0): Eukaryote and prokaryote runs kept separate (see M0).
