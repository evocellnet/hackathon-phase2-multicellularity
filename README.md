(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.


[Project-specific text] The goal of this specific project was to investigate multicellularity across the tree of life. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 24 hours on a computational cluster with 24 cpus, 150GB RAM. if needed, the user will be able to switch the node to a node with a A100 GPU in the future, but at the moment do not rely on it.


Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.


The rest of this repository is entirely AI generated. It includes code for implementation of the predictive model and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

<!-- PI note (AI): the human text above is copied verbatim from PROMPT.md and is left unchanged. Actual effort: round 1 ran from
2026-10-28 12:21 (START.txt) for ~13.6 h on 24 CPUs/150 GB; round 2 ran from 2026-10-29 15:07 (START_ROUND2.txt) for ~17.5 h on
48 CPUs/500 GB and one A100 GPU, after an independent review. The text mentions a "predictive model" and "benchmarking against other
methods", which this project does not have. -->

## Start here

- **Final technical report (both rounds):** [`report/REPORT.md`](report/REPORT.md) (self-contained HTML: `report/report.html`;
  supplement: [`report/REPORT_SUPPLEMENT.md`](report/REPORT_SUPPLEMENT.md)). Round-1 report kept as `report/REPORT_round1.md`.
- **Independent senior review of round 1** that shaped round 2: [`reviews/round2_department_head_review.txt`](reviews/round2_department_head_review.txt)
  (human-readable page: `reviews/review_round1.html`); review of the final report: `reviews/round2_final_report_review.txt`.
- **Research record:** [`RESEARCH_TREE.md`](RESEARCH_TREE.md) (milestones, decisions, every result), [`LESSONS.md`](LESSONS.md),
  pre-registrations [`PREREG.md`](PREREG.md) (round 1, tag `prereg-v1`) and [`DESIGN_ROUND2.md`](DESIGN_ROUND2.md) §6 (round 2, tags
  `prereg-r2`, `prereg-r2-p1`); round-1 checkpoint tag `phase1-final`.

## Round 2 in brief (2026-10-29/30)

An independent reviewer judged that round 1 had asked a question it could not answer (seven very different origins, 13–35% power) and
found three overstated statements (corrected in `results/round1_errata/ERRATA.txt`, no conclusion changed). Round 2 therefore studied
single origins densely and the step from simple to complex multicellularity, with three pre-registered tests.

| Question | Result | Where |
|---|---|---|
| Do complex multicellular genomes carry more **cell-death machinery** than expected for their proteome size? | **Supported** (1.21-fold [1.04, 1.42], one-sided p = 0.0047, Holm 0.014; independently re-derived). The excess is mostly lineage-specific immune / self–non-self receptor machinery (NB-ARC–LRR in land plants, HET/NACHT in fungi, Death/NACHT/CARD in animals); weaker without the ten largest proteomes. | `results/r2/w4`, `results/r2/w4_celldeath` |
| Are **linear-motif-binding domains** (ELM) more abundant in multicellular species? | No (0.98-fold, p = 0.72). Cell-type transcription-factor families: no (exploratory). | `results/r2/w4` |
| Do the four unseen **simple→complex** events carry more new gene families than controls? | Not significant (p = 0.0504 as specified; Holm 0.10). Large per-event ratios are a tip-versus-internal-branch artefact; internal-only estimate 2.70× [1.29, 5.63], exploratory. | `results/r2/panels` |
| Do two orthology methods (Ancestra vs **EukaryOMA**) agree on when stem families arose? | Poorly: 29% of Ancestra's stem gains and <8% of duplications confirmed; per-stem counts are method-dependent. | `results/eukaryoma_xval` |
| Does **protein structure** (ProstT5 3Di on GPU, 0.655 agreement with AlphaFold) help Ancestra at deep nodes, or reveal structural novelty or constraint shifts? | No: structure mode failed its pre-registered rule; the apparent structural novelty was explained by short/disordered proteins; structure-vs-sequence divergence inconclusive (negative control failed). | `results/r2/w5`, `results/r2/w6`, verifications in `results/r2/verification/` |
| Do stem duplications favour subunits of ancient protein complexes (Cox et al. co-fractionation map)? | No consistent preference. | `results/r2/w7` |

Eight dense per-origin Ancestra runs (Holozoa, Fungi, Streptophyta, Volvocales, red algae, brown algae, Amoebozoa, cyanobacteria;
18–40 genome-based proteomes each) are in `runs/r2/` (outputs not committed) with metadata in `data/dense/`.

## About this repository (round 1)

**Question.** Multicellularity evolved many times independently. Do the independent origins share genomic changes — new gene families,
expansions of transcription-factor or adhesion/signalling families, loss of unicellular functions such as the flagellum, or a shift in
the rate of protein evolution? We reconstructed ancestral gene content with **Ancestra** (a recently developed, AI-built, scalable
hierarchical-orthology method; https://github.com/evocellnet/ancestra-mvp-ai, frozen tag `freeze-v1`) for 187 eukaryotic (2.97 M
proteins) and 82 prokaryotic (0.41 M) proteomes and compared the branches on which multicellularity arose with matched unicellular branches.

**Main results (all AI-generated, not yet verified by a human; details and caveats in the report).**
- None of five pre-registered tests was significant after Holm correction across seven independent eukaryotic origins (animals, fungal
  hyphae, dictyostelids, streptophytes, red algae, the brown-algal lineage, volvocine algae): novel genes p = 0.11, transcription-factor
  duplications p = 0.36, adhesion/signalling duplications p = 0.77, flagellar loss p = 0.020 (Holm 0.10), protein-evolution rate p = 0.27.
- The rate test excluded a ≥ 1.5-fold genome-wide rate shift (pooled ratio 1.10, 95% CI 0.96–1.26; the interval is optimistic because all
  seven transitions share the same 10 control branches, and a sensitivity analysis gives 1.16, 1.01–1.34 — still below 1.5).
- Ancestra recovered 11 of 16 pre-registered positive controls; its failures show that deep-node gene origins and duplication timing are
  unreliable for large multi-domain families, and simulations show gene gains are inflated 2.1–2.6-fold (in one calibration draw on our
  tree, more on transition branches than on controls).
- Exploratory signals are weak: more new domain combinations on transition branches (fragile to proteome completeness); a lower
  share of pre-existing (co-opted) domain families in animals and streptophytes than in comparable unicellular clades (Pfam-based only);
  an apparent excess of convergent domain loss disappeared under post-hoc checks.

**Read first:** [`report/REPORT.md`](report/REPORT.md) (the full report, with figures in `figures/`); a self-contained HTML rendering is
`report/report.html` (built by `scripts/report/build_html.py`).

## How the project was run
A PI agent planned, delegated to subagents, reviewed and decided; the process record is part of the result:
| file | content |
|---|---|
| `PROMPT.md` | the task, verbatim |
| `RESEARCH_TREE.md` | status of every analysis, milestone summaries and the decision log (D0–D8) |
| `LITERATURE.md` | literature review (≈9,000 words, tool-verified references, hypotheses H1–H16) |
| `DESIGN.md`, `reviews/` | design (v3.1) and three rounds of independent agent review |
| `PREREG.md`, `PREREG_HASHES.tsv` | pre-registration frozen before any eukaryote result (git tag `prereg-v1`), with dated deviations |
| `LESSONS.md`, `agents/retrospective_M3.md` | process lessons and a fresh-agent retrospective |
| `DATA.md` | every data source with version, date, checksum and licence |
| `results/verification/` | independent re-derivation of the headline numbers (all checks PASS) |

## Repository layout
| path | content |
|---|---|
| `mcell/` | analysis library: HOG → per-branch events (`hogs`), branch statistics and nulls (`branchstats`), gain-evidence ladder (`gainladder`), rates (`rates`, `phylo`, `krkc`), analyses (`analysis`, `explore`, `controls`) |
| `scripts/` | all pipelines: data build (`data/`), annotation (`annotate/`), Ancestra launcher, gain ladder, phylogenomics (`phylo/`), primary and exploratory analyses (`analysis/`), robustness (`robustness/`), structure (`structure/`), figures (`report/`), pre-registration (`prereg/`) |
| `tests/` | 126 tests (toy data with known answers, planted effects, null calibration); `python -m pytest -q -m "not slow" tests` |
| `data/` | small metadata only: species tables, guide trees, transitions, pre-registered control/toolkit lists (no sequences) |
| `results/` | all result tables and summaries (primary, exploratory, robustness, verification, Ancestra run metadata) |
| `figures/` | Figures 1–7 and S1 (PNG + SVG) with `captions.md` |
| `envs/` | conda environment specifications and lock files (round 2 adds `gpu`: torch 2.6 + transformers for ProstT5) |
| `data/dense/` | round 2: eight dense per-origin panels (species tables with uni/simple/complex and D/T state coding, guide trees, events) |
| `data/r2_sets/`, `data/r2_prereg/` | round 2 pre-registered family sets (ELM motif-binding, cell death, cell-type TFs, transposon mask) and hashes |
| `results/r2/` | round 2 results: `panels/` (dense runs, novelty, simple→complex), `w4*/` (family sets, cell death), `w5*/` (structure), `w6/`, `w7/`, `verification/` (independent re-derivations) |
| `results/eukaryoma_xval/`, `results/round1_errata/` | EukaryOMA cross-check; round-1 errata |
| `scripts/gpu/` | ProstT5 3Di prediction on GPU (encoder + CNN; windowing) |
| `scripts/r2/`, `scripts/dense/`, `scripts/eukaryoma/`, `scripts/errata/`, `scripts/report/` | round-2 pipelines |
| `report/` | final report (`REPORT.md`, `report.html`, supplement, `figures_v2/`), claims table, citation and numbers audits |
| `reviews/` | design reviews (both rounds), department-head reviews |

## Reproducing
1. Environments: `envs/ancestra.lock.yml` (Ancestra + pipeline; numpy 1.26.4 / scipy 1.12.0 pinned), plus `annot` (pyhmmer, metapredict)
   and `busco` envs (`envs/*.pip-freeze.txt`). Ancestra: clone https://github.com/evocellnet/ancestra-mvp-ai (commit ffb0bbe) into `external/`.
2. Data: `bash scripts/data/run_all.sh` downloads and preprocesses all proteomes (sources and checksums in `DATA.md`); Pfam-A 38.2 and
   BUSCO odb10 as in `DATA.md`.
3. Ancestra: `scripts/run_ancestra.sh data/proteomes/euk data/guide_tree_euk.nwk runs/euk_fast 12 fast` (≈ 9.3 h, peak 84 GB on
   24 CPUs; see `results/ancestra_runs/`), likewise for `prok`.
4. Analyses: `results/README_run_primary.md` and `results/README_run_exploratory.md` give the exact commands; every analysis script first
   verifies the pre-registration hashes (`python scripts/prereg/check_prereg.py`). Seeds: 11 (default), 23, 37, 41, 53.
5. Round 2: panels `python scripts/dense/build_panel.py` (see `data/dense/STATUS.txt`); Ancestra per panel
   `scripts/run_ancestra_patched.sh data/dense/<R>/proteomes data/dense/<R>/guide_tree.nwk runs/r2/<R> 16 fast 150` (15–46 min each);
   ProstT5 3Di `scripts/gpu/run_prostt5_all.sh` (A100; ~41k residues/s); analyses `scripts/r2/*`; figures `scripts/r2/figures/make_all.sh`;
   report HTML `conda_envs/annot/bin/python scripts/report/render_html.py report/REPORT.md report/report.html`.
Raw proteomes, Ancestra outputs (≈ 8 GB) and AlphaFold models are not committed (size; third-party licences); they are regenerated by the steps above.
