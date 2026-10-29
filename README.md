(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.


[Project-specific text] The goal of this specific project was to investigate multicellularity across the tree of life. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 24 hours on a computational cluster with 24 cpus, 150GB RAM. if needed, the user will be able to switch the node to a node with a A100 GPU in the future, but at the moment do not rely on it.


Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.


The rest of this repository is entirely AI generated. It includes code for implementation of the predictive model and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

<!-- PI note (AI): the human text above is copied verbatim from PROMPT.md. Actual elapsed time: START 2026-10-28 12:21 (START.txt);
the report was finalised on 2026-10-29 (see the last commit time), i.e. within the 24 h stated above if the final commit is before
2026-10-29 12:21. The text mentions a "predictive model" and "benchmarking against other methods", which this project does not have;
it is left unchanged because it is human-written. -->

## About this repository

**Question.** Multicellularity evolved many times independently. Do the independent origins share genomic changes — new gene families,
expansions of transcription-factor or adhesion/signalling families, loss of unicellular functions such as the flagellum, or a shift in
the rate of protein evolution? We reconstructed ancestral gene content with **Ancestra** (a recently developed, AI-built, scalable
hierarchical-orthology method; https://github.com/evocellnet/ancestra-mvp-ai, frozen tag `freeze-v1`) for 187 eukaryotic (2.97 M
proteins) and 82 prokaryotic (0.41 M) proteomes and compared the branches on which multicellularity arose with matched unicellular branches.

**Main results (all AI-generated, not yet verified by a human; details and caveats in the report).**
- None of five pre-registered tests was significant after Holm correction across seven independent eukaryotic origins (animals, fungal
  hyphae, dictyostelids, streptophytes, red algae, the brown-algal lineage, volvocine algae): novel genes p = 0.11, transcription-factor
  duplications p = 0.36, adhesion/signalling duplications p = 0.77, flagellar loss p = 0.020 (Holm 0.10), protein-evolution rate p = 0.27.
- The rate test was sensitive enough to exclude a ≥ 1.5-fold genome-wide rate shift (pooled ratio 1.10, 95% CI 0.96–1.26).
- Ancestra recovered 11 of 16 pre-registered positive controls; its failures show that deep-node gene origins and duplication timing are
  unreliable for large multi-domain families, and simulations show gene gains are inflated (more on transition branches).
- Exploratory signals that survived skeptical checks are weak: more new domain combinations on transition branches (fragile to proteome
  completeness) and lower co-option in animals and streptophytes (Pfam-based only).

**Read first:** [`report/REPORT.md`](report/REPORT.md) (the full report, with figures in `figures/`).

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
| `envs/` | conda environment specifications and lock files |

## Reproducing
1. Environments: `envs/ancestra.lock.yml` (Ancestra + pipeline; numpy 1.26.4 / scipy 1.12.0 pinned), plus `annot` (pyhmmer, metapredict)
   and `busco` envs (`envs/*.pip-freeze.txt`). Ancestra: clone https://github.com/evocellnet/ancestra-mvp-ai (commit ffb0bbe) into `external/`.
2. Data: `bash scripts/data/run_all.sh` downloads and preprocesses all proteomes (sources and checksums in `DATA.md`); Pfam-A 38.2 and
   BUSCO odb10 as in `DATA.md`.
3. Ancestra: `scripts/run_ancestra.sh data/proteomes/euk data/guide_tree_euk.nwk runs/euk_fast 12 fast` (≈ 9.3 h, peak 84 GB on
   24 CPUs; see `results/ancestra_runs/`), likewise for `prok`.
4. Analyses: `results/README_run_primary.md` and `results/README_run_exploratory.md` give the exact commands; every analysis script first
   verifies the pre-registration hashes (`python scripts/prereg/check_prereg.py`). Seeds: 11 (default), 23, 37, 41, 53.
Raw proteomes, Ancestra outputs (≈ 8 GB) and AlphaFold models are not committed (size; third-party licences); they are regenerated by the steps above.
