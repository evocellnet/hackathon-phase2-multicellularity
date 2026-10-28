# Review round 1: data, algorithms and compute (DESIGN.md v1)

Evidence: Ancestra docs and code, `results/scaling_large/summary.md`, `data/species_final.tsv` (196 taxa, 3.76 M raw proteins, 0.94 M from 36 transcriptomes), `results/validation_hogparser/summary.md`, and network probes of AlphaFold DB (AFDB).

**Verdict.** The design is sound in outline, but raw HOG gains cannot carry the gain-based tests (H1/H2, H5, H8, H14, Q6). In simulation, gains are inflated **2.1–2.6×**: 6–11 % of genes are given a too-young origin, mostly on long branches and on the children of the root. Losses and duplications are usable (r ≥ 0.98), except under whole-genome duplication (WGD), where internal losses are 1.5–1.8× too high. Memory is a hard risk for the full run.

## 1. Biases of Ancestra at eukaryote-wide depth

- **Split HOGs produce spurious gains at transition nodes.** A failed merge at P (last unicellular ancestor) yields two gains (M and sister S) plus a loss at M. Transition stems are long branches, so the artefact is confounded with transition status. The log(branch length) covariate fixes this only if the artefact is linear in it.
- **The S0 family-size cap is an absolute number.** `s2_gmax = 7,500` genes was 3 % of N on OrthoBench but is 0.22 % of N at 3.4 M proteins. Kinase, WD40, LRR, ankyrin, homeobox/C2H2 and ABC superfamilies (15–70 k genes each) will be Leiden-split, and these are the toolkit categories of H4–H6. Splitting creates spurious deep gains and hides expansions. Keep the frozen config and diagnose with `families_split_flags.tsv`.
- **Low recall** (TreeFam-A 0.58 vs 0.74 for OrthoFinder) fragments HOGs. Losses on single-species M branches mostly reflect incompleteness; WGD lineages (land plants, ciliates) need independent loss checks (H12).

## 2. Mitigations for gains (pre-register before seeing results)

Classify each HOG gained on a branch into M with an **evidence ladder**:

1. **Family-novel (free).** Re-join the Leiden-split parts to their source component (`split_from_component`) and add the scoring-stage rescue joins. The gain counts as family-novel if that component has no gene outside M. This test is conservative, because S0 families are superfamilies.
2. **Sequence-novel.** Search ≤ 3 representatives per gained HOG, from different species, with DIAMOND `--more-sensitive -e 1e-5 --query-cover 30` against one whole-dataset database. Then filter out hits to species inside M.
   - Scope: internal branches, transition terminals and matched controls (≤ 300 k queries); about 1 h at 24 threads.
3. **Architecture-novel (free once Pfam is done).** The majority ordered Pfam architecture of the HOG is absent outside M.
4. **Split-HOG re-merge (2–3 agent-hours; compute in minutes).** Use `scores.npz` to find an older HOG X that existed at P, is "lost" at M, and is the reciprocal best match of the gained HOG G by mean bitscore. Relabel G as the continuation of X, which removes one gain and one loss. Apply the same rule at the children of the root.
   - **Acceptance test (three existing simulations):** gain ratio ≤ 1.3×, loss/duplication r within 0.01 of current; otherwise report levels 1–2 only.

**Primary gain statistic** (H1/H2/H8): sequence-novel gains; raw HOG gains are an upper bound.

**Calibration (about 1 h).** Simulate on the real guide topology (pruned to about 60 taxa, concatenation branch lengths), and model spurious gains by branch length and depth as a GLM offset or check. Add a `root_child` indicator or drop root children from the controls.

**Orthology-free view.** Run Wagner parsimony per domain on the species × Pfam copy-number matrix (gain penalty 2, as in Count; seconds to run). For **expansions (H4/H5)**, make Pfam counts primary and HOG duplications secondary.

## 3. Transcriptomes

**Redundancy.** *Codosiga* (61 k), *Telonema* (59 k) and *Diphylleia* (41 k) are redundant. `cd-hit -c 0.95` keeps contained fragments, so use `cd-hit -c 0.95 -G 0 -aS 0.9 -n 5`. Flag sets with BUSCO duplicated (D) or fragmented (F) > 15 %; fragments become false in-paralogs.

**Contamination.** The direction of the bias matters:
- Prey contamination in a *unicellular sister* (the predators *Pigoraptor*, *Syssomonas*, *Rhodelphis*, *Telonema*) creates false presence at P. That inflates co-option (Q3) and deflates gains.
- Contamination in M creates false gains and false horizontal gene transfer (H14).

**Screen (about 1 h).** Search the transcriptome proteins with DIAMOND `--fast` against the genome proteomes plus the prokaryote set. Flag a protein if its best hit lies in another supergroup or domain at ≥ 80 % identity (≥ 70 % for prokaryotes), with a bitscore margin ≥ 50 over the best hit in its own group. Prioritise by LukProt's OMArk contamination %; flag rather than delete.

**Masking in the parser (no Ancestra rerun needed):**
- Presence that rests only on transcriptome taxa needs ≥ 2 species.
- Losses are recomputed with transcriptome taxa and BUSCO < 60 leaves set to "unknown", because transcriptomes miss stage-specific genes.
- *Guttulinopsis* counts for presence only.

The sisters of Metazoa (5/8), Dictyostelia (4/4), *Copromyxa* and *Guttulinopsis* are mostly transcriptomes.

## 4. Guide tree

- **Canonicalise Newick child order.** Ancestra's output depends on it (DISCLOSURES #11).
- **Write a re-merge driver.** The input-stage key contains the tree md5, so any tree change reruns S0 and S2 (about 5 h). `scripts/remerge.py` should load the cached `genes.tsv` and `scores.npz` and call `run_merge`. It must be byte-identical on the main tree. Each re-merge then costs about 0.3 h (`fast`).
- **Variants to run:**
  - Root of the eukaryotes: polytomy of Amorphea / Diaphoretickes / Discoba / Metamonada, versus a Discoba–Metamonada root.
  - Contested placements: Pluriformea/Ichthyosporea, *Fonticula*/*Parvularia*, *Copromyxa*, and Telonemia/Haptista/Cryptista.
  - A reversed-child-order run, to give a noise floor.
- **Flag fragile transitions.** In `transitions.tsv`, flag transitions whose P is a polytomy or whose sister group is contested. Call a result tree-sensitive if its sign flips across variants.
- Keep polytomies for Ancestra (k-way merge); resolve them by constrained IQ-TREE search only for RER and branch lengths.

## 5. Compute (24 CPU, 150 GB)

| step | estimate | basis |
|---|---|---|
| Ancestra `fast`, 3.4 M proteins | **5–6.5 h at 12 threads; peak RSS 120–130 GB** | n500: 3.32 h / 68.8 GB; n1000: 10.6 h / 221 GB; local exponents 1.56 (wall), 1.57 (RSS) |
| `accurate` merge from copied `fast` stages | +3–4.5 h at 12 threads | n250: merge 0.87 h (accurate) against 0.05 h (fast) |
| Pfam | 132 CPU-h, i.e. 11 h at 12 threads | measured |
| BUSCO / metapredict | about 5 CPU-h / 2.4 h at 4 threads | |
| Gain DIAMOND + contamination screen | about 2 h at 24 threads | |
| Concatenation (BUSCO single-copy markers, about 150 loci, about 50 k sites) | 1–2 h at 8 threads, about 5 GB | LG+F+G4 with constrained search; C20 needs about 5× the RAM, so skip it |
| RER, ≤ 3,000 clade-level + 300 global genes | 15–25 CPU-h | 20 CPU-s per 40-taxon gene, 1–2 CPU-min per 196-taxon gene; MAFFT, trimAl, IQ-TREE `-te` LG+G4 |

**Memory is binding.** DESIGN's estimate of "10–15 h accurate" ignores memory. The plan:
- Measure peak RSS on the 40-proteome test and fit the scaling exponent.
- If the predicted peak exceeds 110 GB, trim the optional animals and fungi to ≤ 3.0 M proteins (about 95 GB predicted).
- During the scoring stage (S2), run only Pfam alongside Ancestra, and log RSS every 60 s.

**Order:** freeze proteomes (any md5 change reruns all stages); `fast` (12 threads) ‖ Pfam (12 threads); copy `stages/` and run the `accurate` merge; compare presets on the **full** data (per-branch Spearman, transition-effect sign agreement). The 40-proteome run is then only a timing and memory test.

**Checkpoints:** Ancestra stages; per-species Pfam/BUSCO/metapredict and per-gene RER outputs.

## 6. Structure (§3.6): feasible only if narrowed

**Access (tested).**
- The AFDB Google Cloud bucket lists, but anonymous download returns **401** (and *Capsaspora*'s tar is 2.0 GB).
- EBI per-accession files work (`alphafold.ebi.ac.uk/files/AF-<acc>-F1-model_v6.cif`, about 80 KB each). Match the API's sequence md5 to ours and drop mismatches.

**Scope.** Only transitions with UniProt species on both sides qualify: Metazoa, choanoflagellate colonies, *Acrasis*, the four fungal transitions, Volvocales, Rhodophyta and Phaeophyceae. Dictyostelia, Streptophyta, *Copromyxa* and *Guttulinopsis* are excluded.
- **Queries:** sequence-novel gains (§2, level 2) with mean pLDDT ≥ 70, about 200–2,000 per transition.
- **Targets:** the full proteomes of 1–2 sister species.

**Cost.** About 150–200 k CIF files (12–16 GB), 1.5–2 h with 8 parallel downloads; Foldseek `easy-search --exhaustive-search` takes minutes. Homolog = E < 1e-3 and TM-score ≥ 0.5. Controls: sequence-novel gains on matched control branches (negative); genes present at P (positive, about 100 % expected). Keep "structural change" descriptive.

## Ranked changes

**Must**
1. Gain-evidence ladder, with sequence-novel gains as the primary statistic; re-merging of split HOGs gated on the simulator acceptance test (ratio ≤ 1.3×).
2. Simulation calibration on the real tree, plus a root-child covariate.
3. Memory gate from the 40-proteome test (if > 110 GB is predicted, cut to ≤ 3.0 M proteins). Run `fast`, then the `accurate` merge from its checkpoints, and compare them on the full data.
4. Transcriptomes: CD-HIT with `-G 0 -aS 0.9`, the BUSCO D/F gate, the contamination screen, and masking rules in the parser.
5. Canonical child order, `remerge.py`, and 2–3 tree variants; flag transitions with polytomous or contested P.

**Should**
6. Pfam counts as the primary measure of expansion; report the share of genes in Leiden-split families.
7. Concatenation of BUSCO markers to get a resolved tree for RER (≤ 3,300 genes).
8. The narrowed structure plan, using EBI per-accession files and md5 checks.
9. Check H12 losses with curated flagellar and phagotrophy HMMs in genome-only species.

**Could**
10. ProstT5 3Di for non-UniProt sisters (about 1 h per 20 k proteins on CPU).
11. LG+C20-PMSF branch lengths.
