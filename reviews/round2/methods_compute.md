# Round 2: methods implementability and compute (DESIGN v2)

**Evidence.**
- `data/proteomes/euk`: 196 files, **3,295,587 proteins**. BUSCO 76/196 done. No guide tree or `transitions.tsv` yet.
- Ancestra code and `scaling_large`; live DIAMOND/AFDB probes.

**Verdict.** All of v2 can be implemented. Four details would silently bias results or break the schedule.

## 1. Memory gate
- 3.30 M proteins exceeds the 3.0 M cap.
- Scaling from n500→n1000 (RSS exponent 1.57) predicts a **≈125 GB** peak. At 3.0 M: ≈107 GB and about 5.5 h at 12 threads.
- **Cut to ≤3.0 M before the freeze.** Candidates: optional background taxa (Tricvag, Tetrstr, Symbmic, Emilhux, Chrovel, Anaefla: 231 k).
- A 40-proteome test (0.6 M) says nothing about the 3 M regime. Use an RSS watchdog that kills at 145 GB instead.

## 2. `accurate` from `fast` checkpoints (verified)
- The `scores` key excludes the preset and the thread count; only `merge` and `output` rerun.
- The `input` key hashes **every non-dot file** in the proteome directory plus the tree md5. A stray file reruns S0 and S2 (about 5 h). So:
  1. `chmod -R a-w data/proteomes/euk`
  2. `cp -al runs/euk_fast/stages/{input,families,scores} runs/euk_acc/stages/`
  3. Run `--preset accurate`.
  4. Abort unless `scores: up to date`. Runtime 3–4.5 h.
- **Pre-register `fast` as the primary preset.**
- S0 and S2 do not use the tree. A ~40-line `remerge.py` (cached genes and scores → `run_merge` on a variant tree) re-merges each tree variant in about 0.3 h. Must be byte-identical on the main tree.

## 3. Novelty search: one DIAMOND search, filtered by clade afterwards
- **Pitfall.** DIAMOND's default `-k 25` lets in-clade paralogues fill a query's hit list. Gains then look falsely novel. Use `-k 0`.
- **Database.** Eukaryotes plus the 82 prokaryotes, with a fixed `--dbsize`. Pre-register this composition. Flag prokaryote-only hits as HGT candidates.
- **Pass 1.** One representative per gained HOG on internal branches (about 150–250 k queries):
  `--very-sensitive -k 0 --max-hsps 1 -e 1e-5 --query-cover 30 -b4 -c1 -f 6 qseqid sseqid pident evalue bitscore qcovhsp`
- **Reducer** (`mcell/novelty.py`): per query, apply the clade species mask; keep the best out-of-clade hit, out-of-clade species count, a transcriptome-only flag and the top-50 out-of-clade hits.
- **Pass 2.** Two more representatives, only for HOGs still without a hit.
- **Ultra-sensitive mode.** Rerun 2 k pass-2 "novel" HOGs with `--ultra-sensitive`. Switch pass 2 to ultra mode if more than 5% gain hits.
- **Cost.** A probe of 38 queries took 91 s, almost all fixed database-scan overhead. So **pilot 20 k queries first**. Estimate: 1–3 h at 20 threads, 25 GB. Do not overlap with Ancestra S2.
- **Sensitivity check.** Ignore support that comes only from transcriptomes at ≥90% identity.

## 4. Split re-merge cannot use `scores.npz`
`scores.npz` holds only pairs within one S0 family, and Leiden or connected-component splits put G and X in different families. Use the reducer's top hits instead:
- X = the HOG at P (`HogResult.members(P)`) with the highest summed bitscore and no descendant at c.
- Accept the re-merge only if G is X's best-scoring gain at c.
- Run identical code on the simulations for the ≤1.3× gate.

## 5. BUSCO concatenated tree
- Markers: single-copy BUSCOs present in ≥70% of taxa (about 150–200 of 255).
- Commands: `mafft --auto`, then `trimal -automated1`, then `iqtree3 -m LG+F+G4 -g constraint.nwk -T 6`. Expect 2–4 h and 15 GB.

## 6. Gene trees for the relative-rate test
- Per transition: taxa M ∪ S ∪ 2–4 outgroups with genomes.
- Loci: ≤300, from BUSCO plus HOGs single-copy at parent(P).
- Commands: `mafft --localpair`, then `trimal`, then `iqtree3 -te pruned.nwk -m LG+G4 -T 1`.
- Statistic: log(mean P→tip distance in M ÷ mean in S).
- Drop loci with branches >2 substitutions/site.
- Cost about 20 CPU-h. Cap gene-level RER at **1,000 HOGs** (25 CPU-h).

## 7. Pfam and metapredict
- **Now:**
  1. The ~150 toolkit, control and H12 models, with `--cut_ga`. This takes <1 CPU-h and unblocks H4, H5 and H12.
  2. Full Pfam-A (132 CPU-h) at 8–10 threads, checkpointed per species, done by t ≈ 15–17 h. Resolve clan overlaps.
- metapredict (2 CPU-h) after BUSCO.
- v2's plan to run `accurate`, DIAMOND and gene trees in parallel at 10–15 h, alongside Pfam, oversubscribes the 24 cores. Run DIAMOND at 10–13 h, then `accurate` at 13–17 h.

## 8. AFDB
- **Requests with curl's default User-Agent get HTTP 403;** with a User-Agent header (`-A 'mcell-phase2/1.0'`) they succeed.
- Keep files gzipped: 151 KB against 720 KB raw for a 761-aa protein. Expect about 200 k files, 20 GB and 1–1.5 h.
- The API's `sequenceChecksum` is the md5 of our sequence (verified). Use it with `globalMetricValue` ≥ 70 as a filter before downloading.
- Only the 113 UniProt-sourced species qualify.

## Ranked changes
**Must**
1. ≤3.0 M proteins, plus a 145 GB watchdog.
2. Novelty search with `-k 0`, a pre-registered euk+prok database, a clade-mask reducer and two passes.
3. Re-merge from the novelty hits.
4. `fast` as primary; `accurate` from hard-linked stages with a read-only proteome directory.
5. Pfam subset pass plus full Pfam, both started now.

**Should**
6. AFDB fetching with a User-Agent header, gzipped files and the API filter.
7. Clade-restricted gene trees; RER ≤1,000 HOGs.
8. BUSCO concatenation as in §5, and `remerge.py`.
9. Very- vs ultra-sensitive chosen from the calibration.

**Could**
10. Novelty on terminal branches for exploratory tips.
11. ProstT5 3Di for sisters not in UniProt.
