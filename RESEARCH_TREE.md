# Research tree — genomic transitions to multicellularity (Ancestra)

Legend: `[ ]` planned · `[~]` in progress · `[x]` done · `[-]` abandoned/negative · `[!]` key finding

START: 2026-10-28T12:21:22+01:00 (see START.txt). Budget: 48 h allocation; target finish ≤ 36 h after start, leaving margin.

## ROUND 2 (started 2026-10-29 15:07; 48 CPU, 500 GB, A100; soft deadline ~30 h; START_ROUND2.txt)

### R2-M2 — 2026-10-29 22:03 (round-2 t ≈ 6.9 h): confirmatory family closed; one positive, two nulls
Holm over the three confirmatory tests: P3 (cell death) rejected (adjusted 0.014); P1 (more sequence-novel families at four unseen
simple→complex events) not significant (Stouffer p = 0.072, adjusted 0.144); P2 (motif-binding domains) null. P1 was underpowered
by design (3–4 events, 18–41% power) and, in retrospect, should have been labelled estimation (LESSONS rule 5). Its largest per-event
ratios (Volvox carteri, Astrephomene ≈ 13–14×) are a tip-versus-internal-branch artefact: reversal tips score ≈ 21× by the same
comparison, and against tip controls the effect is ≈ 0.6–1.0. A cleaner secondary estimate, internal simple→complex stems only, is
2.70× [1.29, 5.63] (k = 7) versus 0.84 for internal origins — the one gene-content lead that survives in round 2, exploratory.

All 8 dense panels ran in 15–60 min each (vs 10 h for the round-1 eukaryote run). Positive controls recover 14/18 on panels, including
the deep-node T-box control (C07) that failed in round 1; cadherin and integrin controls still fail. The LECA complex map (W7) adds
nothing. Running: independent verification of P1/panels, W5 (Ancestra 3Di structure mode on R1; structural novelty), W6 (3Di vs
amino-acid divergence, exploratory); GPU 3Di 93/322 species. Next: report drafting starts after verification (target t ≈ 10 h).

### R2-M1 — 2026-10-29 18:18 (round-2 t ≈ 3.2 h): first confirmatory result; cross-method check is sobering
Confirmatory R2-P3 is supported: genomes of multicellular species carry more cell-death machinery than expected for their proteome
size (1.21-fold [1.04, 1.42], raw one-sided p = 0.0047 < Holm 0.0167 regardless of P1; mainly complex multicellularity; sign
holds leaving out any origin), independently re-derived with a symmetric inclusion rule and a random-set check of the pipeline
(2.5% false positives); caveat: p = 0.025 without the ten largest proteomes, and the proportion-of-proteome version is not
significant. R2-P2 (motif-binding domains) is clearly null (0.98-fold); cell-type TF families (exploratory) null. Pfam-wide
expansions do not recur across origins beyond chance; contractions do (incl. bacterial-like glycosidases likely reflecting
contamination/bacterivory in sister transcriptomes).

EukaryOMA cross-validation shows that two orthology methods agree poorly on where stem families arise (29% of Ancestra stem
gains confirmed; duplications <8%), so per-stem HOG counts are method-dependent; the sequence-level novelty statistic agrees in
direction. Errata confirm the reviewer's three points without changing any round-1 conclusion. Seven of nine panel runs are done
and verified (15–60 min each); R4 (3 of 4 confirmatory P1 events) and R7 running; ProstT5 3Di validated (0.655 vs AlphaFold) and
R1 3Di nearly complete for the structure-mode test. Next: P1 and per-event analyses, structure mode, report skeleton by t≈10 h.

### R2-M0 — 2026-10-29 15:41 (round-2 t ≈ 0.6 h): checkpoint, external review, new plan
Round 1 checkpointed (tag `phase1-final`). The round-2 instructions are saved verbatim (PROMPT_ROUND2.md). A department-head reviewer
(Fable agent) judged round 1 honest and methodologically valuable but aimed at a question it could not answer (7 heterogeneous stems,
13–35% power), and found three overstated report statements (depth SMD 2.06 not 1.77; power table used placeholder branch lengths;
5 of 11 recovered positive controls never use Ancestra events). The review is public as a page (reviews/review_round1.html;
https://claude.ai/artifact/NfrexbqRBFmS3A97VvTqcR). Scouts: EukaryOMA = OMA instance, 1,328 genomes, HOG OrthoXML downloaded,
122/187 of our species overlap; Cox et al. (Marcotte lab) CF-MS complex map is LECA-centred (no lineage-specific complexes), 16/31 species overlap.

GPU: Foldseek's GPU build fell back to CPU for ProstT5; an own encoder+CNN ProstT5 implementation runs at 48 proteins/s on the A100 and
agrees 97.1% per residue with Foldseek's ProstT5 (5,000 proteins); 3Di prediction for all 269 round-1 proteomes started (~17–20 h).
Plan: DESIGN_ROUND2.md v0 (W0 errata, W1 EukaryOMA cross-validation, W2 eight dense per-event Ancestra panels, W3 simple→complex,
W4 specific families incl. ELM motif-binding domains, cell death, cell-type TFs, W5/W6 structure, W7 complexes) under one review round
(statistics, evolutionary biology, compute/AI). Running: errata agent, EukaryOMA agent, dense-panel data agent, ProstT5.

### Round-2 tree
- [x] R2-0 Checkpoint (tag phase1-final), PROMPT_ROUND2.md, START_ROUND2.txt
- [x] R2-1 Department-head review → reviews/round2_department_head_review.txt; public page
- [x] R2-2 Scouting: EukaryOMA (docs/EUKARYOMA_SCOUTING.txt), Cox et al. CF-MS (docs/MARCOTTE_SCOUTING.txt), ELM interaction domains (200 Pfam)
- [x] R2-3 GPU ProstT5 pilot (Foldseek GPU fell back to CPU: 15 prot/s; HF autoregressive: too slow; encoder+CNN: 48 prot/s, 97.1% agreement)
- [~] R2-4 ProstT5 3Di for all proteomes (data/structure3di/all/); [x] validated: 0.655 per-residue identity vs AFDB-derived 3Di (n=4,981; chance 0.167; published 0.41–0.65), windowed 0.649; integrity PASS (results/r2/w5_3di_validation)
- [x] R2-5 Design round 2 v1 after 3 reviews; pre-registered (tag prereg-r2)
- [x] W4 (results/r2/w4): [!] P3 cell death SUPPORTED raw one-sided p=0.0047 (1.21x [1.04,1.42]; 148 genomes; LOO sign 7/7; mainly complex MC; independently re-derived PASS; symmetric rule p=0.0068, all 22 families p=0.0078; random matched sets 2.5% p<0.05; caveats: weaker without 10 largest proteomes p=0.025, proportion-of-proteome version p=0.14 → "more than expected for proteome size") · [-] P2 motif-binding null (0.98x, p=0.72) · [-] P4 TFs p=0.30 · scans: expansions not recurrent beyond chance (family-shuffle null, exploratory), contractions recurrent (p=0.001)
- [!] W1 EukaryOMA cross-validation (results/eukaryoma_xval): 125 shared species, 87% proteins mapped; methods agree poorly on WHERE stem families arise: 29% of Ancestra stem gains confirmed (controls 11%; depth-adjusted expectation 22%), reverse 9%; duplications <8% either way; disagreements ~8 nodes apart; Metazoa stem dups Ancestra 6,617 vs EukaryOMA 196 → per-stem HOG counts are method-dependent; descriptive H1/H2 T 0.63 (EukaryOMA) vs 0.72 (Ancestra), direction agrees 4/5. C07 passes in Ancestra on the shared tree (round-1 failure was level labelling); C03/C06 fail in both.
- [x] W0 errata (results/round1_errata/ERRATA.txt): SMD 2.06 confirmed (1.77 stale); power table unchanged with real branch lengths (≤0.025); Ancestra-event controls 7/10 (deep-node 0/3); H9 exact over 1,024 orientations p=0.31, ratio 1.06, no ≥1.5-fold shift in any orientation
- [x] W2 runs: all 8 panels + R4_curated done and verified (runs/r2/)
- [x] W2/W3 panel analyses (results/r2/panels): [-] confirmatory P1 NOT significant: Stouffer one-sided p=0.072 (Astrephomene 0.143, V. carteri 0.143, V. africanus 0.571, Fucophycidae 0.167); pooled RR 4.6 [0.58, 37]; no depth-matched controls (events 7–10 nodes deep vs controls 2–6). [!] Tip-vs-internal artefact: reversal tips score ~21× vs internal controls; against tip controls the tip events give 0.56–0.80 (P1 p=0.57). Secondary: internal-only simple→complex RR 2.70 [1.29, 5.63] (k=7) vs internal origins 0.84; positive controls 14/18 (C07 recovered on R1; C03/C06 fail; C10 fails on R3); domain recurrence at complex events = size effect (capped p=0.77); kinases in novel families at 5/7 events (p=0.023, q=0.06, suggestive). Final Holm (P1–P3): P3 rejected (0.014), P1 0.144, P2 n.s.
- [~] W5 structure mode + structural novelty · [ ] W6 3Di/AA divergence
- [-] W7 complexes (results/r2/w7): stem duplications show no consistent preference for LECA-complex subunits (OR 0.79–1.48, none significant after BH); co-retention slightly lower in multicellular species relative to null (paired p=0.008; exploratory, weak)
- [ ] verification · [ ] report update

## ROUND 1 (frozen at tag phase1-final)

## Status summary (latest first)

### M5 — 2026-10-29 01:55 (t ≈ 13.6 h): report final
The report (report/REPORT.md, ~8,900 words of prose, Figures 1–7 + S1) went through a citation audit (51 citations, 48 claim–citation
pairs: 39 supported, 8 partially → reworded, 0 unsupported, 1 uncheckable) and an external critical review (10 major problems, all
addressed — notably: the co-option index is equally high on control branches so it cannot favour co-option; the stem-specific gain bias
rests on one simulation draw; the rate-test interval is widened for shared controls to about 0.94–1.29; 'the only convergent loss is
the flagellum' softened). Revision log: report/revision_notes.txt. README.md carries the human disclaimer verbatim.

Bottom line: across seven independent eukaryotic origins, no shared genomic signature of multicellularity was detected with the
pre-registered tests; a ≥1.5-fold average shift in protein evolutionary rate is excluded; Ancestra's deep-node origins and duplication
placement are the main limitation for gene-content questions at this phylogenetic depth.

### M4 — 2026-10-29 01:18 (t ≈ 13.0 h): all pre-registered tests done; nothing significant after Holm
The eukaryote Ancestra run finished at t = 10.3 h (families 5.3 h, scoring 3.7 h, merge 10 min; verified). None of the five primary
tests is significant after Holm: H1/H2 sequence-novel gains p = 0.113 (Holm 0.45), H4 TF duplications p = 0.36, H5 adhesion/signalling
p = 0.77, H9 rate shift p = 0.27 (all Holm 0.81), H12 flagellar loss p = 0.020 (Holm 0.10; fungi and dictyostelids lose the whole
module as expected, but two stems give little power). Positive controls 11/16 — pre-registered criterion met, but the three
high-confidence failures expose real limits of HOG inference at deep animal stems (domains gained inside older gene lineages are
invisible; duplications placed one node too deep with compensating losses). An independent agent re-derived every headline number
from raw files with its own code (all PASS) and confirmed that the earlier stem-path bug did not touch H9.

Exploratory signals were attacked with post-hoc checks: convergent domain loss disappears (it is the fungi × Dictyostelia flagellum loss);
more new domain architectures on stems survives leave-one-out but not a completeness filter; lower co-option in animals and
streptophytes survives subsampling but is Pfam-only. Structure: ~1/4 of modelled sequence-novel gains are structurally old, equally
on stems and controls. Next: report (results/discussion drafting in progress), figure for exploratory results, final citation audit.

### M3 — 2026-10-28 20:09 (t ≈ 7.8 h): waiting on the eukaryote Ancestra run; primary H9 negative
The eukaryote Ancestra run is the bottleneck. Its family stage took 5.3 h because single-linkage similarity joined 73% of all
eukaryotic proteins into one component that Ancestra must cut into ≤7,500-gene families with Leiden (D8); within-family scoring is
~60% done (26 DIAMOND blocks, ETA ~22:00), then the merge. The eukaryote gain-ladder search starts automatically when Ancestra
finishes (est. 2.5–5 h), after which the primary gene-content tests run with one command. Full Pfam 119/187, disorder 112/187.
The prokaryote run is complete and served as the pipeline dry run (both prokaryotic positive controls recovered).

Results so far: the pre-registered rate test H9 is negative (p = 0.27; pooled rate ratio 1.10 [0.96, 1.26], i.e. no ≥1.5-fold
shift) with heterogeneous per-transition directions; exploratory gene-level RER and Kr/Kc found nothing. Before any eukaryote HOG
output, two method facts were established and turned into pre-committed rules: raw gains are inflated 2.1–2.6× (simulations), and on
our tree stems are gain-inflated 1.7–2.2× more than controls, so a "burst" claim now needs a pooled effect > 2.2 (Dev-4). A fresh-agent
retrospective led to pruned LESSONS and a precise PREREG clarification of what had been seen when each deviation was made.

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

- [x] 0. Setup: environments (ancestra pinned, annot, busco), taxdump, Pfam 38.2, BUSCO odb10
- [x] 1. Literature review → LITERATURE.md (16 hypotheses; 254 refs)
- [x] 2. Data scouting → data/species_candidates*.tsv, docs/DATA_SOURCES_SCOUTING.md
- [x] 3. Design → DESIGN.md v3.1 after 3 review rounds; PREREG frozen (tag prereg-v1); Dev-1–Dev-4
- [x] 4. Dataset: 187 euk (2.97 M proteins) + 82 prok (0.41 M); BUSCO; guide trees (62 edits); transitions.tsv
- [x] 5a. Method validation: HOG parser exact on true HOGs; Ancestra gains 2.1–2.6× inflated (sims); statistics engine calibrated (permute; parametric bootstrap for H1/H2); gain ladder (split repair rejected; novel precision 0.42–0.56); tree calibration (stem gain bias 1.7–2.2×)
- [x] 5b. Ancestra prokaryote run (fast; 1.4 h) — dry run; gain ladder; descriptive stats; P1/P2 controls recovered
- [x] 5c. Ancestra eukaryote run (fast): done 2026-10-28 22:36 (t=10.3 h); families 5.3 h, scores 3.7 h (peak 84 GB), merge 10 min; verify PASS; 2.22 M genes in HOGs, 750k singletons
- [~] 5d. Annotation: Pfam pre-registered subset done (269/269); full Pfam 119/187 euk; disorder 112/187 euk
- [~] 6. PRIMARY tests (H1/H2, H4, H5, H12 on euk; Holm with H9) + positive controls
  - [-] H1/H2 (sequence-novel gains, burst vs co-option): T=0.711, parametric-bootstrap p=0.113, Holm 0.45 — no burst detected (secondary nested complex set p=0.016, outside Holm)
  - [-] H4 (TF-family expansion, compositional): T=0.548, p=0.36, Holm 0.81; orthology-free version agrees
  - [-] H5 (adhesion/ECM+signalling expansion): T=0.41, p=0.77, Holm 0.81; orthology-free version flips sign (direction not robust)
  - [-] H12 (flagellar loss where expected): T=0.738, p=0.020, Holm 0.10 — not significant after correction
  - [x] Determinism: separate and combined runs give identical p-values
  - [!] Positive controls: FINAL 11/16 recovered (C10 GRAS recovered with gain ladder) → pre-registered count criterion met; high-confidence failures C03/C06/C07 explained (below). No code bugs; C03/C06/C07 failures are Ancestra inference artefacts (domain gained within an older HOG is invisible; mixed-domain HOGs; duplications placed one node too deep with compensating losses — 52% of Choanozoa-stem duplication products are animal-only); C08/C12 are expectation problems. Implication: HOG-level origins/duplication placement unreliable at deep animal stems → H4/H5 robustness (orthology-free) + post-hoc re-placement running.
  - [x] H9 (genome-wide rate shift): not detected — p = 0.27, pooled ratio 1.10 [0.96, 1.26] (no ≥1.5-fold shift); heterogeneous: faster in Dictyostelia/volvocines/red algae/fungi, slower in Metazoa/PX; seed-sensitive p 0.19–0.61
- [~] 7. Secondary: nested complex H1/H2 p=0.016 (outside Holm; Dev-4 rule applies); prokaryote descriptive done
- [ ] 8. EXPLORATORY
  - [x] gene-level RER (250 markers): 0 genes at q < 0.10; none consistent across ≥ 4 transitions
  - [x] Kr/Kc radical/conservative: no support for relaxed constraint (p = 0.60 / 0.16)
  - [-] Pfam-wide convergent expansion/gain: not above chance (k≥2 p=0.79; top 'convergent' domains are transposon domains)
  - [-] Pfam-wide convergent LOSS (k≥2: 421 vs 192, p=0.043) DISAPPEARS under post-hoc checks: entirely fungi×Dictyostelia shared flagellar losses (= H12 signal); without Dictyostelia p=0.60, genome-only p=0.77
  - [~] New domain architectures (H11) on stems > controls (T=0.80, p=0.005): survives LOO (≤0.017, incl. without Metazoa) and genome-only; weakens with BUSCO≥80% (p=0.16) → suggestive only
  - [~] Co-option index (Pfam): Metazoa and Streptophyta less co-opted than controls (u=0.024; 0.071 after subsampling M to 3–5 species); volvocine/red/PX ≈100% co-opted; HOG-based version shows no effect; Pfam annotation bias caveat
  - [-] Protein traits (length/disorder on shared single-copy HOGs): within-transition shifts common but not distinct from controls
  - [~] Reversals: odds ratios mixed (0.5–2.0), see results/exploratory/euk/reversals
  - [-] Structure (Q6): 85% of stem vs 79% of control seq-novel gains lack a usable AFDB model (OR 1.51, p=0.025); of modelled ones ~24% vs 21% structurally old (Swiss-Prot ref; OR 0.96, p=0.92) — no stem-specific structural novelty signal; underpowered
- [x] 9a. Independent verification: event counts, novel gains, Holm, T_obs, H9 (unaffected by bug), controls C01/C05/C10, PREREG timeline — all PASS (results/verification/)
- [x] 9b. Citation audit: 51 citations; 48 claim–citation pairs audited (39 supported, 8 partial → reworded, 0 unsupported)
- [x] 10. Report final: report/REPORT.md (+ external review, revision notes, references); Figures 1–7, S1

## Decision log
- D0 (t=0): Eukaryote and prokaryote runs kept separate (see M0).
- D1 (t=0.3h): Dataset build launched before design review completes; taxon changes from review are cheap before the Ancestra run.
- D2 (t=0.4h): "Eukaryoma" not found → proceed with UniProt + EukProt + LukProt + NCBI; question logged for the user.
- D3 (t=1.2h): Eukaryote Ancestra `fast` launched (12 threads, 140 GB watchdog) with the guide tree as built; if the round-3 audit changes the tree, only the merge stage is re-run (S0/S2 cached). Prokaryote run (4 threads) and Pfam (8 threads) in parallel.
- D4 (t=1.3h): Volvocine primary transition coded at the TGV stem (all sampled colonial species descend from it); single- vs two-origin interpretation is a reported caveat; Tetrabaena-tip and Gonium+Volvocaceae rows kept exploratory.
- D5 (t=1.3h): Eligibility rule refined before results: low-quality species (BUSCO C<50%) do not count toward the ≥2 per side, rather than disqualifying a transition; M must be a single clade. Result: 7 primary eukaryotic origins, 7 secondary (3 nested complex, 4 prokaryotic). Script: scripts/prereg/apply_eligibility.py.
- D6 (t=1.6h): Statistics engine calibration on the real tree (before results): literal relabelling null anti-conservative (14.5% at α=0.05) because |U|=30 and matched sets overlap → primary scheme `permute` (4.9%); H1/H2 statistic = NB-GLM residual fitted on controls (matching poor: depth SMD 1.77). Prokaryotes |U|=2 → descriptive only.
- D7 (t=1.7h): Round-3 audit F1/F2 adopted on phylogenetic grounds (Lindsey 2021; Kiss 2019): volvocine and fungal-hyphae stems become 2-edge path stems ending at Gonium_Volvocaceae and Amastigomycota. Consequence (not the reason): fungal stem becomes flagellum-loss-expected, so H12 has 2 loss-expected stems and stays primary.
- D8 (t=5.9h): Eukaryote Ancestra S0 stage-2 graph has ONE connected component of 2.16 M genes (73%; 473k stage-1 reps); stricter coverage does not dissolve it (min-cov 0.5: 1.48 M; 0.7: 0.96 M). The frozen leiden_split rescans all component edges per part and runs serially → many hours. Kept the original run going; launched in parallel an output-identical patched wrapper (parallel independent parts, per-part edge partitioning; identity verified on random graphs and the prokaryote run before use). Whichever finishes first is used; both must give identical families. Scientific note for the report: in eukaryote-wide data single-linkage collapses most of the proteome into one 'family', so Leiden cut points define family boundaries — the mechanism behind spurious HOG gains that the DIAMOND novelty ladder guards against.
- Bug found by an independent agent (t≈7h, before any eukaryote results): exploratory T_volvocales_gv shadowed the primary volvocine path stem in mcell.branchstats; fixed (priority order) — verify H9 unaffected during verification.
- D8 outcome (t=6.3h): the ORIGINAL frozen run finished the family stage first (18,944 s; peak 33.6 GB; 553,686 families incl. 468,618 singletons; max family 7,565; Σn² = 7.1e9). The patched run (verified byte-identical on prokaryotes) was stopped unused at 3% progress; its code is kept in scripts/ancestra_patch/ as a documented, tested speed-up for future runs.
- D9 (R2 t=0.5h): ProstT5 3Di computed for ALL proteomes on the otherwise idle GPU (~17-20 h), multicellular species first, rather than only targeted sets: it enables Ancestra 3Di rescue on any panel and complete structural novelty checks; cost is GPU time only.
- D10 (R2 t=0.9h): ProstT5 windowing (1,000 aa windows, 250 overlap) adopted: 34k vs 21k residues/s; identical for ≤1,000 aa (93% of proteins), 90.3% per-residue agreement with full-length prediction for 1,000–3,000 aa. GPU queue reordered by panel (Holozoa first, prokaryotes last) because the full set (~1.43 G residues) cannot finish in time with certainty.
- D11 (R2 t=0.7h, 15:48): Round-2 pre-registration frozen (tag prereg-r2): DESIGN_ROUND2 v1 + frozen sets (ELM motif-binding 115+4, cell death 22, cell-type TFs 26, transposon mask 223), sha256 in data/r2_prereg/PREREG_R2_HASHES.tsv. Correction: the v1 header says "16:15"; the actual freeze time is 15:48 (header written ahead of time; file left unedited because it is hashed). No round-2 test output existed at freeze. Review changes adopted: confirmatory family = P1–P3 (Holm over 3); P1 only unseen simple→complex events; panels use `fast` (comparability) with patched Leiden, caps and go/no-go at t=8 h; W4 launched on the 187 round-1 species (annotations complete).
- D12 (R2 16:02): P1 confirmatory events frozen (tag prereg-r2-p1): R4 Astrephomene, Volvox carteri, Volvox africanus/reticuliferus; R6 Fucophycidae stem; Pleodorina sensitivity. Panels R2 (Fungi, 40 proteomes, 0.40 M proteins) and R3 (Streptophyta, 21, 0.47 M) launched (`fast`, 16 threads each, patched Leiden wrapper, 150 GB watchdog); no separate pilot because panels are 6–7× smaller than the round-1 run and the wrapper parallelises giant components; the family stage timing serves as the pilot.
- D13 (R2 16:22, t=1.25h): Panel R2 (Fungi) Ancestra `fast` finished in 18 min (peak 5.4 GB) — panels are cheap, so all panels (incl. R7) will run; go/no-go deadline not binding. Launched R6, R1 (16 thr), R5, R8 (8 thr). W2/W3 panel-analysis agent started (reuses round-1 mcell code; novelty DB = panel + all 187/82 round-1 proteomes outside clade).
- D14 (R2 17:28, t=2.3h): Panels R1, R2, R3, R5, R6, R8 done and verified (verify_rc 0; each 15–60 min); R7 launched; R4 waits for de novo AUGUSTUS proteomes (Eudorina; cap 17:45). ProstT5 queue rebuilt: all 37 R1 proteomes first (24 missing incl. new panel species), then other new panel species (53 new in total), then the rest; driver re-pointed at panel sources.
- D15 (R2 17:33, t=2.4h): All 8 panels built (data/dense/STATUS.txt; 18–40 proteomes, 0.19–0.63 M proteins each). R4 includes two AUGUSTUS de novo proteomes (Yamagishiella, Eudorina; BUSCO 71/69%, comparable to annotated Gonium 71%); R4 and the pre-agreed sensitivity R4_curated (without them) launched. Known gaps: no annotated genomes for Physarum, Zygnema, Penium, hornworts, Coleochaete, Nuclearia; R1 has 10 and R7 9 transcriptomes.
