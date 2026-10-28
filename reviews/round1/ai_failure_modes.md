# Round 1 review: AI-agent failure modes

Files read: PROMPT.md, DESIGN.md v1, RESEARCH_TREE.md, LESSONS.md, LITERATURE.md (§7, reference conventions), `data/species_final.tsv` header, `mcell/` and `tests/` listing.

## Assessment
The design is strong on statistical pre-registration (§1b, §3.7), but pre-registration is currently prose only. Nothing makes it tamper-evident, and files it depends on (`toolkit_pfam.tsv`, `positive_controls.tsv`, `transitions.tsv`) do not yet exist. Main risks: (1) post-hoc drift of those files, (2) circularity between literature-derived labels and literature "confirmation", (3) untested analysis code beyond the OrthoXML parser, (4) report numbers and citations drifting from committed results.

## MUST

**M1. Hash-frozen pre-registration, before any Ancestra output is parsed.**
Implementation: after round 3, write `PREREG.md` with sha256 and git commit of DESIGN.md, `transitions.tsv`, `toolkit_pfam.tsv`, `positive_controls.tsv`, primary-test list and seeds; tag `prereg-v1`. `scripts/check_prereg.py` runs at the start of every analysis script and the report build and fails on mismatch. Later changes go in a dated "Deviations" section (what, why, which results had been seen), never in-place edits.

**M2. Single source of truth for numbers and figures.**
Implementation: analyses write only to committed `results/*.tsv`, indexed in `results/MANIFEST.tsv` (file, script, git hash, seed, status). `scripts/build_report.py` fills `{{placeholders}}` in the report template from those tables; figures come from scripts only. A lint step fails on literal numbers in the template outside an allow-list. Report agents edit prose only.

**M3. Citation integrity gate.**
LITERATURE.md already tags [unverified] (e.g. the "at least 25 origins" figure), but the report agent may paraphrase from memory, and §7 rows list 3-6 authors per claim, inviting attribution creep.
Implementation: a script checks that every `[Author Year]` in the report resolves to a reference in LITERATURE.md; [unverified] items are dropped or worded "reportedly". A fresh agent audits at least 20 claim-citation pairs against abstracts (`agents/literature/abstracts.py`) and logs pass/fail in `reviews/citation_audit.md`. Tool and database versions asserted in DESIGN (Pfam-A 38.2, UniProt 2026_03, metapredict v3, EukProt v3, LukProt v1.5.1) are checked against downloaded file headers and recorded in DATA.md, not copied from agent memory.

**M4. Toy-data tests for every analysis module.**
Only `mcell.hogs` is tested. The NB-GLM, matched-control permutation null, co-option index, expansion test, RER normalisation, sister contrasts and architecture counting are new, agent-written code.
Implementation: each `mcell/<x>.py` gets `tests/test_<x>.py` with (a) a planted effect (rate ratio 2 recovered, CI covers), (b) a null simulation (false-positive rate about alpha over 200 reps), (c) a hand-computable 3-species case. No result is cited unless tests pass at the commit recorded in MANIFEST. The reviewer agent differs from the author and receives code plus tests, not the author's summary.

**M5. Break label circularity.**
Transitions, types, toolkit Pfam lists and positive controls all come from LITERATURE.md; §7 then "tests" claims from the same papers. Agreement on cadherin, integrin or Wnt is nearly guaranteed by construction.
Implementation: (i) add equal-sized negative controls with expected null to `positive_controls.tsv` (housekeeping domain classes; random Pfam sets matched on size and abundance) so a method that calls everything "expanded" fails; (ii) the report states that transition classification is an input, not a finding; (iii) leave-one-transition-out and relabelling sensitivity for contested labels (fungal hyphae, Volvocales, Acrasidae), with a frozen confidence column in `transitions.tsv`; (iv) Q7 reports concordance as reproduces / contradicts / not addressed.

## SHOULD

**S1. Analysis registry against cherry-picking.** `results/ANALYSIS_LOG.tsv` is appended automatically by a decorator in the stats module for every test run (id, tier primary/secondary/exploratory, script, hash). The report states the total number of tests; only pre-registered primary tests may be headlined; exploratory p-values are BH-adjusted across the full logged family.

**S2. ID and name reconciliation.** `scripts/check_ids.py`, run after each build, asserts that the sets of species in `species_final.tsv`, FASTA prefixes, guide-tree tips, `transitions.tsv`, OrthoXML and annotation tables are identical. Synonyms are resolved only through a committed alias table. Also assert that each transition node M exists and is monophyletic for its listed species after tree editing (the tree edits in §2 are the likely breakage point).

**S3. Partial-run guard.** Long jobs write `DONE.json` (expected vs produced counts, exit code, wall time, output md5) only on success; MANIFEST status derives from it; analysis scripts refuse inputs without it. Per-branch conservation checks (§3.1) are hard failures. Timed-out or subsampled runs go in a limitations table with the completed fraction and subsampling seeds.

**S4. Handoff template and claims register.** Each agent writes `agents/<name>/HANDOFF.md`: task, input paths and hashes, outputs, commands, tests passed, known failures, proposed claims with file:row evidence. A claims register (claim, evidence, verifier, status) is kept; the PI marks [x]/[!] in RESEARCH_TREE.md only after a verifier entry exists. Tree entries carry commit hash and result path, so the tree rebuilds context after a reset. Agents append a three-line summary to the tree before exiting.

**S5. Independent re-derivation of headline claims.** For about 5-8 claims marked [!], a different agent re-derives the number from raw inputs with its own implementation; more than 5% disagreement blocks citation.

**S6. Novelty wording.** "Consistent with" is the default. "Novel" is allowed only after the literature agent has searched for it (query log saved) and is worded "to our search". Limitations must say Ancestra is new, AI-developed and untested at this scale; absence of prior report is not novelty.

## COULD
- C1. Add LESSONS rules on handoffs, DONE.json and prereg after M1 (current four rules are all environment-related).
- C2. Red-team agent at about t=28 h: for each headline result, propose one artefact explanation (transcriptome annotation bias, taxon sampling density, branch-length covariate) and run one cheap check.
- C3. Hold out about 10% of positive controls when setting thresholds (e.g. the "supported gains" definition).
- C4. A `make verify` target (tests, prereg check, ID check) with a one-line summary the PI reads at each milestone.

## Minor
- H1 and H2 are opposed but share one test; state that the output is an effect size with CI, and what outcome would count as "neither".
- DESIGN §5 targets at most 36 h while PROMPT's public text says 24 h; the report must quote actual wall time from START.txt.
