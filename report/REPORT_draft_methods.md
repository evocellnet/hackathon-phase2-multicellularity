# [TITLE — to be set by the PI]

<!-- Working title suggestion: "Gene content at seven independent origins of multicellularity: an agent-run, pre-registered comparison with Ancestra". -->
<!-- Draft status: sections written BEFORE the main eukaryote results existed. Results, discussion and abstract to be added by the PI. Every number is followed by an HTML comment naming its source file. -->

(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.

[Project-specific text] The goal of this specific project was to investigate multicellularity across the tree of life. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within [TIME] on a computational cluster with 24 cpus, 150GB RAM. if needed, the user will be able to switch the node to a node with a A100 GPU in the future, but at the moment do not rely on it.

Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for implementation of the predictive model and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

<!-- PI-NOTE: the human text said "finished within 24 hours". START.txt: 2026-10-28T12:21:22+01:00, 48 h allocation. Replace [TIME] with the actual elapsed wall time; if it exceeds 24 h, tell the user (their text, their edit). Otherwise copied verbatim from PROMPT.md, including the last paragraph, which mentions a "predictive model" and "benchmarking against other methods" that this project does not have. -->

---

## 1. Why this question

**Multicellularity evolved many times, independently, in very different organisms. Each origin is a natural experiment. If the same genomic changes recur at independent origins, they point to what a multicellular lifestyle requires; if they do not, the origins took different routes.** This project asks what changed in gene content, and in the rate of protein evolution, on the branches of the tree of life where multicellularity appeared, and whether those branches differ from comparable branches where it did not.

### What "multicellular" means

Multicellularity is not one trait but a set of separable properties (LITERATURE.md §1):

- **Clonal vs aggregative.** In *clonal* multicellularity, daughter cells stay attached after division, so all cells descend from one cell and share its genotype (animals, plants, *Volvox*). In *aggregative* multicellularity, free-living cells, which may differ in genotype, come together, as in the fruiting bodies of the social amoeba *Dictyostelium* (Grosberg & Strathmann 2007; Brunet & King 2017). Grosberg & Strathmann (2007) argued that passing through a single cell is the main defence against conflict between cell lineages; Fisher, Cornwallis & West (2013) found that clonal groups were more likely to have become obligately multicellular and to have more cell types.
- **Facultative vs obligate.** A *facultative* multicellular organism can live as single cells and forms groups only under some conditions (starvation in dictyostelids; a bacterial signal in the choanoflagellate *Salpingoeca rosetta*). An *obligate* one is multicellular throughout its life. Márquez-Zacarías et al. (2021) note that all known aggregative organisms are facultative.
- **Simple vs complex.** Knoll (2011) reserves *complex* multicellularity for organisms with differentiated tissues and cell-to-cell channels for nutrients and signals, and restricts it to animals, land plants, florideophyte red algae, laminarialean brown algae and two groups of fungi; Nagy, Kovács & Krizsán (2018) instead count 8–11 origins of complex fruiting bodies within fungi. Multicellular but not complex (filaments, colonies, simple sheets) is *simple*.
- **Grey zones.** Many organisms sit between categories: cells joined by cytoplasmic bridges after incomplete division (Chaigne & Brunet 2022); *coenocytes*, single cells with many nuclei (ichthyosporeans, the plasmodium of the slime mould *Physarum*); facultative colonies such as choanoflagellate rosettes; culture-induced aggregates such as those of *Capsaspora*. Queller & Strassmann (2009) frame the question as a continuum of "organismality" rather than a yes/no trait.

We coded these properties separately for each species and origin, and kept grey-zone cases out of the main tests (Section 3).

### How many origins?

Lamża (2023) counts 45 independent lineages of "simple" multicellular eukaryotes of six types and estimates the true number is likely close to a hundred. Complex multicellularity is much rarer: five or six eukaryotic clades (Knoll 2011), plus the additional fungal origins of Nagy et al. (2018). Several bacterial lineages are multicellular (Claessen et al. 2014); archaeal cases are rare (Rados et al. 2025). Only a few origins have genome data on both sides of the transition with enough relatives to test against: animals, fungi, dictyostelids, volvocine algae, streptophytes (the green-algal lineage that includes land plants), red algae and brown algae (LITERATURE.md §2).

### Competing ideas about the genomic basis

Four ideas from the literature frame our tests (LITERATURE.md §3–4, §7):

1. **Co-option.** Most genes later used for multicellular life already existed in the unicellular ancestor and were re-used. The choanoflagellate *Monosiga brevicollis* has cadherins and tyrosine kinases (King et al. 2008); *Capsaspora* has integrin adhesion machinery (Sebé-Pedrós et al. 2010); about 80% of the proteins essential for dictyostelid development are present in their unicellular relatives (Glöckner et al. 2016); the *Volvox* proteome differs little in coding potential from that of unicellular *Chlamydomonas* (Prochnik et al. 2010); fungal hyphae arose by co-opting ancient genes, without expansions of kinases, receptors or adhesive proteins (Kiss et al. 2019).
2. **A burst of novelty.** Paps & Holland (2018) describe "an unprecedented increase in the extent of genomic novelty" at the origin of animals; Bowles et al. (2020) find two bursts of genomic novelty in the plant lineage; Denoeud et al. (2024) find a marked gain of new gene families at the emergence of brown algae. Such estimates depend on how densely relatives are sampled: adding 19 choanoflagellate transcriptomes moved about 372 supposedly animal-specific families to before the animal–choanoflagellate split (Richter et al. 2018).
3. **Expansion of a "toolkit".** Transcription-factor families (de Mendoza et al. 2013; Lang et al. 2010) and adhesion, extracellular-matrix and signalling families (Knoll 2011; Cock et al. 2010) are proposed to expand repeatedly at transitions; Pereira Lobo et al. (2025) found shared functional themes across 81 eukaryotes, but for complex lineages only. A related idea is that transitions shed unicellular functions such as the flagellum (Merényi et al. 2023: early fungal evolution as a "gradual shedding of protist genes").
4. **Relaxed selection through drift.** Groups of cells that reproduce through bottlenecks have smaller effective population sizes, so selection is weaker relative to chance (genetic drift) and slightly harmful changes accumulate (Lynch & Conery 2003; Lynch 2007). Bingham & Ratcliff (2024) extend this to multicellularity, and snowflake-yeast experiments show more drift in clonal than in aggregative groups (Pentz et al. 2023). The literature review found no genome-wide test of this across independent origins with matched unicellular comparisons (LITERATURE.md §4.3; an absence in our searches, not proof of absence).

What has been missing is one pipeline, run on many origins at once, that compares each origin with comparable non-origin branches (LITERATURE.md, main point 2). Ancestra, a new scalable method for reconstructing ancestral gene content, makes this feasible. The transitions were defined from the literature before any data were analysed; they are inputs, not findings.

---

## 2. How the work was organised

**AI agents ran the project. The design was reviewed three times and frozen before any eukaryote result existed. The reviews changed the design substantially, and two of the most consequential problems were found by the agents' own checks before the results, not after.**

### Agents and the research record

A principal-investigator (PI) agent (Claude Opus 5.5) planned the work, delegated to subagents (Opus 5.5 or Sonnet 5.5), reviewed their output and made the decisions. Subagents wrote the literature review, scouted data, built the dataset, wrote and tested each analysis module, and acted as reviewers. The PI kept a *research tree* (RESEARCH_TREE.md), a status list of every analysis plus a dated decision log (D0–D8) that also served as its memory across context resets. A lessons file recorded process rules learned the hard way, for example capping each agent's CPU use after four agents together pushed the load to 73 on 24 CPUs and stalled the main computation. <!-- src: LESSONS.md rule 6 -->

The literature review came first: about 9,000 words, 249 references each looked up by tool rather than recalled, and 16 testable hypotheses (H1–H16). <!-- src: RESEARCH_TREE.md M1 -->

### Three rounds of design review

<!-- src: reviews/round1-3/*.md; DESIGN.md §9 -->
- **Round 1** (five reviewers: an external PI, an evolutionary biologist, a statistician, a data/algorithms/compute reviewer, and an AI-agent-failure-modes reviewer). The external PI named the most likely misleading headline in advance: "transitions are marked by a burst of new gene families". The data reviewer showed why the pipeline leaned towards it. On simulated genomes with a known history, Ancestra's gene *gains* were 2.1–2.6 times the true number, because it sometimes splits one old family and the younger pieces look like new genes; losses and duplications were accurate (Section 4.2). This led to the gain-evidence ladder: a gain counts as new only if it has no detectable relative outside the clade. The biologist recoded transitions as changes between states. Dictyostelids became simple aggregative, budding yeasts and unicellular zygnematophyte algae became *reversals* to unicellularity, filamentous relatives showed that filaments predate brown algae, and volvocines got two origins (Lindsey et al. 2021). The statistician separated three branch classes (transition stems, inside multicellular clades, inside unicellular lineages). They also showed that RERconverge-style relative rates remove genome-wide shifts by construction, so the drift hypothesis needed a different test. The AI-failure-modes reviewer required tamper-evident pre-registration, one source of truth for report numbers, and toy-data tests for every module. Horizontal gene transfer and site-level convergence were cut for time.
- **Round 2** (three reviewers). The design had promised ≥200 matched controls per transition; only about 40–60 suitable branches exist, so controls were rebuilt around one shared pool (Section 4.3). With about seven transitions, a rank test could never pass the multiple-testing threshold, so each test was given a table of its smallest attainable p-value and detectable effect before results. The expected direction of flagellum loss was wrong in two of three cases (the flagellate *Allomyces* sits inside the hyphal clade; mosses and *Chara* have flagellate sperm). A technical pitfall was caught: DIAMOND reports only 25 hits by default, which close in-clade relatives can fill, so a gene would look falsely new.
- **Round 3** (pre-registration audit). Two transitions were mis-specified on the tree. The node chosen as the start of volvocine multicellularity is inferred to be unicellular in the very study we cited (Lindsey et al. 2021). The fungal-hyphae node rested on one species whose hyphal status Kiss et al. (2019) consider uncertain. Both were fixed with *path stems*, chains of two consecutive branches from the last ancestor in the old state to the first in the new one (D7). <!-- src: reviews/round3/prereg_audit.md F1/F2 -->

Still before any eukaryote output, the PI calibrated the statistics on the real tree. The planned null distribution turned out to be *anti-conservative*: under no true effect it declared significance in 13.8% of simulated datasets at a nominal 5%, because the control sets of different transitions overlap. It was replaced by a calibrated scheme (Section 4.3; D6, Dev-1). <!-- src: results/validation_branchstats/null_calibration.tsv (spec, two-sided, α 0.05: 0.138) -->

### Pre-registration and deviations

Before any eukaryote Ancestra output existed, the following were frozen: the design, transitions, guide trees, species sets, control definitions, toolkit and flagellar gene lists, literature scorecard and power table. The 15 files were hashed with SHA-256 and the repository tagged `prereg-v1`; every analysis script re-checks the hashes and refuses to run on a mismatch. <!-- src: PREREG.md; PREREG_HASHES.tsv (15 entries) --> Later changes are dated deviations, all made before eukaryote results existed:
- **Dev-1:** the gene-gain test (H1/H2) uses a parametric-bootstrap null, because permuting residuals was anti-conservative when transition branches lie outside the controls' covariate range (Section 4.3).
- **Dev-2:** DIAMOND ran at `--very-sensitive` instead of `--ultra-sensitive`. On 500 real proteins it recovered 99.65% of the ultra-sensitive hits at about a third of the cost; the ultra-sensitive setting was estimated at 150–700 CPU-hours. The search covered only the branches the tests use, identically for stems and controls. <!-- src: PREREG.md Dev-2; results/gainladder/diamond_timing.json -->
- **Dev-3:** the planned "split-repair" step failed its pre-registered acceptance test on simulations and was dropped (Section 4.2). This was a test outcome, not a choice.
<!-- PI-NOTE: PREREG.md stamps Dev-2/Dev-3 "~19:00", but the git commit introducing them is 18:19 on 2026-10-28; fix one or the other. -->

### Verification by agents that did not produce the work

Every analysis module was tested on toy data with known answers, planted effects and null simulations (e.g. 38 tests for the statistics engine) before use. <!-- src: docs/branchstats.md --> Independent agents reviewed code and inputs. One found, before any eukaryote result, that an exploratory volvocine transition sharing the same end node silently replaced the primary one in the statistics code; this was fixed. <!-- src: RESEARCH_TREE.md bug entry t≈7h --> Planned: re-derivation of headline numbers by an agent that did not produce them, and an audit of ≥20 claim–citation pairs. <!-- PI-NOTE: update with what verification actually did (RESEARCH_TREE item 7). -->

---

## 3. Data

**We assembled 187 eukaryotic and 82 prokaryotic proteomes (the full set of predicted proteins of a species) and placed them on literature-curated species trees. We coded 52 candidate transitions, of which 7 independent eukaryotic origins met a pre-registered rule for the primary tests (Figure 1).** <!-- src: DATA.md; data/transitions.tsv (52 rows) -->

![Figure 1](../figures/fig1_euk_tree.png)

### Proteomes

| domain | source | species | proteins |
|---|---|---|---|
| eukaryotes | UniProt reference proteomes 2026_03 | 112 | 1,629,372 |
| | EukProt v3 (Richter et al. 2022) | 59 | 1,060,055 |
| | NCBI annotations (incl. Phaeoexplorer brown algae) | 12 | 217,918 |
| | LukProt v1.5.1 (Sobala 2024) | 4 | 61,916 |
| | **total** | **187** | **2,969,261** |
| prokaryotes | UniProt 2026_03 | 75 | 366,790 |
| | NCBI annotations | 7 | 38,874 |
| | **total** | **82** | **405,664** |
<!-- src: DATA.md summary table; data/DATASET_DONE.json -->

Species were chosen to cover both sides of each documented origin. The starting pools were 151 core and 104 optional eukaryotic candidates and 79 core prokaryotic ones. <!-- src: RESEARCH_TREE.md M1 --> All sources are open (UniProt, EukProt, LukProt CC BY 4.0; NCBI/INSDC unrestricted). Restricted JGI data were not used and no data files are committed. DATA.md lists the source, version, date, MD5 and licence of every file, and a script re-creates the dataset. "Eukaryoma" could not be found (the closest match, EUKARYOME, is an rRNA database without proteomes) and was flagged to the user. *Sorogena* and *Sorodiplophrys* have no public proteome. <!-- src: RESEARCH_TREE.md D2; LITERATURE.md §6 -->

**Preparation.** Ancestra expects one protein per gene, so splice variants were collapsed: UniProt canonical sets, and the longest protein per gene for NCBI and genome-based EukProt sets. Thirty-six eukaryotic proteomes come from *transcriptomes* (sequenced messenger RNA), which miss unexpressed genes and contain fragments and near-duplicates. These were reduced to the longest open reading frame per gene and clustered at 95% identity (cd-hit, shorter sequence ≥90% covered). Proteins under 50 amino acids were removed. <!-- src: DATA.md "Processing"; DATASET_DONE.json n_transcriptome 36 -->

**Completeness.** BUSCO 6.1.0 (Manni et al. 2021) measures the share of near-universal single-copy genes (odb10 eukaryota/bacteria/archaea sets) found in a proteome. Median completeness is 84.7% for eukaryotes and 97.6% for prokaryotes. <!-- src: median of busco_C in data/species_final.tsv / species_final_prok.tsv; PI-NOTE: computed ad hoc, recompute in a committed script if kept --> Five eukaryotic proteomes below 50% were kept but flagged as low quality: *Tetrabaena socialis* 32.2%, *Coleochaete scutata* 45.9%, *Porphyra umbilicalis* 43.9%, *Discosporangium mesarthrocarpum* 45.9%, *Guttulinopsis vulgaris* 11.8%. Eight transcriptome sets with >15% duplicated or fragmented BUSCO genes were also flagged. <!-- src: DATA.md "QC flags" -->

**Memory cap.** Ancestra's memory grows faster than linearly with protein number, so the 150 GB node set a cap of 3.0 million eukaryotic proteins. A written rule dropped optional taxa largest-first, never removing a species a transition side depended on. Ten taxa were dropped (e.g. *Trichomonas vaginalis*, *Emiliania huxleyi*), taking the total from about 3.3 million proteins to 2,969,261. <!-- src: DATA.md cap section; data/species_cap_drops.tsv; reviews/round2/methods_compute.md -->

### Guide trees

Ancestra and every test need a fixed rooted species tree. We edited the NCBI taxonomy to current phylogenomic consensus and logged 62 operations on the eukaryote tree (57 applied, 5 skipped because the node was absent), each with rationale and reference. For example, *Gonium* was placed as sister to *Astrephomene* + Volvocaceae following Lindsey et al. (2021). Contested relationships were left unresolved: the tree has 25 nodes with more than two branches (polytomies), including the order of sponges, comb jellies and other animals. The root lies between Amorphea + CRuMs and all other eukaryotes; alternative roots affect only backbone branches, not transition branches. The prokaryote tree needed three edits. <!-- src: data/guide_tree_edits.md -->

Branch lengths (substitutions per site) were estimated once on this fixed topology. We used 250 BUSCO marker genes that are single-copy in at least half of the species: each was aligned (MAFFT) and trimmed (trimAl), giving 64,727 concatenated columns, and the tree was fitted with IQ-TREE 3.1.4 under LG+G4. Branch lengths serve as covariates and for the rate test; they never change which branch is a transition. <!-- src: results/phylo/summary.md; notes.md --> <!-- PI-NOTE: MAFFT, trimAl, IQ-TREE, DIAMOND, cd-hit, pyhmmer are named with versions only; their papers are not in LITERATURE.md. -->

### Transitions

Each transition is a **state change** between unicellular, simple and complex, with type labels (clonal/aggregative; obligate/facultative). **M** is the most recent common ancestor of the sampled species in the new state (the *crown*); **P** is the last sampled ancestor in the prior state. The branch or chain of branches from P to M is the **transition stem**, where the change must have happened. Eukaryotic tip states (Figure 1): 80 unicellular, 36 simple, 39 complex, 22 grey-zone, 10 reversals. <!-- src: data/species_traits.tsv -->

A transition is primary only if all of these hold (D5):
- M is an internal node.
- At least two sampled species are in the new state and at least two are in the sister lineage; low-quality proteomes do not count towards the two.
- The state is not grey-zone. <!-- src: DESIGN.md §2; RESEARCH_TREE.md D5; scripts/prereg/apply_eligibility.py -->

Seven eukaryotic origins qualify (Table 1). Nine further transitions are secondary: three nested simple→complex ones (Pezizomycotina and Agaricomycetes fruiting bodies, land plants) and four prokaryotic ones (myxobacteria, streptomycetes, multicellular magnetotactic bacteria, *Methanosarcina*; Figure S1). Six reversals serve as converse tests, 30 transitions are exploratory, and 2 lack data. <!-- src: data/transitions.tsv analysis_set: euk primary 7, secondary 3, reversal 6, exploratory 22, excluded 2; prok secondary 4, exploratory 8 -->

**Table 1. Primary transitions** (species: sampled in the multicellular clade / its sister lineage)

| # | transition | stem P → M | type | species | notes |
|---|---|---|---|---|---|
| 1 | Animals | Choanozoa → Metazoa | clonal, complex, obligate | 15 / 8 | sisters partly transcriptomes |
| 2 | Fungal hyphae | Eumycota → Blastocladiomycota+Amastigomycota → Amastigomycota | hyphal, simple, obligate | 26 / 4 | path stem (D7) |
| 3 | Dictyostelids | Eumycetozoa → Dictyostelia | aggregative, simple, facultative | 6 / 2 | both sisters transcriptomes |
| 4 | Streptophytes | Streptophyta → Klebsormidiophyceae+Phragmoplastophyta | clonal, simple, obligate | 9 / 2 | contains 3 unicellular reversals |
| 5 | Red algae | Rhodophytina → Eurhodophytina | clonal, simple, obligate | 5 / 2 | M rests on two incomplete Bangiales |
| 6 | Filamentous ochrophytes incl. brown algae | Ochrophyta → PX clade | clonal, simple, obligate | 9 / 6 | P unresolved |
| 7 | Volvocine colonies | Chlamydomonas+TGV → TGV clade → Gonium+Volvocaceae | clonal, simple, obligate | 5 / 3 | path stem (D7) |
<!-- src: data/transitions.tsv (label, stem_path, type, n_mc_species, n_sister_species, flags); red-algae note: reviews/round3/prereg_audit.md §1 -->

---

## 4. Methods

**The core idea is simple: reconstruct which genes existed at every ancestor, count gene gains, losses and duplications on every branch, and ask whether transition stems stand out against comparable branches where multicellularity did not arise. Most of the effort went into making that comparison fair (Figure 2).**

![Figure 2](../figures/fig2_workflow.png)

### 4.1 Ancestral gene content with Ancestra

**Ancestra** (evocellnet/ancestra-mvp-ai, commit ffb0bbe, tag `freeze-v1`) is a recently developed, itself AI-built method for *orthology inference*: deciding which genes in different species descend from the same ancestral gene. It works in three steps:
1. It groups all proteins into families of related sequences with DIAMOND, without a reference database.
2. It scores similarity only within families.
3. It walks the guide tree from tips to root, merging the genes of descendant lineages into ancestral genes by mutual best matches and resolving ambiguous cases with small local gene trees.

Each ancestral gene at a node, with all its descendants, is a **hierarchical orthologous group (HOG)**: one inferred gene of that ancestor. The approach descends from the bottom-up methods SYNERGY and Hieranoid (Schreiber & Sonnhammer 2013) and from the HOG framework of OMA and FastOMA (Altenhoff et al. 2013; Majidian et al. 2025). Its addition is a reference-free family step feeding one merge over all families. <!-- PI-NOTE: add SYNERGY (Wapinski et al. 2007) to LITERATURE.md. --> Its benchmarks were dominated by animals and fungi; behaviour on deep, protist-rich trees like ours was untested (LITERATURE.md §5.1). We ran the pre-registered `fast` preset separately for eukaryotes and prokaryotes (orthology across domains adds cost without informing any transition; D0), with a watchdog that killed the run above 145 GB. <!-- src: DESIGN.md §3.2; RESEARCH_TREE.md D0, D3; scripts/run_ancestra.sh -->

**From HOGs to events.** A parser converts HOGs into per-branch events:
- a **gain**: a HOG first appearing on the branch;
- a **loss**: an ancestral gene present at the upper node and absent at the lower one;
- a **duplication**: one ancestral gene becoming two.

Given the true HOGs of simulated data, the parser reproduced the true events exactly, including with gene transfer, whole-genome duplication and multifurcating trees. <!-- src: results/validation_hogparser/summary.md (a) -->

### 4.2 Why raw gains are inflated, and the gain-evidence ladder

**Ancestra's gains cannot be taken at face value; losses and duplications can, with caveats.** On three simulations with Ancestra's own simulator, inferred losses and duplications tracked the truth closely (per-branch correlation mostly ≥0.98; totals within about 2–25%). The exception was a high-turnover simulation with a whole-genome duplication, where internal-branch losses were inflated 1.5–1.8-fold. Gains were inflated 2.1–2.6-fold in every simulation. The cause: Ancestra sometimes splits one true family into several HOGs at deep nodes. For 6–11% of genes the inferred family origin was too young (almost never too old), and spurious gains concentrated on long branches and near the root. <!-- src: results/validation_hogparser/summary.md (b) and Key points --> This bias favours exactly the misleading headline anticipated in review, because transition stems are long, deep branches.

The eukaryote run itself showed why. Linking proteins by similarity joined 73% of all proteins (2.16 million) into one connected group. Ancestra cuts such groups into families of at most 7,500 genes with the Leiden community-detection algorithm, so many family boundaries are set by cut points rather than biology, and each artificial boundary can create a spurious gain. <!-- src: RESEARCH_TREE.md D8 -->

Every gained HOG was therefore ranked on a pre-registered **gain-evidence ladder**:
1. **HOG gain**: the raw call; an upper bound, never a headline.
2. **Sequence-novel gain** (primary): up to three representative proteins are searched with DIAMOND (`--very-sensitive`, all hits kept with `-k 0`, E ≤ 1e-5, query coverage ≥30%) against all 269 proteomes. The gain is novel only if no representative hits a protein outside the clade. <!-- src: DESIGN.md §3.2; PREREG.md Dev-2; diamond_timing.json -->
3. **Architecture-novel gain** (secondary): the HOG's ordered Pfam domain combination is absent outside the clade.

A planned intermediate step re-joined split HOGs to the older HOG they came from. It had to bring simulated inflation to ≤1.3-fold but left it at 2.12–2.54-fold, because most spurious gains are extra copies of a family that continues into the clade (missed duplications), not gains paired with losses. It was dropped (Dev-3). <!-- src: results/gainladder/validation_sim/summary.md §1 -->

The novelty filter is not a cure. On simulations it kept 91–97% of true gains (excluding transferred genes), but only 42–56% of its "novel" calls were real. The simulator lets about a third of families diverge beyond detection: spurious novel gains were about 4.6–4.9 substitutions per site from their nearest outside relative. Real deep families usually keep conserved domains, so the simulation is probably pessimistic. Still, novel-gain counts are upper bounds. This is the homology-detection failure that Weisman et al. (2020) showed explains many, but not all, lineage-specific genes. Conclusions rest on comparing stems and controls under the same filter, never on absolute counts. <!-- src: results/gainladder/validation_sim/summary.md §2–3 -->

### 4.3 Branch-level statistical design

**Unit and comparison.** The unit is a branch. Each primary stem is compared with a **control pool** U of branches that meet all of these:
- unicellular at both ends and inside unicellular lineages;
- at least two descendant species on each side;
- not a child of the root and not adjacent to any stem.

On our tree |U| = 30. <!-- src: DESIGN.md §4.6; docs/branchstats.md --> Branches inside multicellular clades are excluded, as they would mix being multicellular with becoming so.

**Matching.** Some branch properties change event counts regardless of biology. Longer branches accumulate more events, branches above many species collect more gains, and branches whose relatives are transcriptomes gain spuriously from genes missed there. For each stem, the 20 controls closest in depth, branch length, descendant and sister species numbers, and transcriptome fractions are chosen (Mahalanobis distance on standardised covariates). Imbalance is reported as a standardised mean difference (SMD: difference in means in control standard deviations). Matching is poor: stems are deeper than controls (depth SMD 1.77). <!-- src: RESEARCH_TREE.md D6 -->

**One statistic per test.** Each stem's value is ranked among its matched controls, giving a percentile u between 0 and 1. The test statistic is the mean u over the seven stems, so every transition weighs equally.

| test | question | branch statistic | direction |
|---|---|---|---|
| H1/H2 | co-option vs burst of novelty | sequence-novel gains, as residual from a negative-binomial regression fitted on controls | two-sided |
| H4 | transcription-factor expansion | share of stem-duplicated HOGs carrying a TF DNA-binding domain, vs the share among genes present at the branch start | one-sided |
| H5 | adhesion/ECM and signalling expansion | as H4 for adhesion/ECM + signalling domains | one-sided |
| H9 | rate shift | median over marker genes of log(distance P→M tips / P→sister tips) (4.5) | two-sided |
| H12 | flagellar loss where expected | fraction of four flagellar modules lost, signed by a pre-set direction | one-sided |
<!-- src: PREREG.md primary family; DESIGN.md §4.6 -->

The H4/H5 domain list (110 Pfam families in six categories, including plant, algal and other non-animal adhesion families so as not to favour animals; C2H2 zinc fingers reported separately) and the H12 modules (axoneme, intraflagellar transport, BBSome, transition zone) were frozen in advance. Loss is expected only on the fungal-hyphae and dictyostelid stems. Red algae are excluded from H12 because their flagellum was lost before the stem. <!-- src: data/toolkit_pfam.tsv (110 rows); data/h12_sets.tsv; data/h12_expected_direction_frozen.tsv -->

**Gains as a regression residual.** Because matching is poor, raw gain counts would mostly measure depth. For H1/H2, a *negative-binomial generalised linear model* (NB-GLM; regression for over-dispersed counts) is fitted on controls only. It predicts sequence-novel gains from branch length, descendant and sister species numbers, and transcriptome fraction. A branch's statistic is log(observed/predicted); each control is predicted from a fit that leaves it out. <!-- src: DESIGN.md §4.6 v3.1 (D6) -->

**Null distributions, and how the planned one failed.** A p-value needs a *null distribution*: the statistic's values when nothing special happens on the stems. The design first relabelled random controls as "transitions". On the real tree this rejected 13.8% of null datasets at a nominal 5% (4.5% at 1%), because the seven stems share many of the same 30 controls. Permuting values over stems and controls jointly was calibrated (5.1%, 1.3%) and is the primary scheme for H4, H5, H9, H12. For the gain residual, permutation failed once stems lie outside the controls' covariate range, as ours do: 8.3% false positives at a 1-SD shift, 28.3% at 2.5 SD. A **parametric bootstrap** stayed at 4.0–5.3% and became the H1/H2 null (Dev-1): it simulates counts for all branches from the control fit, refits, and recomputes the statistic. The permutation p is a sensitivity analysis. <!-- src: results/validation_branchstats/null_calibration.tsv; docs/branchstats.md -->

**Multiple testing and effect sizes.** The five primary p-values are Holm-corrected, which controls the chance of any false positive among them; the first threshold is 0.05/5 = 0.01. Permutations use 10,000 draws, seed 11 (pre-registered seeds 11, 23, 37, 41, 53), so the smallest attainable p is 1/10,001. Per transition, we report observed/expected with the control range, plus a random-effects pooled estimate with leave-one-transition-out and leave-one-supergroup-out checks. Non-significant results are worded "not detected at this power" with the compatible interval; "no ≥1.5-fold effect" is used only if the confidence interval lies within 1/1.5–1.5. All other analyses are exploratory, logged in an analysis registry and labelled. <!-- src: DESIGN.md §4.6, §7; power_table.tsv (min_attainable_p 9.999e-05) -->

**Power, stated before results.** With seven stems and 30 controls only large effects are detectable (80% power, α = 0.05; effects planted on the real tree):

| test | minimum detectable effect | power at 1.5-fold |
|---|---|---|
| H1/H2 | 2.4-fold (moderate over-dispersion) to 5.9-fold (high) | 0.35 / 0.13 |
| H4 | odds ratio 2.8 | 0.28 |
| H5 | odds ratio 2.2 | 0.35 |
| H9 | rate ratio 1.26 (low gene noise) to 1.69 (high) | 1.00 / 0.66 |
| H12 | not reached (power ≤0.33 even with complete loss) | — |
<!-- src: results/prereg/power_table.tsv; docs/branchstats.md -->

In plain words, only the rate test is sensitive. The gene-content tests can detect a doubling or more, not the subtler differences the literature debates. H12, with two loss-expected stems, cannot reach 80% power at all; it stayed primary because the pre-registered rule required two loss-expected stems, not adequate power. The simulations assumed equal branch lengths (the branch-length tree did not yet exist). They also showed the gain test slightly over-rejecting (8.5% at nominal 5%) before the switch to the bootstrap. <!-- src: power_table.tsv (placeholder_1.0; size 0.085); PREREG.md Power -->

### 4.4 Annotation

All proteins were annotated with Pfam-A 38.2 domains (Mistry et al. 2021; 30,134 families; pyhmmer 0.12.3 with curated gathering thresholds). They were also annotated with predicted *intrinsic disorder*, the fraction of residues lacking a stable fold, by metapredict 3.1 (Emenecker et al. 2021; score ≥0.5 counts as disordered). A HOG belongs to a domain class if at least half its members carry a class domain. The class tests use a pre-registered subset of Pfam models, run first; full Pfam supports domain-wide and architecture analyses. <!-- src: DATA.md; data/annotations/*/*.DONE.json; results/README_run_primary.md gap 1 -->

### 4.5 Rate of protein evolution (H9)

If multicellular lineages experience more drift, their proteins should evolve faster on average. For each of 250 marker genes, IQ-TREE fitted a gene tree on the guide tree pruned to the species present (taxa covering <50% of the alignment removed). For each transition and gene we took the log ratio of the mean path length from P to multicellular tips over that from P to sister tips. The transition's value is the median over genes, with a 2,000-resample gene bootstrap interval. Comparing two lineages from the same ancestor cancels changes before P. Controls are the 10 control-pool nodes with two descendant clades of ≥2 species, randomly oriented. <!-- src: DESIGN.md §4.6 H9; results/phylo/notes.md; results/h9/summary.md (methods paragraphs only) --> A difference is a **rate shift**, not proof of relaxed selection. Generation time and mutation rate also change rates; claiming relaxed constraint needs a within-protein contrast such as radical vs conservative amino-acid changes, run as an exploratory check. Codon dN/dS was not used because synonymous sites saturate at these depths (LITERATURE.md §4). A gene-by-gene RERconverge-style analysis (Chikina et al. 2016; Kowalczyk et al. 2019) was exploratory.

### 4.6 Structure (exploratory plan)

Sequence searches miss distant relatives that share a fold. The plan applies to sequence-novel gains with an AlphaFold model (Jumper et al. 2021; AlphaFold DB v6, Varadi et al. 2024) of mean per-residue confidence (pLDDT) ≥70. Foldseek (van Kempen et al. 2024) searches them against models of the unicellular P-side species (11 species, 123,882 proteins), or as fallback against 590,183 Swiss-Prot models excluding the clade. A gain is "structurally old" with a hit at E ≤1e-3 and TM-score ≥0.5 (TM >0.5 generally means the same fold). Control-branch gains get the same treatment. In tests, 16/16 proteins were found in a database containing them, and 0/42 shuffled models hit anything. The analysis is limited to UniProt species and by the database's rate limit (~3 requests/s; ~12–13 h for the reference set), and is exploratory. <!-- src: docs/STRUCTURE_DATA.md -->

### 4.7 Controls and reproducibility

Sixteen eukaryotic and two prokaryotic **positive controls** were frozen with explicit pass rules: events the literature firmly expects, such as a Wnt-family gain on the animal stem. The pipeline passes if ≥11 of 16 eukaryotic controls are recovered and every high-confidence failure is explained. <!-- src: data/positive_controls.tsv (C01 PF00110); PREREG.md --> The prokaryote run was the pipeline's dry run. Reproducibility measures:
- code, small results and decisions are in git;
- long jobs write a completion file only on success;
- species identifiers are cross-checked across files;
- seeds are fixed;
- every test is appended to an analysis registry, so the total number of tests is known;
- report numbers are to be traced to committed tables. <!-- src: DESIGN.md §7 -->

---

## 5. Limitations known before the results

- **Few origins, low power.** Seven stems and 30 controls give 13–35% power for a 1.5-fold gene-content difference; flagellar-loss power never reaches 80%. Most negative results will be inconclusive. <!-- src: power_table.tsv -->
- **Gains are noisy.** Simulated gain inflation is 2.1–2.6-fold and novelty-filter precision 0.42–0.56; only stem-vs-control comparisons are interpreted.
- **Stems are unlike controls** (depth SMD 1.77); the regression and bootstrap correct this only as far as the model is right.
- **Uneven data quality.** Transcriptomes dominate some sister lineages (both dictyostelid sisters; part of the choanoflagellates). Two stems end in clades represented through low-quality proteomes (*Porphyra* for red algae; *Tetrabaena* on the volvocine path). Missed genes in sisters inflate stem gains; missed genes in the clade deflate them.
- **The tree is fixed.** One root and literature-based edits; polytomies remain (e.g. the animal base); alternative trees were not re-run.
- **Ancestra is new.** It is AI-developed and benchmarked mainly on animals and fungi. It merged 73% of eukaryotic proteins into one similarity component before splitting it at a size cap, so family boundaries in large superfamilies (kinases, TFs, adhesion domains) are partly arbitrary. <!-- src: RESEARCH_TREE.md D8 -->
- **Rate-test controls are mixed.** Only 10 control nodes qualify, and in 9 of them a descendant clade contains multicellular species. <!-- src: results/phylo/notes.md H9 control caveat -->
- **A quirk of the H4/H5 statistic.** Its continuity correction gives a positive value when a stem has no TF duplications but the baseline share is very low; this matters only for stems with few duplications. <!-- src: results/README_run_primary.md -->
- **Labels come from the literature.** Transitions, toolkit lists and positive controls come from the same literature the results are compared with, so some agreement is built in.
- **Proteins only.** Regulatory change, which much of the literature emphasises, is invisible here.
