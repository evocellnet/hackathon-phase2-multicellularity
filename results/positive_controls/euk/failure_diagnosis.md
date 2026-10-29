# Positive controls (euk): why five controls failed

Independent verifier, 2026-10-28. Inputs: `runs/analysis/euk/*.parquet` (Ancestra fast run, primary guide tree), pre-registered Pfam
subset at GA, `results/positive_controls/euk/controls.tsv` (10/16 recovered). Nothing frozen was changed (`check_prereg.py` → PREREG OK
before and after this work), and the scorer was **not** modified.
Evidence files:
- `failure_diagnosis_relaxed_hmm.tsv`: a single-HMM pyhmmer search with no GA, sequence E ≤ 0.01, for PF01049, PF00362, PF00907,
  PF01414, PF06816 and PF07684 on all 69 Obazoa proteomes. The `passes_GA` column flags hits that reach the gathering cutoff.
- `failure_diagnosis_scripts/`: the analysis scripts (`s1`–`s12`, `relax.py`). Run each with the ancestra env Python. `relax.py`
  needs the annot env.

## Summary

| control | conf. | failing component(s) | verdict category | code bug? |
|---|---|---|---|---|
| C03 classical cadherin PF01049 | high | (a) no PF-HOG with origin on stem(Metazoa) | **(iv)** HOG-level origin cannot see a domain gained inside an older gene lineage. Ancestra calls it a duplication on the Metazoa stem of a Choanozoa cadherin HOG. The family root is also pushed to Eukaryota by 19 stray genes. | no |
| C06 integrin β PF00362 | high | (a), (b), (c) all fail | **(iv)** for (a): the Thecamonas integrin sits only in an over-merged EGF-superfamily HOG (19% integrin), so it is not a PF-HOG. **(i)** for (b) and (c): Fonticula has no integrin, so the loss is on the Holomycota stem, not the Fungi stem. *Didymoeca costata* has a real integrin β, as Richter et al. 2018 also report, so there is no loss on the Choanoflagellata stem. | no |
| C07 T-box PF00907 | high | (c) no duplication on stem(Metazoa) | **(iv)** Ancestra puts the 9 T-box duplications on the **Choanozoa stem** (Filozoa→Choanozoa). Every child HOG is metazoan-only, and the duplications are paid for with losses on the Choanoflagellata stem. This placement is systematic (52% of Choanozoa-stem duplication products are metazoan-only). | no |
| C08 DSL + NOD/NODP | moderate | (a) and (b) | **(i)** for (b): the cited source (Richter 2018) states that NOD/NODP are absent from *every* choanoflagellate and sponge, so (b) could never pass. **(i)/(ii)** for (a): no DSL at GA in our 8 choanoflagellates, only sub-GA DSL-like EGF proteins in 4 of them, which sit in non-DSL HOGs. Richter's "strong" Delta needs MNNL, which is not in the pre-registered Pfam subset. | no |
| C12 pherophorin PF12499 | moderate | ratio ≥ 2 | **(i)**, an operational threshold above the published effect. Observed 138 vs 71 = 1.94× (Δ = 67). Prochnik et al. 2010 report 49 vs 27 = **1.81×**. The "expansion" is also not multicellularity-specific (unicellular *Edaphochlamys* has 135). | no |

**No code bug.** For every failed control I compared the scorer's logic with the pass_rule text and the `prereg_definitions.tsv` terms,
and they match. One point deserves explicit confirmation. `Lineage.origin_row` walks up to the root HOG of the lineage, and this
equals origin(H) as defined. All 825,157 gain events in `events.parquet` are either root HOGs (72,429) or species-specific genes
(752,728 singletons on leaf branches). No sub-HOG ever carries a gain event. The C12 counts were re-checked against the raw
`.pfam.tsv.gz` files (138 / 71). I am therefore proposing no patch.

If these are accepted as explanations, the pass count stays at 10/16 (C10 is still undetermined). The pre-registered criterion
needs ≥ 11 recovered, so it is **not met** whatever the explanations. The C03, C06 and C07 explanations are offered for the
"every high-confidence failure explained" clause.

---

## C03: classical cadherin cytoplasmic domain (PF01049, CADH_Y-type_LIR); expected gain on stem(Metazoa)

**Verdict: (iv) Ancestra / operationalisation granularity.** The biology holds (the domain is animal-only), but the domain arises
inside a gene lineage that is older than Metazoa.

The hits are clean. At GA: 44 proteins in 10 metazoan species and **0 outside Metazoa**, so (b) passes. A relaxed search finds only
noise outside Metazoa (scores 14.7–17.1, E ≈ 0.005–0.01 in Thectra, Acanspe, Didycos and five fungi). Sponges Amphque and Sycocil,
*Mnemiopsis* and *Trichoplax* have no hit even when relaxed. Nichols et al. 2012 detected sponge classical-cadherin cytoplasmic
domains by other means, so this is a minor sensitivity limit of the Pfam model (ii). It has no bearing on (a).

Where the 44 proteins sit:

| lineage root (family) | origin | genes in root | PF01049 proteins | PF-HOGs inside lineage (levels) |
|---|---|---|---|---|
| Eukaryota\|HOG0058630 | Eukaryota | 2,496 (2,467 Holozoa; 19 non-opisthokont genes from 11 species; 10 Holomycota) | 40 | 12 (Metazoa, ParaHoxozoa, Bilateria, Cnidaria, …) |
| Obazoa\|HOG0012365 | Obazoa | 20 | 3 (Branbel) | 2 (Bilateria, Deuterostomia) |
| Eukaryota\|HOG0064100 | Eukaryota | 94 | 1 (Nemavec) | 1 (Cnidaria) |

All Metazoa-level PF-HOGs have Choanozoa-level parents that **contain choanoflagellate cadherins**:

| Metazoa-level PF-HOG [PF01049 / genes] | parent at Choanozoa [PF01049 / genes; members] |
|---|---|
| Metazoa\|HOG0058630.294 [9/18] | Choanozoa\|HOG0058630.288 [9/27; 21 Metazoa + 6 Choanoflagellata]. **dup event on stem(Metazoa)** → .294 + .310 |
| Metazoa\|HOG0058630.386 [22/48] | Choanozoa\|HOG0058630.380 [22/60; 48 Metazoa + 12 Choanoflagellata] |
| Metazoa\|HOG0058630.251 / .277 | Choanozoa\|HOG0058630.246 [7/53; 48 Metazoa + 5 Choanoflagellata] |

All of these continue up through Filozoa\|HOG0058630.244 (356 genes, 287 carrying cadherin repeats PF00028) to the Eukaryota root.
The cytoplasmic domain is therefore a metazoan architectural novelty grafted onto a cadherin lineage that is shared with
choanoflagellates. Ancestra represents it as a **duplication on the Metazoa stem**, not as a gain. That is a defensible HOG call; a
HOG-level origin simply cannot detect a domain gain. The extra depth of the root, which is Eukaryota rather than Choanozoa or
Holozoa, comes from 19 of 2,496 genes (0.8%) in non-opisthokonts. That is a typical stray-member over-deepening.

**Implication:** HOG-level origins (ladder tier 1) do not date domain or architecture novelties. Many "animal innovations" in the
literature (classical cadherin, RTK architectures, Notch/Delta) are defined by architecture. They will appear as duplications or
inherited lineages on the stem, not as gains. Deep roots of large multidomain families are unreliable because a handful of stray
members pushes them to Eukaryota. Gain counts are less affected (one spurious deep root removes at most one stem gain), but origins
of large families should not be quoted at face value. The architecture-novel tier of the gain ladder is the appropriate test for
this kind of control.

## C06: integrin β (PF00362); expected present with Thecamonas, lost on Fungi and Choanoflagellata stems

**Verdict:** (a) fails from **(iv)**, over-merging combined with the 50% PF-HOG rule. (b) and (c) fail from **(i)**: the expectation
does not hold in our sample, and for (c) it is contradicted by Richter et al. 2018.

GA hits by group, with relaxed-search counts in parentheses:

| group | species with PF00362 |
|---|---|
| Apusozoa | Thectra 1 (1), Pygsbif 1 (1) |
| Metazoa | all 15 (1–40 each) |
| Choanoflagellata | **Didycos 2** (5 relaxed); the other 7 have 0 (0) |
| Filasterea | Capsowc 4, Minivib 3, Pigochi 6 |
| Ichthyosporea / Pluriformea | Creofra 1, Sphaarc 1; Coralim 1, Syssmul 1 |
| Holomycota (35 spp. incl. Fonticula) | **0 at GA**; relaxed only noise (Coprcin 16.7, Crypneo 15.6) |

The HOG lineage of the Thecamonas integrin (Thectra_A0A0L0DEK2), with [PF00362 proteins / genes]:
Eukaryota\|HOG0060409 [129/830] → Podiata .6 → Amorphea .7 → Obazoa .40 [129/755] → **Apusomonadida_Opisthokonta\|HOG0060409.42
[126/657 = 19%]**. The Thectra gene sits directly in this HOG → Opisthokonta .43 [125/656] → Holozoa .138 [70/140] and Holozoa .252
[55/99]. These are the first PF-HOGs.

- **(a)** The only HOGs containing the Thecamonas integrin are 15–19% integrin, so none is a PF-HOG and (a) = 0. The other 531
  members of the Opisthokonta-level HOG carry EGF (PF00008), EGF_CA (PF07645), PK_Tyr (PF07714) and laminin-EGF domains. This is an
  EGF-superfamily over-merge driven by integrin β's I-EGF repeats. The integrin lineage itself is inferred at the right depth (root
  at Eukaryota, Thectra included), but the frozen PF-HOG ≥ 50% rule cannot see it.
- **(b)** Fonticula alba (genome, BUSCO C 72%) has no integrin hit even when relaxed. The loss therefore maps to the Holomycota stem,
  as the pi_resolution caveat anticipated. A test on stem(Fungi) cannot pass. In Ancestra the over-merged lineage is even "present"
  in Holomycota, via Holomycota\|HOG0060409.44, which holds 6 *Parvularia* laminin-EGF proteins with no integrin domain.
- **(c)** The *Didymoeca costata* proteins P010559 and P019539 (GA scores 84.6 / ≈80; 1,670–1,690 aa) sit in the integrin-β PF-HOG
  Filozoa\|HOG0060409.140 [62/67]. That HOG is therefore present_at(Choanoflagellata), and there is no loss on the Choanoflagellata
  stem. These are not contaminants: DIAMOND best hits against all Holozoa are only 24–30% identical (Oscarella, Corallochytrium,
  Capitella; E ≈ 1e-17 to 1e-20). Richter et al. 2018 (eLife, doi:10.7554/eLife.34226) found "'strong' evidence for integrin β (in the
  form of the diagnostic Integrin_beta Pfam domain) in only a single species of choanoflagellate". The C06 expectation comes from
  Sebé-Pedrós et al. 2010, which had only *Monosiga* and *Salpingoeca*. The newer source, which we sample, contradicts it.

**Implication:** the presence pattern is right, but losses scored on specific stems are fragile. One transcriptome species that
retains the gene, or one unsampled sister lineage (Fonticula), moves or removes the loss. The over-merging of EGF-rich proteins into
600–2,500-gene deep HOGs (also seen for C03 and C08) affects the class-composition statistics (H5 adhesion/ECM ∪ signalling). Such
families are a single lineage at deep levels and fall below the 50% class threshold.

## C07: T-box (PF00907); expected present in Holozoa, lost in choanoflagellates, duplicated on the Metazoa stem

**Verdict: (iv) Ancestra duplication placement.** (a) and (b) pass. (c) fails because the duplications sit one branch too deep.

The biology holds in our sample:

| group | T-box at GA (relaxed) |
|---|---|
| Metazoa | 5–21 per species |
| Choanoflagellata | 0 in all 8 (0 relaxed) |
| Filasterea, Ichthyosporea | 2–6 per species |
| Fungi | Basimer 2, Morthyg 4, Pirofin 1, Rozeall 1, Spizpun 1 (early-diverging, as in Sebé-Pedrós 2013) |

171 of 175 T-box proteins fall in one family, Opisthokonta\|HOG0012132 (173 genes). Its events:

| branch | event | HOG → children |
|---|---|---|
| Holozoa→Filozoa | dup ×4 | Holozoa\|.7 → Filozoa\|.9, .100, Pigochi singletons |
| **Filozoa→Choanozoa** | **dup ×9** | Filozoa\|.9 → Choanozoa\|.11, .28, .44, .58, .69, .93, .95 (+2 singletons). **Every child is metazoan-only** (22, 36, 18, 24, 32, 3, 3 genes) |
| Choanozoa→Choanoflagellata | losses | the 9 Choanozoa-level copies |
| Choanozoa→Metazoa | **no dup** | each Choanozoa HOG passes 1:1 into a Metazoa HOG |

Because no choanoflagellate carries a T-box, gene presence cannot separate "duplicate on the Choanozoa stem, then lose 9 copies in
choanoflagellates" from "lose 1 copy in choanoflagellates, duplicate on the Metazoa stem". Ancestra picks the former, which is less
parsimonious. The scorer reads the events table correctly (`branch_child == "Metazoa"`, parent or child in the PF-HOG set), and there
is simply no such event.

**This placement is systematic.** Of 5,728 child HOGs produced by duplications on stem(Choanozoa), 2,994 (52%) are metazoan-only and
1,571 are choanoflagellate-only. Only 1,163 (20%) contain both. All 7,070 children of duplications on stem(Metazoa) are, as expected,
metazoan-only. The pattern also inflates Choanoflagellata-stem losses (14,531 loss events). It is the same split/ordering artefact
that the failed split-repair acceptance test (Dev-3) points to.

**Implication, the most consequential one here for the primary tests:** duplications of paralogs whose sister lineage lacks the
family are dated one branch too early. This shifts duplications from the Metazoa stem onto the Choanozoa branch, which is P(Metazoa)
→ a control-pool or neighbouring branch. H4/H5 use stem duplications on primary stems, so the Metazoa stem duplication count d, and
the TF-DBD share that C07 is about, are probably underestimated. The bias is not specific to Metazoa: any stem whose sister lacks
the family is affected, and stems are compared against control branches that may suffer it less or more. I recommend that the PI
quantify the same metric (share of duplication products on branch P→X that lack any member of X's sister clade) for every primary
stem and its matched controls, and report it next to H4/H5. That would be a sensitivity analysis, not a rule change.

## C08: Delta (DSL PF01414) and Notch (NOD PF06816 / NODP PF07684) present at Choanozoa with choanoflagellate members

**Verdict:** (b) fails from **(i)**: the pre-registered operationalisation contradicts the cited source. (a) fails from **(i)/(ii)**:
the sampled species and annotation sensitivity.

- **(b) NOD/NODP.** There are zero choanoflagellate hits, at GA and relaxed (E ≤ 0.01). Capsaspora, all sponges and Mnemiopsis also have
  none. Nematostella and Hydra have only sub-GA NOD hits (3 and 2). Richter et al. 2018 state explicitly: "In our data, the NOD and NODP domains were not present in any
  sponge or choanoflagellate species", and they identified Notch by the architecture EGF + Notch(LNR) + TM + Ank instead. Component
  (b) could not have passed with these Pfams on any choanoflagellate sample, and this is a design error in the frozen rule. It is
  not a scorer or Ancestra issue.
- **(a) DSL.** No choanoflagellate protein reaches DSL GA (24.2). Sub-GA DSL-like hits appear in Choafle (20.3), Mylnflu (20.6),
  Codohol (18.1) and Acanspe (17.2), E ≈ 4e-4 to 6e-3. They sit in EGF-family HOGs with no GA-DSL member: Choanozoa\|HOG0067099.184
  (58 genes, 53 metazoan), Choanozoa\|HOG0058660.3 and Choanozoa\|HOG0073219.22. These could not become DSL PF-HOGs. Richter et al.
  2018 report a Delta (MNNL + DSL) in one choanoflagellate and one of the two domains in four others. The species names are lost in
  the PMC text, so I cannot confirm whether the Delta species is among our 8. MNNL (PF07657) is not in the pre-registered Pfam subset,
  so Richter's "strong" criterion cannot be checked here. Metazoan DSL proteins are themselves scattered across at least 6 families,
  several of them EGF-superfamily HOGs rooted at Eukaryota or Podiata (e.g. Metazoa\|HOG0060929.446, Filozoa\|HOG0019749.186). DSL
  outside Metazoa occurs in Capsowc, Pigochi and Creofra, and in non-opisthokonts such as Physpol, Guilthe and Vermver.

**Implication:** this failure says little about Ancestra. It reflects a mis-specified control, which the moderate confidence label
already allowed for. It does reinforce the C03 point that architecture-defined genes are poorly served by HOG + single-Pfam rules.

## C12: pherophorin (PF12499) copy number, Volvox ≥ 2 × Chlamydomonas and Δ ≥ 5

**Verdict: (i)**, an operational threshold set above the published effect size. This is not a data or code problem.

| species | state | PF12499 proteins (GA) | proteome |
|---|---|---|---|
| Volvox carteri | multicellular | **138** | UniProt, 14,332 |
| Chlamydomonas reinhardtii | unicellular | **71** | UniProt, 17,611 |
| Volvox africanus | multicellular | 113 | |
| Pleodorina starrii | multicellular | 43 | |
| Gonium pectorale | multicellular | 39 (BUSCO 71%) | |
| Tetrabaena socialis | multicellular | 63 (BUSCO 32%) | |
| Edaphochlamys debaryana | unicellular | **135** | |
| C. incerta / C. schloesseri | unicellular | 93 / 87 | |

The ratio is 1.94 and Δ = 67, so the rule fails on the ratio by 4 proteins (Volvcar would need ≥ 142). The counts were verified
against the raw annotation files, and both proteomes are UniProt canonical with isoforms collapsed (BUSCO D 0.8%). The scorer's
counts are exact. Prochnik et al. 2010 (Science, doi:10.1126/science.1188800) report "pherophorins (49 versus 27 members)", a ratio
of **1.81**. The pre-registered 2× is therefore above the literature effect, and our observed ratio is stronger than the published
one. Across species, high pherophorin counts are not specific to multicellular lineages (Edaphochlamys 135, C. incerta 93), so
"Volvox-lineage expansion" is weak support for a multicellularity-linked expansion in any case.

**Implication:** none for Ancestra (this control uses no HOGs). It is a near-miss on an operational threshold.

---

## What this means for trusting Ancestra's origin calls

1. **Gain origins of compact, single-domain families are trustworthy.** C01 (Wnt), C02 (TGF-β), C10 and C13 recovered. Where a HOG
   failed here, it was never because the right genes were absent or split into too-young HOGs. For C06 the Thecamonas integrin is in
   the correct, deep lineage. No case showed the "split into young HOGs" failure mode (category iv-a).
2. **The failure modes are over-merging and over-deepening.** Large multidomain families (cadherin, EGF/integrin, Delta) form deep
   HOGs of 600–2,500 genes whose roots are pushed to Eukaryota by a few stray members, and whose class or Pfam fraction falls below
   50%. Origins of such families and their class membership (H5) should be read with care.
3. **Duplication timing is biased one branch deeper when the sister lineage lacks the family (C07).** 52% of duplication products on
   the Choanozoa stem are metazoan-only. This is the result most relevant to H4/H5 and deserves a sensitivity check (see C07).
4. **Domain or architecture novelties inside older lineages are invisible at HOG tier 1 (C03, and partly C08).** They need the
   architecture-novel tier.
5. Three of the five failures (C06 b/c, C08 b, C12) come from the expectations or thresholds themselves, and two of those
   contradict the more recent or cited literature. They are not evidence against the pipeline.


**Update 2026-10-29:** after the eukaryote gain ladder became available, C10 (GRAS) was re-scored as recovered, so the final count is **11/16** (results/positive_controls/euk/summary.md). The diagnoses above are unchanged.
