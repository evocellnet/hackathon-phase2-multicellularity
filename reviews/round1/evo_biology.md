# Round 1 review — evolution of multicellularity / comparative genomics (prior-work assessor)

Reviewer: evo-biology agent. Inputs: PROMPT.md, DESIGN.md v1, LITERATURE.md, `data/species_candidates.tsv`, `data/species_final.tsv` (196 species), `docs/DATA_SOURCES_SCOUTING.md`.

Conventions:
- All Pfam accessions were checked against the local `data/raw/pfam/Pfam-A.clans.tsv.gz` (Pfam 38.2).
- [abstract-checked]: confirmed in a PubMed abstract during this review.
- [verify]: from reviewer knowledge; the PI should check it before freezing.

**Verdict.** In spirit the design follows the literature: co-option vs novelty, matched control branches, outgroup density. The transition model has one systematic flaw. It assumes every transition goes from a last *unicellular* ancestor P to a multicellular M. In fact:
- Several key transitions are nested (simple → complex).
- Several "closest unicellular sisters" in the species table are reversals, or are multicellular themselves.

Two §3.7 positive controls contradict the cited literature (Hox, Notch/Delta). A third cannot be measured with Pfam: RTK via PF07714. All of this can be fixed before compute.

## 1. Transition definitions

**Cross-cutting (must).** Code each transition in `data/transitions.tsv` as a state change at one of three levels:
- **U→S**: unicellular to simple multicellular.
- **S→X**: simple to complex (Knoll 2011: tissues plus cell–cell channels; or ≥3 cell types).
- **M→U**: reversal to unicellularity.

P is then the last ancestor in the *prior* state, not the last unicellular ancestor. Without this, nested transitions (fruiting bodies, kelps, Florideophyceae, embryophytes) are compared with P nodes that were already multicellular. That dilutes the transition effect in the NB-GLM. Run every test under a strict and a permissive coding (LITERATURE §1.3).

| Transition | M / P | Status | Rationale |
|---|---|---|---|
| Metazoa | Metazoa crown / Choanozoa | **primary U(F)→X** | P state: "facultatively colonial?". Run with both the ctenophore-first and the sponge-first root. |
| Choanoflagellate colonies | no internal M | **grey, tip-level, exploratory** | Coloniality is scattered across Craspedida, facultative and partly aggregative (*C. flexa*; Brunet et al. 2019), so no origin node is defensible. The table gives the same species both the "unicellular sister" and the "multicellular" role. Resolve this. |
| *Capsaspora* aggregation | tip / Filasterea | **grey, exploratory** | One species; the aggregative stage is a culture phenotype. |
| Ichthyosporea | — | **exclude from primary** | Coenocytic growth is likely ancestral to Teretosporea (*Corallochytrium*; Kożyczkowska et al. 2021 [verify]), so P and M have the same state. |
| Fungal hyphae | Blastocladio+Zoopago+Mucoro+Dikarya / Fungi minus *Rozella* | **primary U→S** (coenocytic, so grey under strict coding) | Kiss et al. 2019: early origin by co-option, with *no* kinase, receptor or adhesion expansion [abstract-checked]. *Gonapodya*: separate tip origin [verify]. |
| Pezizomycotina fruiting bodies | Pezizomycotina / Saccharomyceta | **primary S→X** | Saccharomycotina are a **reversal** (Nagy et al. 2014), not a unicellular sister. Recode them as M→U. |
| Agaricomycete fruiting bodies | Agaricomycetes+Dacrymycetes / Agaricomycotina | **primary S→X** | Tremellomycetes: possible separate origin (Nagy et al. 2018); use as a sensitivity coding. ***Wallemia* is filamentous, not a yeast.** |
| *Neolecta* | tip / Taphrinomycotina | secondary S→X | BUSCO 71. |
| Yeast reversals | Saccharomycotina, *Schizosaccharomyces*, *Saitoella*, *Cryptococcus*, *Malassezia*, *Rhodotorula* | **M→U class** | Converse test for H1 and H12. |
| Dictyostelia | Dictyostelia / Eumycetozoa | **primary U→S (aggregative)** | The table codes Dictyostelia as "complex", but Knoll 2011 and LITERATURE §2 say simple. **Fix this, because it drives H4 and H8.** *Physarum* is the coenocytic (grey) sister clade, not a "unicellular" one. |
| *Fonticula*, *Acrasis*, *Copromyxa* | tips | secondary aggregative tips | Acrasidae fall within Tetramitia (Brown et al. 2012b), so *Neovahlkampfia* is a better sister than *Naegleria*. Verify *Copromyxa*'s sister (Brown et al. 2011). |
| *Guttulinopsis* | — | **presence-only** | BUSCO 11. |
| *Sorogena*, *Sorodiplophrys* | — | **drop** | No data. |
| Volvocales | **two** M nodes | **primary U→S ×2** | Lindsey et al. 2021: multicellularity arose "at least twice" and Goniaceae is not monophyletic [abstract-checked]. Set the guide tree first. *Volvox* is polyphyletic (Herron & Michod 2008). Code *V. carteri* and *Pleodorina* as "differentiated", not "complex". |
| Streptophyta | M1 = Klebsormidiophyceae+Phragmoplastophyta / Streptophyta (U→S); M2 = Embryophyta / Zygnematophyceae+Embryophyta (S→X); *Chara* tip S→X | **primary** | ***Spirogloea*, *Mesotaenium* and *Closterium* sit inside M1** (Zygnematophyceae are sister to land plants), so they are reversals, not unicellular sisters. *Chlorokybus* packets make P grey. |
| Rhodophyta | M1 = Eurhodophytina / Proteo+Eurhodophytina (U→S); M2 = Florideophyceae (S→X) | **primary** | Stylonematophyceae and Compsopogonophyceae filaments are independent origins (Muñoz-Gómez et al. 2017 [verify]). The flagellum was lost before M. Bangiales have low BUSCO, so flag inferred losses there. |
| Phaeophyceae | U→S within the PX clade (P includes Raphidophyceae); Phaeophyceae = S→X | **primary S→X** | Filaments **predate** Phaeophyceae: *Schizocladia*, *Phaeothamnion* and *Tribonema* are filamentous. The table labels them "unicellular sister" yet codes them "simple; obligate". Fix this. *Ectocarpus* and *Discosporangium* are filamentous, so code them simple. |
| Ulvophyceae | no single M (Del Cortona et al. 2020) | *Ulva* tip U→S; Bryopsidales grey | Add unicellular ulvophytes if data exist [verify]. |
| Oomycetes, *Vaucheria*, Myxogastria, *Plasmodiophora*, *Labyrinthula* | — | **grey (coenocytic or plasmodial); permissive coding only** | Lamża 2023, types (i) and (ii). |

## 2. Species choice

- **Missing (add if data exist):**
  - *Paraphelidium* (Aphelida, the sister of Fungi; Torruella et al. 2018). With *Rozella* alone, gains on the fungal stem are inflated (the lesson of Richter et al. 2018).
  - *Hyphochytrium* or *Developayella* as oomycete sisters.
  - *Yamagishiella* or *Eudorina* to break the long Volvocaceae stem [verify].
  - *Interfilum* (Ma et al. 2026) [verify].
- **Roles:** make `role` specific to each transition *and* level, and give every tip a state code (U, S, X or grey). Keep *Thecamonas* in the core set, because the integrin control needs it.

## 3. Positive controls in §3.7

1. **Hox at Metazoa is wrong.** Sponges and ctenophores lack Hox and ParaHox genes (Ryan et al. 2010 [abstract-checked]; Larroux et al. 2008). PF00046 is pan-eukaryotic and cannot identify Hox. Drop this control.
2. **Notch/Delta "complete at Metazoa" contradicts LITERATURE §3.1.** Richter et al. 2018 found homologues in choanoflagellates. Re-specify it (A1 #8).
3. **RTK via PF07714.** In Pfam 38.2, PF07714 is "PK_Tyr_Ser-Thr": it matches Ser/Thr and TKL kinases, including plant RLKs. Tyrosine kinases are also already present in *Capsaspora* (Suga et al. 2012). Use the SH2+kinase architecture or the catalytic-loop motif (HRDLAARN for TKs vs HRDLKPEN for Ser/Thr). The expected node is Filozoa.
4. **Cadherin repeat is pre-choanozoan.** *Capsaspora* has a cadherin (Nichols et al. 2012 [abstract-checked]). The classical cadherin (PF01049) is the true Metazoa gain.
5. **Integrin β is correct but incomplete.** It predates Opisthokonta and was lost independently in fungi and choanoflagellates (Sebé-Pedrós et al. 2010 [abstract-checked]). That makes it a good *loss* control.
6. **Kinase expansion in Dictyostelia has no transition-specific support.** Replace it with PKS expansion (low-to-moderate confidence).
7. **Pherophorin (PF12499) is correct**, but the expansion is at the *Volvox* tip. VARL/regA has no Pfam model, so anchor a HOG on *V. carteri* regA instead.

**Replacement:** Appendix A1 lists 16 eukaryote controls plus 2 prokaryote controls. They cover four kinds of expectation: gain, co-option, loss and expected null. Freeze them in `data/positive_controls.tsv` with a pass criterion fixed in advance. Suggested criterion: ≥11/16 recovered, and every high-confidence control either recovered or explained. Report every failure.

## 4. Toolkit and H12 sets

- **Toolkit (Appendix A2):**
  - Exclude LRR proteins that carry NB-ARC, NACHT or TIR domains from "adhesion". Immune LRR expansions in plants and *Ectocarpus* would otherwise create false convergence for H5.
  - Report C2H2 (PF00096) outside the primary H4 test, because lineage-specific repeat expansions dominate it.
- **H12 is not testable as stated.** Loss of flagella and phagotrophy tracks trophic mode and terrestriality, not multicellularity:
  - Animals, volvocines and brown algae keep their flagella.
  - Red algae lost them before M.
  - Aggregative amoebae stay phagotrophic.
  - Hyphae *co-opted* phagocytosis genes (Kiss et al. 2019).

  Pre-register an expected direction for each transition (A3) and run one-sided tests only where a loss is predicted.
  - Flagellar set: axoneme, IFT and transition-zone domains only. Exclude dynein heavy-chain and kinesin domains, which also cover cytoplasmic paralogues.
  - Score the BBSome separately: it is the module lost first (van Dam et al. 2013 [abstract-checked]).
  - Phagotrophy set: map the gene set of Burns et al. 2018 to HOGs. Arp2/3 and WASP are shared with endocytosis and tip growth, so they do not work as markers.

## 5. Hypotheses: attribution and testability

- **H4:** Nagy et al. 2014 concerns Zn-cluster TFs at *yeast reversals*. Move it to the reversal analysis. With about 6 S→X transitions, the complex×transition interaction has low power; state this in advance.
- **H8:** Grosberg & Strathmann and Fisher et al. predict *organismal* differences, so label H8 a derived prediction. As designed, it is confounded: every aggregative case except Dictyostelia is a single-species terminal branch, often from a transcriptome, while the clonal cases are internal branches. Restrict H8 to contrasts matched for level and branch type, or make it exploratory.
- **H9:** Use sister contrasts within one level (U→S separately from S→X). Protein length and disorder are strongly lineage-structured (Schad et al. 2011).
- **H6:** Needs kinase classification beyond PF07714 counts (§3, item 3). Hyphae are the predicted null (Kiss et al. 2019).
- **H12:** See §4. The reversal half (Yang et al. 2026) can only be tested in the prokaryote run.
- **H13:** Retinoblastoma (RB) is a *co-option* control in volvocines (Hanschen et al. 2016), not a duplication prediction.
- **H14:** Use GRAS gained at Zygnematophyceae+Embryophyta as an HGT positive control (Cheng et al. 2019 [abstract-checked]). Compare *Acrasis* against matched tips.
- **H1–H3, H5, H7, H10, H11, H15, H16:** attributions match LITERATURE.md.

## 6. Ranked changes

**Must**
1. Recode the transitions as U→S, S→X or M→U with an explicit P state (§1). Fix the contradictory roles: Saccharomycotina, Zygnematophyceae, the PX-clade filaments, *Wallemia*, *Physarum*, the choanoflagellates, and Dictyostelia (complex → simple).
2. Fix the guide tree for Volvocales (two origins; Lindsey et al. 2021) and for the streptophytes before deriving M nodes.
3. Replace the §3.7 controls with A1, frozen before results. Drop Hox and "Notch at Metazoa"; stop treating PF07714 as a TK marker.
4. Re-scope H12 with expected directions per transition and axoneme-specific sets (A3).
5. Remove *Guttulinopsis*, *Sorogena*, *Sorodiplophrys* and Ichthyosporea from the primary tests. Treat choanoflagellate colonies and *Capsaspora* as grey-zone and exploratory.

**Should**
6. Add *Paraphelidium*, *Hyphochytrium*/*Developayella*, *Yamagishiella*/*Eudorina* and *Interfilum* if data exist.
7. Adopt the A2 toolkit list, with immune LRRs and C2H2 excluded as in §4.
8. Fix the H4 and H8 attributions and restrict H8 to matched contrasts.

**Could**
9. Add a reversal analysis (yeasts, Zygnematophyceae, *Interfilum*) as the converse of H1 and H12.
10. Add sensitivity codings for Tremellomycetes and *Astrephomene*.

---

## Appendix A1 — proposed pre-registered positive controls (eukaryote run unless noted)

| # | Feature (Pfam / HOG) | Expected node | Direction | Source | Conf. |
|---|---|---|---|---|---|
| 1 | Wnt PF00110 | Metazoa stem | gain (absent outside) | King et al. 2008; Richter et al. 2018 | high |
| 2 | TGF-β PF00019 (+ PF00688 propeptide) | Metazoa stem | gain | Richter et al. 2018; Huminiecki et al. 2009 [verify] | high |
| 3 | Classical cadherin cytoplasmic PF01049 | Metazoa stem | gain | Nichols et al. 2012 | high |
| 4 | Cadherin repeat PF00028 | Filozoa (present at P of Metazoa) | co-option | Nichols et al. 2012; Abedin & King 2008 | high |
| 5 | SH2 PF00017 + PF07714 architecture; PTP PF00102 | Filozoa or earlier | co-option | Manning et al. 2008; Pincus et al. 2008 [verify]; Suga et al. 2012 | high |
| 6 | Integrin β PF00362 | pre-Opisthokonta (present with *Thecamonas*) | co-option + **losses** on Fungi and Choanoflagellata stems | Sebé-Pedrós et al. 2010 | high |
| 7 | T-box PF00907 | Opisthokonta/Holozoa; lost on Choanoflagellata stem; diversified (duplications) at Metazoa | co-option + loss + duplication | Sebé-Pedrós et al. 2013 [abstract-checked: in ichthyosporeans, filastereans, fungi] | high |
| 8 | DSL PF01414; Notch NOD PF06816 / NODP PF07684 | Choanozoa (in some choanoflagellates) | co-option | Richter et al. 2018; Gazave et al. 2009 | moderate |
| 9 | TALE homeobox KN PF05920 | present at P of Streptophyta, Phaeophyceae and Rhodophyta transitions | co-option (convergent recruitment) | Lee et al. 2008 [verify]; Arun et al. 2019 [abstract-checked]; Petroll et al. 2025 | high |
| 10 | GRAS PF03514 | Zygnematophyceae+Embryophyta | gain by HGT (best hit bacterial) | Cheng et al. 2019 | high |
| 11 | Retinoblastoma RB_A PF01858 / RB_B PF01857 | present in *Chlamydomonas* (P of volvocine M) | co-option, no gain | Hanschen et al. 2016 | high |
| 12 | Pherophorin PF12499 | *V. carteri* vs *C. reinhardtii* copy number | expansion (tip/Volvocaceae) | Prochnik et al. 2010 | moderate |
| 13 | Hex1 (eIF5A-paralog HOG; PF01287) | Pezizomycotina stem | duplication/gain (Woronin body) | Jedd & Chua 2000 [verify] | moderate |
| 14 | Hydrophobin PF01185 (class I) | Agaricomycetes | expansion; absent in yeasts | Wösten 2001 [verify]; Krizsán et al. 2019 | moderate |
| 15 | Kinase/receptor/adhesion classes (A2) | fungal hyphal node | **no** expansion (expected null) | Kiss et al. 2019 | high |
| 16 | Polyketide synthase (PF00109 + PF00698 architecture) | Dictyostelia | expansion | Eichinger et al. 2005 | low–mod |
| P1 | HetR_C PF18460; HetP PF28069; FraC PF24301 (prokaryote run) | filamentous cyanobacteria | presence vs unicellular | Herrero et al. 2016; Hammerschmidt et al. 2021 | moderate |
| P2 | Rodlin PF25848; BldD-like_C PF21179 (prokaryote run) | *Streptomyces* | presence | McCormick & Flärdh 2012 | moderate |

## Appendix A2 — toolkit categories (Pfam 38.2 accessions, verified locally)

- **Adhesion / ECM:** PF00028 Cadherin, PF01049 cadherin cytoplasmic, PF00362 Integrin_beta, PF01839 FG-GAP, PF08441 Integrin_A_Ig_1, PF01044 Vinculin, PF00041 fn3, PF00047 ig, PF07679 I-set, PF07686 V-set, PF00008 EGF, PF07645 EGF_CA, PF00053 EGF_laminin, PF00055 Laminin_N, PF00052 Laminin_B, PF02210/PF00054 Laminin_G, PF01391 Collagen, PF00090 TSP_1, PF00084 Sushi, PF00431 CUB, PF00754 F5_F8_type_C, PF00059 Lectin_C, PF00092 VWA, PF00100 Zona_pellucida, PF01390 SEA, PF02494 HYR, PF02469 Fasciclin, PF12499 Pherophorin, PF00057 Ldl_recept_a. LRR (PF00560, PF13855) only if the protein lacks NB-ARC PF00931, NACHT PF05729 and TIR PF01582.
- **Signalling (ligands/receptors):** PF00110 wnt, PF01534 Frizzled, PF01392 Fz, PF00019 TGF_beta, PF03166 MH2 (SMAD), PF01414 DSL, PF00066 Notch LNR, PF06816 NOD, PF07684 NODP, PF01085 HH_signal, PF02460 Patched, PF00001/PF00002/PF00003 7tm_1/2/3, PF05462 Dicty_CAR, PF01825 GPS, PF00029 Connexin, PF00876 Innexin.
- **Kinases / pTyr:** PF00069 Pkinase, PF07714 PK_Tyr_Ser-Thr (subclassify TK vs S/T by catalytic-loop motif), PF00017 SH2, PF00102 Y_phosphatase, PF03637 Mob1_phocein (Hippo). Receptor kinase = kinase + TM + an extracellular domain.
- **TF DNA-binding domains:** PF00046 Homeodomain, PF05920 Homeobox_KN, PF00157 Pou, PF00292 PAX, PF00010 HLH, PF00170 bZIP_1, PF00907 T-box, PF00853 Runt, PF00250 Forkhead, PF00178 Ets, PF00505 HMG_box, PF00104/PF00105 nuclear receptor, PF00870 P53, PF00554 RHD_DNA_bind, PF01285 TEA, PF09271 LAG1-DNAbind, PF02864 STAT_bind, PF00319 SRF-TF (MADS), PF00249 Myb, PF00847 AP2, PF03106 WRKY, PF02365 NAM, PF00320 GATA, PF03110 SBP, PF03634 TCP, PF03195 LOB, PF02362 B3, PF03514 GRAS, PF02701 Zn_ribbon_Dof, PF00172 Zn_clus, PF02319 E2F/DP, PF00447 HSF, PF02376 CUT, PF01388 ARID, PF00751 DM, PF02257 RFX. C2H2 PF00096 is reported separately.
- **Programmed cell death:** PF00656 Peptidase_C14 (caspase/metacaspase), PF00619 CARD, PF00531 Death, PF01335 DED, PF00452 Bcl-2, PF00653 BIR, PF01823 MACPF.
- **Cell cycle / cytokinesis / bridges:** PF01858/PF01857 RB, PF02319 E2F, PF00134/PF02984 Cyclin, PF00735 Septin, PF03357 Snf7 (ESCRT-III).

## Appendix A3 — H12 sets and expected direction

**Flagellum (score modules separately):**
- IFT-B: PF10498 IFT57, PF14931 IFT20, PF12317 IFT46_B_C, PF18383 IFT81_CH, PF10234 Cluap1.
- BBSome: PF14779 BBS1, PF14781/PF14783/PF14782 BBS2, PF14777 BBIP10.
- Axoneme: PF04712 Radial_spoke, PF06098 Radial_spoke_3, PF03148 Tektin, PF10211 Ax_dynein_light, PF14738 CFAP91, PF11069 CFAP298, PF10629 CMI2B-like.
- Transition zone: PF07162 B9-C2, PF09773 Meckelin, PF07773 TCTN.
- Exclude PF03028/PF00225, which include cytoplasmic paralogues.
- Optionally anchor HOGs on the *C. reinhardtii* flagellar proteome (Pazour et al. 2005 [verify]).

**Phagotrophy:** Burns et al. 2018 model genes mapped to HOGs; supplementary domains PF04727 ELMO_CED12 and PF06920 DHR-2 (DOCK), labelled exploratory.

**Expected direction:**
- *Loss expected:* fungal hyphal node (flagellum, phagotrophy); Dictyostelia stem (flagellum only); Zygnematophyceae+Embryophyta lineage (flagellum lost in Zygnematophyceae).
- *No loss expected:* Metazoa, Volvocales, Phaeophyceae, Ulva, *Acrasis*/*Fonticula*/*Copromyxa* (phagotrophy retained).
- *Not applicable:* Rhodophyta (lost before M).
