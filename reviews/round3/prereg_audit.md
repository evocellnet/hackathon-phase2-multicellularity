# Round 3: pre-registration input audit

Leaf sets were checked with ete3. Quality data come from `species_final.tsv` (low quality = BUSCO C < 50%). PubMed full text: Lindsey 2021 (doi:10.1186/s12915-021-01087-0), Kiss 2019 (doi:10.1038/s41467-019-12085-w).

## 1. Primary transitions

| transition | M / P / sister | verdict |
|---|---|---|
| T_metazoa | Metazoa / Choanozoa / Choanoflagellata (8) | OK |
| T_dictyostelia | Dictyostelia / Eumycetozoa / Physpol, Cerafru | OK. Both sisters are plasmodial transcriptomes, so the P state rests on Variosea and *Mastigamoeba*. Cerafru is coded `unicellular` but `coenocytic`. |
| T_streptophyta | Klebsormidiophyceae_Phragmoplastophyta / Streptophyta / Mesovir, Chloatm | OK. Only one sister is in the prior state (Chloatm is grey). |
| T_rhodophyta | Eurhodophytina / Rhodophytina / Porppur, Rhodmar | OK. *Rhodosorus* and *Porphyridium* are unicellular. **Bias:** M exists only through two incomplete Bangiales (Porpumb C = 44%, Pyroyez C = 50%). Gains missed in both shift to the nested stem, which deflates this stem's gain count. |
| T_px_filaments | PX_clade / Ochrophyta polytomy / 6 unicells | Defensible. The filamentous ancestor is inferred, and leaving Raphidophyceae in the polytomy is conservative. `analysis_tier` contradicts `analysis_set`. |
| **T_volvocine_tgv** | TGV_clade / Chlamydomonas_TGV / Chlamydomonas | **Error (F1)** |
| **T_fungi_hyphae** | Blastocladiomycota_Amastigomycota / Eumycota / chytridiomyceta | **Error (F2)** |

**D4.** The citation is correct. Lindsey et al. 2021 find:
- multicellularity "independently evolved at least twice";
- Tetrabaenaceae sister to a unicell (MLBS = 100);
- monophyly of the colonial species rejected by an AU test (p = 2.8e−38);
- Astrephomene sister to Volvocaceae, which matches edits #39–41.

D4 nevertheless treats the single-origin model (Herron & Michod 2008) as an equal alternative. Under the cited phylogeny TGV_clade is unicellular, so the "transition stem" is unicellular→unicellular. M also exists only through Tetrsoc (C = 32%). TGV-stem gains count only if detected in that proteome; otherwise they fall to the Gonium_Volvocaceae stem.

**T_fungi_hyphae.** Among sampled taxa, M contains no chytrid-like unicell (only yeast reversals). But M exists only through Allomac. Kiss et al. 2019 (our own source):
- describe Blastocladiomycota as mostly "unicellular thalli", with likely convergent hypha-like forms;
- give a hyphal posterior of 0.53 at Blastocladio + higher fungi versus 0.92 at Zoopagomycota+;
- place Blastocladiomycota *before* Chytridiomycota. Li et al. 2021 report basal polytomies.

The stem therefore rests on one species of contested placement and contested homology.

## 2. Required fixes

**F1.** Set M = Gonium_Volvocaceae and P = Chlamydomonas_TGV, and make the stem the *path* P→TGV_clade→M, with Tetrsoc as a side branch (the D5 logic). The transition lies on this path under either model.
- Fallback: demote the volvocines to secondary.

**F2.** Use the path Eumycota→Blastocladiomycota_Amastigomycota→Amastigomycota with M = Amastigomycota. Allomac becomes an exploratory tip origin, like Gonapodya. The sister (3 prior-state chytrids) is unchanged.
- Fallback: keep the current M, with the path as a pre-registered sensitivity analysis.

**Path stems (add to §4.6):**
- Gains: origin on either edge. Losses: present at P and absent at M. Duplications: summed over both edges.
- Covariates: BL = sum of the edges, depth = depth(M), n_desc = tips(M), sister = sister of the path.
- H9: M tips exclude the side-branch species.
- U: excludes all edges adjacent to the path.

**F3 (stale or inconsistent fields):**
- T_volvocine_tgv flags `sister_contains_multicellular(Tetrsoc)`, and its M_kind is `internal`.
- T_px_filaments tier vs set.
- `guide_tree_edits.md` still says "PX secondary, Phaeophyceae primary".
- §4.6 places Gonapodya "inside the hyphal crown" (it is in the sister clade).
- C11 says "two Volvocales primary transitions".
- `h12_expected_direction.tsv` is superseded.

**F4.** Five cited works are missing from LITERATURE.md, which violates §7: Herron & Michod 2008, van Dam et al. 2013, Muñoz-Gómez et al. 2017, James et al. 2006, Wickett et al. 2014.

## 3. H12 (`data/h12_expected_direction_frozen.tsv`)

The rule is "present at P, absent at M".
- **Current coding:** loss is expected only for Dictyostelia. Metazoa, Streptophyta, PX, the volvocines and fungi (*Allomyces* is flagellate and inside M) are no_loss. Rhodophyta is not_applicable. With 1 loss-expected stem, **H12 is exploratory** and Holm covers 4 tests.
- **With F2:** all 26 Amastigomycota are non-flagellate, so the fungal path becomes loss-expected. That gives **2 stems → H12 primary**.
- **Caution:** decide F2 on the phylogenetic grounds above, not to rescue H12, and log it as D6 before any output is parsed. The F2 row is marked conditional.

## 4. Positive controls (`data/positive_controls_nodes.tsv`)

All expected nodes exist, except:
- **C14:** no "Agaricomycetes" node. The 3 sampled species are Agaricales, so use node `Agaricales`.
- **C06:** the integrin loss may map to the Holomycota rather than the Fungi stem; report both.
- **C11:** wording (single P).
- **C15:** follows F2.
- **C09:** the Phaeophyceae crown depends on low-quality Discmes.

## 5. Secondary set (flags only)

- **Sister not in the prior state:** Pezizomycotina and Embryophyta (the sisters are reversals; Embryophyta's only prior-state reference, Colescu, is low quality).
- **Demote to exploratory:**
  - Agaricomycetes: its sister contains Tremella, which has the new state.
  - Magnetotactic MMB: M = Desulfobacterales is half unicellular (metagenomic drafts).
- **Uncertain P state:** Streptomyces.

## 6. Should

- Add min child-clade BUSCO C to the §4.6 covariates; Rhodophyta and Phaeophyceae are otherwise mismatched.
- State in PREREG that eligibility counts grey sister species (Chloatm, Physpol) toward "≥2 sister species".

No further error invalidates a primary test. Without F1 and F2, two primary stems are mis-specified and biased toward deficits.
