"""D17 recurrence ladder, hypothesis level. Rules: results/r2/recurrence_ladder/RULES.txt (fixed before computing).
Reads only existing per-event outputs; no new per-event test. Outputs:
  hypothesis_events.tsv  one row per (ladder, hypothesis, event) with the per-event p and significance flags
  hypothesis_ladder.tsv  one row per (ladder, hypothesis, criterion, k) with observed, expected, P(X>=k), Holm/BH per tier
"""
import math
import sys
import os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import R, R1_ORIGINS, R1_NESTED, LABEL, out, poibin_sf, holm, bh  # noqa: E402

ALPHAS = (0.05, 0.10)
rows = []


def add_rank(ladder, hyp, desc, event, unit, u=None, k=None, p=None, source=""):
    """Rank/permutation statistic. Either mid-rank u against k controls, or a permutation p with k controls."""
    if k is None or (isinstance(k, float) and math.isnan(k)):
        return
    k = int(k)
    if p is None:
        if u is None or (isinstance(u, float) and math.isnan(u)):
            return
        p = 1 - float(u) + 0.5 / (k + 1)
    p = float(min(1.0, max(p, 1.0 / (k + 1))))
    r = dict(ladder=ladder, hypothesis=hyp, description=desc, event=event, unit=unit, test_type="rank (calibrated)",
             p_event=p, k_controls=k, p_floor=1.0 / (k + 1), source=source)
    for a in ALPHAS:
        r[f"sig_{a}"] = int(p <= a + 1e-9)
        r[f"pi_{a}"] = math.floor(a * (k + 1) + 1e-9) / (k + 1)
    r["sig_top"] = int(p <= 1.0 / (k + 1) + 1e-9)
    r["pi_top"] = 1.0 / (k + 1)
    rows.append(r)


def add_nominal(ladder, hyp, desc, event, unit, p, source=""):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return
    r = dict(ladder=ladder, hypothesis=hyp, description=desc, event=event, unit=unit, test_type="nominal p (uncalibrated)",
             p_event=float(p), k_controls=np.nan, p_floor=0.0, source=source)
    for a in ALPHAS:
        r[f"sig_{a}"] = int(p <= a)
        r[f"pi_{a}"] = a
    r["sig_top"] = np.nan
    r["pi_top"] = np.nan
    rows.append(r)


def add_sign(ladder, hyp, desc, event, unit, stat, source=""):
    if stat is None or (isinstance(stat, float) and math.isnan(stat)):
        return
    r = dict(ladder=ladder, hypothesis=hyp, description=desc, event=event, unit=unit, test_type="sign only (no per-event test)",
             p_event=np.nan, k_controls=np.nan, p_floor=np.nan, source=source, sig_sign=int(stat > 0), pi_sign=0.5, stat=stat)
    rows.append(r)


def r1_ladder(t):
    if t in R1_ORIGINS:
        return "R1 origins", t
    if t in R1_NESTED:
        return "R1 nested", R1_NESTED[t]
    return None, None


# ---------------- Round 1, 187 eukaryotes, 20 matched controls ----------------
pt = f"{R}/primary/euk"
H = {"H1H2": ("H1H2", "more sequence-novel gene families on the stem (round 1)", "per_transition_H1H2.tsv"),
     "H4": ("H4", "more transcription-factor duplications on the stem", "per_transition_H4.tsv"),
     "H5": ("H5", "more adhesion/ECM and signalling gains+duplications on the stem", "per_transition_H5.tsv"),
     "H12": ("H12", "loss of flagellar modules as predicted per lineage", "per_transition_H12.tsv")}
for key, (hyp, desc, fn) in H.items():
    d = pd.read_csv(f"{pt}/{fn}", sep="\t")
    for _, x in d.iterrows():
        lad, unit = r1_ladder(x.transition_id)
        if lad:
            add_rank(lad, hyp, desc, x.transition_id, unit, u=x.u, k=x.k, source=f"results/primary/euk/{fn}")

# H9 rate of protein evolution (10 controls, floor 1/11 > 0.05)
d = pd.read_csv(f"{R}/h9/h9_transitions.tsv", sep="\t")
for _, x in d.iterrows():
    lad, unit = r1_ladder(x.transition_id)
    if lad:
        add_rank(lad, "H9", "faster protein evolution in the multicellular clade than its sister (r)", x.transition_id, unit,
                 u=x.u_percentile, k=10, source="results/h9/h9_transitions.tsv")

# H11 new domain architectures (20 controls)
d = pd.read_csv(f"{R}/exploratory/euk/architectures/per_transition.tsv", sep="\t")
for _, x in d.iterrows():
    lad, unit = r1_ladder(x.transition_id)
    if lad:
        add_rank(lad, f"H11 {x.statistic}", "more new Pfam architectures on the stem", x.transition_id, unit, u=x.u, k=x.k,
                 source="results/exploratory/euk/architectures/per_transition.tsv")

# H1 co-option index (descriptive; one-sided upper = more co-option than controls)
d = pd.read_csv(f"{R}/exploratory/euk/cooption/per_transition.tsv", sep="\t")
for _, x in d.iterrows():
    lad, unit = r1_ladder(x.transition)
    if lad:
        add_rank(lad, f"H1 co-option {x.version}", "higher co-option index (expanded families pre-existing) than controls",
                 x.transition, unit, u=x.u, k=x.k_controls, source="results/exploratory/euk/cooption/per_transition.tsv")

# H10 disorder / length of new genes vs controls (u_vs_controls; controls per transition not stored -> k from u grid)
d = pd.read_csv(f"{R}/exploratory/euk/protein_traits/tests.tsv", sep="\t")
d = d[d.test.isin(["disorder_sc", "length_sc"]) & (d.subset == "all")]
for _, x in d.iterrows():
    lad, unit = r1_ladder(x.transition)
    if lad and pd.notna(x.u_vs_controls):
        add_rank(lad, f"H10 {x.test}", "new genes more disordered / longer than at control branches", x.transition, unit,
                 u=x.u_vs_controls, k=x.k_controls, source="results/exploratory/euk/protein_traits/tests.tsv")
        add_nominal(lad + " [nominal]", f"H10 {x.test} (within-stem Wilcoxon)",
                    "new genes more disordered / longer than parental genes, within stem, uncalibrated",
                    x.transition, unit, x.p / 2 if x.hl > 0 else 1 - x.p / 2,
                    source="results/exploratory/euk/protein_traits/tests.tsv")

# EukaryOMA novelty per stem (cross-method; u against 16-17 controls)
d = pd.read_csv(f"{R}/eukaryoma_xval/h1h2_descriptive.tsv", sep="\t")
for _, x in d.iterrows():
    lad, unit = r1_ladder(x.transition_id)
    if lad and pd.notna(x.u):
        add_rank(lad, f"H1H2 novelty ({x.method}, 125 shared species)", "more new families on stem, labels stripped",
                 x.transition_id, unit, u=x.u, k=x.k_controls, source="results/eukaryoma_xval/h1h2_descriptive.tsv")

# W7 complexes: duplicated genes enriched in LECA-complex subunits (Fisher vs unicellular branches; nominal)
w7map = {"Metazoa": "T_metazoa", "Fungal hyphae": "T_fungi_hyphae", "Dictyostelia": "T_dictyostelia",
         "Streptophyta": "T_streptophyta", "PX clade": "T_px_filaments", "Volvocine": "T_volvocine_tgv"}
d = pd.read_csv(f"{R}/r2/w7/dup_fraction_stem_vs_control_fisher.tsv", sep="\t")
d = d[d.control == "U"]
for _, x in d.iterrows():
    t = w7map[x.stem]
    # file p is two-sided Fisher; one-sided (enrichment) read-out: halve if OR > 1, else 1 - p/2
    p1 = x.p / 2 if x.OR > 1 else 1 - x.p / 2
    add_nominal("R1 origins", "W7 complex subunits among stem duplications", "duplicates enriched in ancient-complex subunits",
                t, t, p1, source="results/r2/w7/dup_fraction_stem_vs_control_fisher.tsv (two-sided p halved)")

# Round-2 family sets: PGLS sister-clade contrasts (no per-origin test; sign only)
d = pd.read_csv(f"{R}/r2/w4/pgls/sister_contrasts.tsv", sep="\t")
d = d[d.species_set == "primary148"]
setname = {"P2_elm": "R2-P2 motif-binding set", "P3_celldeath": "R2-P3 cell-death set", "P4_tf": "R2-P4 cell-type TF set"}
for _, x in d.iterrows():
    add_sign("R1 origins [sign]", setname[x.set], "multicellular clade above its sister (mean standardised log copy)",
             x.origin, x.origin, x["diff"], source="results/r2/w4/pgls/sister_contrasts.tsv")

# Round-1 reversals: genes gained/duplicated at the preceding transition preferentially lost (Fisher, nominal)
d = pd.read_csv(f"{R}/exploratory/euk/reversals/reversal_tests.tsv", sep="\t")
d = d[(d.group == "gained_or_dup") & (d.role == "reversal") & (~d.pseudo.astype(bool))]
for _, x in d.iterrows():
    add_nominal("R1 reversals", "H12 reversal loses transition genes", "transition-gained/duplicated families lost more often",
                x.reversal, x.reversal, x.p_greater, source="results/exploratory/euk/reversals/reversal_tests.tsv")

# Round-1 prokaryotes (4 origins; descriptive exploratory; 13 or 9 controls)
for fn, col_t, col_k, col_s in (("architectures", "transition_id", "k", "statistic"), ("cooption", "transition", "k_controls",
                                                                                      "version")):
    d = pd.read_csv(f"{R}/exploratory/prok/{fn}/per_transition.tsv", sep="\t")
    for _, x in d.iterrows():
        nm = f"H11 {x[col_s]}" if fn == "architectures" else f"H1 co-option {x[col_s]}"
        add_rank("R1 prokaryotes", nm, "as eukaryotes", x[col_t], x[col_t], u=x.u, k=x[col_k],
                 source=f"results/exploratory/prok/{fn}/per_transition.tsv")

# ---------------- Round 2 dense panels ----------------
ev = pd.read_csv(f"{R}/r2/panels/pooled/final_8panels/event_residuals.tsv", sep="\t")
ev = ev[ev.kind != "sensitivity"].copy()


def panel_ladder(x):
    if x.panel == "R8":
        return "Panels prokaryotes (internal)" if not x.is_tip else "Panels prokaryotes (tips)"
    if x.kind == "reversal":
        return "Panels reversals"
    if x.is_tip:
        return "Panels tips"
    if x.kind in ("origin", "origin_complex"):
        return "Panels origins (internal)"
    if x.kind == "nested":
        return "Panels nested (internal)"
    return None


ev["ladder"] = ev.apply(panel_ladder, axis=1)
evk = ev.set_index("event_id")
SRC = "results/r2/panels/pooled/final_8panels/event_residuals.tsv"
for _, x in ev.iterrows():
    if x.ladder is None:
        continue
    if not x.is_tip:
        add_rank(x.ladder, "Novelty (panels)", "more sequence-novel families than the panel's control branches", x.event_id,
                 x.event_id, p=x.p_perm, k=x.n_matched, source=SRC)
    elif pd.notna(x.tip_p_perm):
        add_rank(x.ladder, "Novelty (panels, tip vs tip controls)", "tip event vs unicellular tip controls", x.event_id,
                 x.event_id, p=x.tip_p_perm, k=x.tip_n_controls, source=SRC)

# cell-death families among gained / duplicated families at the event (Fisher vs control branches; nominal)
d = pd.read_csv(f"{R}/r2/w4_celldeath/panels_event_celldeath.tsv", sep="\t")
for _, x in d.iterrows():
    if x.event_id not in evk.index or evk.at[x.event_id, "ladder"] is None:
        continue
    add_nominal(evk.at[x.event_id, "ladder"], f"Cell-death families among {x.kind} families",
                "cell-death share above control branches", x.event_id, x.event_id, x.p_one_sided,
                source="results/r2/w4_celldeath/panels_event_celldeath.tsv")

# W5 structural novelty: share of sequence-novel families without a structural match, ranked against panel controls
b = pd.read_csv(f"{R}/r2/w5/novelty_summary/per_branch.tsv", sep="\t")
b = b[b.kind == "novel"]
for pan, g in b.groupby("panel"):
    ctrl = g[g.role == "control"].frac_new.to_numpy()
    for _, x in g[g.role == "event"].iterrows():  # tips excluded: controls are internal branches (tip artefact)
        eid = str(x.event).split(":")[0]
        if eid not in evk.index or evk.at[eid, "ladder"] is None:
            continue
        p = ((ctrl >= x.frac_new).sum() + 1) / (len(ctrl) + 1)
        add_rank(evk.at[eid, "ladder"], "Structural novelty (unadjusted)",
                 "more sequence-novel families without a structural match than panel controls", eid, eid, p=p, k=len(ctrl),
                 source="results/r2/w5/novelty_summary/per_branch.tsv (rank read-out)")

# W6 structure-per-sequence change at the event branch (real data, and the shuffled-letter negative control)
d = pd.read_csv(f"{R}/r2/w6/stem_tests.tsv", sep="\t")
d = d[d.analysis == "primary"]
for pipe, nm in (("real", "Structure/sequence change ratio"), ("shuffled3di", "NEG CONTROL: ratio with shuffled 3Di")):
    for _, x in d[d.pipeline == pipe].iterrows():
        if x.event not in evk.index or evk.at[x.event, "ladder"] is None:
            continue
        add_rank(evk.at[x.event, "ladder"], nm, "more 3Di change per amino-acid change at the event than controls",
                 x.event, x.event, p=x.p_one_sided, k=x.n_ctrl, source="results/r2/w6/stem_tests.tsv")

E = pd.DataFrame(rows)

# ---------------- origin-collapsed ladder (round 1) ----------------
coll = []
r1 = E[E.ladder.isin(["R1 origins", "R1 nested"])]
for hyp, g in r1.groupby("hypothesis"):
    if not (g.ladder == "R1 nested").any():
        continue
    for unit, gu in g.groupby("unit"):
        if unit not in R1_ORIGINS:
            continue
        r = dict(ladder="R1 origins collapsed (origin OR nested)", hypothesis=hyp, description=gu.description.iloc[0],
                 event=";".join(gu.event), unit=unit, test_type=gu.test_type.iloc[0], p_event=gu.p_event.min(),
                 k_controls=gu.k_controls.min(), p_floor=gu.p_floor.min(), source=gu.source.iloc[0])
        for c in ("0.05", "0.1", "top"):
            r[f"sig_{c}"] = int(gu[f"sig_{c}"].max())
            r[f"pi_{c}"] = 1 - np.prod(1 - gu[f"pi_{c}"].to_numpy())
        coll.append(r)
E = pd.concat([E, pd.DataFrame(coll)], ignore_index=True)
E.to_csv(out("hypothesis_events.tsv"), sep="\t", index=False)

# ---------------- ladders ----------------
L = []
for (lad, hyp), g in E.groupby(["ladder", "hypothesis"], sort=False):
    crits = [("sign", "sig_sign", "pi_sign")] if g.test_type.iloc[0].startswith("sign") else \
        [("alpha=0.05", "sig_0.05", "pi_0.05"), ("alpha=0.10", "sig_0.1", "pi_0.1")] + \
        ([("top rank", "sig_top", "pi_top")] if g.test_type.iloc[0].startswith("rank") else [])
    for cname, sc, pc in crits:
        pis = g[pc].to_numpy(float)
        sig = g[sc].to_numpy(float)
        n = len(g)
        n_reach = int((pis > 0).sum())
        kobs = int(np.nansum(sig))
        evs = ";".join(g.event[sig == 1])
        for k in range(1, n + 1):
            L.append(dict(ladder=lad, hypothesis=hyp, criterion=cname, test_type=g.test_type.iloc[0], n_events=n,
                          n_events_can_reach=n_reach, observed_k=kobs, expected=float(pis.sum()), tier_k=k,
                          tier_reached=int(kobs >= k), p_ge_k_null=poibin_sf(pis, k),
                          p_observed=poibin_sf(pis, kobs) if kobs > 0 else 1.0, significant_events=evs,
                          description=g.description.iloc[0]))
L = pd.DataFrame(L)
# multiple testing per tier across hypotheses within the ladder family (ladder x criterion)
L["p_tier"] = np.where(L.tier_reached == 1, L.p_ge_k_null, 1.0)
L["holm_tier"] = np.nan
L["bh_tier"] = np.nan
for (lad, cr, k), g in L.groupby(["ladder", "criterion", "tier_k"]):
    L.loc[g.index, "holm_tier"] = holm(g.p_tier.to_numpy())
    L.loc[g.index, "bh_tier"] = bh(g.p_tier.to_numpy())
L.to_csv(out("hypothesis_ladder.tsv"), sep="\t", index=False)

# compact console view: one row per ladder x hypothesis x criterion at the observed tier
C = L[(L.tier_k == L.observed_k.clip(lower=1))].copy()
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 500)
pd.set_option("display.max_colwidth", 60)
print(C[["ladder", "hypothesis", "criterion", "n_events", "n_events_can_reach", "observed_k", "expected", "p_observed",
         "holm_tier", "bh_tier", "significant_events"]].to_string())
