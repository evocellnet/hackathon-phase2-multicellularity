#!/usr/bin/env python
"""Apply the DESIGN.md v3 §2 eligibility rule to data/transitions.tsv and assign analysis sets.

Rule (DESIGN §2): a transition is eligible for confirmatory testing only if M is an internal node with >= 2 sampled
species in the new state, >= 2 sampled species in the unicellular/prior-state sister lineage(s), the state change is not
grey-zone or a reversal. Species flagged low quality (BUSCO C < 50%) do not count towards the >= 2 on either side
(refinement D5, made before any results: one low-quality species must not disqualify an otherwise well-sampled transition).

Analysis sets:
  primary        eligible eukaryotic origins (unicellular -> simple/complex) — the Holm family of DESIGN §4.6
  secondary      eligible nested simple -> complex eukaryotic transitions, and eligible prokaryotic transitions
                 (same statistics, reported separately, not in the Holm family)
  exploratory    everything else with data (single-species, grey-zone, low quality, one sister species)
  reversal       reversals to unicellularity (exploratory reversal analysis)
  excluded       no data

PI decision D4 (RESEARCH_TREE.md): the volvocine transition is coded at the TGV stem (M = TGV_clade,
P = Chlamydomonas_TGV), because every sampled colonial species descends from TGV_clade; single-origin
(Herron & Michod 2008) vs two-origin (Lindsey et al. 2021) interpretations are reported as a caveat. The two-origin rows
(Tetrabaena tip; Gonium_Volvocaceae) are kept as exploratory.
"""
from __future__ import annotations

import csv
import sys

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"


def main() -> int:
    path = f"{P}/data/transitions.tsv"
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    fields = list(rows[0].keys())
    traits = {}
    for fn in ("species_final.tsv", "species_final_prok.tsv"):
        for r in csv.DictReader(open(f"{P}/data/{fn}"), delimiter="\t"):
            traits[r["code"]] = r

    # PI decision D4: add the volvocine TGV-stem transition if absent
    if not any(r["transition_id"] == "T_volvocine_tgv" for r in rows):
        gv = next(r for r in rows if r["transition_id"] == "T_volvocales_gv")
        tb = next(r for r in rows if r["transition_id"] == "T_volvocales_tetrabaena")
        new = dict(gv)
        mc = sorted(set(gv["mc_species"].split(",")) | set(tb["mc_species"].split(",")))
        new.update(transition_id="T_volvocine_tgv", label="Volvocine colonial lineage (TGV stem)", M="TGV_clade",
                   M_kind="internal", P="Chlamydomonas_TGV", sister_clade_nodes="Chlamydomonas",
                   mc_species=",".join(mc), n_mc_species=str(len(mc)),
                   sister_species=gv["prior_state_ref_species"], n_sister_species=gv["n_prior_state_ref"],
                   notes="PI decision D4: TGV stem; shared changes of all sampled colonial volvocines. Under the two-origin "
                         "model (Lindsey et al. 2021) this stem is pre-transition (predisposition); under the single-origin "
                         "model (Herron & Michod 2008) it is the transition stem.")
        rows.append(new)

    # PI decision D7 (round-3 audit F1, F2): path stems. A path stem is the chain of consecutive edges from P to M;
    # events and branch lengths are summed over its edges (DESIGN §4.6 v3.1).
    if "stem_path" not in fields:
        fields += ["stem_path"]
    for r in rows:
        r.setdefault("stem_path", "")
        r["stem_path"] = f"{r['P']}>{r['M']}" if r["M"] and r["P"] and "," not in r["M"] else ""
        if r["transition_id"] == "T_volvocine_tgv":
            gv = next(x for x in rows if x["transition_id"] == "T_volvocales_gv")
            r.update(M="Gonium_Volvocaceae", P="Chlamydomonas_TGV", stem_path="Chlamydomonas_TGV>TGV_clade>Gonium_Volvocaceae",
                     mc_species=gv["mc_species"], n_mc_species=gv["n_mc_species"], flags="path_stem",
                     label="Volvocine colonial lineage (Gonium+Volvocaceae; path stem from Chlamydomonas_TGV)",
                     notes="D7 (audit F1): Lindsey et al. 2021 find colonial volvocines non-monophyletic, so TGV_clade is inferred "
                           "unicellular and exists only through Tetrabaena (BUSCO 32%); the stem is the 2-edge path "
                           "Chlamydomonas_TGV>TGV_clade>Gonium_Volvocaceae. Tetrabaena = separate exploratory origin.")
        if r["transition_id"] == "T_fungi_hyphae":
            r.update(M="Amastigomycota", stem_path="Eumycota>Blastocladiomycota_Amastigomycota>Amastigomycota",
                     flags="path_stem",
                     notes="D7 (audit F2): Kiss et al. 2019 place hyphae at Zoopagomycota+Dikarya (=Amastigomycota here) and treat "
                           "Blastocladiomycota as mostly unicellular; the stem is the 2-edge path from Eumycota. Allomyces "
                           "(hyphal blastoclad) is a separate exploratory origin.")
            mc = r["mc_species"].split(",")
            r["mc_species"] = ",".join(c for c in mc if c not in ("Allomac", "Allomyc", "Blasemer"))
            r["n_mc_species"] = str(len(r["mc_species"].split(",")))
    if "analysis_set" not in fields:
        fields += ["analysis_set", "eligibility_reason"]
    for r in rows:
        reasons = []
        n_mc = int(r["n_mc_species"] or 0)
        n_s = int(r["n_sister_species"] or 0)
        if n_mc == 0:
            r["analysis_set"], r["primary_eligible"] = "excluded", "0"
            r["eligibility_reason"] = "no data"
            continue
        if r["to_state"] == "reversal":
            r["analysis_set"], r["primary_eligible"] = "reversal", "0"
            r["eligibility_reason"] = "reversal"
            continue
        if r["to_state"] == "grey":
            reasons.append("grey-zone state")
        if "," in r["M"]:
            reasons.append("M is not a single clade")
        if r.get("M_kind") == "leaf" or n_mc < 2:
            reasons.append("single multicellular species")
        if n_s < 2:
            reasons.append("<2 sister species")
        def good(field):
            return [c for c in r[field].split(",") if c and traits.get(c, {}).get("low_quality", "0") != "1"]
        if n_mc >= 2 and len(good("mc_species")) < 2:
            reasons.append("<2 good-quality MC species")
        if n_s >= 2 and len(good("sister_species")) < 2:
            reasons.append("<2 good-quality sister species")
        if r["transition_id"] in ("T_volvocales_gv", "T_volvocales_tetrabaena"):
            reasons.append("superseded by T_volvocine_tgv (D4)")
        if reasons:
            r["analysis_set"], r["primary_eligible"] = "exploratory", "0"
        elif r["kingdom"] == "euk" and r["from_state"] == "unicellular":
            r["analysis_set"], r["primary_eligible"] = "primary", "1"
        else:
            r["analysis_set"], r["primary_eligible"] = "secondary", "0"
        r["eligibility_reason"] = "; ".join(reasons) if reasons else "meets DESIGN §2 rule"
        if r["transition_id"] == "T_px_filaments" and not reasons:
            r["eligibility_reason"] += "; caveat: filamentous state of the PX ancestor is inferred (evo_biology.md review)"

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for s in ("primary", "secondary", "exploratory", "reversal", "excluded"):
        ids = [r["transition_id"] for r in rows if r["analysis_set"] == s]
        print(f"{s} ({len(ids)}): {', '.join(ids)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
