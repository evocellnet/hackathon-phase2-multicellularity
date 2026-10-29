#!/usr/bin/env python
"""Errata item 3: which of the 16 eukaryotic positive controls test Ancestra's reconstruction?

Each control's frozen pass rule (data/positive_controls.tsv) and its scoring function
(mcell/controls.py) were read; the class below is recorded per control with the clause of the rule
that determines it. Classes:
  event_placement   the rule requires Ancestra to place an event on a specific branch: a HOG origin
                    on a stem (gain), an origin at/above a node (i.e. NOT gained on the stem), a loss
                    or a duplication on a stem.  (C15 is the one null control: "no expansion" of
                    stem duplications, which a method that under-calls duplications also passes.)
  hog_presence      the rule only needs some Ancestra HOG carrying the domain to be present at an
                    ancestral node with members on both sides; no event placement.
  species_counts    per-species Pfam/architecture counts on the proteomes; Ancestra is not used.
A regex cross-check on the pass-rule text guards against transcription errors.
Output: results/round1_errata/positive_controls_classified.tsv (+ printed tallies)
"""
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CLASS = {
    "C01": ("event_placement", "gain on stem", "origin of PF00110 HOG on stem(Metazoa)"),
    "C02": ("event_placement", "gain on stem", "origin of PF00019 HOG on stem(Metazoa)"),
    "C03": ("event_placement", "gain on stem (deep node)", "origin of PF01049 HOG on stem(Metazoa)"),
    "C04": ("event_placement", "origin at/above P (deep node)", "PF00028 HOG origin at or above Choanozoa, present at Choanozoa"),
    "C05": ("hog_presence", "presence at P + species counts", "(a) SH2+PK_Tyr architecture per species; (b) PF00102 HOG present at Choanozoa"),
    "C06": ("event_placement", "origin pre-Opisthokonta + 2 stem losses (deep node)", "PF00362 HOG origin >= Apusomonadida+Opisthokonta; lost on Fungi and Choanoflagellata stems"),
    "C07": ("event_placement", "presence + stem loss + stem duplication (deep node)", "PF00907 HOG present at Holozoa; lost on Choanoflagellata stem; duplicated on Metazoa stem"),
    "C08": ("hog_presence", "presence at P", "DSL and NOD/NODP HOGs present at Choanozoa with metazoan + choanoflagellate members"),
    "C09": ("hog_presence", "presence at P", "PF05920 HOG present at each of three P nodes with a member in M"),
    "C10": ("event_placement", "gain on stem", "origin of PF03514 HOG on stem(Anydrophyta) + outside-hit check"),
    "C11": ("event_placement", "origin at/above P (no gain on stem)", "PF01858/PF01857 HOG origin at/above P, present at P"),
    "C12": ("species_counts", "per-species counts", "PF12499 proteins Volvox vs Chlamydomonas"),
    "C13": ("event_placement", "origin at crown", "Hex1 HOG origin = crown(Pezizomycotina) (branch i)"),
    "C14": ("species_counts", "per-species counts", "PF01185 median per species, Agaricomycetes vs other Dikarya"),
    "C15": ("event_placement", "null: no expansion of stem duplications", "share of stem-duplicated HOGs in toolkit classes on the fungal-hyphae stem"),
    "C16": ("species_counts", "per-species counts", "PF00109+PF00698 architecture per species, Dictyostelia vs other Amoebozoa"),
}
DEEP_HOG_PLACEMENT = {"C03", "C06", "C07"}   # high-confidence failures diagnosed as Ancestra inference errors
REVIEWER_NON_ANCESTRA = {"C05", "C09", "C12", "C14", "C16"}


def main():
    rules = pd.read_csv(os.path.join(ROOT, "data/positive_controls.tsv"), sep="\t")
    rules = rules[rules.run == "euk"].set_index("control_id")
    res = pd.read_csv(os.path.join(ROOT, "results/positive_controls/euk/controls.tsv"), sep="\t").set_index("control")
    rows = []
    for cid in res.index:
        cls, kind, what = CLASS[cid]
        rule = str(rules.loc[cid, "pass_rule"])
        uses_hog = bool(re.search(r"HOG", rule))
        uses_event = bool(re.search(r"origin|loss|duplicat", rule))
        # cross-check: species_counts <=> rule never mentions a HOG; event_placement => event words present
        assert (cls == "species_counts") == (not uses_hog), cid
        if cls == "event_placement":
            assert uses_event or cid == "C15", cid
        rows.append(dict(control=cid, confidence=res.loc[cid, "confidence"], recovered=res.loc[cid, "recovered"],
                         uses_ancestra=cls != "species_counts", test_class=cls, what_is_tested=kind,
                         clause=what, deep_node_hog_placement=cid in DEEP_HOG_PLACEMENT,
                         in_reviewer_non_ancestra_list=cid in REVIEWER_NON_ANCESTRA,
                         pass_rule=rule))
    df = pd.DataFrame(rows)
    out = os.path.join(ROOT, "results/round1_errata/positive_controls_classified.tsv")
    df.to_csv(out, sep="\t", index=False)
    rec = df.recovered.eq("yes")
    print(df[["control", "confidence", "recovered", "test_class", "what_is_tested"]].to_string(index=False))
    print("\nrecovered / total by class:")
    print(df.assign(rec=rec).groupby("test_class").rec.agg(["sum", "count"]).to_string())
    a = df.uses_ancestra
    print(f"\nuses Ancestra at all: {int((rec & a).sum())}/{int(a.sum())}")
    ev = df.test_class.eq("event_placement")
    ev_aff = ev & df.control.ne("C15")
    print(f"event placement: {int((rec & ev).sum())}/{int(ev.sum())}; excluding the null control C15: "
          f"{int((rec & ev_aff).sum())}/{int(ev_aff.sum())}")
    hc = df.confidence.eq("high")
    print(f"high-confidence event placement (excl. C15): {int((rec & ev_aff & hc).sum())}/{int((ev_aff & hc).sum())}")
    print(f"deep-node HOG placement C03/C06/C07: {int((rec & df.deep_node_hog_placement).sum())}/3")
    rv = df.in_reviewer_non_ancestra_list
    print("reviewer's list:", df.loc[rv, ["control", "recovered", "test_class"]].to_dict("records"))
    print(f"not using Ancestra: {int((rec & ~a).sum())}/{int((~a).sum())}")


if __name__ == "__main__":
    main()
