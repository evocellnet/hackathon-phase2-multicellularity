"""Apply the P3_rule.txt lineage-restriction rule to the frozen sets and write inputs for the PGLS script.
Outputs (results/r2/w4/pgls_input/): set_families.tsv, counts_sets.tsv, species.tsv
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import P, OUT, out, load_transitions, primary_origins  # noqa: E402

SETS = {"P2_elm": "elm_motif_binding.tsv", "P3_celldeath": "cell_death.tsv", "P4_tf": "tf_celltype.tsv"}


def main():
    counts = pd.read_csv(f"{OUT}/matrix/counts.tsv.gz", sep="\t", index_col=0)
    sp = pd.read_csv(f"{OUT}/matrix/species.tsv", sep="\t")
    sp["primary"] = ((sp.data_type == "genome") & (sp.busco_C >= 50)).astype(int)
    assert sp.primary.sum() == 148
    tr = load_transitions()
    po = primary_origins(tr)
    assert len(po) == 7, po.transition_id.tolist()
    # complex D-or-T coding (P3_rule item 6)
    sp["complex_T"] = (sp.mc_state == "complex").astype(int)
    sp["complex_DT"] = ((sp.mc_state == "complex") | sp.code.isin(["Volvcar", "Volvafr", "Astrgub", "Pleosta"])).astype(int)
    for _, r in po.iterrows():
        sp[f"mc_{r.transition_id}"] = sp.code.isin(r.mc_list).astype(int)
        sp[f"sis_{r.transition_id}"] = sp.code.isin(r.sis_list).astype(int)
    prim_codes = sp.code[sp.primary == 1]
    rows = []
    for setname, fn in SETS.items():
        s = pd.read_csv(f"{P}/data/r2_sets/{fn}", sep="\t")
        if "set" not in s.columns:
            s["set"] = "primary"
        for _, r in s.iterrows():
            f = r.pfam
            inm = f in counts.columns
            c = counts[f] if inm else pd.Series(0, index=counts.index)
            n_orig = sum(int((c.loc[[x for x in o.mc_list if x in c.index]] > 0).any()) for _, o in po.iterrows())
            frac_prim = float((c.loc[prim_codes] > 0).mean())
            rows.append(dict(set=setname, pfam=f, name=r["name"], subset=r["set"], in_matrix=int(inm),
                             n_species_present_all=int((c > 0).sum()), frac_present_primary=round(frac_prim, 4),
                             n_primary_origins_mc_present=n_orig,
                             widespread=int((n_orig >= 3) or (frac_prim >= 0.5)),
                             absent_primary=int(frac_prim == 0)))
    fam = pd.DataFrame(rows)
    allf = sorted(set(fam.pfam))
    cs = counts.reindex(columns=allf, fill_value=0)
    fam.to_csv(out("pgls_input/set_families.tsv"), sep="\t", index=False)
    cs.to_csv(out("pgls_input/counts_sets.tsv"), sep="\t")
    sp.to_csv(out("pgls_input/species.tsv"), sep="\t", index=False)
    print(fam.groupby(["set", "subset"]).agg(n=("pfam", "size"), widespread=("widespread", "sum"),
                                             absent=("absent_primary", "sum"), in_matrix=("in_matrix", "sum")))


if __name__ == "__main__":
    main()
