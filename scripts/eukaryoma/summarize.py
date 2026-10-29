#!/usr/bin/env python
"""One row per transition stem combining counts (both methods) and family-level agreement
(both directions, gains and duplications). Output: results/eukaryoma_xval/per_transition_summary.tsv
and control_baseline_summary.tsv (same columns, summed over exact control branches)."""
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results/eukaryoma_xval")


def main():
    C = pd.read_csv(os.path.join(OUT, "counts_per_branch.tsv"), sep="\t", keep_default_na=False)
    A = pd.read_csv(os.path.join(OUT, "agreement_per_branch.tsv"), sep="\t", keep_default_na=False)
    D = pd.read_csv(os.path.join(OUT, "agreement_stems_depth_adjusted.tsv"), sep="\t")
    bm = pd.read_csv(os.path.join(OUT, "branch_mapping.tsv"), sep="\t", keep_default_na=False)
    rows = []

    def agree(sub, d, e):
        s = sub[(sub["direction"] == d) & (sub["event"] == e)]
        if s.empty:
            return {}
        t = s[["same", "deeper", "shallower", "other", "unmapped"]].sum()
        nm = t.sum() - t["unmapped"]
        tag = ("A2E" if d.startswith("ancestra") else "E2A") + "_" + e
        return {f"{tag}_n": int(t.sum()), f"{tag}_same": int(t["same"]), f"{tag}_deeper": int(t["deeper"]),
                f"{tag}_shallower": int(t["shallower"]), f"{tag}_other": int(t["other"]),
                f"{tag}_unmapped": int(t["unmapped"]), f"{tag}_frac_same": round(t["same"] / nm, 3) if nm else None}

    for role, ids in (("stem", None), ("control", None)):
        sub_b = bm[bm["role"] == role]
        for _, b in sub_b.iterrows():
            if b["status"] == "missing":
                if role == "stem":
                    rows.append(dict(branch_id=b["branch_id"], role=role, tier=b["tier"], status="missing",
                                     note=b["reason"]))
                continue
            r = dict(branch_id=b["branch_id"], role=role, tier=b["tier"], analysis_set=b["analysis_set"],
                     status=b["status"], pruned_stem=b["pruned_edges"], extra_orig_edges=b["extra_orig_edges"],
                     note=b["reason"])
            c = C[(C["role"] == role) & (C["branch_id"] == b["branch_id"])].set_index("method")
            for m, tag in (("ancestra", "A"), ("eukaryoma", "E"), ("ancestra_hints", "Ahints")):
                if m in c.index:
                    for k in ("gains_supported", "losses", "dups"):
                        r[f"{tag}_{k}"] = int(c.loc[m, k])
            s = A[(A["role"] == role) & (A["branch_id"] == b["branch_id"])]
            for d in ("ancestra->eukaryoma", "eukaryoma->ancestra"):
                for e in ("gain", "dup"):
                    r.update(agree(s, d, e))
            if role == "stem":
                for _, x in D[D["branch_id"] == b["branch_id"]].iterrows():
                    tag = ("A2E" if x["direction"].startswith("ancestra") else "E2A") + "_" + x["event"]
                    r[f"{tag}_exp_frac_same_depth_adj"] = x["expected_frac_same"]
            rows.append(r)
    T = pd.DataFrame(rows)
    T[T["role"] == "stem"].drop(columns=["role"]).to_csv(os.path.join(OUT, "per_transition_summary.tsv"), sep="\t", index=False)
    T[T["role"] == "control"].drop(columns=["role", "tier", "analysis_set"]).to_csv(
        os.path.join(OUT, "per_control_summary.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    cols = ["branch_id", "status", "A_gains_supported", "E_gains_supported", "A_losses", "E_losses", "A_dups", "E_dups",
            "A2E_gain_frac_same", "E2A_gain_frac_same", "A2E_dup_frac_same", "E2A_dup_frac_same"]
    print(T[T["role"] == "stem"][[c for c in cols if c in T.columns]].to_string())


if __name__ == "__main__":
    main()
