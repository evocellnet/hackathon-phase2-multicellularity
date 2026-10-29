#!/usr/bin/env python
"""Map every round-1 branch of interest onto the pruned (EukaryOMA-overlap) tree.

Rows: the primary and secondary eukaryotic transition stems (data/transitions.tsv,
analysis_tier in {primary, secondary}; path stems via stem_path), the 30 round-1 control branches
(results/primary/euk/branch_values.tsv, in_U), and all other non-root branches of the original
tree (role 'other', used as a broader baseline). Status per xval_lib.branch_map:
exact = the pruned branch/path covers exactly the original edge(s); approx = it covers extra
original edges (e.g. the crown's first split or the parent's other lineages are unsampled);
missing = no sampled species below the crown.
Output: results/eukaryoma_xval/branch_mapping.tsv
"""
import os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd  # noqa: E402
import xval_lib as xl  # noqa: E402


def main():
    orig = xl.load_tree(os.path.join(ROOT, "runs/euk_fast/stages/input/guide_tree.nwk"))
    pr = xl.load_tree(os.path.join(ROOT, "runs/eukaryoma_xval/tree_pruned.nwk"))
    tr = pd.read_csv(os.path.join(ROOT, "data/transitions.tsv"), sep="\t", keep_default_na=False)
    tr = tr[(tr["kingdom"] == "euk") & tr["analysis_tier"].isin(["primary", "secondary"])]
    bv = pd.read_csv(os.path.join(ROOT, "results/primary/euk/branch_values.tsv"), sep="\t", keep_default_na=False)
    rows = []
    for _, r in tr.iterrows():
        path = r["stem_path"].split(">") if r["stem_path"] else [r["P"], r["M"]]
        top, bottom = path[0], path[-1]
        m = xl.branch_map(orig, pr, top, bottom)
        rows.append(dict(branch_id=r["transition_id"], role="stem", tier=r["analysis_tier"],
                         analysis_set=r["analysis_set"], orig_top=top, orig_bottom=bottom,
                         n_orig_edges=len(path) - 1, **m))
    ctrl = set(bv.loc[bv["in_U"].astype(str) == "True", "branch"])
    for nd in orig.traverse("preorder"):
        if nd.up is None:
            continue
        role = "control" if nd.name in ctrl else "other"
        m = xl.branch_map(orig, pr, nd.up.name, nd.name)
        rows.append(dict(branch_id=nd.name, role=role, tier="", analysis_set="", orig_top=nd.up.name,
                         orig_bottom=nd.name, n_orig_edges=1, is_leaf=nd.is_leaf(), **m))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(ROOT, "results/eukaryoma_xval/branch_mapping.tsv"), sep="\t", index=False)
    s = out[out["role"] != "other"]
    print(s.groupby(["role", "status"]).size().to_string())
    print(out[out["role"] == "stem"][["branch_id", "tier", "status", "top", "bottom", "extra_orig_edges", "reason"]]
          .to_string())


if __name__ == "__main__":
    main()
