#!/usr/bin/env python
"""Score the HOG-based opisthokont positive controls (C01-C04, C06, C07; focus C03/C06/C07, which
failed in round 1 because of deep-node placement) with the unchanged round-1 rules
(mcell/controls.py) on (i) EukaryOMA HOGs and (ii) Ancestra HOGs, both projected onto the
125-species pruned tree with the same parser and level rule.

Pfam: the pre-registered subset (data/annotations/pfam_prereg/euk), attached to our protein ids;
EukaryOMA genes carry our ids where mapped (MD5 / RBH). For EukaryOMA the PF-HOG >= 50 % rule is
computed over mapped (annotatable) members only (unmapped members have no annotation and are
neither counted for nor against). Diagnostics: for C03 and C07, the branches on which the
PF01049 / PF00907 lineages originate and duplicate in each method.
Outputs: results/eukaryoma_xval/positive_controls.tsv, positive_controls_diagnostics.tsv
"""
import os, sys, warnings

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import analysis as an  # noqa: E402
from mcell import controls as pc  # noqa: E402
from mcell.hogs import load_species_tree  # noqa: E402
from xval_analysis import W, OUT, load_rename  # noqa: E402

CONTROLS = ["C01", "C02", "C03", "C04", "C06", "C07"]


class LiteCache:
    def __init__(self, tree, hogs, genes, events, hits):
        self.tree, self.hogs, self.genes, self.events, self.hits = tree, hogs, genes, events, hits
        self._lineage = None

    @property
    def lineage(self):
        if self._lineage is None:
            self._lineage = an.Lineage(self)
        return self._lineage


def load(name, tree, hits, rename=None):
    d = os.path.join(W, "parsed", name)
    hogs = pd.read_parquet(os.path.join(d, "hogs.parquet"))
    genes = pd.read_parquet(os.path.join(d, "genes.parquet"))
    events = pd.read_parquet(os.path.join(d, "events.parquet"))
    if rename is not None:
        genes["gene"] = genes["gene"].map(lambda g: rename.get(g, g))
        mapped = genes["gene"].isin(set(rename.values())).to_numpy()
        hogs = hogs.copy()
        hogs["n_genes_all"] = hogs["n_genes"]
        hogs["n_genes"] = an.accumulate(hogs, genes, mapped.astype(float), tree).round().astype(np.int64)
    return LiteCache(tree, hogs, genes, events, hits)


def diag(name, c, pfam):
    """origin and duplication branches of the lineages that contain PF-HOG(pfam) rows."""
    L = c.lineage
    rows = an.pf_hogs(c, [pfam])
    fam = set(c.hogs["family_id"].to_numpy()[rows])
    ev = c.events[c.events["family_id"].isin(fam)]
    dup = ev[ev["event"] == "dup"].groupby("branch_child")["n_events"].sum().sort_values(ascending=False)
    gains = ev[ev["event"] == "gain"]["branch_child"].value_counts()
    origins = pd.Series([L.origin_level(int(h)) for h in rows]).value_counts()
    return dict(method=name, pfam=pfam, n_pf_hogs=len(rows), n_families=len(fam),
                pf_hog_origin_levels="; ".join(f"{k}:{v}" for k, v in origins.head(8).items()),
                family_gain_branches="; ".join(f"{k}:{v}" for k, v in gains.head(8).items()),
                dup_branches="; ".join(f"{k}:{v}" for k, v in dup.head(10).items()),
                dups_on_Metazoa=int(dup.get("Metazoa", 0)), dups_on_Choanozoa=int(dup.get("Choanozoa", 0)))


def main():
    sm = pd.read_csv(os.path.join(OUT, "species_map.tsv"), sep="\t")
    tree = load_species_tree(os.path.join(W, "tree_pruned.nwk"), sm["code"].tolist())
    p = an.run_paths("euk")
    hits = an.read_pfam_dir(p["pfam"], sm["code"].tolist())
    species = pd.read_csv(p["species"], sep="\t")
    tr = pd.read_csv(p["transitions"], sep="\t")
    ctrl = pd.read_csv(os.path.join(ROOT, "data/positive_controls.tsv"), sep="\t")
    ctrl = ctrl[(ctrl["run"] == "euk") & ctrl["control_id"].isin(CONTROLS)]
    r1 = pd.read_csv(os.path.join(ROOT, "results/positive_controls/euk/controls.tsv"), sep="\t").set_index("control")
    ren = load_rename()
    out, dg = [], []
    for name, rename in (("eukaryoma", ren), ("ancestra", None)):
        c = load(name, tree, hits, rename)
        ctx = pc.Ctx(cache=c, species=species, transitions=tr, paths=p)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = pc.score_all(ctx, ctrl, pc.EUK_RULES)
        res.insert(0, "method", f"{name} (125 spp, pruned tree)")
        out.append(res)
        for pf in ("PF01049", "PF00907", "PF00362"):
            dg.append(diag(name, c, pf))
    r = r1.loc[CONTROLS].reset_index()[["control", "expected", "confidence", "recovered", "observed", "components"]]
    r.insert(0, "method", "ancestra round 1 (187 spp, full tree)")
    res = pd.concat([r] + out, ignore_index=True)
    res.to_csv(os.path.join(OUT, "positive_controls.tsv"), sep="\t", index=False)
    pd.DataFrame(dg).to_csv(os.path.join(OUT, "positive_controls_diagnostics.tsv"), sep="\t", index=False)
    print(res[["method", "control", "recovered", "components"]].to_string())
    print(pd.DataFrame(dg).to_string())


if __name__ == "__main__":
    main()
