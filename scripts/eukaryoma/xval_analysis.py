#!/usr/bin/env python
"""Event-by-event cross-validation of Ancestra stem events against EukaryOMA (task steps 3a/3b).

Inputs (prepare_inputs.py, branch_mapping.py, map_proteins.py, rbh_map.py):
  runs/eukaryoma_xval/parsed/{ancestra,eukaryoma}/  mcell.hogs tables on the pruned tree
  results/eukaryoma_xval/branch_mapping.tsv          stems / controls on the pruned tree
  runs/eukaryoma_xval/protein_map.tsv.gz (+ protein_map_rbh.tsv.gz if present)

Per branch (stem path or control edge, as represented on the pruned tree; path stems sum their
edges) and method:
  gains_all, gains_supported (gain events whose family has >= 2 species), losses, dups (k-1 summed
  over duplication events inside existing lineages; duplications at a family's origin are not
  counted), dup_events.
Family-level agreement, both directions (src -> dst):
  gain  : the gained family's member genes (mapped to dst by protein id) -> dst ancestral gene at
          the gain node c (majority over mapped members); its dst origin o:
          same       o is on the stem path (strictly below the path top T)
          deeper     o is at/above T: dst says the family is older (for src = Ancestra: Ancestra
                     split an older family -> spurious gain)
          shallower  no mapped member has a dst ancestral gene at c (dst origin inside the clade or
                     species-specific)
          unmapped   no member maps to dst
  dup   : each child copy (ancestral gene at c) -> majority dst ancestral gene at c;
          same       >= 2 copies map to distinct dst genes at c that share one dst ancestor at T
          deeper     >= 2 distinct dst genes at c, all with different dst ancestors at T (paralogs
                     already separate before the stem)
          shallower  all copies map to one dst gene at c, or none is present at c (dst: one gene at c)
          other      >= 2 distinct dst genes at c, some originating on the path (dst calls a gain)
          unmapped
Outputs in results/eukaryoma_xval/: counts_per_branch.tsv, agreement_per_branch.tsv,
agreement_pooled.tsv, agreement_events.tsv.gz
"""
import os, sys, time
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

from mcell.hogs import HogResult, load_species_tree  # noqa: E402

W = os.path.join(ROOT, "runs/eukaryoma_xval")
OUT = os.path.join(ROOT, "results/eukaryoma_xval")
METHODS = ("ancestra", "eukaryoma")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


class M:
    """One parsed method on the pruned tree with fast lineage navigation."""

    def __init__(self, name, tree, rename=None, d=None):
        d = d or os.path.join(W, "parsed", name)
        self.name = name
        self.T = tree
        self.hogs = pd.read_parquet(os.path.join(d, "hogs.parquet"))
        self.genes = pd.read_parquet(os.path.join(d, "genes.parquet"))
        self.events = pd.read_parquet(os.path.join(d, "events.parquet"))
        self.branches = pd.read_csv(os.path.join(d, "branch_summary.tsv"), sep="\t", keep_default_na=False)
        self.gene_row_orig = None
        if rename is not None:      # EukaryOMA ids -> our protein ids where mapped
            self.gene_row_orig = pd.Series(np.arange(len(self.genes)), index=self.genes["gene"].to_numpy())
            self.genes["gene"] = self.genes["gene"].map(lambda g: rename.get(g, g))
        self.res = HogResult(tree=tree, branches=self.branches, hogs=self.hogs, genes=self.genes, events=self.events)
        self.level = self.hogs["level_idx"].to_numpy()
        self.parent = self.hogs["parent_idx"].to_numpy()
        self.hog_row = pd.Series(np.arange(len(self.hogs)), index=self.hogs["hog_id"].to_numpy())
        self.gene_row = pd.Series(np.arange(len(self.genes)), index=self.genes["gene"].to_numpy())
        self.gene_parent = self.genes["parent_idx"].to_numpy()
        self.gene_sp = self.genes["species"].astype(str).to_numpy()
        # origin row of every HOG row
        o = np.arange(len(self.hogs))
        for _ in range(max(tree.depth) + 2):
            p = self.parent[o]
            m = p >= 0
            if not m.any():
                break
            o = np.where(m, p, o)
        self.origin = o
        self.fam_genes = self.genes.groupby("family_id").indices
        self._at = {}

    def anc_at(self, X):
        """per gene row: HOG row of its ancestral gene at internal node X (-1 if none)."""
        if X not in self._at:
            self._at[X] = self.res._gene_hog_at(self.T.names[X])
        return self._at[X]

    def anc_key(self, X, rows):
        """ancestral-gene key at node X for gene rows: HOG row (internal X) or -(gene_row+2) (leaf X)."""
        rows = np.asarray(rows, dtype=np.int64)
        if self.T.is_leaf[X]:
            return np.where(self.gene_sp[rows] == self.T.names[X], -(rows + 2), -1)
        return self.anc_at(X)[rows]

    def key_ancestor_at(self, key, T):
        """HOG row of the ancestor at node T of an ancestral-gene key (None if the lineage is not at T)."""
        h = self.gene_parent[-key - 2] if key < -1 else key
        dT = self.T.depth[T]
        while h >= 0:
            L = self.level[h]
            if L == T:
                return int(h)
            if self.T.depth[L] < dT:
                return None
            h = self.parent[h]
        return None

    def key_origin_level(self, key):
        if key < -1:
            g = -key - 2
            p = self.gene_parent[g]
            return self.T.index[self.gene_sp[g]] if p < 0 else int(self.level[self.origin[p]])
        return int(self.level[self.origin[key]])

    def members_at(self, X, hog_ids):
        """gene rows of each ancestral gene (hog id at internal node X, or a gene id at leaf X)."""
        if self.T.is_leaf[X]:
            gr = self.gene_row_orig if self.gene_row_orig is not None else self.gene_row
            return {h: [int(gr[h])] if h in gr.index else [] for h in hog_ids}
        if ("mem", X) not in self._at:
            a = self.anc_at(X)
            idx = np.flatnonzero(a >= 0)
            self._at[("mem", X)] = pd.Series(idx).groupby(a[idx]).indices, idx
        grp, idx = self._at[("mem", X)]
        out = {}
        for h in hog_ids:
            r = self.hog_row.get(h)
            out[h] = idx[grp[r]] if r is not None and r in grp else []
        return out


def mapped_rows(src, dst, rows):
    names = src.genes["gene"].to_numpy()[rows]
    r = dst.gene_row.reindex(names).dropna().astype(np.int64).to_numpy()
    return r


def classify_gain(src, dst, fam, c, path_nodes, top):
    rows = src.fam_genes.get(fam, [])
    mr = mapped_rows(src, dst, rows)
    out = dict(n_members=len(rows), n_mapped=len(mr))
    if not len(mr):
        return dict(out, cls="unmapped")
    keys = dst.anc_key(c, mr)
    v = keys[keys != -1]
    if not len(v):
        un = float(dst.genes["species_specific"].to_numpy()[mr].mean())
        return dict(out, cls="shallower", purity=np.nan, dst_unassigned_frac=un)
    e, n = Counter(v.tolist()).most_common(1)[0]
    o = dst.key_origin_level(e)
    purity = n / len(mr)
    T = src.T
    if o in path_nodes:
        return dict(out, cls="same", purity=purity)
    if T.is_anc(o, top):
        return dict(out, cls="deeper", purity=purity, nodes_deeper=T.depth[top] - T.depth[o] + 1)
    return dict(out, cls="shallower", purity=purity)   # cannot happen (e is at c); kept for safety


def classify_dup(src, dst, child_ids, c, top):
    mem = src.members_at(c, child_ids)
    es, n_mapped = [], 0
    for h in child_ids:
        mr = mapped_rows(src, dst, mem[h])
        n_mapped += len(mr)
        if not len(mr):
            continue
        keys = dst.anc_key(c, mr)
        v = keys[keys != -1]
        if len(v):
            es.append(Counter(v.tolist()).most_common(1)[0][0])
    out = dict(n_copies=len(child_ids), n_mapped=n_mapped)
    if n_mapped == 0:
        return dict(out, cls="unmapped")
    D = sorted(set(es))
    out["n_dst_genes_at_c"] = len(D)
    if len(D) <= 1:
        return dict(out, cls="shallower")
    anc = [dst.key_ancestor_at(e, top) for e in D]
    a = [x for x in anc if x is not None]
    if len(a) > len(set(a)):
        return dict(out, cls="same")
    if len(a) == len(anc):
        return dict(out, cls="deeper")
    return dict(out, cls="other")


def load_rename():
    pm = pd.read_csv(os.path.join(W, "protein_map.tsv.gz"), sep="\t")
    ren = dict(zip(pm["oma_id"], pm["our_id"]))
    p2 = os.path.join(W, "protein_map_rbh.tsv.gz")
    if os.path.exists(p2):
        rb = pd.read_csv(p2, sep="\t")
        used = set(ren.values())
        for o, u in zip(rb["oma_id"], rb["our_id"]):
            if o not in ren and u not in used:
                ren[o] = u
                used.add(u)
    return ren


def main():
    sm = pd.read_csv(os.path.join(OUT, "species_map.tsv"), sep="\t")
    tree = load_species_tree(os.path.join(W, "tree_pruned.nwk"), sm["code"].tolist())
    bm = pd.read_csv(os.path.join(OUT, "branch_mapping.tsv"), sep="\t", keep_default_na=False)
    ren = load_rename()
    log(f"protein map: {len(ren):,} EukaryOMA ids renamed to our ids")
    meth = {"ancestra": M("ancestra", tree), "eukaryoma": M("eukaryoma", tree, ren)}
    cmeth = dict(meth)
    if os.path.exists(os.path.join(W, "parsed", "ancestra_hints")):
        cmeth["ancestra_hints"] = M("ancestra_hints", tree)
    for m in meth.values():
        log(f"{m.name}: {len(m.genes):,} genes, {len(m.hogs):,} HOGs, {len(m.events):,} events")
    # branches to evaluate: stems + controls (exact/approx) + other exact internal branches (baseline 2)
    sel = bm[(bm["status"].isin(["exact", "approx"])) &
             ((bm["role"] != "other") | ((bm["status"] == "exact") & (bm["is_leaf"].astype(str) == "False")))]
    stem_bottoms = set(bm.loc[bm["role"] == "stem", "orig_bottom"]) | set(bm.loc[bm["role"] == "stem", "bottom"])
    sel = sel[~((sel["role"] == "other") & sel["branch_id"].isin(stem_bottoms))]
    sel = sel.drop_duplicates(subset=["role", "branch_id"])
    counts, ev_rows = [], []
    for _, b in sel.iterrows():
        edges = [tuple(e.split(">")) for e in b["pruned_edges"].split(";")]
        top = tree.index[b["top"]]
        path_nodes = {tree.index[c] for _, c in edges}
        for mn, m in cmeth.items():
            br = m.branches.set_index("branch_child")
            ev = m.events
            kids = [c for _, c in edges]
            e = ev[ev["branch_child"].isin(kids)]
            g = e[e["event"] == "gain"]
            d = e[(e["event"] == "dup") & ~e["at_origin"]]
            counts.append(dict(branch_id=b["branch_id"], role=b["role"], tier=b["tier"], status=b["status"],
                               method=mn, n_edges=len(edges), top=b["top"], bottom=b["bottom"],
                               gains_all=int(br.loc[kids, "gains"].sum()),
                               gains_supported=int((g["n_species"] >= 2).sum()),
                               losses=int(br.loc[kids, "losses"].sum()), dups=int(d["n_events"].sum()),
                               dup_events=len(d), depth=tree.depth[tree.index[b["bottom"]]],
                               n_desc=int(tree.clade_mask[tree.index[b["bottom"]]].bit_count()), n_parent=int(br.loc[kids[0], "n_parent"]),
                               n_child=int(br.loc[kids[-1], "n_child"])))
        for sn, dn in (("ancestra", "eukaryoma"), ("eukaryoma", "ancestra")):
            src, dst = meth[sn], meth[dn]
            for _, c in edges:
                ci = tree.index[c]
                e = src.events[src.events["branch_child"] == c]
                for r in e[(e["event"] == "gain") & (e["n_species"] >= 2)].itertuples():
                    res = classify_gain(src, dst, r.family_id, ci, path_nodes, top)
                    ev_rows.append(dict(branch_id=b["branch_id"], role=b["role"], tier=b["tier"], status=b["status"],
                                        direction=f"{sn}->{dn}", event="gain", edge_child=c, hog_id=r.hog_id,
                                        family_id=r.family_id, n_species=r.n_species, **res))
                for r in e[(e["event"] == "dup") & ~e["at_origin"]].itertuples():
                    ids = [x for x in r.child_hog_ids.split(";") if x]
                    res = classify_dup(src, dst, ids, ci, top)
                    ev_rows.append(dict(branch_id=b["branch_id"], role=b["role"], tier=b["tier"], status=b["status"],
                                        direction=f"{sn}->{dn}", event="dup", edge_child=c, hog_id=r.hog_id,
                                        family_id=r.family_id, n_species=r.n_species, **res))
        log(f"{b['role']:8s} {b['branch_id']}: done ({len(ev_rows):,} classified events so far)")
    C = pd.DataFrame(counts)
    C.to_csv(os.path.join(OUT, "counts_per_branch.tsv"), sep="\t", index=False)
    E = pd.DataFrame(ev_rows)
    if "dst_unassigned_frac" not in E.columns:
        E["dst_unassigned_frac"] = np.nan
    E.to_csv(os.path.join(OUT, "agreement_events.tsv.gz"), sep="\t", index=False)
    classes = ["same", "deeper", "shallower", "other", "unmapped"]
    keys = ["branch_id", "role", "tier", "status", "direction", "event"]
    A = E.groupby(keys)["cls"].value_counts().unstack(fill_value=0).reindex(columns=classes, fill_value=0)
    A["n"] = A.sum(axis=1)
    A["n_mapped"] = A["n"] - A["unmapped"]
    for k in ("same", "deeper", "shallower"):
        A[f"frac_{k}"] = (A[k] / A["n_mapped"]).round(4)
    A = A.reset_index()
    A.to_csv(os.path.join(OUT, "agreement_per_branch.tsv"), sep="\t", index=False)
    # pooled: stems (exact; primary+secondary tiers) vs controls (exact) vs other exact internal branches
    rows = []
    for (direc, evt), sub in A.groupby(["direction", "event"]):
        grp = {"stems_exact": sub[(sub["role"] == "stem") & (sub["status"] == "exact")],
               "stems_approx": sub[(sub["role"] == "stem") & (sub["status"] == "approx")],
               "controls_exact": sub[(sub["role"] == "control") & (sub["status"] == "exact")],
               "other_internal_exact": sub[(sub["role"] == "other")]}
        for gname, s in grp.items():
            tot = s[classes].sum()
            nm = tot.sum() - tot["unmapped"]
            per = s["frac_same"].dropna()
            rows.append(dict(direction=direc, event=evt, group=gname, n_branches=len(s), n_events=int(tot.sum()),
                             n_mapped=int(nm), **{k: int(tot[k]) for k in classes},
                             frac_same=round(tot["same"] / nm, 4) if nm else np.nan,
                             frac_deeper=round(tot["deeper"] / nm, 4) if nm else np.nan,
                             frac_shallower=round(tot["shallower"] / nm, 4) if nm else np.nan,
                             median_branch_frac_same=round(per.median(), 4) if len(per) else np.nan))
        st, ct = grp["stems_exact"][classes].sum(), grp["controls_exact"][classes].sum()
        if st.sum() and ct.sum():
            tab = [[st["same"], st[["deeper", "shallower", "other"]].sum()],
                   [ct["same"], ct[["deeper", "shallower", "other"]].sum()]]
            orr, p = fisher_exact(tab)
            rows.append(dict(direction=direc, event=evt, group="stems_vs_controls_fisher", odds_ratio_same=orr,
                             fisher_p=p))
    # depth/size-adjusted baseline: binomial GLM of per-branch 'same' fraction on log(n_desc) and depth,
    # fitted on non-stem exact internal branches (controls + other); stems scored as observed/expected.
    import statsmodels.api as sm
    cov = pd.DataFrame(counts).query("method == 'ancestra'").drop_duplicates(["role", "branch_id"]).set_index(
        ["role", "branch_id"])[["depth", "n_desc"]]
    adj = []
    for (direc, evt), sub in A.groupby(["direction", "event"]):
        sub = sub.join(cov, on=["role", "branch_id"])
        sub = sub[sub["n_mapped"] > 0]
        base = sub[(sub["role"] != "stem") & (sub["status"] == "exact")]
        X = lambda d: sm.add_constant(np.column_stack([np.log(d["n_desc"]), d["depth"]]), has_constant="add")
        y = np.column_stack([base["same"], base["n_mapped"] - base["same"]])
        fit = sm.GLM(y, X(base), family=sm.families.Binomial()).fit(scale="X2")
        st = sub[sub["role"] == "stem"].copy()
        st["exp_frac_same"] = fit.predict(X(st))
        res_b = base["frac_same"] - fit.predict(X(base))
        for _, r in st.iterrows():
            rr = r["frac_same"] - r["exp_frac_same"]
            adj.append(dict(direction=direc, event=evt, branch_id=r["branch_id"], status=r["status"],
                            n_mapped=int(r["n_mapped"]), frac_same=r["frac_same"],
                            expected_frac_same=round(r["exp_frac_same"], 4),
                            obs_over_exp=round(r["frac_same"] / r["exp_frac_same"], 3),
                            percentile_among_baseline_residuals=round(float((res_b < rr).mean()), 3)))
        ex = st[st["status"] == "exact"]
        if len(ex):
            rows.append(dict(direction=direc, event=evt, group="stems_exact_depth_adjusted",
                             n_branches=len(ex), frac_same=round(ex["same"].sum() / ex["n_mapped"].sum(), 4),
                             expected_frac_same=round((ex["exp_frac_same"] * ex["n_mapped"]).sum() / ex["n_mapped"].sum(), 4),
                             n_baseline_branches=len(base)))
    pd.DataFrame(adj).to_csv(os.path.join(OUT, "agreement_stems_depth_adjusted.tsv"), sep="\t", index=False)
    P = pd.DataFrame(rows)
    P.to_csv(os.path.join(OUT, "agreement_pooled.tsv"), sep="\t", index=False)
    log("written counts_per_branch.tsv, agreement_per_branch.tsv, agreement_pooled.tsv, agreement_events.tsv.gz")
    print(P.to_string())


if __name__ == "__main__":
    main()
