"""Evolutionary-rate statistics on fixed-topology gene trees (DESIGN.md §4.5, §4.6 H9).

Inputs are *gene trees*: the guide tree pruned to the taxa present in a marker
(:func:`mcell.phylo.prune_named`; every kept node carries the guide name of the
MRCA of its present tips) with branch lengths from IQ-TREE ``-te`` (LG+G4),
as written by ``scripts/phylo/gene_trees.py`` (``<marker>.bl.nwk``) or as the
long table ``branch_lengths.tsv`` (gene, node, parent, length). Loaders:
:func:`load_gene_trees`, :func:`gene_trees_from_table`.

(a) H9 relative-rate statistic (reviews/round2/statistics.md M6, binding)
-------------------------------------------------------------------------
For a unit (P, M, S) -- a transition stem (P = parent of M, S = the sister
nodes of M, from ``transitions.tsv``) or a control node (P = v, M/S = two child
clades with >= 2 tips) -- and a gene g::

    r_g = log( L̄(P→M tips) / L̄(P→S tips) )

with L̄ the *equal-split weighted mean path length*: weight 1 starts at the
first node below P on that side and is split equally among the children at
every node, so L̄(X) = (1/k) Σ_children c (ℓ_c + L̄(c)) with L̄(tip) = 0, and

    L̄(P→M tips) = ℓ(P⇝M) + L̄(M),     L̄(P→S tips) = mean_{s ∈ S present} (ℓ_s + L̄(s)),

where ℓ(P⇝M) is the gene-tree path length from P to M (one edge, or the
summed edges of a D7 path stem P > X > M).

Gene inclusion (M6): ≥ 1 outgroup tip (outside P, so P is identifiable in the
unrooted fit), ≥ 2 M tips whose MRCA is M (M is a kept internal node of the
pruned tree), ≥ 1 S tip, ≥ 50 % site coverage (applied per taxon upstream in
gene_trees.py: taxa covering < 50 % of the trimmed columns are not in the gene
tree). r_i = median_g r_g with a 2,000× gene-bootstrap percentile CI
(seed 11): :func:`h9_gene_table`, :func:`h9_summary`, :func:`h9_unit`.

Controls: :func:`control_triples` turns the control nodes of
``mcell.branchstats.h9_control_branches`` (nodes in U with two child clades of
≥ 2 tips) into (node, child1, child2) triples in a FIXED orientation (child1 =
alphabetically first of the two largest qualifying children); the random
orientation of §4.6 is applied once, by ``branchstats.h9_values(orient=True)``
(seed 11). Do not flip signs here as well.

(b) RERconverge-style relative evolutionary rates (exploratory, §4.5)
---------------------------------------------------------------------
Python re-implementation following Chikina et al. 2016 / Partha et al. 2019:

1. *Paths.* A gene-tree edge (parent p → child c) spans the guide edges from
   c up to p; it is stored in the column of its child node c (RERconverge's
   tree2Paths convention). The two root-child edges of every gene tree are
   dropped (their split is not identifiable from an unrooted fit).
2. *Genome-wide average branch-length vector* m (the "master tree"; robust version):
   iteratively, s_g = Σ_e ℓ_ge / Σ_e E_ge (gene scale) and
   m_j = mean_g(ℓ_gj / s_g) over genes in which guide edge j is observed as a
   single (non-merged) edge; E_ge = Σ_{j ∈ path(e)} m_j. Initial m = a
   reference tree (e.g. the concatenated species tree) if given, else the
   per-edge means of the raw lengths. Unobserved guide edges keep the
   (rescaled) reference length.
3. *Residuals.* y = √ℓ_ge, x = √E_ge; per gene a least-squares fit
   y = a_g + b_g x (projection on the expected-length vector), weighted by
   1/σ²(x), where σ²(x) = c·x^k is a power law fitted to the binned squared
   OLS residuals pooled over genes (the mean–variance correction of
   RERconverge ``weighted=TRUE``). RER_ge = (y − ŷ)/σ(x) (standardised).
4. *Association.* Kendall's τ_b between RER and the binary foreground
   indicator (foreground = stems + multicellular crowns: every edge whose child
   is M or a descendant of M), per gene, with ≥ ``min_sp`` tips, ≥ 2
   foreground and ≥ 2 background edges.
5. *Permulation null* (Saputra et al. 2021, binary-trait analogue): each
   permutation places the observed foreground clades at random non-overlapping
   positions of the guide tree with matching tip counts (within a factor
   ``size_tol``; the nearest available size otherwise), preserving the number
   and the phylogenetic clustering of foreground edges; τ_b is recomputed for
   every gene against the same permuted foregrounds. Two-sided per-gene
   p = (1 + #{|τ*| ≥ |τ|}) / (1 + n valid permutations); BH-FDR across genes.

τ_b with a binary variable has a closed form in the ranks, so all genes and
permutations are evaluated with two matrix products (:func:`_tau_b_binary`).
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

SEED = 11
N_BOOT = 2000

__all__ = [
    "GeneTree", "gene_tree_from_ete", "load_gene_trees", "gene_trees_from_table",
    "equal_split_mean_path", "h9_gene_stat", "h9_gene_table", "h9_summary", "h9_unit",
    "stem_triples", "control_triples", "drop_tips",
    "master_lengths", "relative_rates", "foreground_nodes", "rer_correlation",
    "permuted_foregrounds", "rer_permulation_test", "bh",
]


# =============================================================================
# gene trees
# =============================================================================

@dataclass
class GeneTree:
    """Rooted tree with guide node names: parent/length/children maps."""
    name: str
    root: str
    parent: dict
    length: dict
    children: dict = field(default_factory=dict)
    _tips: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        if not self.children:
            ch = {}
            for c, p in self.parent.items():
                ch.setdefault(p, []).append(c)
            self.children = ch
        self._tips = {}

    @property
    def nodes(self):
        return set(self.parent) | {self.root}

    def is_leaf(self, n):
        return n not in self.children

    def tips(self, n=None) -> frozenset:
        n = self.root if n is None else n
        if n not in self._tips:
            if self.is_leaf(n):
                self._tips[n] = frozenset([n])
            else:
                self._tips[n] = frozenset().union(*(self.tips(c) for c in self.children[n]))
        return self._tips[n]


def gene_tree_from_ete(t, name="") -> GeneTree:
    parent, length = {}, {}
    for n in t.traverse():
        if n.up is not None:
            parent[n.name] = n.up.name
            length[n.name] = float(n.dist)
    return GeneTree(name=name, root=t.get_tree_root().name, parent=parent, length=length)


def load_gene_trees(gdir: str, genes: Iterable[str] | None = None) -> dict:
    """{gene: GeneTree} from ``<gdir>/<gene>.bl.nwk``."""
    import ete3
    if genes is None:
        genes = sorted(f[:-7] for f in os.listdir(gdir) if f.endswith(".bl.nwk"))
    out = {}
    for g in genes:
        p = os.path.join(gdir, f"{g}.bl.nwk")
        if os.path.exists(p):
            out[g] = gene_tree_from_ete(ete3.Tree(open(p).read().strip(), format=1), g)
    return out


def gene_trees_from_table(df: pd.DataFrame) -> dict:
    """{gene: GeneTree} from a long table with columns gene, node, parent, length."""
    out = {}
    for g, d in df.groupby("gene", sort=True):
        parent = dict(zip(d["node"], d["parent"]))
        roots = set(parent.values()) - set(parent)
        if len(roots) != 1:
            raise ValueError(f"{g}: {len(roots)} roots")
        out[g] = GeneTree(name=g, root=roots.pop(), parent=parent,
                          length=dict(zip(d["node"], d["length"].astype(float))))
    return out


def drop_tips(gt: GeneTree, drop: Iterable[str]) -> GeneTree | None:
    """Gene tree without the tips in ``drop`` (unary nodes removed, lengths summed;
    the lower node keeps its name, as in mcell.phylo.prune_named). None if < 2 tips remain."""
    drop = set(drop) & gt.tips()
    if not drop:
        return gt
    parent, length = {}, {}

    def rec(n):
        if gt.is_leaf(n):
            return n if n not in drop else None
        kids = [k for k in (rec(c) for c in gt.children[n]) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            length[kids[0]] = length[kids[0]] + gt.length.get(n, 0.0)
            return kids[0]
        for k in kids:
            parent[k] = n
        length[n] = gt.length.get(n, 0.0)
        return n

    for leaf in gt.tips():
        length[leaf] = gt.length.get(leaf, 0.0)
    for n in gt.children:
        length.setdefault(n, gt.length.get(n, 0.0))
    root = rec(gt.root)
    if root is None or gt.is_leaf(root):
        return None
    return GeneTree(name=gt.name, root=root, parent=parent, length={k: length[k] for k in parent})


def _guide_index(guide):
    """name -> (parent name, frozenset tips, children names) for an ete3 guide tree."""
    idx = {}
    for n in guide.traverse("postorder"):
        tips = frozenset([n.name]) if n.is_leaf() else frozenset().union(*(idx[c.name][1] for c in n.children))
        idx[n.name] = (n.up.name if n.up is not None else None, tips, [c.name for c in n.children])
    return idx


# =============================================================================
# (a) H9
# =============================================================================

def equal_split_mean_path(gt: GeneTree, node: str) -> float:
    """L̄ from ``node`` down to its tips with equal-split weights (0 for a tip)."""
    if gt.is_leaf(node):
        return 0.0
    kids = gt.children[node]
    return sum(gt.length[c] + equal_split_mean_path(gt, c) for c in kids) / len(kids)


def h9_gene_stat(gt: GeneTree, P: str, M: str, S: Sequence[str], guide_idx: Mapping) -> dict:
    """r_g for one gene and one unit (P, M, S) with the M6 inclusion rules.

    ``guide_idx`` = :func:`_guide_index` of the full guide tree (to know which
    present tips belong to M, S and the outgroup).
    """
    tipsP, tipsM = guide_idx[P][1], guide_idx[M][1]
    tipsS = frozenset().union(*(guide_idx[s][1] for s in S))
    present = gt.tips()
    n_m, n_s = len(present & tipsM), len(present & tipsS)
    n_out = len(present - tipsP)
    base = dict(gene=gt.name, n_M=n_m, n_S=n_s, n_out=n_out, L_M=np.nan, L_S=np.nan, r=np.nan)
    if n_out < 1:
        return {**base, "included": False, "reason": "no_outgroup"}
    if n_m < 2:
        return {**base, "included": False, "reason": "lt2_M_tips"}
    if M not in gt.parent or gt.is_leaf(M):
        return {**base, "included": False, "reason": "MRCA_of_M_tips_not_M"}
    if n_s < 1:
        return {**base, "included": False, "reason": "no_S_tip"}
    stem, n = 0.0, M          # path stem P > ... > M (D7): sum the gene-tree edges from M up to P
    while n != P:
        if n not in gt.parent:
            raise AssertionError(f"{gt.name}: {P} is not an ancestor of {M} in the gene tree")
        stem += gt.length[n]
        n = gt.parent[n]
    s_kids = [c for c in gt.children[P] if gt.tips(c) <= tipsS]
    if not s_kids:
        return {**base, "included": False, "reason": "no_S_tip"}
    L_M = stem + equal_split_mean_path(gt, M)
    L_S = sum(gt.length[c] + equal_split_mean_path(gt, c) for c in s_kids) / len(s_kids)
    if not (L_M > 0 and L_S > 0):
        return {**base, "L_M": L_M, "L_S": L_S, "included": False, "reason": "zero_length"}
    return {**base, "L_M": L_M, "L_S": L_S, "r": math.log(L_M / L_S), "included": True, "reason": ""}


def stem_triples(transitions: pd.DataFrame, guide, sets: Sequence[str] | None = ("primary",),
                 kingdom: str | None = None) -> list:
    """[(transition_id, P, M, [S...])] for transitions whose M is an internal guide node.

    S = ``sister_clade_nodes`` (';'-separated), each of which must be a child of
    P; P must be the parent of M, or ``stem_path`` ("P>X>M", DESIGN v3.1 D7)
    must be a parent->child chain from P to M. For a path stem, L̄(P→M tips)
    includes every edge of the path, and tips below the intermediate nodes
    but outside M (e.g. Allomyces under Blastocladiomycota_Amastigomycota)
    belong to neither side. Rows with M a tip or a tip list are skipped
    (they cannot meet the ≥ 2 M tips rule).
    """
    idx = _guide_index(guide)
    out = []
    for _, r in transitions.iterrows():
        if kingdom and r.get("kingdom") != kingdom:
            continue
        if sets is not None and r.get("analysis_set") not in sets:
            continue
        M, P = str(r["M"]), str(r["P"])
        if M not in idx or P not in idx or not idx[M][2]:
            continue
        path = [x for x in str(r.get("stem_path", "") or "").split(">") if x and x != "nan"]
        if path:
            if path[0] != P or path[-1] != M or any(idx[b][0] != a for a, b in zip(path, path[1:])):
                raise ValueError(f"{r['transition_id']}: stem_path {path} is not a chain P={P} > ... > M={M}")
        elif idx[M][0] != P:
            raise ValueError(f"{r['transition_id']}: P={P} is not the parent of M={M} and no stem_path")
        S = [s for s in str(r["sister_clade_nodes"]).split(";") if s]
        bad = [s for s in S if idx.get(s, (None,))[0] != P]
        if bad:
            raise ValueError(f"{r['transition_id']}: sister nodes {bad} are not children of P={P}")
        out.append((r["transition_id"], P, M, S))
    return out


def control_triples(guide, nodes: Iterable[str], min_tips: int = 2) -> list:
    """[(node, node, child1, [child2])] for control nodes with ≥ 2 children of ≥ ``min_tips`` tips.

    Fixed orientation: the two largest qualifying children (ties by name), child1
    = alphabetically first. The unit id is the node name (= branch id in
    branchstats). Random orientation is applied later by
    ``branchstats.h9_values(orient=True)``.
    """
    idx = _guide_index(guide)
    out = []
    for v in nodes:
        big = [c for c in idx[v][2] if len(idx[c][1]) >= min_tips]
        if len(big) < 2:
            continue
        big = sorted(big, key=lambda c: (-len(idx[c][1]), c))[:2]
        c1, c2 = sorted(big)
        out.append((v, v, c1, [c2]))
    return out


def h9_gene_table(gene_trees: Mapping[str, GeneTree], guide, units: Sequence[tuple],
                  max_branch: float | None = None, unit_exclude: Mapping[str, Iterable[str]] | None = None
                  ) -> pd.DataFrame:
    """Long table unit × gene with L_M, L_S, r, included, reason.

    ``units``: (unit_id, P, M, [S...]) tuples from :func:`stem_triples` /
    :func:`control_triples`. ``max_branch``: optional exclusion of genes whose
    tree has any branch longer than this (methods review §6 suggests 2 subst./site;
    not part of the binding M6 rules, so off by default). ``unit_exclude``: optional
    unit -> tips removed from every gene tree before computing that unit (e.g. the
    multicellular tips inside the child clades of a control node; off by default,
    the binding M6 rule uses all tips).
    """
    idx = _guide_index(guide)
    long_genes = set()
    if max_branch is not None:
        long_genes = {g for g, gt in gene_trees.items() if max(gt.length.values()) > max_branch}
    rows = []
    for uid, P, M, S in units:
        excl = set((unit_exclude or {}).get(uid, ()))
        for g in sorted(gene_trees):
            gt = drop_tips(gene_trees[g], excl) if excl else gene_trees[g]
            if gt is None:
                rows.append({"unit": uid, "P": P, "M": M, "S": ";".join(S), "gene": g, "included": False,
                             "reason": "lt2_tips_after_exclusion"})
                continue
            d = h9_gene_stat(gt, P, M, S, idx)
            if d["included"] and g in long_genes:
                d.update(included=False, reason="branch_gt_max")
            rows.append({"unit": uid, "P": P, "M": M, "S": ";".join(S), **d})
    return pd.DataFrame(rows)


def _boot_median_ci(r: np.ndarray, n_boot: int, seed: int, level: float = 0.95):
    rng = np.random.default_rng(seed)
    meds = np.median(r[rng.integers(0, len(r), size=(n_boot, len(r)))], axis=1)
    a = (1 - level) / 2
    return float(np.quantile(meds, a)), float(np.quantile(meds, 1 - a))


def h9_summary(table: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED, min_genes: int = 1) -> pd.DataFrame:
    """Per unit: n_genes, r = median r_g, 95 % gene-bootstrap CI, exp(r), exclusion counts."""
    rows = []
    for uid, d in table.groupby("unit", sort=False):
        inc = d.loc[d["included"], "r"].to_numpy(float)
        row = dict(unit=uid, P=d["P"].iat[0], M=d["M"].iat[0], S=d["S"].iat[0], n_genes=len(inc),
                   r=np.nan, ci_lo=np.nan, ci_hi=np.nan, rate_ratio=np.nan, frac_pos=np.nan,
                   n_boot=n_boot, seed=seed)
        row.update({f"excl_{k}": int(v) for k, v in d.loc[~d["included"], "reason"].value_counts().items()})
        if len(inc) >= min_genes and len(inc) > 0:
            lo, hi = _boot_median_ci(inc, n_boot, seed)
            row.update(r=float(np.median(inc)), ci_lo=lo, ci_hi=hi, rate_ratio=float(np.exp(np.median(inc))),
                       frac_pos=float(np.mean(inc > 0)))
        rows.append(row)
    return pd.DataFrame(rows).set_index("unit")


def h9_unit(gene_trees, guide, P, M, S, n_boot=N_BOOT, seed=SEED, unit="unit") -> dict:
    """Convenience: one unit -> summary dict."""
    tab = h9_gene_table(gene_trees, guide, [(unit, P, M, list(S))])
    return h9_summary(tab, n_boot, seed).loc[unit].to_dict()


# =============================================================================
# (b) RERconverge-style relative rates
# =============================================================================

def _edge_table(gene_trees: Mapping[str, GeneTree], guide_idx) -> pd.DataFrame:
    rows = []
    for g, gt in gene_trees.items():
        for c, p in gt.parent.items():
            path = []
            n = c
            while n != p:
                path.append(n)
                n = guide_idx[n][0]
                if n is None:
                    raise ValueError(f"{g}: {p} is not a guide ancestor of {c}")
            rows.append((g, c, p, gt.length[c], tuple(path), p == gt.root, len(gt.tips())))
    return pd.DataFrame(rows, columns=["gene", "node", "parent", "length", "path", "root_child", "n_tips"])


def master_lengths(gene_trees: Mapping[str, GeneTree], guide, reference: Mapping[str, float] | None = None,
                   n_iter: int = 5, stat: str = "median") -> tuple[pd.Series, pd.Series]:
    """Genome-wide average branch-length vector m (guide node -> length) and gene scales s_g.

    ``stat``: 'median' (default; robust to a minority of genes with lineage-specific
    shifts, which would otherwise bias every other gene's RER in the opposite
    direction) or 'mean' (RERconverge's averaging)."""
    idx = _guide_index(guide)
    et = _edge_table(gene_trees, idx)
    et = et[~et["root_child"]]
    nodes = [n for n in idx if idx[n][0] is not None]
    single = et[et["path"].map(len) == 1]
    if reference is not None:
        m = pd.Series({n: float(reference[n]) for n in nodes})
    else:
        m = single.groupby("node")["length"].agg(stat).reindex(nodes)
        m = m.fillna(m.median())
    m = m.clip(lower=1e-6)
    ref0 = m.copy()
    s = pd.Series(1.0, index=sorted(gene_trees))
    for _ in range(n_iter):
        E = et["path"].map(lambda pth: sum(m[j] for j in pth))
        s = (et["length"] / E).groupby(et["gene"]).median() if stat == "median" else \
            et.groupby("gene")["length"].sum() / E.groupby(et["gene"]).sum()
        sg = single.assign(rel=single["length"] / single["gene"].map(s))
        new = sg.groupby("node")["rel"].agg(stat)
        obs = new.index
        k = float(np.median(new / m.reindex(obs))) if len(obs) else 1.0
        m = (ref0 * k).where(~m.index.isin(obs), new.reindex(m.index)).clip(lower=1e-6)
        ref0 = ref0 * k
    return m, s


def relative_rates(gene_trees: Mapping[str, GeneTree], guide, master: pd.Series | None = None,
                   reference: Mapping[str, float] | None = None, weighted: bool = True,
                   min_sp: int = 10, n_bins: int = 20) -> pd.DataFrame:
    """Genes × guide-node matrix of relative evolutionary rates (NaN = absent / root-child / too few taxa)."""
    idx = _guide_index(guide)
    if master is None:
        master, _ = master_lengths(gene_trees, guide, reference)
    et = _edge_table(gene_trees, idx)
    et = et[~et["root_child"] & (et["n_tips"] >= min_sp)].copy()
    et["x"] = np.sqrt(et["path"].map(lambda pth: sum(master[j] for j in pth)))
    et["y"] = np.sqrt(et["length"].clip(lower=0))

    def fit(d, w=None):
        X = np.column_stack([np.ones(len(d)), d["x"].to_numpy()])
        y = d["y"].to_numpy()
        ww = np.ones(len(d)) if w is None else w
        beta, *_ = np.linalg.lstsq(X * np.sqrt(ww)[:, None], y * np.sqrt(ww), rcond=None)
        return y - X @ beta

    et["res"] = (pd.concat([pd.Series(fit(d), index=d.index) for _, d in et.groupby("gene", sort=False)])
                 if len(et) else pd.Series(dtype=float))
    if weighted and len(et) > n_bins:
        # power-law mean-variance relationship of squared residuals vs x (pooled over genes)
        q = pd.qcut(et["x"], n_bins, duplicates="drop")
        b = et.groupby(q, observed=True).agg(x=("x", "mean"), v=("res", lambda r: np.mean(r ** 2)))
        b = b[(b["x"] > 0) & (b["v"] > 0)]
        k, logc = np.polyfit(np.log(b["x"]), np.log(b["v"]), 1)
        var = np.exp(logc) * np.clip(et["x"].to_numpy(), b["x"].min(), b["x"].max()) ** k
    else:
        var = np.full(len(et), float(np.mean(et["res"] ** 2)) if len(et) else 1.0)
    et["var"] = var
    out = []
    for g, d in et.groupby("gene", sort=False):
        r = fit(d, 1.0 / d["var"].to_numpy())
        out.append(pd.Series(r / np.sqrt(d["var"].to_numpy()), index=d["node"].to_numpy(), name=g))
    nodes = [n for n in idx if idx[n][0] is not None]
    return pd.DataFrame(out).reindex(columns=nodes)


def foreground_nodes(guide, M_nodes: Iterable[str]) -> set:
    """Stems + crowns: every guide node equal to or below one of ``M_nodes``."""
    idx = _guide_index(guide)
    fg = set()
    for M in M_nodes:
        fg.update(n for n in idx if _is_desc(idx, n, M))
    return fg


def _is_desc(idx, n, anc):
    while n is not None:
        if n == anc:
            return True
        n = idx[n][0]
    return False


def _tau_b_binary(R: np.ndarray, mask: np.ndarray, tx: np.ndarray, F: np.ndarray):
    """Kendall τ_b between per-gene values (ranks R, NaN-free, 0 where absent) and binary rows of F.

    R, mask: genes × edges; tx: per-gene tie term Σ t(t−1)/2 of the values;
    F: perms × edges (0/1). Returns (tau, n1, n0) of shape genes × perms.
    """
    sumR = R @ F.T
    n1 = mask @ F.T
    n = mask.sum(1, keepdims=True)
    n0 = n - n1
    S = 2.0 * (sumR - n1 * (n1 + 1) / 2.0) - n1 * n0
    N0 = n * (n - 1) / 2.0
    ty = n1 * (n1 - 1) / 2.0 + n0 * (n0 - 1) / 2.0
    with np.errstate(invalid="ignore", divide="ignore"):
        tau = S / np.sqrt((N0 - tx[:, None]) * (N0 - ty))
    return tau, n1, n0


def _rank_matrix(rer: pd.DataFrame):
    from scipy.stats import rankdata
    X = rer.to_numpy(float)
    mask = ~np.isnan(X)
    R = np.zeros_like(X)
    tx = np.zeros(len(X))
    for i in range(len(X)):
        v = X[i, mask[i]]
        R[i, mask[i]] = rankdata(v)
        _, cnt = np.unique(v, return_counts=True)
        tx[i] = float(np.sum(cnt * (cnt - 1) / 2.0))
    return R, mask.astype(float), tx


def rer_correlation(rer: pd.DataFrame, fg: set, min_fg: int = 2) -> pd.DataFrame:
    """Per gene: Kendall τ_b(RER, foreground), n edges, n foreground edges."""
    F = np.array([[1.0 if c in fg else 0.0 for c in rer.columns]])
    R, mask, tx = _rank_matrix(rer)
    tau, n1, n0 = _tau_b_binary(R, mask, tx, F)
    out = pd.DataFrame({"tau": tau[:, 0], "n_edges": mask.sum(1), "n_fg": n1[:, 0]}, index=rer.index)
    out.loc[(out["n_fg"] < min_fg) | (out["n_edges"] - out["n_fg"] < 2), "tau"] = np.nan
    return out


def _clade_roots(idx, fg_M: Iterable[str]) -> list:
    roots = sorted(set(fg_M))
    return [r for r in roots if not any(o != r and _is_desc(idx, r, o) for o in roots)]


def permuted_foregrounds(guide, M_nodes: Iterable[str], n_perm: int = 1000, seed: int = SEED,
                         size_tol: float = 1.5, max_tries: int = 200) -> list:
    """Random clade placements matching the observed foreground clades (see module doc).

    Returns a list of foreground node sets (each = union of chosen clades).
    """
    idx = _guide_index(guide)
    roots = _clade_roots(idx, M_nodes)
    sizes = sorted((len(idx[r][1]) for r in roots), reverse=True)
    cand = [n for n in idx if idx[n][0] is not None]
    csize = np.array([len(idx[n][1]) for n in cand])
    rng = np.random.default_rng(seed)
    pools = []
    for k in sizes:
        ratio = np.abs(np.log(csize / k))
        ok = np.where(ratio <= math.log(size_tol))[0]
        if len(ok) == 0:
            ok = np.where(ratio == ratio.min())[0]
        pools.append(ok)
    out = []
    for _ in range(n_perm):
        for _try in range(max_tries):
            chosen = []
            fail = False
            for pool in pools:
                opts = [i for i in pool if all(not _is_desc(idx, cand[i], c) and not _is_desc(idx, c, cand[i])
                                               for c in chosen)]
                if not opts:
                    fail = True
                    break
                chosen.append(cand[opts[rng.integers(len(opts))]])
            if not fail:
                break
        if fail:
            continue
        out.append(foreground_nodes(guide, chosen))
    return out


def bh(p) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    adj = np.minimum.accumulate((pv[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    res = np.empty(n)
    res[o] = np.minimum(adj, 1.0)
    q[ok] = res
    return q


def rer_permulation_test(rer: pd.DataFrame, guide, M_nodes: Iterable[str], n_perm: int = 1000,
                         seed: int = SEED, min_fg: int = 2, alternative: str = "two-sided",
                         size_tol: float = 1.5) -> pd.DataFrame:
    """Per gene τ_b, permulation p and BH q (foreground = stems + crowns of ``M_nodes``)."""
    M_nodes = list(M_nodes)
    fg = foreground_nodes(guide, M_nodes)
    obs = rer_correlation(rer, fg, min_fg)
    perms = permuted_foregrounds(guide, M_nodes, n_perm, seed, size_tol)
    F = np.array([[1.0 if c in s else 0.0 for c in rer.columns] for s in perms])
    R, mask, tx = _rank_matrix(rer)
    tau_p, n1, n0 = _tau_b_binary(R, mask, tx, F)
    valid = (n1 >= min_fg) & (n0 >= 2) & ~np.isnan(tau_p)
    t = obs["tau"].to_numpy()[:, None]
    if alternative == "two-sided":
        hit = np.abs(tau_p) >= np.abs(t) - 1e-12
    elif alternative == "greater":
        hit = tau_p >= t - 1e-12
    else:
        hit = tau_p <= t + 1e-12
    nv = valid.sum(1)
    p = (1 + (hit & valid).sum(1)) / (1 + nv)
    p = np.where(np.isnan(obs["tau"].to_numpy()) | (nv == 0), np.nan, p)
    out = obs.copy()
    out["p_perm"] = p
    out["n_perm_valid"] = nv
    out["q_bh"] = bh(p)
    out.attrs.update(n_perm=len(perms), seed=seed, M_nodes=M_nodes, alternative=alternative)
    return out
