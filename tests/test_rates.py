"""Tests for mcell.rates (H9 relative-rate statistic and RER) and mcell.phylo.

Simulation: an ultrametric random tree (equal root-to-tip time), per-gene
branch lengths = gene rate × time × lognormal noise, random missing taxa
(pruned with mcell.phylo.prune_named, so merged edges are exercised). No
sequence simulation (branch lengths are simulated directly).
"""
import math

import ete3
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mcell import rates
from mcell.phylo import prune_named, transfer_branch_lengths, read_tree


# ----------------------------------------------------------------------------- helpers

def ultrametric_tree(n_tips, seed):
    """Random rooted binary ultrametric tree (height 1) with named nodes n0.., tips t0..."""
    rng = np.random.default_rng(seed)
    root = ete3.TreeNode(name="n0")
    root.add_feature("h", 1.0)
    leaves = [root]
    ni = 1
    while len(leaves) < n_tips:
        # split the leaf with the largest height (keeps splits roughly time-ordered)
        leaves.sort(key=lambda x: -x.h)
        v = leaves.pop(0)
        for _ in range(2):
            c = v.add_child(name=f"n{ni}")
            ni += 1
            c.add_feature("h", v.h * rng.uniform(0.35, 0.85))
            leaves.append(c)
    k = 0
    for lf in root.get_leaves():
        lf.h = 0.0
        lf.name = f"t{k}"
        k += 1
    for n in root.traverse():
        if n.up is not None:
            n.dist = n.up.h - n.h
    return root


def simulate_genes(tree, n_genes, seed, planted_nodes=(), planted_genes=(), factor=2.0,
                   sd_noise=0.3, p_missing=0.1, keep_always=()):
    """{gene: GeneTree}; edges into planted_nodes (and their descendants) × factor in planted_genes."""
    rng = np.random.default_rng(seed)
    fg = rates.foreground_nodes(tree, planted_nodes) if planted_nodes else set()
    tips = tree.get_leaf_names()
    out = {}
    for g in range(n_genes):
        name = f"g{g:04d}"
        gr = rng.lognormal(0, 0.4)
        t = tree.copy()
        for n in t.traverse():
            if n.up is not None:
                mult = factor if (name in planted_genes and n.name in fg) else 1.0
                n.dist = n.dist * gr * mult * rng.lognormal(-sd_noise ** 2 / 2, sd_noise)
        keep = [x for x in tips if x in keep_always or rng.random() > p_missing]
        out[name] = rates.gene_tree_from_ete(prune_named(t, keep), name)
    return out


def pick_clade(tree, lo, hi, exclude=()):
    """A non-root internal node with lo..hi tips whose parent is not the root."""
    for n in tree.traverse("levelorder"):
        if n.up is not None and n.up.up is not None and not n.is_leaf() and lo <= len(n) <= hi \
                and n.name not in exclude:
            return n
    raise RuntimeError("no clade")


# ----------------------------------------------------------------------------- hand-computable L̄

HAND = "(((A:1,B:3)M:2,(C:1,(D:2,E:4)X:1)S:0.5)P:1,O:5)R;"


def test_equal_split_hand():
    t = read_tree(HAND)
    gt = rates.gene_tree_from_ete(t, "hand")
    # L̄(M) = (1+3)/2 = 2 ; L̄(S) = ½(1 + (1 + (2+4)/2)) = 2.5
    assert rates.equal_split_mean_path(gt, "M") == pytest.approx(2.0)
    assert rates.equal_split_mean_path(gt, "S") == pytest.approx(2.5)
    idx = rates._guide_index(t)
    d = rates.h9_gene_stat(gt, "P", "M", ["S"], idx)
    assert d["included"]
    assert d["L_M"] == pytest.approx(4.0)      # 2 + 2
    assert d["L_S"] == pytest.approx(3.0)      # 0.5 + 2.5 (unweighted tip mean would be 3.5)
    assert d["r"] == pytest.approx(math.log(4 / 3))


def test_inclusion_rules_hand():
    guide = read_tree(HAND)
    idx = rates._guide_index(guide)

    def stat(keep):
        gt = rates.gene_tree_from_ete(prune_named(guide, keep), "g")
        return rates.h9_gene_stat(gt, "P", "M", ["S"], idx)

    assert stat(["A", "B", "C", "O"])["included"]
    assert stat(["A", "B", "C", "D", "E"])["reason"] == "no_outgroup"
    assert stat(["A", "C", "D", "O"])["reason"] == "lt2_M_tips"
    assert stat(["A", "B", "O"])["reason"] == "no_S_tip"
    # pruning D,E,C... S represented by X only: S path merges S->X edge (0.5+1) then X's tips
    d = stat(["A", "B", "D", "E", "O"])
    assert d["included"] and d["L_S"] == pytest.approx(1.5 + 3.0)


def test_inclusion_mrca_not_M():
    guide = read_tree("((((A:1,B:1)Y:1,C:1)M:1,D:1)P:1,O:1)R;")
    idx = rates._guide_index(guide)
    gt = rates.gene_tree_from_ete(prune_named(guide, ["A", "B", "D", "O"]), "g")
    d = rates.h9_gene_stat(gt, "P", "M", ["D"], idx)
    assert not d["included"] and d["reason"] == "MRCA_of_M_tips_not_M"


def test_path_stem_hand():
    # D7 path stem P > X > M; Z (under X, outside M) belongs to neither side
    guide = read_tree("((((A:1,B:3)M:2,Z:7)X:0.5,(C:1,D:1)S:1)P:1,O:5)R;")
    idx = rates._guide_index(guide)
    gt = rates.gene_tree_from_ete(guide, "g")
    d = rates.h9_gene_stat(gt, "P", "M", ["S"], idx)
    assert d["L_M"] == pytest.approx(0.5 + 2 + 2) and d["L_S"] == pytest.approx(2.0)
    gt2 = rates.gene_tree_from_ete(prune_named(guide, ["A", "B", "C", "O"]), "g2")   # Z absent: X pruned
    assert rates.h9_gene_stat(gt2, "P", "M", ["S"], idx)["L_M"] == pytest.approx(4.5)
    tr = pd.DataFrame([dict(transition_id="T1", analysis_set="primary", M="M", P="P", sister_clade_nodes="S",
                            stem_path="P>X>M")])
    assert rates.stem_triples(tr, guide) == [("T1", "P", "M", ["S"])]
    with pytest.raises(ValueError):
        rates.stem_triples(tr.assign(stem_path=""), guide)


def test_drop_tips_matches_prune():
    guide = read_tree("((((A:1,B:3)M:2,Z:7)X:0.5,(C:1,D:1)S:1)P:1,O:5)R;")
    gt = rates.gene_tree_from_ete(guide, "g")
    a = rates.drop_tips(gt, {"Z", "D"})
    b = rates.gene_tree_from_ete(prune_named(guide, ["A", "B", "C", "O"]), "g")
    assert a.parent == b.parent and a.length == pytest.approx(b.length) and a.root == b.root


# ----------------------------------------------------------------------------- tree utilities

def test_transfer_and_prune():
    guide = read_tree("((a,b,c)A,(d,(e,f,g)C)B)R;")
    fitted = read_tree("(a:0.1,b:0.2,c:0.3,(d:0.4,(e:0.5,f:0.6,g:0.7)C:0.05)B:0.8)A;")
    L = transfer_branch_lengths(guide, fitted)
    assert L["A"] == pytest.approx(0.4) and L["B"] == pytest.approx(0.4)   # bifurcating root edge split
    assert L["C"] == pytest.approx(0.05) and L["g"] == pytest.approx(0.7)
    bad = read_tree("(a:0.1,b:0.2,c:0.3,(d:0.4,((e:0.5,f:0.6):0.1,g:0.7)C:0.05)B:0.8)A;")
    with pytest.raises(ValueError):
        transfer_branch_lengths(guide, bad)
    p = prune_named(read_tree("((a:1,b:1)A:1,(d:1,(e:1,f:1)C:2)B:1)R;"), ["a", "b", "e"])
    names = {n.name: n.dist for n in p.traverse()}
    assert "B" not in names and "C" not in names and names["e"] == pytest.approx(4.0)


# ----------------------------------------------------------------------------- H9 simulations

def _h9_setup(seed, planted):
    tree = ultrametric_tree(40, seed)
    M = pick_clade(tree, 6, 12)
    P = M.up
    S = [c.name for c in P.children if c is not M]
    genes = simulate_genes(tree, 150, seed + 1, planted_nodes=[M.name] if planted else (),
                           planted_genes={f"g{g:04d}" for g in range(150)}, factor=2.0)
    return tree, P.name, M.name, S, genes


def test_h9_planted_detected():
    tree, P, M, S, genes = _h9_setup(3, planted=True)
    res = rates.h9_unit(genes, tree, P, M, S)
    assert res["n_genes"] > 50
    assert res["ci_lo"] > 0.3                           # 2x shift: r ≈ log 2 = 0.69
    assert 0.45 < res["r"] < 0.95


def test_h9_null_near_zero():
    for seed in (5, 7):
        tree, P, M, S, genes = _h9_setup(seed, planted=False)
        res = rates.h9_unit(genes, tree, P, M, S)
        assert abs(res["r"]) < 0.12
        assert res["ci_lo"] < 0 < res["ci_hi"] or abs(res["r"]) < 0.05


def test_h9_controls_and_table():
    tree, P, M, S, genes = _h9_setup(9, planted=False)
    nodes = [n.name for n in tree.traverse() if not n.is_leaf() and n.up is not None]
    trip = rates.control_triples(tree, nodes)
    assert trip and all(len(t[3]) == 1 for t in trip)
    tab = rates.h9_gene_table(genes, tree, trip)
    summ = rates.h9_summary(tab, n_boot=200)
    ok = summ[summ["n_genes"] >= 30]
    assert len(ok) >= 3
    assert ok["r"].abs().median() < 0.15            # null controls centred on 0


# ----------------------------------------------------------------------------- RER simulations

def _rer_setup(seed, n_planted, factor=3.0, n_genes=200):
    tree = ultrametric_tree(48, seed)
    c1 = pick_clade(tree, 4, 8)
    c2 = pick_clade(tree, 4, 8, exclude={n.name for n in c1.traverse()} | {a.name for a in c1.get_ancestors()})
    planted = {f"g{g:04d}" for g in range(n_planted)}
    genes = simulate_genes(tree, n_genes, seed + 1, planted_nodes=[c1.name, c2.name], planted_genes=planted,
                           factor=factor, sd_noise=0.35, p_missing=0.08)
    return tree, [c1.name, c2.name], genes, planted


def test_rer_detects_planted_genes():
    tree, Ms, genes, planted = _rer_setup(21, n_planted=30)
    rer = rates.relative_rates(genes, tree, min_sp=10)
    res = rates.rer_permulation_test(rer, tree, Ms, n_perm=400, seed=11)
    pl = res.loc[sorted(planted)]
    nul = res.drop(index=sorted(planted))
    assert (pl["tau"] > 0).mean() > 0.95
    assert (pl["p_perm"] < 0.05).mean() >= 0.8
    # 15 % of genes share one 3x shift: the genome-wide master is pulled towards them and the
    # other genes look slightly slower in the foreground (inherent to RER normalisation; the
    # robust median master keeps this small)
    assert (nul["p_perm"] < 0.05).mean() < 0.15
    # top hits are the planted genes
    top = res.sort_values("p_perm").index[:20]
    assert len(set(top) & planted) >= 17


def test_rer_null_uniform():
    tree, Ms, genes, _ = _rer_setup(33, n_planted=0)
    rer = rates.relative_rates(genes, tree, min_sp=10)
    # RER is centred within genes (no residual dependence on the master length)
    assert abs(np.nanmean(rer.to_numpy())) < 0.05
    res = rates.rer_permulation_test(rer, tree, Ms, n_perm=400, seed=11)
    p = res["p_perm"].dropna().to_numpy()
    assert len(p) > 150
    assert stats.kstest(p, "uniform").pvalue > 0.01
    assert 0.0 < np.mean(p < 0.05) < 0.11


def test_tau_closed_form_matches_scipy():
    rng = np.random.default_rng(1)
    X = pd.DataFrame(rng.normal(size=(5, 30)))
    X.iloc[0, :3] = np.nan
    X.iloc[1, 4] = X.iloc[1, 5]                   # a tie
    fg = {0, 3, 7, 8, 20, 21}
    got = rates.rer_correlation(X, fg)["tau"].to_numpy()
    for i in range(5):
        v = X.iloc[i].to_numpy()
        ok = ~np.isnan(v)
        y = np.array([1 if c in fg else 0 for c in X.columns])[ok]
        assert got[i] == pytest.approx(stats.kendalltau(v[ok], y).statistic, abs=1e-12)


def test_master_lengths_recovers_times():
    tree = ultrametric_tree(30, 4)
    genes = simulate_genes(tree, 300, 5, sd_noise=0.2, p_missing=0.1)
    m, s = rates.master_lengths(genes, tree)
    true = pd.Series({n.name: n.dist for n in tree.traverse() if n.up is not None and n.up.up is not None})
    est = m.reindex(true.index)
    assert stats.spearmanr(est, true).statistic > 0.95
