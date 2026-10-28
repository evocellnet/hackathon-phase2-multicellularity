"""Tests for mcell.krkc: toy alignments with hand-countable radical/conservative changes, and a
planted simulation (more radical replacements on the M side) vs a null."""
import numpy as np
import pytest

from mcell import krkc, rates
from mcell.phylo import read_tree

TREE = "(((A:1,B:1)M:1,(C:1,D:1)S:1)P:1,O:1)R;"


def gt_of(nwk):
    return rates.gene_tree_from_ete(read_tree(nwk), "g")


def test_classification():
    R = krkc.radical_matrix("charge")
    i = krkc.AA_IDX
    assert R[i["K"], i["D"]] and not R[i["K"], i["R"]] and not R[i["A"], i["L"]] and R[i["A"], i["E"]]
    Rp = krkc.radical_matrix("polarity")
    assert Rp[i["A"], i["S"]] and not Rp[i["K"], i["E"]] and not Rp[i["L"], i["I"]]


def test_toy_counts():
    # site 1: M side changes K->D (radical in charge; both polar -> conservative in polarity)
    # site 2: S side changes K->R (conservative in both)
    # site 3: tip A changes L->S (charge conservative, polarity radical)
    # site 4: gap in B: never counted as a change
    seqs = {"A": "DKSL", "B": "DKL-", "C": "KRLL", "D": "KRLL", "O": "KKLL"}
    gt = gt_of(TREE)
    ec = krkc.edge_counts(gt, seqs)
    assert ec["n_changes"].sum() == 3
    assert ec.loc["M", "n_changes"] == 1 and ec.loc["S", "n_changes"] == 1 and ec.loc["A", "n_changes"] == 1
    assert ec.loc["M", "R_charge"] == 1 and ec.loc["M", "C_polarity"] == 1
    assert ec.loc["S", "C_charge"] == 1 and ec.loc["S", "C_polarity"] == 1
    assert ec.loc["A", "C_charge"] == 1 and ec.loc["A", "R_polarity"] == 1
    idx = rates._guide_index(read_tree(TREE))
    tab = krkc.unit_table({"g": ec}, {"g": gt}, idx, [("T", "P", "M", ["S"])])
    r = tab.iloc[0]
    assert (r.RM_charge, r.CM_charge, r.RS_charge, r.CS_charge) == (1, 1, 0, 1)
    st = krkc.unit_table({"g": ec}, {"g": gt}, idx, [("T", "P", "M", ["S"])], scope="stem").iloc[0]
    assert (st.RM_charge, st.CM_charge, st.RS_charge, st.CS_charge) == (1, 0, 0, 1)


def test_polytomy_fitch():
    # polytomy (A,B,C) with K,K,R and outgroup K: exactly one change, on C
    gt = gt_of("((A:1,B:1,C:1)X:1,O:1)R;")
    ec = krkc.edge_counts(gt, {"A": "K", "B": "K", "C": "R", "O": "K"})
    assert ec["n_changes"].sum() == 1 and ec.loc["C", "n_changes"] == 1


def _simulate(tree_nwk, n_genes, L, planted, seed, rad_boost=4.0, rate=0.15):
    """Sequences evolved down the tree; on M-side edges radical (charge) targets are up-weighted."""
    rng = np.random.default_rng(seed)
    t = read_tree(tree_nwk)
    R = krkc.radical_matrix("charge")
    m_side = {n.name for n in t.search_nodes(name="M")[0].traverse()}
    genes = {}
    for g in range(n_genes):
        seq = {t.name: rng.integers(0, 20, L)}
        for n in t.traverse("preorder"):
            if n.up is None:
                continue
            s = seq[n.up.name].copy()
            nchg = rng.random(L) < rate * n.dist
            for i in np.where(nchg)[0]:
                w = np.ones(20)
                w[s[i]] = 0
                if planted and n.name in m_side:
                    w[R[s[i]]] *= rad_boost
                s[i] = rng.choice(20, p=w / w.sum())
            seq[n.name] = s
        genes[f"g{g}"] = {lf: "".join(krkc.AA[x] for x in seq[lf]) for lf in t.get_leaf_names()}
    return t, genes


SIM_TREE = ("((((A1:1,A2:1)a:1,(A3:1,A4:1)b:1)M:1,((B1:1,B2:1)c:1,(B3:1,B4:1)d:1)S:1)P:1,"
            "(O1:1,O2:1)O:1)R;")


@pytest.mark.parametrize("planted", [True, False])
def test_simulated_lr(planted):
    t, genes = _simulate(SIM_TREE, 40, 150, planted, seed=3 if planted else 4)
    gt = rates.gene_tree_from_ete(t, "g")
    idx = rates._guide_index(t)
    counts = {g: krkc.edge_counts(gt, s) for g, s in genes.items()}
    gts = {g: rates.GeneTree(name=g, root=gt.root, parent=gt.parent, length=gt.length) for g in genes}
    tab = krkc.unit_table(counts, gts, idx, [("T", "P", "M", ["S"])])
    res = krkc.summarise(tab, n_boot=500).loc["T"]
    if planted:
        assert res.LR_charge > 0.3 and res.ci_lo_charge > 0
    else:
        assert abs(res.LR_charge) < 0.15 and res.ci_lo_charge < 0 < res.ci_hi_charge
