"""Planted-effect tests for the analysis layer (mcell.analysis, mcell.controls, scripts/analysis/*).

A synthetic balanced species tree (32 tips, named internal nodes) with one primary transition
(P = n01, M = n010) and a synthetic OrthoXML whose families carry planted events:

* every internal non-root branch: 12 single-duplication families, ~10 % of them TF-DBD (Homeodomain PF00046);
  the stem n01 -> n010 instead gets 12 duplications of which 10 are TF-DBD   -> H4 u(stem) at its maximum;
* 20 flagellar-module families (one per h12_sets.tsv member Pfam of the 4 modules) present everywhere and lost on the stem
                                                                            -> H12 f(stem) = 1, sign by direction;
* a Wnt (PF00110) family gained on the stem (present only in M)              -> control rule 'gain' recovered;
  a TGF-beta (PF00019) family gained on the stem but with hits in 3 outside species -> leak, not recovered;
* one control branch with a 4-fold duplication excess                       -> WGD/duplication-excess flag.
Expected values are hand-derivable (see each test).
"""
from __future__ import annotations

import gzip
import io
import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts", "analysis"))

from mcell import analysis as an  # noqa: E402
from mcell import branchstats as bs  # noqa: E402
from mcell import controls as pc  # noqa: E402

DEPTH = 5
STEM_P, STEM_M = "n01", "n010"
TF = "PF00046"


# -----------------------------------------------------------------------------
# synthetic data
# -----------------------------------------------------------------------------

def make_tree(depth=DEPTH):
    """Balanced binary tree; internal nodes 'n' + path bits, leaves 's' + path bits."""
    def rec(path):
        if len(path) == depth:
            return f"s{path}"
        return f"({rec(path + '0')},{rec(path + '1')})n{path}"
    return rec("") + ";"


class Tree:
    def __init__(self, depth=DEPTH):
        self.depth = depth
        self.nodes = []
        self.children = {}

        def rec(path):
            nm = f"s{path}" if len(path) == depth else f"n{path}"
            self.nodes.append(nm)
            if len(path) < depth:
                self.children[nm] = [rec(path + "0"), rec(path + "1")]
            else:
                self.children[nm] = []
            return nm
        self.root = rec("")

    def leaves(self, nm):
        if not self.children[nm]:
            return [nm]
        return [x for c in self.children[nm] for x in self.leaves(c)]


class Family:
    """One gene family: origin node, duplication nodes (dup on the branch into the node), loss nodes."""

    def __init__(self, fid, origin, dups=(), losses=(), pfam=None):
        self.fid, self.origin, self.dups, self.losses, self.pfam = fid, origin, set(dups), set(losses), pfam


def emit(T: Tree, fam: Family, genes: list):
    gid = [0]

    def rec(nm, suf, allow_dup=True):
        if nm in fam.losses:
            return None
        if allow_dup and nm in fam.dups:
            a, b = rec(nm, suf + "a", False), rec(nm, suf + "b", False)
            kids = [x for x in (a, b) if x]
            if len(kids) == 2:
                return f"<paralogGroup>{kids[0]}{kids[1]}</paralogGroup>"
            return kids[0] if kids else None
        if not T.children[nm]:
            g = f"{nm}_{fam.fid}{suf}"
            genes.append((nm, g, fam.pfam))
            return f'<geneRef id="{g}"/>'
        kids = [x for x in (rec(c, suf) for c in T.children[nm]) if x]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        gid[0] += 1
        return f'<orthologGroup id="{fam.fid}_{gid[0]}">{"".join(kids)}</orthologGroup>'

    x = rec(fam.origin, "")
    if x is None or x.startswith("<geneRef"):
        return None                      # a single gene stays ungrouped (species-specific singleton)
    return x


def build_synthetic(tmp, seed=11):
    rng = np.random.default_rng(seed)
    T = Tree()
    internal = [n for n in T.nodes if T.children[n] and n != T.root]
    fams = []
    k = 0
    # duplication families: 12 per internal non-root branch (stem: 10 TF), one branch with 48 (4x excess)
    for b in internal:
        nd = 48 if b == "n100" else 12
        for i in range(nd):
            if b == STEM_M:
                tf = i < 10
            else:
                tf = rng.random() < 0.1
            k += 1
            fams.append(Family(f"D{k}", T.root, dups=[b], pfam=TF if tf else None))
    # background core families (TF share ~10 %)
    for i in range(150):
        k += 1
        fams.append(Family(f"B{k}", T.root, pfam=TF if rng.random() < 0.1 else None))
    # flagellar modules: one family per component, lost on the stem
    h12 = pd.read_csv(os.path.join(ROOT, "data", "h12_sets.tsv"), sep="\t")
    h12 = h12[(h12["role"] == "member") & h12["module"].isin(an.H12_MODULES)]
    for _, r in h12.iterrows():
        k += 1
        fams.append(Family(f"F{k}", T.root, losses=[STEM_M], pfam=r["pfam_acc"]))
    # gains on the stem: Wnt (clean) and TGF-beta (leaky)
    fams.append(Family("WNT", STEM_M, pfam="PF00110"))
    fams.append(Family("TGF", STEM_M, pfam="PF00019"))
    genes, groups = [], []
    for f in fams:
        x = emit(T, f, genes)
        if x:
            groups.append(x)
    # leaky TGF-beta hits: species-specific genes outside M
    outside = [s for s in T.leaves(T.root) if s not in T.leaves(STEM_M)]
    for s in outside[:3]:
        genes.append((s, f"{s}_TGFleak", "PF00019"))
    by_sp = {}
    for sp, g, _ in genes:
        by_sp.setdefault(sp, []).append(g)
    xml = ['<?xml version="1.0" encoding="utf-8"?>',
           '<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="synthetic" originVersion="1">']
    for sp in sorted(by_sp):
        xml.append(f'<species name="{sp}" NCBITaxId="0"><database name="syn" version="1"><genes>')
        xml += [f'<gene id="{g}" protId="{g}"/>' for g in by_sp[sp]]
        xml.append("</genes></database></species>")
    xml.append("<groups>" + "".join(groups) + "</groups></orthoXML>")
    d = tmp
    os.makedirs(os.path.join(d, "run", "stages", "input"), exist_ok=True)
    open(os.path.join(d, "run", "hogs.orthoxml"), "w").write("\n".join(xml))
    nwk = make_tree()
    open(os.path.join(d, "run", "stages", "input", "guide_tree.nwk"), "w").write(nwk)
    # Pfam tables
    pdir = os.path.join(d, "pfam")
    os.makedirs(pdir, exist_ok=True)
    hits = pd.DataFrame([(g, pf, sp) for sp, g, pf in genes if pf], columns=["protein", "pfam_acc", "species"])
    for sp in by_sp:
        t = hits[hits["species"] == sp].drop(columns="species").assign(env_from=1, dom_score=50.0)
        with gzip.open(os.path.join(pdir, f"{sp}.pfam.tsv.gz"), "wt") as fh:
            t.to_csv(fh, sep="\t", index=False)
    leaves = T.leaves(T.root)
    tr = pd.DataFrame([dict(transition_id="T_syn", kingdom="syn", label="syn", type="clonal", from_state="unicellular",
                            to_state="simple", analysis_set="primary", primary_eligible=1, M=STEM_M, P=STEM_P,
                            stem_path=f"{STEM_P}>{STEM_M}", mc_species=",".join(T.leaves(STEM_M)))])
    tr.to_csv(os.path.join(d, "transitions.tsv"), sep="\t", index=False)
    traits = pd.DataFrame(dict(code=leaves, data_type="genome", busco_C=95.0))
    traits.to_csv(os.path.join(d, "traits.tsv"), sep="\t", index=False)
    paths = an.run_paths("syn")
    paths.update(run="syn", run_dir=os.path.join(d, "run"), orthoxml=os.path.join(d, "run", "hogs.orthoxml"),
                 tree=os.path.join(d, "run", "stages", "input", "guide_tree.nwk"), cache=os.path.join(d, "cache"),
                 pfam=pdir, transitions=os.path.join(d, "transitions.tsv"), traits=os.path.join(d, "traits.tsv"))
    return T, paths, tr, traits, nwk


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("syn"))
    T, paths, tr, traits, nwk = build_synthetic(d)
    an.build_cache("syn", paths, levels="all", log=lambda m: None)
    cache = an.load_cache("syn", paths)
    ete = bs.load_tree(nwk)
    rng = np.random.default_rng(3)
    for n in ete.traverse():
        n.dist = float(rng.uniform(0.05, 0.5))
    bt = bs.build_branch_table(ete, traits, tr, branch_lengths="tree")
    return dict(T=T, paths=paths, cache=cache, bt=bt, tr=tr, traits=traits)


# -----------------------------------------------------------------------------
# cache / annotation
# -----------------------------------------------------------------------------

def test_cache_events_match_plan(syn):
    br = syn["cache"].branches.set_index("branch_child")
    assert br.at[STEM_M, "dups"] == 12
    assert br.at["n100", "dups"] == 48
    assert br.at["n11", "dups"] == 12
    assert br.at[STEM_M, "gains"] == 2                    # WNT + TGF
    h12 = pd.read_csv(os.path.join(ROOT, "data", "h12_sets.tsv"), sep="\t")
    n_flag = int(((h12["role"] == "member") & h12["module"].isin(an.H12_MODULES)).sum())
    assert br.at[STEM_M, "losses"] == n_flag              # one flagellar family per module Pfam
    assert (br["conservation_ok"]).all()


def test_class_annotation_and_accumulation(syn):
    c = syn["cache"]
    hc = c.hog_classes.set_index("hog_id")
    assert (hc["n_genes"].to_numpy() == c.hogs.set_index("hog_id")["n_genes"].reindex(hc.index).to_numpy()).all()
    # TF families: every member carries PF00046 -> every HOG of a TF family is in tf_dbd
    fam = c.hogs.set_index("hog_id")["family_id"]
    tf_genes = set(c.hits.loc[c.hits["pfam_acc"] == TF, "protein"])
    tf_fams = set(c.genes.loc[c.genes["gene"].isin(tf_genes), "family_id"])
    in_tf = hc["in_tf_dbd"]
    exp = fam.reindex(hc.index).isin(tf_fams)
    assert (in_tf == exp).all()
    assert hc["mod_IFT_B"].sum() > 0 and hc["mod_BBSome"].sum() > 0


def test_lrr_exclusion_rule():
    spec = an.load_class_spec(os.path.join(ROOT, "data", "toolkit_pfam.tsv"), os.path.join(ROOT, "data", "h12_sets.tsv"))
    hits = pd.DataFrame({"protein": ["p1", "p2", "p2", "p3", "p3", "p4"],
                         "pfam_acc": ["PF00560", "PF00560", "PF01582", "PF00028", "PF01582", "PF00096"]})
    f = an.protein_flags(hits, spec)
    assert f.at["p1", "adhesion_ecm"] and not f.at["p2", "adhesion_ecm"]   # LRR alone vs LRR + TIR
    assert f.at["p3", "adhesion_ecm"]                                       # cadherin + TIR stays adhesion
    assert f.at["p4", "c2h2"] and not f.at["p4", "tf_dbd"]                  # C2H2 not in tf_dbd
    assert f.at["p3", "adh_sig"] and not f.at["p4", "adh_sig"]


# -----------------------------------------------------------------------------
# branch statistics
# -----------------------------------------------------------------------------

def test_h4_planted_tf_burst_gives_high_u(syn):
    c, bt = syn["cache"], syn["bt"]
    comp = bs.compositional_statistic(c.events, c.hogs, c.class_hogs("tf_dbd"), bt, min_d=10)
    assert comp.at[STEM_M, "d"] == 12 and comp.at[STEM_M, "a"] == 10
    su = an.run_relabel_suite("H4", comp["s"], bt, {"T_syn": STEM_M}, "greater", scale="difference",
                              extra=comp[["log_d"]], n_draws=2000)
    k = su.match.k
    assert su.relabel.u["T_syn"] == pytest.approx((k + 0.5) / (k + 1))    # stem above every control
    assert su.effects.at["T_syn", "est"] > 2
    assert su.relabel.p < 0.1


def test_wgd_flag_detects_duplication_excess(syn):
    w = an.wgd_flags(syn["cache"], syn["bt"], None)
    assert w.at["n100", "flag_dup_excess"]
    assert w.loc[w["flag_dup_excess"]].index.tolist() == ["n100"]
    comp = bs.compositional_statistic(syn["cache"].events, syn["cache"].hogs, syn["cache"].class_hogs("tf_dbd"),
                                      syn["bt"], min_d=10, exclude=set(w.index[w["wgd_flag"]]))
    assert comp.at["n100", "reason"] == "excluded"


@pytest.mark.parametrize("direction,expect_high", [("loss", True), ("no_loss", False)])
def test_h12_planted_flagellar_loss_sign(syn, direction, expect_high):
    c, bt = syn["cache"], syn["bt"]
    mh, cov = an.h12_module_hogs(c)
    assert set(cov.loc[cov["present"], "module"]) == set(an.H12_MODULES)
    f = bs.module_loss_statistic(c.events, mh, {m: 0 for m in an.H12_MODULES}, bt, hogs=c.hogs)
    assert f.at[STEM_M, "f"] == pytest.approx(1.0)
    assert f.loc[bt.index[bt["in_U"]], "f"].max() == 0
    signs = bs.h12_signs({"T_syn": direction})
    su = an.run_relabel_suite("H12", f["f"], bt, {"T_syn": STEM_M}, "greater", scale="difference", signs=signs,
                              n_draws=2000)
    k = su.match.k
    if expect_high:
        assert su.relabel.u["T_syn"] == pytest.approx((k + 0.5) / (k + 1))
        assert su.effects.at["T_syn", "est"] == pytest.approx(1.0)
    else:
        assert su.relabel.u["T_syn"] == pytest.approx(0.5 / (k + 1))
        assert su.effects.at["T_syn", "est"] == pytest.approx(-1.0)


def test_h1h2_from_ladder_and_bootstrap(syn, tmp_path):
    bt = syn["bt"]
    rng = np.random.default_rng(5)
    internal = bt.index[~bt["is_terminal"]]
    lad = pd.DataFrame({"branch_child": internal, "seq_novel": rng.poisson(5, len(internal))})
    lad.loc[lad["branch_child"] == STEM_M, "seq_novel"] = 60
    lad.loc[lad["branch_child"] == STEM_P, "seq_novel"] = 3      # P edge (not part of the stem)
    p = tmp_path / "ladder.tsv"
    lad.to_csv(p, sep="\t", index=False)
    y = an.h1h2_counts(str(p), bt)
    assert y[STEM_M] == 60                                         # single-edge stem: not summed with P
    rs = bs.count_residual_statistic(y, bt)
    su = an.run_relabel_suite("H1H2", rs["s"], bt, {"T_syn": STEM_M}, "two-sided", effect_values=y,
                              scale="log_ratio", n_draws=2000)
    assert su.relabel.u["T_syn"] > 0.9
    assert su.effects.at["T_syn", "est"] > 1.5
    boot = bs.residual_bootstrap_test(y, bt, su.match, n_boot=19)
    assert 0 < boot["p"] <= 1 and np.isfinite(boot["T_obs"])


def test_path_stem_sums_edges(syn):
    tr = syn["tr"].copy()
    tr["stem_path"] = "n0>n01>n010"            # 2-edge path stem
    ete = bs.load_tree(make_tree())
    bt2 = bs.build_branch_table(ete, syn["traits"], tr)
    y = pd.Series(1.0, index=bt2.index)
    y["n01"], y[STEM_M] = 4.0, 6.0
    assert bs.path_sum(y, bt2)[STEM_M] == 10.0


# -----------------------------------------------------------------------------
# positive-control primitives
# -----------------------------------------------------------------------------

def _ctx(syn):
    sp = pd.DataFrame(dict(code=syn["traits"]["code"], data_type="genome", low_quality=0))
    return pc.Ctx(cache=syn["cache"], species=sp, transitions=syn["tr"], paths=syn["paths"], bt=syn["bt"])


def test_control_gain_rule_recovers_planted_gain(syn):
    ctx = _ctx(syn)
    s = pc.rule_gain(ctx, "Cx", "PF00110", STEM_M)
    assert s.recovered == "yes", s.observed
    assert "gain ladder not available" in s.observed


def test_control_gain_rule_leak_fails(syn):
    ctx = _ctx(syn)
    s = pc.rule_gain(ctx, "Cy", "PF00019", STEM_M)
    assert s.components == {"a": True, "b": False}
    assert s.recovered == "no"
    # wrong node: no gain on stem(n01)
    assert pc.rule_gain(ctx, "Cz", "PF00110", STEM_P).recovered == "no"


def test_control_coption_and_loss_primitives(syn):
    ctx = _ctx(syn)
    L = ctx.L
    # a TF background family is present at P with members on both sides of the stem
    rows = pc.co_option_hogs(ctx, [TF], STEM_P, [ctx.clade(STEM_M), ctx.clade("n011")], origin_at_or_above=STEM_P)
    assert len(rows) > 0
    assert all(L.origin_level(h) == "n" for h in rows)
    # flagellar HOGs: lost on stem(M), present at P
    fl = an.pf_hogs(ctx.cache, ["PF10498"])
    top = [h for h in fl if L.level[h] == 0]
    assert len(top) == 1
    assert L.lost_on_stem(top[0], ctx.node(STEM_M)) and not L.lost_on_stem(top[0], ctx.node("n011"))
    assert L.present_at(top[0], ctx.node(STEM_P)) and not L.present_at(top[0], ctx.node(STEM_M))


# -----------------------------------------------------------------------------
# Holm helper
# -----------------------------------------------------------------------------

def test_holm_primary_compute(tmp_path):
    import holm_primary
    pt = pd.DataFrame(dict(hypothesis=["H1H2", "H4", "H5", "H12"], set="primary", analysis="primary",
                           p=[0.01, 0.04, 0.03, 0.5], p_mcse=0.001, tier=["primary"] * 3 + ["exploratory"],
                           n=7, k=20, scheme="permute"))
    h9 = pd.DataFrame(dict(analysis=["primary"], label=["primary"], scheme=["permute"], n=[7], k=[10], T_obs=[0.6],
                           p=[0.2], p_mcse=[0.004]))
    pt.to_csv(tmp_path / "pt.tsv", sep="\t", index=False)
    h9.to_csv(tmp_path / "h9.tsv", sep="\t", index=False)
    out, missing = holm_primary.compute(str(tmp_path / "pt.tsv"), str(tmp_path / "h9.tsv"))
    assert missing == []
    o = out.set_index("hypothesis")
    assert np.isnan(o.at["H12", "holm_p"])                       # demoted -> not in the family
    assert o.at["H1H2", "holm_p"] == pytest.approx(0.04)         # 4 tests: 4 x 0.01
    assert o.at["H5", "holm_p"] == pytest.approx(0.09)
    assert o.at["H4", "holm_p"] == pytest.approx(0.09)
    assert o.at["H9", "holm_p"] == pytest.approx(0.2)


# -----------------------------------------------------------------------------
# end-to-end: scripts/analysis/primary_tests.py (euk code path) on the synthetic run
# -----------------------------------------------------------------------------

def test_primary_tests_end_to_end(syn, tmp_path, monkeypatch):
    import argparse
    import holm_primary
    import primary_tests
    d = tmp_path
    tr = pd.concat([syn["tr"], pd.DataFrame([
        dict(transition_id="T_syn2", kingdom="syn", label="syn2", type="clonal", from_state="unicellular",
             to_state="simple", analysis_set="primary", primary_eligible=1, M="n110", P="n1",
             stem_path="n1>n11>n110", mc_species=",".join(syn["T"].leaves("n110"))),
        dict(transition_id="T_nested", kingdom="syn", label="nested", type="clonal", from_state="simple",
             to_state="complex", analysis_set="secondary", primary_eligible=0, M="n0100", P=STEM_M,
             stem_path=f"{STEM_M}>n0100", mc_species=",".join(syn["T"].leaves("n0100")))])], ignore_index=True)
    ete = bs.load_tree(make_tree())
    rng = np.random.default_rng(3)
    for n in ete.traverse():
        n.dist = float(rng.uniform(0.05, 0.5))
    bt = bs.build_branch_table(ete, syn["traits"], tr, branch_lengths="tree")
    assert bt.at["n110", "stem_path"] == "n11;n110" and not bt.at["n11", "in_U"]
    internal = bt.index[~bt["is_terminal"]]
    lad = pd.DataFrame({"branch_child": internal, "seq_novel": rng.poisson(5, len(internal))})
    lad.loc[lad["branch_child"] == STEM_M, "seq_novel"] = 60
    lad.to_csv(d / "ladder.tsv", sep="\t", index=False)
    sp = pd.DataFrame(dict(code=syn["traits"]["code"],
                           lineage_group=["GroupA-x" if c.startswith("s0") else "GroupB-x" for c in syn["traits"]["code"]],
                           data_type="genome", low_quality=0))
    sp.to_csv(d / "species.tsv", sep="\t", index=False)
    pd.DataFrame(dict(transition_id=["T_syn", "T_syn2"], expected=["loss", "no_loss"])).to_csv(
        d / "h12dir.tsv", sep="\t", index=False)
    p = dict(syn["paths"], ladder=str(d / "ladder.tsv"), species=str(d / "species.tsv"), h12_dir=str(d / "h12dir.tsv"),
             qc_flags=str(d / "none.tsv"))
    h9 = pd.DataFrame(dict(analysis=["primary"], label=["primary"], scheme=["permute"], n=[2], k=[10], T_obs=[0.5],
                           p=[0.9], p_mcse=[0.003]))
    h9.to_csv(d / "h9.tsv", sep="\t", index=False)
    monkeypatch.setattr(holm_primary, "H9", str(d / "h9.tsv"))
    a = argparse.Namespace(out=str(d / "out"), allow_missing_ladder=False, tests="all", n_draws=500, seed=11,
                           n_boot=19, n_boot_secondary=9, registry=str(d / "reg.tsv"))
    primary_tests.run_euk(a, p, syn["cache"], tr, syn["traits"], bt)
    t = pd.read_csv(d / "out" / "primary_tests.tsv", sep="\t")
    prim = t[(t["set"] == "primary") & (t["analysis"] == "primary")].set_index("hypothesis")
    assert set(prim.index) == {"H1H2", "H4", "H5", "H12"}
    assert prim.at["H1H2", "scheme"] == "parametric_bootstrap"
    assert prim.at["H12", "tier"] == "exploratory"            # only one loss-expected stem -> demoted
    assert prim["holm_p"].notna().sum() == 3                   # H1H2, H4, H5 + H9 (from file) = family of 4
    holm = pd.read_csv(d / "out" / "holm.tsv", sep="\t")
    assert set(holm.loc[holm["holm_p"].notna(), "hypothesis"]) == {"H1H2", "H4", "H5", "H9"}
    assert {"secondary"} <= set(t["set"])
    pt = pd.read_csv(d / "out" / "per_transition.tsv", sep="\t")
    h4 = pt[(pt["hypothesis"] == "H4") & (pt["set"] == "primary")].set_index("transition_id")
    assert h4.at["T_syn", "u"] > 0.9
    h12 = pt[(pt["hypothesis"] == "H12") & (pt["set"] == "primary")].set_index("transition_id")
    assert h12.at["T_syn", "y"] == pytest.approx(1.0) and h12.at["T_syn2", "y"] == pytest.approx(0.0)
    loo = pd.read_csv(d / "out" / "loo.tsv", sep="\t")
    assert {"transition", "supergroup"} <= set(loo["kind"])
    reg = pd.read_csv(d / "reg.tsv", sep="\t")
    assert len(reg) == len(t) and (reg.loc[reg["tier"] == "primary", "hypothesis"].isin(["H1H2", "H4", "H5"])).all()


# -----------------------------------------------------------------------------
# smoke test: every eukaryote control rule on a synthetic run over the real frozen guide tree
# -----------------------------------------------------------------------------

class EteTree(Tree):
    def __init__(self, newick_path):
        import ete3
        t = ete3.Tree(open(newick_path).read().strip(), format=1)
        self.nodes, self.children = [], {}
        for n in t.traverse("preorder"):
            self.nodes.append(n.name)
            self.children[n.name] = [c.name for c in n.children]
        self.root = t.name


def test_euk_control_rules_smoke(tmp_path):
    """Planted: Wnt (PF00110) gained on Metazoa -> C01 yes; no cadherin-Y (PF01049) anywhere -> C03 no;
    pherophorin (PF12499) 12 copies in Volvcar vs 2 in Chlarei -> C12 yes; T-box (PF00907) present at the root and
    lost on Choanoflagellata -> C07 (a),(b) true; every other rule must run without error."""
    T = EteTree(os.path.join(ROOT, "data", "guide_tree_euk.nwk"))
    fams = [Family("WNT", "Metazoa", pfam="PF00110"),
            Family("TBX", T.root, losses=["Choanoflagellata"], dups=["Metazoa"], pfam="PF00907"),
            Family("CAD", T.root, pfam="PF00028"), Family("PTP", T.root, pfam="PF00102"),
            Family("RB", T.root, pfam="PF01858"), Family("KN", T.root, pfam="PF05920"),
            Family("ITG", T.root, losses=["Fungi", "Choanoflagellata"], pfam="PF00362"),
            Family("GRAS", "Anydrophyta", pfam="PF03514"), Family("EIF", T.root, pfam="PF01287")]
    fams += [Family(f"BG{i}", T.root) for i in range(5)]
    genes, groups = [], []
    for f in fams:
        x = emit(T, f, genes)
        if x:
            groups.append(x)
    genes += [("Volvcar", f"Volvcar_PH{i}", "PF12499") for i in range(12)]
    genes += [("Chlarei", f"Chlarei_PH{i}", "PF12499") for i in range(2)]
    genes += [("Neurcra", "Neurcra_P87252", "PF01287")]
    by_sp = {}
    for sp, g, _ in genes:
        by_sp.setdefault(sp, []).append(g)
    leaves = [n for n in T.nodes if not T.children[n]]
    for sp in leaves:
        by_sp.setdefault(sp, [])
    xml = ['<?xml version="1.0" encoding="utf-8"?>',
           '<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="synthetic" originVersion="1">']
    for sp in sorted(by_sp):
        if not by_sp[sp]:
            continue
        xml.append(f'<species name="{sp}" NCBITaxId="0"><database name="syn" version="1"><genes>')
        xml += [f'<gene id="{g}" protId="{g}"/>' for g in by_sp[sp]]
        xml.append("</genes></database></species>")
    xml.append("<groups>" + "".join(groups) + "</groups></orthoXML>")
    d = tmp_path
    os.makedirs(d / "run" / "stages" / "input")
    (d / "run" / "hogs.orthoxml").write_text("\n".join(xml))
    import shutil
    shutil.copy(os.path.join(ROOT, "data", "guide_tree_euk.nwk"), d / "run" / "stages" / "input" / "guide_tree.nwk")
    os.makedirs(d / "pfam")
    hits = pd.DataFrame([(g, pf, sp) for sp, g, pf in genes if pf], columns=["protein", "pfam_acc", "species"])
    for sp in by_sp:
        t = hits[hits["species"] == sp].drop(columns="species").assign(env_from=1, dom_score=50.0)
        with gzip.open(d / "pfam" / f"{sp}.pfam.tsv.gz", "wt") as fh:
            t.to_csv(fh, sep="\t", index=False)
    paths = an.run_paths("euk")
    paths.update(run_dir=str(d / "run"), orthoxml=str(d / "run" / "hogs.orthoxml"),
                 tree=str(d / "run" / "stages" / "input" / "guide_tree.nwk"), cache=str(d / "cache"), pfam=str(d / "pfam"),
                 qc_flags=str(d / "none.tsv"))
    an.build_cache("euk", paths, levels="auto", log=lambda m: None)
    cache = an.load_cache("euk", paths)
    tr, traits, tree, bt = an.load_setup("euk", paths)
    ctx = pc.Ctx(cache=cache, species=pd.read_csv(paths["species"], sep="\t"), transitions=tr, paths=paths, bt=bt)
    ctrl = pd.read_csv(os.path.join(ROOT, "data", "positive_controls.tsv"), sep="\t")
    res = pc.score_all(ctx, ctrl[ctrl["run"] == "euk"], pc.EUK_RULES).set_index("control")
    res.to_csv(d / "res.tsv", sep="\t")
    assert len(res) == 16
    assert not res["observed"].str.contains("missing node").any(), res["observed"].to_dict()
    assert res.at["C01", "recovered"] == "yes"
    assert res.at["C03", "recovered"] == "no"
    assert res.at["C12", "recovered"] == "yes"
    assert res.at["C07", "components"].startswith("a=True;b=True")
    assert res.at["C06", "components"] == "a=True;b=True;c=True" and res.at["C06", "recovered"] == "yes"
    assert res.at["C04", "recovered"] == "yes" and res.at["C11", "recovered"] == "yes"
    assert res.at["C10", "recovered"] == "undetermined"         # needs the gain ladder
    assert res.at["C13", "recovered"] in ("yes", "no")            # anchor found
