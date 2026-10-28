"""Planted toy datasets for the exploratory analyses (mcell/explore.py; scripts/analysis/{pfam_wide,cooption,
architectures,protein_traits,reversals}.py): known expansion / gain / loss / convergence, known co-option fraction,
known new architecture, known trait shift, known preferential loss after a reversal, plus hand-computable cases."""
import math

import numpy as np
import pandas as pd
import pytest

from mcell import branchstats as bs
from mcell import explore as ex


# ----------------------------------------------------------------------------- helpers
def _contrast(M, S, O, P="P"):
    return ex.Contrast("T1", "transition", "M", P, list(M), list(S), list(S) + list(O), list(O), ["M"])


def test_bh_and_hodges_lehmann_hand_computed():
    q = ex.bh([0.01, 0.04, 0.03, np.nan])
    assert np.allclose(q[:3], [0.03, 0.04, 0.04]) and np.isnan(q[3])
    est, lo, hi = ex.hl_paired([1.0, 2.0, 3.0])          # Walsh averages 1, 1.5, 2, 2, 2.5, 3
    assert est == 2.0 and np.isnan(lo)                   # n = 3: no 95 % CI attainable
    est, lo, hi = ex.hl_two_sample([1.0, 2.0], [0.0])    # differences 1, 2
    assert est == 1.5
    d = np.random.default_rng(1).normal(0.5, 1, 400)
    e1 = ex.hl_paired(d)
    e2 = ex.hl_paired(d, max_pairs=20000)                 # sampled Walsh averages
    assert abs(e1[0] - e2[0]) < 0.05 and e1[1] < 0.5 < e1[2]


def test_mwu_exact_enumeration():
    XM = np.array([[3.0], [4.0]])
    XS = np.array([[1.0], [2.0]])
    p = ex.mwu_greater(XM, XS)
    assert p[0] == pytest.approx(1 / 6)                  # C(4,2) = 6 labellings, observed is the extreme one
    # constant column untestable
    assert np.isnan(ex.mwu_greater(np.ones((2, 1)), np.ones((2, 1)))[0])
    # large-sample normal approximation agrees with scipy
    from scipy import stats
    rng = np.random.default_rng(2)
    a, b = rng.poisson(5, (30, 1)).astype(float), rng.poisson(3, (30, 1)).astype(float)
    p2 = ex.mwu_greater(a, b, max_exact=10)[0]
    assert p2 == pytest.approx(stats.mannwhitneyu(a[:, 0], b[:, 0], alternative="greater").pvalue, rel=0.02)


def _toy_counts():
    M = [f"m{i}" for i in range(5)]
    S = [f"s{i}" for i in range(5)]
    O = [f"o{i}" for i in range(4)]
    sp = M + S + O
    X = pd.DataFrame(0, index=sp, columns=["EXP", "GAIN", "LOST", "FLAT", "FOLD", "OUTONLY"])
    X.loc[M, "EXP"] = [6, 7, 5, 6, 8]
    X.loc[S, "EXP"] = [1, 2, 1, 1, 2]
    X.loc[O, "EXP"] = 1
    X.loc[M[:3], "GAIN"] = 2                              # 60 % of M, nowhere else
    X.loc[S, "LOST"] = 1
    X.loc[O, "LOST"] = 1
    X.loc[sp, "FLAT"] = 3
    X.loc[M, "FOLD"] = [4, 0, 4, 4, 0]                    # mean 2.4 vs 0.2: log2 fold 1.5, 60 % of M
    X.loc[S[:1], "FOLD"] = 1
    X.loc[O, "OUTONLY"] = 2
    return X, M, S, O


def test_domain_contrast_planted_expansion_gain_loss():
    X, M, S, O = _toy_counts()
    t = ex.domain_contrast(X, _contrast(M, S, O)).set_index("pfam_acc")
    assert t.at["EXP", "status"] == "expanded" and t.at["EXP", "expanded_by_test"]
    assert t.at["EXP", "p_mwu"] == pytest.approx(1 / math.comb(10, 5))
    assert t.at["GAIN", "status"] == "gained"
    assert t.at["LOST", "status"] == "lost"
    assert t.at["FLAT", "status"] == ""
    assert t.at["FOLD", "status"] == "expanded" and t.at["FOLD", "expanded_by_fold"]
    assert t.at["OUTONLY", "status"] == ""
    # missing (unannotated) species are dropped, not counted as zeros
    t2 = ex.domain_contrast(X.drop(index=["o0", "o1", "o2", "o3"]), _contrast(M, S, O + ["o9"]))
    assert int(t2["n_out"].iloc[0]) == len(S)


def test_convergence_planted_domain():
    rng = np.random.default_rng(3)
    doms = [f"D{i}" for i in range(200)]
    tids = [f"T{i}" for i in range(6)]
    ctrl = [f"c{i}" for i in range(30)]
    obs = pd.DataFrame(rng.random((200, 6)) < 0.05, index=doms, columns=tids)
    obs.loc["D0"] = True                                   # planted: flagged in all 6 transitions
    C = pd.DataFrame(rng.random((200, 30)) < 0.05, index=doms, columns=ctrl)
    sets = {t: ctrl[i * 5:(i + 1) * 5] + ctrl[:3] for i, t in enumerate(tids)}
    per, glob, info = ex.convergence_test(obs, C, sets, n_perm=2000, seed=11)
    per = per.set_index("pfam_acc")
    assert info["scheme"] == "permute" and info["matched"] == "matched"
    assert per.at["D0", "k_obs"] == 6 and per.at["D0", "p_perm"] < 0.01
    noise = per.drop(index="D0")
    assert (noise["p_perm"].dropna() <= 0.05).mean() < 0.08          # permute null: calibrated per domain
    # null calibration of the global statistic: no planted effect
    obs2 = pd.DataFrame(rng.random((200, 6)) < 0.05, index=doms, columns=tids)
    _, g2, _ = ex.convergence_test(obs2, C, sets, n_perm=2000, seed=11)
    assert g2.loc[g2["scheme"] == "permute", "p_perm"].min() > 0.01


def test_coption_index_known_fraction_and_transcriptome_rule():
    M = ["m0", "m1", "m2", "m3"]
    S = ["s0", "s1"]
    O = ["o0", "o1", "o2"]
    sp = M + S + O
    X = pd.DataFrame(0, index=sp, columns=["A", "B", "C", "D", "E"])
    for d in "ABCD":
        X.loc[M, d] = 4                                   # all four expanded in M
    X.loc["s0", "A"] = 1                                  # A present in a sister genome
    X.loc["o0", "B"] = 1                                  # B in an outgroup genome
    X.loc["o1", "C"] = 1                                  # C in one transcriptome only -> not present at P
    X.loc[["o1", "o2"], "D"] = 1                          # D in two transcriptomes -> present
    X.loc[sp, "E"] = 1                                    # E not expanded
    dtype = {"o1": "transcriptome", "o2": "transcriptome"}
    c = _contrast(M, S, O)
    r = ex.coption_index(X, c, dtype)
    assert sorted(r["F"]) == ["A", "B", "C", "D"]
    assert r["n_P"] == 3 and r["C"] == pytest.approx(0.75)
    # present at P's parent: only species outside P (O)
    assert r["n_Pparent"] == 2 and r["C_Pparent"] == pytest.approx(0.5)
    # split-repaired families count as present
    r2 = ex.coption_index(X, c, dtype, extra_present=["C"])
    assert r2["C"] == pytest.approx(1.0)
    # toolkit 'present' rule restricted to a family list
    r3 = ex.coption_index(X, c, dtype, families=["A", "E"], rule="present")
    assert sorted(r3["F"]) == ["A", "E"] and r3["C"] == pytest.approx(1.0)


def _toy_tree():
    return bs.load_tree("(((a,b)X,(c,d)Y)XY,(e,f)Z)R;")


def test_architecture_new_on_branch_planted():
    tree = _toy_tree()
    ann = list("abcdef")
    rows = []
    for s in ann:
        rows += [(f"{s}_1", "PF1", 1), (f"{s}_2", "PF2", 1)]          # old single domains everywhere
    rows += [("a_3", "PF1-PF2", 2), ("b_3", "PF1-PF2", 2)]            # new combination of old domains on X
    rows += [("a_4", "PF9-PF1", 2), ("b_4", "PF9-PF1", 2)]            # contains a domain new on X
    rows += [("c_5", "PF2-PF1", 2)]                                   # single species: no branch
    rows += [("a_6", "PF1-PF2-PF1", 2), ("e_6", "PF1-PF2-PF1", 2)]    # spans the root split -> origin R
    arch = pd.DataFrame(rows, columns=["protein", "arch", "n_domains"])
    arch["species"] = arch["protein"].str[0]
    org = ex.architecture_origins(arch, tree, ann).set_index("arch")
    assert org.at["PF1-PF2", "origin"] == "X" and org.at["PF1-PF2", "recombined"]
    assert org.at["PF9-PF1", "origin"] == "X" and not org.at["PF9-PF1", "recombined"]
    assert org.at["PF2-PF1", "origin"] == ""
    assert org.at["PF1-PF2-PF1", "origin"] == "R"
    bt = bs.branch_table(tree)
    br = ex.branch_architecture_counts(org.reset_index(), bt, tree, ann)
    assert br.at["X", "new_arch_multi"] == 2 and br.at["X", "new_arch_recombined"] == 1
    assert br.at["Y", "new_arch_all"] == 0
    # unannotated species (e, f) are ignored: PF1-PF2-PF1 then has one annotated carrier -> no origin
    org2 = ex.architecture_origins(arch, tree, list("abcd")).set_index("arch")
    assert org2.at["PF1-PF2-PF1", "origin"] == "" and org2.at["PF1-PF2", "origin"] == "X"
    assert br.at["X", "cov_outside"] == 1.0


def test_architectures_from_hits_order_overlap_and_repeats():
    h = pd.DataFrame({"protein": ["p"] * 4, "pfam_acc": ["B", "A", "A", "C"],
                      "ali_from": [300, 10, 110, 305], "ali_to": [400, 100, 200, 395],
                      "dom_score": [50.0, 30.0, 30.0, 20.0]})       # C overlaps B (lower score) -> dropped
    a = ex.architectures_from_hits(h).set_index("protein")
    assert a.at["p", "arch"] == "A-B" and a.at["p", "n_domains"] == 2


def test_paired_trait_shift_planted():
    rng = np.random.default_rng(5)
    M = ["m0", "m1", "m2"]
    S = ["s0", "s1"]
    rows, length, dis = [], {}, {}
    for h in range(60):
        base = rng.uniform(200, 800)
        for s in M + S:
            g = f"{s}_{h}"
            rows.append((h, g, s))
            f = 1.5 if s in M else 1.0
            length[g] = base * f * rng.uniform(0.98, 1.02)
            dis[g] = 0.2 + (0.1 if s in M else 0.0) + rng.normal(0, 0.01)
    rows.append((0, "m0_extra", "m0"))                   # HOG 0 is two-copy in m0 -> excluded
    length["m0_extra"] = 100.0
    memb = pd.DataFrame(rows, columns=["hog_row", "gene", "species"])
    sc = ex.single_copy_hogs(memb, M, S)
    assert 0 not in set(sc["hog_row"]) and sc["hog_row"].nunique() == 59
    med = ex.per_hog_side_medians(sc, pd.Series(length))
    t = ex.paired_test(np.log2(med["M"] / med["S"]).to_numpy())
    assert t["n"] == 59 and t["hl"] == pytest.approx(math.log2(1.5), abs=0.02) and t["p"] < 1e-8
    assert t["hl_lo"] < math.log2(1.5) < t["hl_hi"]
    med = ex.per_hog_side_medians(sc, pd.Series(dis))
    t = ex.paired_test((med["M"] - med["S"]).to_numpy())
    assert t["hl"] == pytest.approx(0.1, abs=0.01)
    # null: no shift -> HL near 0, CI covers 0
    t0 = ex.paired_test(rng.normal(0, 1, 200))
    assert t0["hl_lo"] < 0 < t0["hl_hi"]
    # presence rule: a HOG missing from all sisters is not used
    memb2 = memb[~((memb["hog_row"] == 5) & memb["species"].isin(S))]
    assert 5 not in set(ex.single_copy_hogs(memb2, M, S)["hog_row"])
    # two-sample whole-proteome version
    w = ex.two_sample_test([3.0, 3.1, 3.2], [2.0, 2.1])
    assert w["hl"] == pytest.approx(1.05, abs=1e-9)     # differences 0.9,1,1,1.1,1.1,1.2


def _toy_hog_tables(n_bg=100, n_gain=20, lost_gain=15, lost_bg=10, n_dup=5, lost_dup=4):
    """Species tree R -> A -> B -> {C, D}; R -> E. Stem edge = A. Reversal edge B -> C."""
    names = ["R", "A", "B", "C", "D", "E"]
    is_leaf = [False, False, False, True, True, True]
    hogs, genes, events = [], [], []

    def add_hog(hid, lvl, parent, n_genes):
        hogs.append(dict(hog_id=hid, level_idx=lvl, parent_idx=parent, n_genes=n_genes))
        return len(hogs) - 1

    for k in range(n_bg + n_gain):
        gain = k >= n_bg
        dup = (not gain) and k < n_dup
        lose = (k - n_bg < lost_gain) if gain else (k < lost_dup if dup else (n_dup <= k < n_dup + lost_bg))
        nC = 0 if lose else 1
        ng = nC + 1 + (0 if gain else 1)
        if gain:
            a = add_hog(f"A|g{k}", 1, -1, ng)
            events.append(dict(branch_child="A", event="gain", hog_id=f"A|g{k}", at_origin=True))
        else:
            r = add_hog(f"R|h{k}", 0, -1, ng)
            a = add_hog(f"A|h{k}", 1, r, ng)
            if dup:
                events.append(dict(branch_child="A", event="dup", hog_id=f"R|h{k}", at_origin=False))
            genes.append(dict(gene=f"E_{k}", species="E", parent_idx=r))
        b = add_hog(f"B|{k}", 2, a, ng)
        genes.append(dict(gene=f"D_{k}", species="D", parent_idx=b))
        if nC:
            genes.append(dict(gene=f"C_{k}", species="C", parent_idx=b))
    return names, is_leaf, pd.DataFrame(hogs), pd.DataFrame(genes), pd.DataFrame(events)


def test_reversal_preferential_loss_planted():
    names, is_leaf, hogs, genes, events = _toy_hog_tables()
    r = ex.reversal_test(hogs, genes, events, names, is_leaf, ["A"], "B", "C")
    assert r["n_hogs_parent"] == 120 and r["n_stem_gained"] == 20 and r["n_stem_dup"] == 5
    g = r["gained"]
    assert (g["a_lost"], g["b_kept"], g["c_bg_lost"], g["d_bg_kept"]) == (15, 5, 14, 86)
    assert g["odds_ratio"] > 5 and g["p_greater"] < 1e-4
    d = r["dup"]
    assert (d["a_lost"], d["b_kept"]) == (4, 1)
    # sister control D keeps everything: no preferential loss
    s = ex.reversal_test(hogs, genes, events, names, is_leaf, ["A"], "B", "D")
    assert s["n_lost"] == 0 and np.isnan(s["gained"]["p_greater"]) or s["gained"]["p_greater"] == 1.0
    # null: losses independent of stem status
    names, is_leaf, hogs, genes, events = _toy_hog_tables(n_bg=100, n_gain=20, lost_gain=2, lost_bg=10,
                                                         n_dup=0, lost_dup=0)
    r0 = ex.reversal_test(hogs, genes, events, names, is_leaf, ["A"], "B", "C")
    assert r0["gained"]["p_greater"] > 0.2


def test_hog_domain_events_and_members():
    names, is_leaf, hogs, genes, events = _toy_hog_tables(n_bg=4, n_gain=2, lost_gain=0, lost_bg=0, n_dup=1,
                                                         lost_dup=0)
    pdm = pd.DataFrame({"protein": ["D_4", "C_4", "D_0", "C_0", "E_0", "D_1"],
                        "pfam_acc": ["PF7", "PF7", "PF8", "PF8", "PF8", "PF8"]})
    gdom = ex.gene_domain_codes(genes, pdm)
    row_of = pd.Series(np.arange(len(hogs)), index=hogs["hog_id"].to_numpy())
    t = ex.hog_domain_events(events, hogs, genes, gdom, ["A"], row_of)
    assert t.at["PF7", "gain"] == 1 and t.at["PF7", "dup"] == 0     # gained HOG A|g4: both genes carry PF7
    assert t.at["PF8", "dup"] == 1                                    # R|h0 duplicated on A: 3/3 genes PF8
    m = ex.members_of(hogs["parent_idx"].to_numpy(), genes["parent_idx"].to_numpy(), [row_of["R|h1"]])
    assert set(genes["gene"].to_numpy()[m["gene_row"]]) == {"E_1", "D_1", "C_1"}


def test_gene_hog_at_and_present_rule():
    names, is_leaf, hogs, genes, events = _toy_hog_tables(n_bg=3, n_gain=1, lost_gain=0, lost_bg=0, n_dup=0,
                                                         lost_dup=0)
    r = ex.gene_hog_at(hogs, genes, 1)
    assert (r[genes["species"].isin(["C", "D"]).to_numpy()] >= 0).all()
    assert (r[(genes["species"] == "E").to_numpy()] == -1).all()
    X = pd.DataFrame({"F": [1, 0, 1]}, index=["g1", "t1", "t2"])
    assert not ex.present_rule(X, ["t1", "t2"], {"t1": "transcriptome", "t2": "transcriptome"})["F"]
    assert ex.present_rule(X, ["g1"], {})["F"]


def test_control_contrasts_and_pool_on_toy_tree():
    tree = _toy_tree()
    tr = pd.DataFrame([dict(transition_id="T_x", kingdom="toy", from_state="unicellular", to_state="simple",
                            M="X", P="XY", analysis_set="primary", sister_species="c,d", stem_path="")])
    bt = bs.build_branch_table(tree, None, tr)
    cons = ex.transition_contrasts(tr, tree, bt, ["T_x"])
    assert cons[0].M == ["a", "b"] and cons[0].S == ["c", "d"] and cons[0].outside_P == ["e", "f"]
    pool, label = ex.control_pool(bt, min_size=10)
    assert label.startswith("relaxed")
    cc = ex.control_contrasts(tree, bt, ["Y"])
    assert cc[0].M == ["c", "d"] and cc[0].S == ["a", "b"] and cc[0].P == "XY"


def test_own_stem_path_wins_when_transitions_share_M():
    """When transitions share M, build_branch_table uses the highest-priority (primary) row's stem_path (fixed
    2026-10-28; previously the first row shadowed it); the exploratory layer's local patch is then a no-op."""
    tree = bs.load_tree("((((a:1,b:1)GV:1,c:1)TGV:2,d:1)CT:1,(e:1,f:1)Z:1)R;")
    tr = pd.DataFrame([
        dict(transition_id="T_first", kingdom="toy", from_state="unicellular", to_state="simple", M="GV", P="TGV",
             analysis_set="exploratory", sister_species="c", stem_path="TGV>GV"),
        dict(transition_id="T_path", kingdom="toy", from_state="unicellular", to_state="simple", M="GV", P="CT",
             analysis_set="primary", sister_species="d", stem_path="CT>TGV>GV")])
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bt = bs.build_branch_table(tree, None, tr, branch_lengths="tree")
    assert bt.at["GV", "stem_path"] == "TGV;GV"                   # primary row's path wins
    cons = ex.transition_contrasts(tr, tree, bt, ["T_path"])
    assert cons[0].edges == ["TGV", "GV"] and cons[0].P == "CT"
    bt2 = ex.patch_path_stems(bt, cons)
    assert bs.path_edges(bt2, "GV") == ["TGV", "GV"] and bt2.at["GV", "branch_length"] == 3.0
    assert bt.at["GV", "branch_length"] == 3.0                    # path length summed in branchstats itself
