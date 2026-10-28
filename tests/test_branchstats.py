"""Tests for mcell.branchstats (DESIGN.md §4.6).

(a) hand-computable u and T on a 3-transition toy (+ exact enumeration of the null)
(b) null simulations: p ~ Uniform (KS over >= 200 replicate datasets), FPR ~ alpha
(c) planted 2-fold effect on stems of a 10-transition toy tree is detected
(d) control-pool exclusion rules on a hand-made tree
(e) Holm and REML/DL against statsmodels / metafor values
plus statistic helpers, matching, power/MDE helpers and the registry.
"""
import itertools
import math
import random
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mcell import branchstats as bs


# -----------------------------------------------------------------------------
# (a) hand-computable toy
# -----------------------------------------------------------------------------
TOY_VALUES = pd.Series({"A": 5, "B": 1, "C": 3,
                        "c1": 1, "c2": 2, "c3": 5, "c4": 1, "c5": 0}, dtype=float)
TOY_STEMS = {"tA": "A", "tB": "B", "tC": "C"}
TOY_SETS = {"tA": ["c1", "c2", "c3"], "tB": ["c3", "c4", "c5"], "tC": ["c1", "c4", "c5"]}


def test_u_score_hand():
    # A=5 vs {1,2,5}: 2 below, 1 tie -> (2 + .5 + .5)/4
    assert bs.u_score(5, [1, 2, 5]) == pytest.approx(3 / 4)
    assert bs.u_score(1, [5, 1, 0]) == pytest.approx(2 / 4)
    assert bs.u_score(3, [1, 1, 0]) == pytest.approx(3.5 / 4)
    assert bs.u_score(0, []) == pytest.approx(0.5)


def _toy_exact_null(alternative):
    """Enumerate all relabellings with distinct controls (uniform), spec scheme."""
    v = TOY_VALUES
    Ts = []
    for tup in itertools.product(*[TOY_SETS[t] for t in TOY_STEMS]):
        if len(set(tup)) < len(tup):
            continue
        us = []
        for t, c in zip(TOY_STEMS, tup):
            others = [v[x] for x in TOY_SETS[t] if x != c]
            us.append(bs.u_score(v[c], others))
        Ts.append(np.mean(us))
    return np.array(Ts)


def test_toy_T_obs_and_null():
    res = bs.relabel_test(TOY_VALUES, TOY_SETS, stems=TOY_STEMS, alternative="greater",
                          n_draws=10_000, seed=11, scheme="spec")
    assert res.u["tA"] == pytest.approx(0.75)
    assert res.u["tB"] == pytest.approx(0.5)
    assert res.u["tC"] == pytest.approx(0.875)
    assert res.T_obs == pytest.approx((0.75 + 0.5 + 0.875) / 3)
    # one relabelling by hand: (c3, c4, c5) -> u = 5/6, 1/2, 1/6 -> T = 1/2
    Ts = _toy_exact_null("greater")
    assert len(Ts) == bs.count_relabellings(list(TOY_SETS.values()))
    assert np.any(np.isclose(Ts, 0.5))
    exact = np.mean(Ts >= res.T_obs - 1e-12)
    # MC p agrees with the enumerated null within 4 MC SE
    assert abs(res.p - exact) < 4 * math.sqrt(exact * (1 - exact) / 10_000) + 2e-4
    assert res.p_mcse == pytest.approx(math.sqrt(res.p * (1 - res.p) / 10_000))
    # duplicate draws are rejected: every accepted T is one of the enumerated values
    assert np.all(np.isclose(res.T_null[:, None], np.unique(Ts)[None, :]).any(axis=1))
    # deterministic with the default seed
    res2 = bs.relabel_test(TOY_VALUES, TOY_SETS, stems=TOY_STEMS, alternative="greater", n_draws=10_000,
                           scheme="spec")
    assert res2.p == res.p and res.seed == 11
    # two-sided uses |T - 1/2|
    r2 = bs.relabel_test(TOY_VALUES, TOY_SETS, stems=TOY_STEMS, alternative="two-sided", n_draws=10_000,
                         scheme="spec")
    ex2 = np.mean(np.abs(Ts - 0.5) >= abs(res.T_obs - 0.5) - 1e-12)
    assert abs(r2.p - ex2) < 4 * math.sqrt(ex2 * (1 - ex2) / 10_000) + 2e-4


def test_toy_signs_and_less():
    # w = -1 for every stratum mirrors the statistic: 'greater' with signs == 'less' without
    signs = {t: -1.0 for t in TOY_STEMS}
    a = bs.relabel_test(TOY_VALUES, TOY_SETS, stems=TOY_STEMS, alternative="greater", signs=signs)
    b = bs.relabel_test(-TOY_VALUES, TOY_SETS, stems=TOY_STEMS, alternative="greater")
    assert a.T_obs == pytest.approx(b.T_obs) and a.p == b.p
    # u(-x) = 1 - u(x) exactly, so 'less' on x == 'greater' on -x (same draws)
    c = bs.relabel_test(TOY_VALUES, TOY_SETS, stems=TOY_STEMS, alternative="less")
    assert c.T_obs == pytest.approx(1 - b.T_obs) and c.p == b.p


def test_min_attainable_p():
    f = bs.min_attainable_p([20, 20], n_draws=10_000)
    assert f["n_relabellings"] == 400 and f["p_floor_exact"] == pytest.approx(1 / 401)
    assert f["min_attainable_p"] == pytest.approx(1 / 401)
    f = bs.min_attainable_p([["a", "b"], ["a", "b"]], n_draws=10_000)
    assert f["n_relabellings"] == 2               # (a,b), (b,a)
    f = bs.min_attainable_p([list("abcdefghij")] * 6, n_draws=10_000)
    assert f["n_relabellings"] == 10 * 9 * 8 * 7 * 6 * 5
    assert f["min_attainable_p"] == pytest.approx(1 / 10_001)


def test_nan_handling():
    v = TOY_VALUES.copy()
    v["B"] = np.nan
    v["c1"] = np.nan
    r = bs.relabel_test(v, TOY_SETS, stems=TOY_STEMS, n_draws=500)
    assert r.n == 2 and r.dropped["tB"] == "stem_value_missing"
    assert r.k["tA"] == 2


# -----------------------------------------------------------------------------
# synthetic configurations for (b)
# -----------------------------------------------------------------------------

def _synthetic_match(rng, n_pool=50, n_stems=7, k=20, overlap="nearest"):
    pool = [f"u{i}" for i in range(n_pool)]
    stems = {f"T{i}": f"s{i}" for i in range(n_stems)}
    if overlap == "nearest":
        X = pd.DataFrame(rng.normal(size=(n_pool + n_stems, 3)), columns=["x1", "x2", "x3"],
                         index=pool + list(stems.values()))
        X["in_U"] = X.index.isin(pool)
        m = bs.match_controls(X, stems, covariates=["x1", "x2", "x3"], k=k)
        return m.sets, stems
    sets = {}
    for i, t in enumerate(stems):
        sets[t] = pool[i * k:(i + 1) * k]       # disjoint
    return sets, stems


def _null_pvalues(scheme, overlap, reps, n_stems, n_pool, seed, discrete=False, n_draws=1000,
                  alternative="two-sided"):
    rng = np.random.default_rng(seed)
    ps = []
    for r in range(reps):
        sets, stems = _synthetic_match(rng, n_pool=n_pool, n_stems=n_stems, overlap=overlap)
        ids = sorted({c for s in sets.values() for c in s}) + list(stems.values())
        y = rng.poisson(5, len(ids)) if discrete else rng.normal(size=len(ids))
        v = pd.Series(y, index=ids, dtype=float)
        ps.append(bs.relabel_test(v, sets, stems=stems, alternative=alternative, n_draws=n_draws,
                                  seed=int(rng.integers(2**31)), scheme=scheme).p)
    return np.array(ps)


@pytest.mark.parametrize("alternative,seed", [("two-sided", 101), ("greater", 102)])
def test_null_uniform_permute_overlapping_sets(alternative, seed):
    """(b) realistic overlap (|U| = 50, k = 20, n = 7, Mahalanobis-nearest sets)."""
    ps = _null_pvalues("permute", "nearest", reps=300, n_stems=7, n_pool=50, seed=seed,
                       alternative=alternative)
    assert stats.kstest(ps, "uniform").pvalue > 0.01
    for a in (0.05, 0.10):
        se = math.sqrt(a * (1 - a) / len(ps))
        assert abs(np.mean(ps <= a) - a) <= 2.5 * se


def test_null_uniform_spec_disjoint_sets():
    """(b) the literal §4.6 scheme is calibrated when the N_i do not overlap."""
    ps = _null_pvalues("spec", "disjoint", reps=300, n_stems=3, n_pool=60, seed=202)
    assert stats.kstest(ps, "uniform").pvalue > 0.01
    for a in (0.05, 0.10):
        se = math.sqrt(a * (1 - a) / len(ps))
        assert abs(np.mean(ps <= a) - a) <= 2.5 * se


def test_null_discrete_counts_fpr():
    """(b) NB/Poisson-like counts with ties: FPR at alpha not inflated (permute)."""
    ps = _null_pvalues("permute", "nearest", reps=250, n_stems=7, n_pool=50, seed=303, discrete=True)
    se = math.sqrt(0.05 * 0.95 / len(ps))
    assert np.mean(ps <= 0.05) <= 0.05 + 2.5 * se


def test_spec_scheme_anticonservative_with_overlap():
    """Documents a property of the literal §4.6 null (see docs/branchstats.md):
    with heavily overlapping N_i (|U| = 30, k = 20, n = 7, as on the eukaryote
    guide tree) the without-replacement draws make T_d less variable than
    T_obs, so the 'spec' scheme rejects too often; 'permute' does not."""
    ps_spec = _null_pvalues("spec", "nearest", reps=250, n_stems=7, n_pool=30, seed=404)
    ps_perm = _null_pvalues("permute", "nearest", reps=250, n_stems=7, n_pool=30, seed=404)
    se = math.sqrt(0.05 * 0.95 / 250)
    assert np.mean(ps_spec <= 0.05) > 0.05 + 2.5 * se
    assert abs(np.mean(ps_perm <= 0.05) - 0.05) <= 2.5 * se


# -----------------------------------------------------------------------------
# toy tree with 10 transitions for (c)
# -----------------------------------------------------------------------------

def _random_tree(n_tips, seed):
    import ete3
    random.seed(seed)
    t = ete3.Tree()
    t.populate(n_tips, random_branches=False)
    for i, n in enumerate(t.traverse("preorder")):
        n.name = (f"sp{i}" if n.is_leaf() else f"n{i}")
    return t


def _toy_tree_with_transitions(n_tr=10, n_tips=220, seed=5):
    t = _random_tree(n_tips, seed)
    rng = random.Random(seed)
    cand = [n for n in t.traverse() if not n.is_leaf() and n.up is not None and n.up.up is not None
            and 3 <= len(n) <= 8 and len(n.up) - len(n) >= 2]
    rng.shuffle(cand)
    chosen, blocked = [], set()
    for n in cand:
        # M clade, its parent and the sister clade must not touch earlier choices
        zone = {x.name for x in n.up.traverse()} | {n.up.up.name}
        if zone & blocked:
            continue
        chosen.append(n)
        blocked |= zone
        if len(chosen) == n_tr:
            break
    assert len(chosen) == n_tr
    rows, state = [], {}
    for i, n in enumerate(chosen):
        rows.append(dict(transition_id=f"T{i}", kingdom="euk", from_state="unicellular",
                         to_state="simple", primary_eligible=1, M=n.name, P=n.up.name))
        for x in n.get_leaf_names():
            state[x] = "simple"
    tr = pd.DataFrame(rows)
    traits = pd.DataFrame({"code": t.get_leaf_names()})
    traits["mc_state"] = traits["code"].map(state).fillna("unicellular")
    traits["data_type"] = ["transcriptome" if i % 5 == 0 else "genome" for i in range(len(traits))]
    traits["busco_C"] = 90.0
    return t, traits, tr


def test_planted_effect_detected_10_transitions():
    """(c) a 2-fold multiplicative effect on 10 stems is detected with high power."""
    t, traits, tr = _toy_tree_with_transitions()
    bt = bs.build_branch_table(t, traits, tr)
    st = bs.stem_branches(bt)
    assert len(st) == 10
    assert bt["in_U"].sum() >= 40
    m = bs.match_controls(bt, st)
    assert m.k == 20
    rng = np.random.default_rng(7)
    mu, disp = 20.0, 0.1
    y = pd.Series(rng.negative_binomial(1 / disp, (1 / disp) / (1 / disp + mu), len(bt)),
                  index=bt.index, dtype=float)
    pw = bs.simulate_power(y, m, effects=(1.0, 2.0), n_sim=40, n_draws=1000, seed=11,
                           alternative="greater", scheme="permute")
    p1 = pw.set_index("effect")
    assert p1.loc[2.0, "power_0.01"] >= 0.9
    assert p1.loc[2.0, "power_0.05"] >= 0.9
    assert p1.loc[1.0, "power_0.05"] <= 0.2
    # literal spec scheme also detects it
    pw2 = bs.simulate_power(y, m, effects=(2.0,), n_sim=20, n_draws=1000, seed=12, scheme="spec")
    assert bs.relabel_test(y, m, n_draws=100).scheme == "permute"          # default (D6)
    assert pw2["power_0.01"].iloc[0] >= 0.9
    d = bs.mde(pd.DataFrame({"effect": [1.0, 1.5, 2.0], "power_0.05": [0.05, 0.6, 1.0]}))
    assert 1.5 < d["alpha_0.05"] < 2.0
    # effects table: log O/E of the planted stems is about log 2
    vals = y.copy()
    for tid, b in st.items():
        vals[b] = 2 * y.reindex(m.sets[tid]).mean()
    eff = bs.effect_table(vals, m)
    assert np.allclose(eff["est"], np.log((eff["y"] + 0.5) / (eff["expected"] + 0.5)))
    pooled = bs.random_effects(eff["est"], eff["var"])
    assert pooled["ci_low"] < math.log(2) < pooled["ci_high"]


# -----------------------------------------------------------------------------
# (d) control pool exclusion on a hand-made tree
# -----------------------------------------------------------------------------
#
#  R
#  |- A                       root child
#  |  |- A1 (A1a(t1,t2), A1b(t3,t4))
#  |  `- A2                   P of the grey transition G
#  |     |- A2a (t5,t6)       sister of the grey crown  -> adjacent
#  |     `- A2b (g1,g2)       grey crown (grey stem)
#  `- B                       root child
#     |- B1                   P of the multicellular transition -> adjacent (edge into P)
#     |  |- M (Ma(m1,m2), Mb(m3,m4))   multicellular crown (stem), Ma/Mb within-mc
#     |  `- S (s1,s2)         sister -> adjacent
#     `- B2
#        |- B2a (u1,u2)
#        `- B2b (u3, X(u4,u5))  X has a 1-tip sister -> excluded
HAND = ("((((t1,t2)A1a,(t3,t4)A1b)A1,((t5,t6)A2a,(g1,g2)A2b)A2)A,"
        "((((m1,m2)Ma,(m3,m4)Mb)M,(s1,s2)S)B1,((u1,u2)B2a,(u3,(u4,u5)X)B2b)B2)B)R;")


def _hand():
    t = bs.load_tree(HAND)
    tips = t.get_leaf_names()
    st = {x: ("complex" if x.startswith("m") else "grey" if x.startswith("g") else "unicellular")
          for x in tips}
    traits = pd.DataFrame({"code": tips, "mc_state": [st[x] for x in tips],
                           "data_type": ["transcriptome" if x in ("t1", "u5") else "genome" for x in tips],
                           "busco_C": 95.0})
    tr = pd.DataFrame([
        dict(transition_id="T_mc", from_state="unicellular", to_state="complex", primary_eligible=1,
             M="M", P="B1"),
        dict(transition_id="T_grey", from_state="unicellular", to_state="grey", primary_eligible=0,
             M="A2b", P="A2"),
    ])
    return t, traits, tr


def test_branch_table_structure():
    t, traits, tr = _hand()
    bt = bs.build_branch_table(t, traits, tr)
    assert len(bt) == len(list(t.traverse())) - 1
    assert bt.loc["A1", "n_desc"] == 4 and bt.loc["A1", "n_sister"] == 4
    assert bt.loc["X", "n_sister"] == 1 and bt.loc["X", "depth"] == 4
    assert bt.loc["A", "root_child"] and not bt.loc["A1", "root_child"]
    assert bt.loc["t1", "is_terminal"] and not bt.loc["A1a", "is_terminal"]
    assert bt.loc["A1a", "frac_tx_desc"] == pytest.approx(0.5)
    assert bt.loc["A1b", "frac_tx_sister"] == pytest.approx(0.5)
    assert bt.loc["B2b", "frac_tx_desc"] == pytest.approx(1 / 3)
    assert bt.loc["A1", "mean_busco_desc"] == pytest.approx(95.0)
    assert bt.loc["A1", "parent"] == "A"
    assert bt["log_bl"].isna().all()           # guide tree without lengths
    bt2 = bs.branch_table(t, traits, branch_lengths="tree")
    assert np.allclose(bt2["branch_length"], 1.0)


def test_control_pool_exclusions():
    t, traits, tr = _hand()
    bt = bs.build_branch_table(t, traits, tr)
    cls = bt["branch_class"]
    assert cls["M"] == "stem" and bt.loc["M", "primary"] and bt.loc["M", "transition_id"] == "T_mc"
    assert cls["Ma"] == "within_multicellular" and cls["m1"] == "within_multicellular"
    assert cls["A2b"] == "grey_zone" and cls["g1"] == "grey_zone"
    assert cls["B1"] == "within_unicellular" and cls["S"] == "within_unicellular"
    U = set(bt.index[bt["in_U"]])
    assert U == {"A1", "A1a", "A1b", "B2", "B2a", "B2b"}
    ex = bt["U_exclusion"]
    assert ex["A"] == "root_child" and ex["B"] == "root_child"
    assert ex["B1"] == "adjacent_to_stem"       # edge into P
    assert ex["S"] == "adjacent_to_stem"        # P -> sister
    assert ex["A2a"] == "adjacent_to_stem"      # sister of the grey crown
    assert ex["M"] == "class:stem"
    assert ex["Ma"] == "class:within_multicellular" and ex["Mb"] == "class:within_multicellular"
    assert ex["A2b"] == "class:grey_zone"
    assert ex["X"] == "sister_tips<2"
    assert ex["t1"] == "terminal"


def test_control_pool_edge_into_grey_P():
    t, traits, tr = _hand()
    bt = bs.build_branch_table(t, traits, tr)
    # A2 is P of the grey transition: the edge A -> A2 is "the edge into P"
    assert bt.loc["A2", "U_exclusion"] == "adjacent_to_stem"


def test_node_states_reversal_and_nested():
    t = bs.load_tree("((((a1,a2)AA,((r1,r2)RR,(b1,b2)BB)RB)MC,(c1,c2)C)P,(d1,(e1,e2)E)D)R;")
    tr = pd.DataFrame([
        dict(transition_id="T1", from_state="unicellular", to_state="simple", primary_eligible=1, M="MC", P="P"),
        dict(transition_id="T2", from_state="simple", to_state="complex", primary_eligible=1, M="AA", P="MC"),
        dict(transition_id="Trev", from_state="simple", to_state="reversal", primary_eligible=0, M="RR", P="RB"),
    ])
    bt = bs.build_branch_table(t, None, tr)
    assert bt.loc["MC", "branch_class"] == "stem"
    assert bt.loc["AA", "branch_class"] == "nested_stem" and bt.loc["AA", "primary"]
    assert bt.loc["RR", "branch_class"] == "other" and bt.loc["r1", "branch_class"] == "other"
    assert bt.loc["BB", "branch_class"] == "within_multicellular"
    st = bs.stem_branches(bt)
    assert dict(st) == {"T1": "MC", "T2": "AA"}


# -----------------------------------------------------------------------------
# matching
# -----------------------------------------------------------------------------

def test_matching_mahalanobis_and_smd():
    rng = np.random.default_rng(3)
    n = 40
    X = pd.DataFrame({"x1": rng.normal(size=n), "x2": rng.normal(size=n) * 10,
                      "c": np.ones(n)}, index=[f"u{i}" for i in range(n)])
    X.loc["s0"] = [0.0, 0.0, 1.0]
    X["in_U"] = X.index != "s0"
    m = bs.match_controls(X, {"T": "s0"}, covariates=["x1", "x2", "c"], k=5)
    assert m.dropped_covariates == {"c": "constant_in_pool"}
    Z = X.loc[X["in_U"], ["x1", "x2"]]
    Zs = (Z - Z.mean()) / Z.std(ddof=1)
    z0 = (X.loc["s0", ["x1", "x2"]] - Z.mean()) / Z.std(ddof=1)
    Si = np.linalg.pinv(np.cov(Zs.to_numpy(), rowvar=False))
    d = Zs.to_numpy() - z0.to_numpy()
    d2 = np.einsum("ij,jk,ik->i", d, Si, d)
    assert m.sets["T"] == list(Zs.index[np.argsort(d2)[:5]])
    row = m.smd[(m.smd.transition_id == "T") & (m.smd.covariate == "x2")].iloc[0]
    assert row["smd"] == pytest.approx((0 - X.loc[m.sets["T"], "x2"].mean()) / Z["x2"].std(ddof=1))
    assert m.k == 5
    m2 = bs.match_controls(X, {"T": "s0"}, covariates=["x1", "x2"], k=100)
    assert m2.k == n


# -----------------------------------------------------------------------------
# (e) Holm, REML, DL
# -----------------------------------------------------------------------------

def test_holm_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(0)
    for _ in range(20):
        p = rng.uniform(0, 0.2, size=rng.integers(1, 8))
        assert np.allclose(bs.holm(p), multipletests(p, method="holm")[1])
    assert np.allclose(bs.holm([0.01, 0.04, 0.03, 0.005]), [0.03, 0.06, 0.06, 0.02])
    out = bs.holm([0.01, np.nan, 0.02])
    assert np.isnan(out[1]) and np.allclose(out[[0, 2]], [0.02, 0.02])


BCG = [(4, 119, 11, 128), (6, 300, 29, 274), (3, 228, 11, 209), (62, 13536, 248, 12619),
       (33, 5036, 47, 5761), (180, 1361, 372, 1079), (8, 2537, 10, 619), (505, 87886, 499, 87892),
       (29, 7470, 45, 7232), (17, 1699, 65, 1600), (186, 50448, 141, 27197), (5, 2493, 3, 2338),
       (27, 16886, 29, 17825)]


def _bcg():
    y = np.array([math.log((a / (a + b)) / (c / (c + d))) for a, b, c, d in BCG])
    v = np.array([1 / a - 1 / (a + b) + 1 / c - 1 / (c + d) for a, b, c, d in BCG])
    return y, v


def test_reml_against_metafor_bcg():
    # metafor: rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg) (REML)
    y, v = _bcg()
    r = bs.random_effects(y, v, method="REML")
    assert r["method"] == "REML"
    assert r["tau2"] == pytest.approx(0.3132, abs=1e-4)
    assert r["mu"] == pytest.approx(-0.7145, abs=1e-4)
    assert r["se"] == pytest.approx(0.1798, abs=1e-4)
    assert r["ci_low"] == pytest.approx(-1.0669, abs=1e-4)
    assert r["ci_high"] == pytest.approx(-0.3622, abs=1e-4)
    # REML optimum agrees with a brute-force grid of the restricted likelihood
    grid = np.linspace(0, 2, 200001)
    ll = [bs._reml_negll(t2, y, v) for t2 in grid[::100]]
    assert abs(grid[::100][int(np.argmin(ll))] - r["tau2"]) < 2e-3


def test_dl_against_statsmodels():
    from statsmodels.stats.meta_analysis import combine_effects
    y, v = _bcg()
    r = bs.random_effects(y, v, method="DL")
    ref = combine_effects(y, v, method_re="dl")
    assert r["tau2"] == pytest.approx(ref.tau2)
    assert r["mu"] == pytest.approx(ref.mean_effect_re)
    assert r["se"] == pytest.approx(ref.sd_eff_w_re)
    assert r["tau2"] == pytest.approx(0.3088, abs=1e-4)      # metafor method="DL"


def test_reml_homogeneous_and_hksj():
    y = np.array([0.1, 0.12, 0.09, 0.11])
    v = np.array([0.04, 0.05, 0.03, 0.04])
    r = bs.random_effects(y, v)
    assert r["tau2"] == 0.0
    w = 1 / v
    assert r["mu"] == pytest.approx(np.sum(w * y) / np.sum(w))
    h = bs.random_effects(y, v, ci_method="hksj")
    assert h["ci_high"] - h["ci_low"] > 0 and h["ci_method"] == "hksj"


def test_leave_one_out():
    y, v = _bcg()
    ids = [f"t{i}" for i in range(len(y))]
    eff = pd.DataFrame({"est": y, "var": v}, index=ids)
    fn = lambda keep: bs.random_effects(eff.loc[keep, "est"], eff.loc[keep, "var"])
    lo = bs.leave_one_out(ids, fn)
    assert len(lo) == 14 and lo.iloc[0]["left_out"] == "none" and lo.iloc[1]["n_kept"] == 12
    assert lo.iloc[0]["mu"] == pytest.approx(-0.7145, abs=1e-4)
    groups = {t: ("g1" if i < 5 else "g2") for i, t in enumerate(ids)}
    lg = bs.leave_one_out(ids, fn, groups=groups)
    assert list(lg["left_out"]) == ["none", "g1", "g2"] and list(lg["n_kept"]) == [13, 8, 5]


# -----------------------------------------------------------------------------
# statistic helpers
# -----------------------------------------------------------------------------

def _events():
    rows = []
    # branch u->b1: 12 dup HOGs (4 in class C), 2 gains (1 seq-novel), 3 losses
    for i in range(12):
        rows.append(("u", "b1", "dup", 1, False, f"u|h{i}"))
    rows.append(("u", "b1", "dup", 2, True, "b1|new0"))     # at-origin: ignored
    rows += [("u", "b1", "gain", 1, True, "b1|g0"), ("u", "b1", "gain", 1, True, "b1|g1")]
    rows += [("u", "b1", "loss", 1, False, f"u|h{i}") for i in (20, 21, 22)]
    # branch u->b2: 5 dup HOGs only (d < 10)
    rows += [("u", "b2", "dup", 1, False, f"u|h{i}") for i in range(5)]
    return pd.DataFrame(rows, columns=["branch_parent", "branch_child", "event", "n_events",
                                       "at_origin", "hog_id"])


def _hogs_at_u():
    # 40 HOGs at u, 10 of them in C (h0..h3 and h30..h35)
    ids = [f"u|h{i}" for i in range(40)]
    return pd.DataFrame({"hog_id": ids, "level": "u"})


def test_count_statistic():
    e = _events()
    s = bs.count_statistic(e, "gain", branches=["b1", "b2", "b3"])
    assert s.to_dict() == {"b1": 2.0, "b2": 0.0, "b3": 0.0}
    s = bs.count_statistic(e, "gain", hogs={"b1|g1"})
    assert s["b1"] == 1
    d = bs.count_statistic(e, "dup")
    assert d["b1"] == 12 and d["b2"] == 5


def test_compositional_statistic():
    e, h = _events(), _hogs_at_u()
    C = {f"u|h{i}" for i in list(range(4)) + list(range(30, 36))}
    bt = pd.DataFrame({"parent": ["u", "u"]}, index=["b1", "b2"])
    r = bs.compositional_statistic(e, h, C, bt=bt)
    pi = 10 / 40
    expect = math.log((4.5 / 13) / (1 - 4.5 / 13)) - math.log(pi / (1 - pi))
    assert r.loc["b1", "d"] == 12 and r.loc["b1", "a"] == 4
    assert r.loc["b1", "s"] == pytest.approx(expect)
    assert np.isnan(r.loc["b2", "s"]) and r.loc["b2", "reason"] == "d<10"
    r2 = bs.compositional_statistic(e, h, C, bt=bt, exclude={"b1"})
    assert r2.loc["b1", "reason"] == "excluded"
    mem = pd.DataFrame({"hog_id": ["x", "x", "x", "y", "y"], "gene": ["g1", "g2", "g3", "g4", "g5"]})
    assert bs.class_hogs_from_members(mem, {"g1", "g2", "g4"}) == {"x", "y"}
    assert bs.class_hogs_from_members(mem, {"g1"}) == set()


def test_module_loss_statistic_and_signs():
    e = _events()
    mh = pd.DataFrame({"level": "u", "hog_id": [f"u|h{i}" for i in (20, 21, 23, 24, 22, 25)],
                       "module": ["IFT", "IFT", "IFT", "IFT", "BBS", "BBS"]})
    bt = pd.DataFrame({"parent": ["u", "u"]}, index=["b1", "b2"])
    r = bs.module_loss_statistic(e, mh, {"IFT": 4, "BBS": 2, "axo": 10}, bt=bt)
    # IFT: 2/4 lost, BBS: 1/2 lost, axoneme absent -> skipped
    assert r.loc["b1", "f"] == pytest.approx(0.5) and r.loc["b1", "n_modules"] == 2
    assert r.loc["b2", "f"] == 0
    r = bs.module_loss_statistic(e, mh, {"IFT": 10, "BBS": 2}, bt=bt)   # IFT 4 < 0.5*10 -> skipped
    assert r.loc["b1", "f"] == pytest.approx(0.5) and r.loc["b1", "n_modules"] == 1
    assert bs.h12_signs({"a": "loss", "b": "no_loss"}) == {"a": 1.0, "b": -1.0}


def test_h9_values():
    v = bs.h9_values({"T1": 0.3}, {"T1": "M1"}, {"c1": 0.2, "c2": -0.1, "c3": 0.5})
    assert v["M1"] == 0.3
    assert set(np.abs(v[["c1", "c2", "c3"]]).round(6)) == {0.2, 0.1, 0.5}
    v2 = bs.h9_values({"T1": 0.3}, {"T1": "M1"}, {"c3": 0.5, "c1": 0.2, "c2": -0.1})
    assert v.sort_index().equals(v2.sort_index())                       # order-independent
    t, traits, tr = _hand()
    bt = bs.build_branch_table(t, traits, tr)
    assert set(bs.h9_control_branches(t, bt)) == {"A1", "B2"}


def test_effect_table_difference_scale():
    v = pd.Series({"S": 2.0, "a": 0.0, "b": 1.0, "c": -1.0, "d": 0.0})
    e = bs.effect_table(v, {"T": ["a", "b", "c", "d"]}, stems={"T": "S"}, scale="difference")
    assert e.loc["T", "est"] == pytest.approx(2.0)
    assert e.loc["T", "var"] == pytest.approx(np.var([0, 1, -1, 0], ddof=1) * 1.25)
    assert e.loc["T", "percentile"] == pytest.approx(4.5 / 5)
    e = bs.effect_table(v, {"T": ["a", "b", "c", "d"]}, stems={"T": "S"}, scale="difference",
                        signs={"T": -1})
    assert e.loc["T", "est"] == pytest.approx(-2.0)


# -----------------------------------------------------------------------------
# registry
# -----------------------------------------------------------------------------

def test_append_registry(tmp_path):
    f = tmp_path / "results" / "analysis_registry.tsv"
    inp = tmp_path / "x.tsv"
    inp.write_text("a\tb\n")
    r = bs.append_registry(str(f), dict(hypothesis="H1/H2", statistic="T_mean_u", n=7, k=20,
                                        estimate=0.61, ci_low=0.1, ci_high=0.9, p=0.0123,
                                        p_mcse=0.001, holm_p=0.05, tier="primary",
                                        input_md5s=[str(inp)], seed=11, n_draws=10000))
    bs.append_registry(str(f), dict(hypothesis="H8", statistic="x", tier="exploratory", p=float("nan")))
    df = pd.read_csv(f, sep="\t")
    assert list(df.columns) == bs.REGISTRY_COLUMNS and len(df) == 2
    assert df.loc[0, "input_md5s"].startswith("x.tsv:") and len(df.loc[0, "input_md5s"]) == 6 + 32
    assert df.loc[0, "git_commit"] and r["timestamp"]
    assert np.isnan(df.loc[1, "p"])
    with pytest.raises(ValueError):
        bs.append_registry(str(f), dict(hypothesis="H1", tier="primary", bogus=1))
    with pytest.raises(ValueError):
        bs.append_registry(str(f), dict(hypothesis="H1", tier="maybe"))


def test_toy_permute_exact_small():
    """'permute' p equals the full-enumeration permutation p on a tiny case."""
    v = pd.Series({"S": 3.0, "a": 0.0, "b": 1.0, "c": 2.0})
    r = bs.relabel_test(v, {"T": ["a", "b", "c"]}, stems={"T": "S"}, alternative="greater",
                        n_draws=20_000, scheme="permute")
    assert r.T_obs == pytest.approx(3.5 / 4)
    # stem value uniform over 4 positions; T >= T_obs only when the max sits at S
    assert r.p == pytest.approx(0.25, abs=0.015)


# -----------------------------------------------------------------------------
# path stems (DESIGN v3.1, D7)
# -----------------------------------------------------------------------------
PATH_TREE = "(((((a1,a2)M,(b1,b2)B)X,(c1,c2)C)P,(h1,h2)H)Q,((d1,d2)D,((e1,e2)E,(f1,f2)F)EF)G)R;"


def _path_setup():
    t = bs.load_tree(PATH_TREE)
    tips = t.get_leaf_names()
    traits = pd.DataFrame({"code": tips, "mc_state": ["simple" if x.startswith("a") else "unicellular"
                                                      for x in tips], "data_type": "genome"})
    tr = pd.DataFrame([dict(transition_id="T", from_state="unicellular", to_state="simple",
                            analysis_set="primary", primary_eligible=1, M="M", P="P",
                            stem_path="P>X>M")])
    bl = {n.name: 1.0 for n in t.traverse() if n.up is not None}
    bl.update(X=2.0, M=3.0)
    return t, traits, tr, bl


def test_path_stem_branch_table():
    t, traits, tr, bl = _path_setup()
    bt = bs.build_branch_table(t, traits, tr, branch_lengths=bl)
    assert bt.loc["M", "stem_path"] == "X;M" and bs.path_edges(bt, "M") == ["X", "M"]
    assert bt.loc["M", "branch_length"] == 5.0 and bt.loc["M", "edge_branch_length"] == 3.0
    assert bt.loc["M", "log_bl"] == pytest.approx(math.log(5))
    assert bt.loc["X", "branch_class"] == "stem" and bt.loc["X", "role"] == "path_edge"
    assert bt.loc["X", "transition_id"] == "T" and not bt.loc["X", "primary"]
    assert dict(bs.stem_branches(bt)) == {"T": "M"}
    ex = bt["U_exclusion"]
    assert ex["P"] == "adjacent_to_stem"         # edge into P
    assert ex["C"] == "adjacent_to_stem"         # P -> sister
    assert ex["B"] == "adjacent_to_stem"         # intermediate node -> side branch
    assert ex["X"] == "class:stem" and ex["M"] == "class:stem"
    assert set(bt.index[bt["in_U"]]) == {"H", "D", "EF", "E", "F"}
    # covariates of the path = those of the last edge (except length)
    assert bt.loc["M", "n_desc"] == 2 and bt.loc["M", "n_sister"] == 2


def _path_events():
    rows = [("P", "X", "gain", True, "X|g1"), ("P", "X", "gain", True, "X|g2"),
            ("X", "M", "gain", True, "M|g3"), ("X", "M", "gain", True, "M|g4"), ("X", "M", "gain", True, "M|g5"),
            ("P", "X", "dup", False, "P|h1"), ("P", "X", "dup", False, "P|h2"),
            ("X", "M", "dup", False, "X|h1"), ("X", "M", "dup", False, "X|h3"), ("X", "M", "dup", False, "X|new"),
            ("P", "X", "loss", False, "P|h4"), ("X", "M", "loss", False, "X|h3")]
    ev = pd.DataFrame(rows, columns=["branch_parent", "branch_child", "event", "at_origin", "hog_id"])
    hogs = pd.DataFrame({"hog_id": [f"P|h{i}" for i in range(10)] + ["X|h1", "X|h3", "X|new"],
                         "level": ["P"] * 10 + ["X"] * 3,
                         "parent_hog_id": [None] * 10 + ["P|h1", "P|h3", None]})
    return ev, hogs


def test_path_stem_statistics():
    t, traits, tr, bl = _path_setup()
    bt = bs.build_branch_table(t, traits, tr, branch_lengths=bl)
    ev, hogs = _path_events()
    g = bs.count_statistic(ev, "gain", branches=bt)
    assert g["M"] == 5 and g["X"] == 2
    comp = bs.compositional_statistic(ev, hogs, {"P|h1", "P|h5"}, bt=bt, min_d=1)
    # union lifted to P: {P|h1, P|h2, P|h3, X|new}; a = 1 (P|h1); pi_u at P = 2/10
    assert comp.loc["M", "parent"] == "P" and comp.loc["M", "d"] == 4 and comp.loc["M", "a"] == 1
    assert comp.loc["M", "s"] == pytest.approx(math.log((1.5 / 5) / (3.5 / 5)) - math.log(0.2 / 0.8))
    mh = pd.DataFrame({"level": "P", "hog_id": ["P|h3", "P|h4", "P|h6", "P|h7"], "module": "IFT"})
    ml = bs.module_loss_statistic(ev, mh, {"IFT": 4}, bt=bt, hogs=hogs)
    assert ml.loc["M", "f"] == pytest.approx(0.5)       # P|h4 (edge X) + X|h3 -> P|h3 (edge M)
    assert ml.loc["X", "f"] == pytest.approx(0.25)
    with pytest.raises(ValueError):
        bs.module_loss_statistic(ev, mh, {"IFT": 4}, bt=bt)


def test_bad_path_falls_back():
    t, traits, tr, bl = _path_setup()
    tr.loc[0, "stem_path"] = "C>X>M"
    with pytest.warns(UserWarning, match="stem_path"):
        bt = bs.build_branch_table(t, traits, tr, branch_lengths=bl)
    assert bt.loc["M", "stem_path"] == "M"


# -----------------------------------------------------------------------------
# H1/H2 residual statistic (DESIGN v3.1, D6)
# -----------------------------------------------------------------------------

def test_fit_count_glm_matches_statsmodels():
    import statsmodels.api as sm
    rng = np.random.default_rng(1)
    X = np.column_stack([np.ones(40), rng.normal(size=(40, 3))])
    mu = np.exp(2.5 + X[:, 1:] @ [0.5, 0.2, -0.3])
    y = rng.negative_binomial(1 / 0.3, (1 / 0.3) / (1 / 0.3 + mu)).astype(float)
    f = bs.fit_count_glm(X, y)
    ref = sm.NegativeBinomial(y, X).fit(disp=0, maxiter=500)
    assert f["method"] == "nb"
    assert np.allclose(f["beta"], ref.params[:-1], atol=1e-4)
    assert f["alpha"] == pytest.approx(ref.params[-1], rel=1e-3)
    p = bs.fit_count_glm(X, y, family="poisson")
    refp = sm.GLM(y, X, family=sm.families.Poisson()).fit()
    assert np.allclose(p["beta"], refp.params, atol=1e-6) and p["method"] == "poisson"


def test_fit_count_glm_fallback(monkeypatch):
    rng = np.random.default_rng(2)
    X = np.column_stack([np.ones(30), rng.normal(size=30)])
    y = rng.poisson(5, 30).astype(float)

    def boom(*a, **k):
        raise FloatingPointError("forced")
    monkeypatch.setattr(bs, "_nb_ll", boom)
    with pytest.warns(UserWarning, match="Poisson with quasi-dispersion"):
        f = bs.fit_count_glm(X, y)
    assert f["method"] == "poisson_quasi" and f["fallback_reason"].startswith("nb_error")
    assert f["phi"] > 0


def _shift_bt(rng, n_pool=30, n_stems=7, shift=1.0, family_mu=True):
    """Stems shifted by ``shift`` pool SDs towards longer branches and larger
    clades (deep stems vs shallow controls); counts depend on both."""
    pool = [f"u{i}" for i in range(n_pool)]
    stems = {f"T{i}": f"s{i}" for i in range(n_stems)}
    idx = pool + list(stems.values())
    isst = np.r_[np.zeros(n_pool, bool), np.ones(n_stems, bool)]
    bt = pd.DataFrame(index=idx)
    bt["log_bl"] = rng.normal(0, 0.4, len(idx)) + isst * shift * 0.4
    bt["log_n_desc"] = rng.normal(1.5, 0.5, len(idx)) + isst * shift * 0.5
    bt["log_n_sister"] = rng.normal(2.0, 0.5, len(idx))
    bt["frac_tx_desc"] = rng.uniform(0, 0.3, len(idx))
    bt["in_U"] = ~isst
    mu = np.exp(1.5 + 0.8 * bt["log_bl"] + 0.5 * bt["log_n_desc"])
    m = bs.match_controls(bt, stems, covariates=list(bs.RESID_COVARIATES), k=20)
    return bt, stems, m, mu


def _nb(rng, mu, a=0.3):
    return rng.negative_binomial(1 / a, (1 / a) / (1 / a + np.asarray(mu)))


def _shift_pvalues(effect, reps, seed, shift=1.0, raw=False, boot=0, family="nb", disp=0.3):
    rng = np.random.default_rng(seed)
    ps = []
    for _ in range(reps):
        bt, stems, m, mu = _shift_bt(rng, shift=shift)
        mu = mu.copy()
        mu[list(stems.values())] *= effect
        y = pd.Series(_nb(rng, mu, disp) if family == "nb" else rng.poisson(mu), index=bt.index, dtype=float)
        if boot:
            ps.append(bs.residual_bootstrap_test(y, bt, m, n_boot=boot, family=family,
                                                 seed=int(rng.integers(2**31)))["p"])
            continue
        s = y if raw else bs.count_residual_statistic(y, bt, family=family)["s"]
        ps.append(bs.relabel_test(s, m, alternative="two-sided", n_draws=500,
                                  seed=int(rng.integers(2**31))).p)
    return np.array(ps)


def test_residual_null_calibrated_under_covariate_shift():
    """Stems 1 pool-SD deeper (longer branch, larger clade; real tree: SMD ~1.1)
    and counts depend on both: raw counts reject far too often, residuals with
    the permutation scheme stay calibrated."""
    ps = _shift_pvalues(1.0, reps=150, seed=21, shift=1.0)
    se = math.sqrt(0.05 * 0.95 / len(ps))
    assert np.mean(ps <= 0.05) <= 0.05 + 2.5 * se
    assert stats.kstest(ps, "uniform").pvalue > 0.001
    raw = _shift_pvalues(1.0, reps=60, seed=22, shift=1.0, raw=True)
    assert np.mean(raw <= 0.05) > 0.3


def test_residual_bootstrap_calibrated_under_strong_shift():
    """Under strong extrapolation (2.5 SD) the stems share the GLM's prediction
    error, so permuting residuals is anti-conservative; the parametric bootstrap
    (refit + LOO in every replicate) is not. Poisson family for speed."""
    boot = _shift_pvalues(1.0, reps=40, seed=31, shift=2.5, boot=19, family="poisson")
    perm = _shift_pvalues(1.0, reps=40, seed=31, shift=2.5, family="poisson")
    se = math.sqrt(0.1 * 0.9 / 40)
    assert np.mean(boot <= 0.1) <= 0.1 + 2.5 * se
    assert np.mean(perm <= 0.1) > np.mean(boot <= 0.1)


def test_residual_detects_planted_effect_under_covariate_shift():
    """2-fold planted effect on 7 deep stems (shift 1 SD), NB dispersion 0.1.
    (At dispersion 0.3 the MDE on the real tree is ~2.4-fold: results/prereg.)"""
    ps = _shift_pvalues(2.0, reps=40, seed=23, shift=1.0, disp=0.1)
    assert np.mean(ps <= 0.05) >= 0.8


def test_residual_statistic_details():
    rng = np.random.default_rng(4)
    bt, stems, m, mu = _shift_bt(rng)
    bt["log_bl"] = 0.0                                  # branch-length placeholder 1.0
    y = pd.Series(_nb(rng, mu), index=bt.index, dtype=float)
    r = bs.count_residual_statistic(y, bt)
    assert r.attrs["dropped_covariates"] == {"log_bl": "constant_in_pool"}
    assert (r.loc[bt.index[bt.in_U], "fit"] == "loo").all() and (r.loc["s0", "fit"] == "full")
    # LOO prediction of one control equals a fit without it
    keep = [b for b in bt.index[bt.in_U] if b != "u0"]
    cols = r.attrs["covariates"]
    Z = (bt[cols] - bt.loc[bt.in_U, cols].mean()) / bt.loc[bt.in_U, cols].std(ddof=1)
    D = np.column_stack([np.ones(len(Z)), Z.to_numpy()])
    f = bs.fit_count_glm(D[[bt.index.get_loc(k) for k in keep]], y[keep].to_numpy())
    assert r.loc["u0", "yhat"] == pytest.approx(math.exp(D[bt.index.get_loc("u0")] @ f["beta"]))
    assert r.loc["u0", "s"] == pytest.approx(math.log((y["u0"] + 0.5) / (r.loc["u0", "yhat"] + 0.5)))


def test_branch_lengths_from_tree():
    t = bs.load_tree(HAND)
    import ete3
    bl = ete3.Tree(HAND.replace("A1a", "A1a:0.7").replace(")A1,", ")A1:0.3,"), format=1)
    for n in bl.traverse():
        if n.name not in ("A1a", "A1"):
            n.dist = 0.1
    d = bs.branch_lengths_from_tree(t, bl)
    assert d["A1a"] == pytest.approx(0.7) and d["A1"] == pytest.approx(0.3) and d["t1"] == pytest.approx(0.1)
    with pytest.warns(UserWarning, match="placeholder"):
        d = bs.branch_lengths_from_tree(t, "/nonexistent/species_tree_bl_euk.nwk")
    assert set(d.values()) == {1.0}
    bt = bs.branch_table(t, None, branch_lengths=bl)
    assert bt.loc["A1a", "branch_length"] == pytest.approx(0.7)
