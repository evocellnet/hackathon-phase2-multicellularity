#!/usr/bin/env python
"""W5 task B statistics: structurally-new fraction of sampled sequence-novel gained families on event branches
vs in-run control branches, per panel and pooled; positive (sequence-old) and negative (shuffled) controls.

Event branches: M branches of events in results/r2/panels/<R>/event_table.tsv with kind in
{origin, origin_complex, nested}; internal M (M_is_tip == 0) are the primary comparison, tip M are reported
descriptively (tip gains are species-specific and not comparable with internal control branches).
Controls: branch_table.tsv control == True (panel agent's definition, DESIGN §6 / METHODS_DEFINITIONS item 4).
Tests (exploratory): family-level 2x2 (structurally new vs old; event vs control) Fisher exact with conditional
OR and 95 % CI per panel; pooled Mantel-Haenszel OR (statsmodels StratifiedTable) with 95 % CI;
branch-level permutation: statistic = mean(frac_new | event) - mean(frac_new | control), event labels permuted
among event + control branches within panel (exact enumeration per panel; pooled: 10,000 stratified Monte Carlo
permutations, seed 11), one- and two-sided p.
Output: results/r2/w5/novelty_summary/{per_branch.tsv, per_panel.tsv, pooled.tsv, controls.tsv}
Usage: novelty_stats.py R2 R3 ...   (panels whose results/r2/w5/novelty/<R>/DONE.json exists)
"""
import itertools, os, sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
NOV = os.path.join(ROOT, "results/r2/w5/novelty")
OUT = os.path.join(ROOT, "results/r2/w5/novelty_summary")
EVENT_KINDS = ("origin", "origin_complex", "nested")


def labels(R):
    bt = pd.read_csv(os.path.join(ROOT, "results/r2/panels", R, "branch_table.tsv"), sep="\t")
    et = pd.read_csv(os.path.join(ROOT, "results/r2/panels", R, "event_table.tsv"), sep="\t")
    et = et[et["kind"].isin(EVENT_KINDS)]
    ev = dict(zip(et["M"], et["event_id"] + ":" + et["kind"]))
    tip = dict(zip(et["M"], et["M_is_tip"].astype(int)))
    lab = {}
    for r in bt.itertuples():
        b = r.branch_child
        if b in ev:
            lab[b] = ("event_tip" if tip[b] else "event", ev[b])
        elif bool(r.control):
            lab[b] = ("control", "")
    return lab


def n_outside(R):
    """number of panel proteins (with 3Di) outside each branch's clade = size of the search target set"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import struct_novelty as sn
    lk = pd.read_csv(os.path.join(ROOT, "results/r2/w5/work/novelty", R, "db/panel.lookup"), sep="\t", header=None,
                     usecols=[1], names=["id"])
    per_sp = lk["id"].str.split("_").str[0].value_counts()
    tot = int(per_sp.sum())
    return {b: tot - int(per_sp.reindex(list(c)).fillna(0).sum()) for b, c in sn.clades(R).items()}


def perm_exact(fe, fc):
    x = np.r_[fe, fc]
    k = len(fe)
    obs = np.mean(fe) - np.mean(fc)
    ds = []
    for idx in itertools.combinations(range(len(x)), k):
        m = np.zeros(len(x), bool)
        m[list(idx)] = True
        ds.append(x[m].mean() - x[~m].mean())
    ds = np.array(ds)
    return obs, float((ds >= obs - 1e-12).mean()), float((np.abs(ds) >= abs(obs) - 1e-12).mean()), len(ds)


def main(panels):
    os.makedirs(OUT, exist_ok=True)
    for f in ("per_branch.tsv", "per_panel.tsv", "pooled.tsv", "controls.tsv"):
        if os.path.exists(os.path.join(OUT, f)):
            sys.exit(f"refusing to overwrite {OUT}/{f}")
    F = []
    for R in panels:
        if not os.path.exists(os.path.join(NOV, R, "DONE.json")):
            sys.exit(f"{R}: novelty not done")
        f = pd.read_csv(os.path.join(NOV, R, "families.tsv"), sep="\t")
        lab = labels(R)
        f["role"] = f["branch_child"].map(lambda b: lab.get(b, ("other", ""))[0])
        f["event"] = f["branch_child"].map(lambda b: lab.get(b, ("other", ""))[1])
        f["n_outside_targets"] = f["branch_child"].map(n_outside(R))
        F.append(f)
    F = pd.concat(F, ignore_index=True)
    F["new"] = (F["struct_class"] == "structurally_new").astype(int)
    B = (F.groupby(["panel", "branch_child", "role", "event", "kind"])
         .agg(n=("new", "size"), n_new=("new", "sum"), n_outside_targets=("n_outside_targets", "first")).reset_index())
    B["frac_new"] = (B["n_new"] / B["n"]).round(4)
    B.to_csv(os.path.join(OUT, "per_branch.tsv"), sep="\t", index=False)
    # controls of the search
    CT = F.groupby(["panel", "kind"]).agg(n_families=("new", "size"), n_with_struct_hit=("struct_old", "sum")).reset_index()
    CT["frac_with_hit"] = (CT["n_with_struct_hit"] / CT["n_families"]).round(4)
    CT.to_csv(os.path.join(OUT, "controls.tsv"), sep="\t", index=False)
    nov = B[B["kind"] == "novel"]
    rows, tabs, strata = [], [], []
    for R, s in nov.groupby("panel"):
        e, c = s[s["role"] == "event"], s[s["role"] == "control"]
        r = dict(panel=R, n_event_branches=len(e), n_control_branches=len(c),
                 event_new=int(e["n_new"].sum()), event_n=int(e["n"].sum()),
                 control_new=int(c["n_new"].sum()), control_n=int(c["n"].sum()))
        r["event_frac_new"] = round(r["event_new"] / r["event_n"], 4) if r["event_n"] else np.nan
        r["control_frac_new"] = round(r["control_new"] / r["control_n"], 4) if r["control_n"] else np.nan
        tip = s[s["role"] == "event_tip"]
        r["tip_events_frac_new"] = round(tip["n_new"].sum() / tip["n"].sum(), 4) if len(tip) else np.nan
        if r["event_n"] and r["control_n"]:
            t = np.array([[r["event_new"], r["event_n"] - r["event_new"]],
                          [r["control_new"], r["control_n"] - r["control_new"]]])
            o = odds_ratio(t)
            ci = o.confidence_interval(0.95)
            r.update(OR=round(o.statistic, 3), OR_lo=round(ci.low, 3), OR_hi=round(ci.high, 3),
                     fisher_p_two_sided=fisher_exact(t)[1], fisher_p_greater=fisher_exact(t, "greater")[1])
            tabs.append(t)
            ob, p1, p2, nperm = perm_exact(e["frac_new"].to_numpy(), c["frac_new"].to_numpy())
            r.update(branch_mean_diff=round(ob, 4), perm_p_greater=round(p1, 4), perm_p_two_sided=round(p2, 4),
                     n_permutations=nperm)
            strata.append((e["frac_new"].to_numpy(), c["frac_new"].to_numpy()))
        rows.append(r)
    P = pd.DataFrame(rows)
    P.to_csv(os.path.join(OUT, "per_panel.tsv"), sep="\t", index=False)
    pooled = []
    if tabs:
        from statsmodels.stats.contingency_tables import StratifiedTable
        st = StratifiedTable([t.astype(float) + (0.5 if (t == 0).any() else 0) for t in tabs])
        lo, hi = st.oddsratio_pooled_confint()
        tot = np.sum(tabs, axis=0)
        # stratified permutation of branch labels
        rng = np.random.default_rng(11)
        def stat(sts):
            fe = np.concatenate([a for a, _ in sts]); fc = np.concatenate([b for _, b in sts])
            return fe.mean() - fc.mean()
        obs = stat(strata)
        null = []
        for _ in range(10000):
            ps = []
            for a, b in strata:
                x = rng.permutation(np.r_[a, b])
                ps.append((x[:len(a)], x[len(a):]))
            null.append(stat(ps))
        null = np.array(null)
        pooled.append(dict(n_panels=len(tabs), event_new=int(tot[0, 0]), event_n=int(tot[0].sum()),
                           control_new=int(tot[1, 0]), control_n=int(tot[1].sum()),
                           event_frac_new=round(tot[0, 0] / tot[0].sum(), 4),
                           control_frac_new=round(tot[1, 0] / tot[1].sum(), 4),
                           MH_OR=round(st.oddsratio_pooled, 3), MH_OR_lo=round(lo, 3), MH_OR_hi=round(hi, 3),
                           MH_test_p=st.test_null_odds().pvalue, breslow_day_homogeneity_p=st.test_equal_odds().pvalue,
                           branch_mean_diff=round(obs, 4), perm_p_greater=round(float((null >= obs - 1e-12).mean()), 4),
                           perm_p_two_sided=round(float((np.abs(null) >= abs(obs) - 1e-12).mean()), 4)))
    # adjustment for search-space size: family-level logistic GLM new ~ event + log(outside targets) + panel,
    # cluster-robust SEs by branch (events + controls, sequence-novel sample)
    import statsmodels.formula.api as smf
    import statsmodels.api as sm_api
    d = F[(F["kind"] == "novel") & F["role"].isin(["event", "control"])].copy()
    d["is_event"] = (d["role"] == "event").astype(int)
    d["log_targets"] = np.log(d["n_outside_targets"])
    d["cl"] = pd.factorize(d["panel"] + ":" + d["branch_child"])[0]
    form = "new ~ is_event + log_targets" + (" + C(panel)" if d["panel"].nunique() > 1 else "")
    m = smf.glm(form, d, family=sm_api.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": d["cl"]})
    ci = m.conf_int().loc["is_event"]
    adj = dict(model=form + " (logit; SE clustered by branch)", n_families=len(d), n_branches=int(d["cl"].nunique()),
               OR_event_adj=round(float(np.exp(m.params["is_event"])), 3), OR_lo=round(float(np.exp(ci[0])), 3),
               OR_hi=round(float(np.exp(ci[1])), 3), p=float(m.pvalues["is_event"]),
               OR_per_log_target=round(float(np.exp(m.params["log_targets"])), 3),
               p_log_target=float(m.pvalues["log_targets"]))
    if pooled:
        pooled[0].update({f"adj_{k}": v for k, v in adj.items()})
    else:
        pooled.append({f"adj_{k}": v for k, v in adj.items()})
    pd.DataFrame(pooled).to_csv(os.path.join(OUT, "pooled.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(CT.to_string()); print(P.to_string()); print(pd.DataFrame(pooled).T.to_string())


if __name__ == "__main__":
    main(sys.argv[1:])
