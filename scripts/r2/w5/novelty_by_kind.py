#!/usr/bin/env python
"""Secondary (exploratory, post hoc): structurally-new fraction by event kind (origin / origin_complex / nested)
vs controls, pooled over panels (Mantel-Haenszel over panels that have that kind; family-level logistic GLM with
log outside-target count and panel, SE clustered by branch). Input: results/r2/w5/novelty_summary/per_branch.tsv.
Output: results/r2/w5/novelty_summary/by_kind.tsv"""
import os, sys
import numpy as np, pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import StratifiedTable
D = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/w5/novelty_summary"
out = os.path.join(D, "by_kind.tsv")
if os.path.exists(out):
    sys.exit("refusing to overwrite")
B = pd.read_csv(os.path.join(D, "per_branch.tsv"), sep="\t")
B = B[(B["kind"] == "novel") & B["role"].isin(["event", "control"])].copy()
B["ekind"] = np.where(B["role"] == "event", B["event"].fillna("").str.split(":").str[1], "control")
B["ekind"] = B["ekind"].replace({"origin_complex": "origin"})   # animals: uni->complex, grouped with origins
rows = []
for k in ("origin", "nested"):
    tabs = []
    for R, s in B.groupby("panel"):
        e, c = s[s["ekind"] == k], s[s["ekind"] == "control"]
        if len(e) and len(c):
            tabs.append(np.array([[e.n_new.sum(), e.n.sum() - e.n_new.sum()], [c.n_new.sum(), c.n.sum() - c.n_new.sum()]], float))
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint()
    sub = B[B["ekind"].isin([k, "control"]) & B["panel"].isin(B.loc[B["ekind"] == k, "panel"])]
    # expand to families
    fam = sub.loc[sub.index.repeat(sub["n"])].copy()
    fam["new"] = np.concatenate([np.r_[np.ones(a), np.zeros(n - a)] for a, n in zip(sub["n_new"], sub["n"])])
    fam["is_event"] = (fam["ekind"] == k).astype(int)
    fam["log_targets"] = np.log(fam["n_outside_targets"])
    fam["cl"] = pd.factorize(fam["panel"] + ":" + fam["branch_child"])[0]
    m = smf.glm("new ~ is_event + log_targets + C(panel)", fam, family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": fam["cl"]})
    ci = np.exp(m.conf_int().loc["is_event"])
    e = B[B["ekind"] == k]
    rows.append(dict(event_kind=k, n_events=len(e), n_panels=len(tabs), events_new=int(e.n_new.sum()), events_n=int(e.n.sum()),
                     events_frac_new=round(e.n_new.sum() / e.n.sum(), 4), events=";".join(e["event"].str.split(":").str[0]),
                     MH_OR=round(st.oddsratio_pooled, 3), MH_lo=round(lo, 3), MH_hi=round(hi, 3), MH_p=st.test_null_odds().pvalue,
                     adj_OR=round(float(np.exp(m.params["is_event"])), 3), adj_lo=round(float(ci[0]), 3), adj_hi=round(float(ci[1]), 3),
                     adj_p=float(m.pvalues["is_event"])))
c = B[B["ekind"] == "control"]
rows.append(dict(event_kind="control", n_events=len(c), events_new=int(c.n_new.sum()), events_n=int(c.n.sum()),
                 events_frac_new=round(c.n_new.sum() / c.n.sum(), 4)))
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
print(pd.DataFrame(rows).drop(columns=["events"]).to_string())
