#!/usr/bin/env python3
"""Approximate H9 pooled-effect CI that accounts for the 7 stems sharing the same 10 controls.

POST HOC, reporting-only check (requested by the external report review, item 6). It does not change the
pre-registered test or its p-value. The pre-registered pooled interval (results/h9/h9_tests.tsv) gives each
stem the variance var(controls)*(1 + 1/k) and averages over stems, which treats the control-mean uncertainty
(var/k) as independent across stems although all stems use the same 10 controls. Here:
    var(pooled) = var(controls)/n_stems + var(controls)/k
with the fixed-effect mean of (r_i - mean control r) as the point estimate (tau^2 is ~0 in h9_tests.tsv).
Output: results/h9/pooled_ci_shared_controls.tsv
"""
import os, math
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
tr = pd.read_csv(os.path.join(ROOT, "results/h9/h9_transitions.tsv"), sep="\t")
ct = pd.read_csv(os.path.join(ROOT, "results/h9/h9_controls.tsv"), sep="\t")
tests = pd.read_csv(os.path.join(ROOT, "results/h9/h9_tests.tsv"), sep="\t")
st = tr[tr["set"] == "primary"]
rc = ct["r_oriented"].astype(float)
n, k = len(st), len(rc)
d = st["r"].astype(float) - rc.mean()
v = rc.var(ddof=1)
est = d.mean()
se = math.sqrt(v / n + v / k)
pre = tests[tests["analysis"] == "primary"].iloc[0]
out = pd.DataFrame([dict(
    n_stems=n, k_controls=k, var_controls=v, est_mean_diff=est, se_shared=se,
    ci_low=est - 1.96 * se, ci_high=est + 1.96 * se,
    rate_ratio=math.exp(est), rr_ci_low=math.exp(est - 1.96 * se), rr_ci_high=math.exp(est + 1.96 * se),
    prereg_pooled_diff=pre["pooled_diff"], prereg_rr_ci_low=math.exp(pre["pooled_ci_low"]),
    prereg_rr_ci_high=math.exp(pre["pooled_ci_high"]),
    note="post hoc reporting check; shared-control variance counted once; not part of the pre-registered analysis")])
p = os.path.join(ROOT, "results/h9/pooled_ci_shared_controls.tsv")
out.to_csv(p, sep="\t", index=False, float_format="%.4f")
print(out.T.to_string())
