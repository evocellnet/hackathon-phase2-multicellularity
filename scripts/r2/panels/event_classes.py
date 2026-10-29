#!/usr/bin/env python
"""Descriptive comparators (exploratory, post hoc): pooled rate ratios (modified HK) by event class and branch type
from a pooled tag's event_residuals.tsv -> pooled/<tag>/event_class_comparators.tsv (refuses to overwrite)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, guard  # noqa: E402
from pooled import meta_mhk, rr  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

tag = sys.argv[1]
d = os.path.join(OUT, "pooled", tag)
dst = os.path.join(d, "event_class_comparators.tsv")
guard(dst)
E = pd.read_csv(os.path.join(d, "event_residuals.tsv"), sep="\t")
E["cls"] = np.where((E["kind"] == "nested") & (E["from_state"] == "simple") & (E["to_state"] == "complex"), "simple->complex",
                    np.where(E["kind"].isin(["origin", "origin_complex"]), "origin",
                             np.where(E["kind"] == "reversal", "reversal", "other")))
rows = []
for (c, tip), g in E[E["cls"] != "other"].groupby(["cls", "is_tip"]):
    m = rr(meta_mhk(g["log_rr"], g["var"]))
    rows.append(dict(event_class=c, branch="tip" if tip else "internal", model="primary (internal controls)", k=m["k"],
                     RR=m["RR"], RR_low=m["RR_low"], RR_high=m["RR_high"], tau2=m["tau2"], events=";".join(g["event_id"])))
    if tip and "tip_log_rr" in g.columns:
        h = g[g["tip_log_rr"].notna()]
        if len(h):
            m = rr(meta_mhk(h["tip_log_rr"], h["tip_var"]))
            rows.append(dict(event_class=c, branch="tip", model="tip controls (sensitivity)", k=m["k"], RR=m["RR"],
                             RR_low=m["RR_low"], RR_high=m["RR_high"], tau2=m["tau2"], events=";".join(h["event_id"])))
pd.DataFrame(rows).to_csv(dst, sep="\t", index=False)
print(pd.DataFrame(rows).drop(columns="events").round(3).to_string())
