"""Exploratory conditional null for recurrence (DEVIATIONS D-W4-2): per event keep the observed number of significant
families per direction, assign them at random among the families tested in that event; events independent. B=1000, seed 11."""
import os
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import OUT, out, SEED  # noqa: E402

B = 1000
obs = pd.read_csv(f"{OUT}/scan/recurrence_observed.tsv", sep="\t")
ev = pd.read_csv(f"{OUT}/scan/events_used.tsv", sep="\t")
prim = ev.event[ev.tier == "primary"].tolist()
rng = np.random.default_rng(SEED)
F = len(obs)
rows, famp = [], pd.DataFrame(dict(pfam=obs.pfam))
for d, sgn in (("expanded", 1), ("contracted", -1)):
    nullrec = np.zeros((B, F), np.int8)
    for e in prim:
        tested = np.where(obs[f"tested_{e}"].to_numpy() == 1)[0]
        k = int((obs[f"sig_{e}"] == sgn).sum())
        for b in range(B):
            nullrec[b, rng.choice(tested, k, replace=False)] += 1
    o = obs[f"rec_primary_{d}"].to_numpy()
    for kk in range(1, 8):
        nl = (nullrec >= kk).sum(1); ob = int((o >= kk).sum())
        rows.append(dict(direction=d, k=kk, observed_n_families=ob, null_mean=nl.mean(), null_q025=np.quantile(nl, .025),
                         null_q975=np.quantile(nl, .975), p_perm=((nl >= ob).sum() + 1) / (B + 1)))
    famp[f"cond_p_{d}"] = ((nullrec >= o[None, :]).sum(0) + 1) / (B + 1)
r = pd.DataFrame(rows)
r.to_csv(out("scan/recurrence_conditional_null_summary.tsv"), sep="\t", index=False)
famp.to_csv(out("scan/recurrence_conditional_family_p.tsv.gz"), sep="\t", index=False, compression="gzip")
print(r.to_string())
