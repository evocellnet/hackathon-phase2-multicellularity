"""D17 family ladder, post hoc sensitivity analyses (added after seeing 02's output; labelled as such in summary.txt).
 (a) leave-one-origin-out ladders (6 origins);
 (b) origins whose sister clade is genome-based only (drops animals: 5/8 sister proteomes from transcriptomes, and
     dictyostelids: both sisters transcriptomes) -> 5 origins;
 (c) detectability-stratified null: within each event, hits are shuffled only among tested families of the same decile of
     total copy number in that event (preserves the number of hits per decile), because large families have more power
     in every event and could recur for that reason alone.
Same scan inputs, BH q < 0.1 within event; B = 10,000; seed 11. Output: family_ladder_sensitivity.tsv
"""
import os
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import R, R1_ORIGINS, LABEL, out  # noqa: E402

B, SEED, QCUT = 10000, 11, 0.1
SCAN = f"{R}/r2/w4/scan"
D = {}
for e in R1_ORIGINS:
    d = pd.read_csv(f"{SCAN}/event_{e}.tsv", sep="\t")
    d["s"] = np.where(d.q_bh < QCUT, np.sign(d.direction), 0).astype(int)
    d["size"] = d.mean_count_mc * d.n_mc + d.mean_count_sister * d.n_sister
    d["dec"] = pd.qcut(d["size"].rank(method="first"), 10, labels=False)
    D[e] = d.set_index("pfam")
allf = sorted(set().union(*[set(d.index) for d in D.values()]))
fi = {f: i for i, f in enumerate(allf)}
F = len(allf)
rng = np.random.default_rng(SEED)


def null_hits(e, sgn, stratified):
    d = D[e]
    M = np.zeros((B, F), bool)
    groups = d.groupby("dec") if stratified else [(0, d)]
    for _, g in groups:
        idx = np.array([fi[f] for f in g.index])
        k = int((g.s == sgn).sum())
        if k == 0:
            continue
        for b0 in range(0, B, 1000):
            nb = min(1000, B - b0)
            pick = np.argpartition(rng.random((nb, len(idx))), k - 1, axis=1)[:, :k]
            M[np.repeat(np.arange(b0, b0 + nb), k), idx[pick.ravel()]] = True
    return M


rows = []
for dname, sgn in (("expanded", 1), ("contracted", -1)):
    obs = {e: np.isin(np.arange(F), [fi[f] for f in D[e].index[D[e].s == sgn]]) for e in R1_ORIGINS}
    for strat in (False, True):
        NH = {e: null_hits(e, sgn, strat) for e in R1_ORIGINS}
        sets = [("all 7 origins", R1_ORIGINS)]
        if not strat:
            sets += [(f"leave out {LABEL[x]}", [e for e in R1_ORIGINS if e != x]) for x in R1_ORIGINS]
        sets += [("genome-based sisters only (no animals, no dictyostelids)",
                  [e for e in R1_ORIGINS if e not in ("T_metazoa", "T_dictyostelia")])]
        for lname, evs in sets:
            ob = sum(obs[e].astype(int) for e in evs)
            nl = sum(NH[e].astype(np.int8) for e in evs)
            for k in range(1, len(evs) + 1):
                o = int((ob >= k).sum())
                n_null = (nl >= k).sum(1)
                if o == 0 and n_null.mean() < 0.01:
                    continue
                rows.append(dict(analysis=lname, null="detectability-stratified (size deciles)" if strat else
                                 "uniform family shuffle", direction=dname, tier_k=k, n_units=len(evs),
                                 observed_n_families=o, null_mean=n_null.mean(), null_q025=np.quantile(n_null, .025),
                                 null_q975=np.quantile(n_null, .975), p_perm=((n_null >= o).sum() + 1) / (B + 1)))
        del NH
S = pd.DataFrame(rows)
S.to_csv(out("family_ladder_sensitivity.tsv"), sep="\t", index=False)
pd.set_option("display.width", 220)
pd.set_option("display.max_rows", 500)
print(S[S.tier_k.between(2, 4)].to_string(index=False))

