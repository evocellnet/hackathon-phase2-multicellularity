"""D17 post hoc sensitivity: frozen-set recurrence (>= k origins) under the detectability-stratified null
(hits shuffled within deciles of total copy number per event), for the 7-origin and the origin-collapsed ladders.
B = 10,000, seed 11. Output: set_families_ladder_stratified.tsv"""
import os
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import P, R, R1_ORIGINS, R1_NESTED, out  # noqa: E402

B, SEED, QCUT = 10000, 11, 0.1
EVENTS = R1_ORIGINS + list(R1_NESTED)
UNIT = {e: e for e in R1_ORIGINS}
UNIT.update(R1_NESTED)
D = {}
for e in EVENTS:
    d = pd.read_csv(f"{R}/r2/w4/scan/event_{e}.tsv", sep="\t")
    d["s"] = np.where(d.q_bh < QCUT, np.sign(d.direction), 0).astype(int)
    d["size"] = d.mean_count_mc * d.n_mc + d.mean_count_sister * d.n_sister
    d["dec"] = pd.qcut(d["size"].rank(method="first"), 10, labels=False)
    D[e] = d.set_index("pfam")
sets = {}
for fn, nm in (("elm_motif_binding.tsv", "motif-binding"), ("cell_death.tsv", "cell-death"),
               ("tf_celltype.tsv", "cell-type TF")):
    s = pd.read_csv(f"{P}/data/r2_sets/{fn}", sep="\t")
    if "set" in s.columns:
        s = s[s.set == "primary"]
    sets[nm] = set(s.pfam)
allset = sorted(set().union(*sets.values()))
fi = {f: i for i, f in enumerate(allset)}
rng = np.random.default_rng(SEED)
rows = []
for dname, sgn in (("expanded", 1), ("contracted", -1)):
    NH, OB = {}, {}
    for e in EVENTS:
        d = D[e]
        M = np.zeros((B, len(allset)), bool)
        for _, g in d.groupby("dec"):
            k = int((g.s == sgn).sum())
            if k == 0:
                continue
            inset = np.array([f in fi for f in g.index])
            if not inset.any():
                continue
            cols = np.array([fi[f] for f in g.index[inset]])
            pos = np.where(inset)[0]
            for b0 in range(0, B, 2000):
                nb = min(2000, B - b0)
                pick = np.argpartition(rng.random((nb, len(g))), k - 1, axis=1)[:, :k]
                hit = np.isin(pick, pos)
                rr, cc = np.nonzero(hit)
                # map picked positions to set columns
                posmap = {p: c for p, c in zip(pos, cols)}
                M[b0 + rr, [posmap[p] for p in pick[rr, cc]]] = True
        NH[e] = M
        OB[e] = np.array([f in d.index and d.at[f, "s"] == sgn for f in allset])
    for lname in ("7 origins", "7 origins collapsed (origin OR nested)"):
        if lname == "7 origins":
            nl = sum(NH[e].astype(np.int8) for e in R1_ORIGINS)
            ob = sum(OB[e].astype(int) for e in R1_ORIGINS)
        else:
            nl = np.zeros((B, len(allset)), np.int8)
            ob = np.zeros(len(allset), int)
            for u in R1_ORIGINS:
                es = [e for e in EVENTS if UNIT[e] == u]
                nl += np.logical_or.reduce([NH[e] for e in es])
                ob += np.logical_or.reduce([OB[e] for e in es])
        for nm, mem_set in sets.items():
            mem = np.array([f in mem_set for f in allset])
            for k in (1, 2, 3):
                o = int(((ob >= k) & mem).sum())
                n_null = (nl[:, mem] >= k).sum(1)
                rows.append(dict(set=nm, direction=dname, ladder=lname, tier_k=k, observed=o, null_mean=n_null.mean(),
                                 p_perm=((n_null >= o).sum() + 1) / (B + 1),
                                 null="detectability-stratified (size deciles)"))
S = pd.DataFrame(rows)
S.to_csv(out("set_families_ladder_stratified.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200)
print(S[S.tier_k >= 2].to_string(index=False))
