import numpy as np, pandas as pd, sys
P = sys.argv[1]
bv = pd.read_csv(f"{P}/results/primary/euk/branch_values.tsv", sep="\t", index_col=0)
ms = pd.read_csv(f"{P}/results/primary/euk/matched_sets.tsv", sep="\t")
pt = pd.read_csv(f"{P}/results/primary/euk/primary_tests.tsv", sep="\t")
dirn = pd.read_csv(f"{P}/data/h12_expected_direction_frozen.tsv", sep="\t").set_index("transition_id")["expected"]
def u(sb, others):
    o = np.asarray(others, float)
    return ((o < sb).sum() + 0.5 * (o == sb).sum() + 0.5) / (o.size + 1)
def run(hyp, col, alt, signs=None, ndraw=10000, seed=20261029):
    m = ms[(ms.hypothesis == hyp) & (ms.set == "primary")]
    sets = {t: list(g.sort_values("rank").control) for t, g in m.groupby("transition_id")}
    stems = dict(m.groupby("transition_id").stem.first())
    tids = sorted(sets)
    v = bv[col]
    for t in tids:
        sets[t] = [c for c in sets[t] if pd.notna(v.get(c))]
    tids = [t for t in tids if pd.notna(v.get(stems[t]))]
    w = {t: (signs[t] if signs else 1.0) for t in tids}
    uo = {t: u(w[t] * v[stems[t]], [w[t] * v[c] for c in sets[t]]) for t in tids}
    T = np.mean(list(uo.values()))
    ids = sorted({c for t in tids for c in sets[t]}) + [stems[t] for t in tids]
    gid = {x: i for i, x in enumerate(ids)}
    vals = v.reindex(ids).to_numpy(float)
    rng = np.random.default_rng(seed)
    Tn = np.empty(ndraw)
    for d in range(ndraw):
        pv = vals[rng.permutation(len(vals))]
        Tn[d] = np.mean([u(w[t] * pv[gid[stems[t]]], w[t] * pv[[gid[c] for c in sets[t]]]) for t in tids])
    if alt == "greater":
        hits = (Tn >= T - 1e-12).sum()
    else:
        hits = (np.abs(Tn - .5) >= abs(T - .5) - 1e-12).sum()
    p = (1 + hits) / (ndraw + 1)
    rep = pt[(pt.hypothesis == hyp) & (pt.set == "primary") & (pt.analysis == "primary")].iloc[0]
    print(f"{hyp}: n={len(tids)} k={[len(sets[t]) for t in tids]} T_obs mine={T:.6f} reported={rep.T_obs:.6f} | "
          f"p mine (own RNG, permute, {ndraw})={p:.4f} +- {np.sqrt(p*(1-p)/ndraw):.4f} reported={rep.p:.4f}")
    print("   u:", {t: round(x, 4) for t, x in uo.items()}, " w:", w)
    return T, p
appl = [t for t in dirn.index if dirn[t] != "not_applicable"]
signs = {t: (1.0 if dirn[t] == "loss" else -1.0) for t in appl}
print("H12 frozen directions:", dict(dirn))
import json
res = {}
res["H12"] = run("H12", "H12_f", "greater", signs)
res["H4"] = run("H4", "H4_s", "greater")
res["H5"] = run("H5", "H5_s", "greater")
