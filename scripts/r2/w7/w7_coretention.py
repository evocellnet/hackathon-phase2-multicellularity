"""W7 step 3: co-retention of complex subunits across the 187 species (presence = >=1 gene in the OG's Ancestra family)."""
import os, pickle, collections
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"; OUT = f"{ROOT}/results/r2/w7"; SCR = os.environ["SCR"]; CUT = os.environ.get("CUT", "cut_2014")
NPERM = 1000; rng = np.random.default_rng(11)
M = pickle.load(open(f"{SCR}/w7_map.pkl", "rb")); og_fams = M["og_fams"]
X = pickle.load(open(f"{SCR}/w7_cx_{CUT}.pkl", "rb")); names, members = X["names"], X["members"]
genes = pd.read_parquet(f"{ROOT}/runs/analysis/euk/genes.parquet", columns=["species", "family_id"])
genes = genes[genes.family_id.str.startswith("HOG")].drop_duplicates()
tr = pd.read_csv(f"{ROOT}/data/species_traits.tsv", sep="\t"); tr = tr[tr.kingdom == "euk"].set_index("code")
sp = sorted(genes.species.unique()); sidx = {s: i for i, s in enumerate(sp)}
fam_sp = genes.groupby("family_id").species.apply(lambda s: [sidx[x] for x in s])
def vec(fams):
    v = np.zeros(len(sp), bool)
    for f in fams:
        if f in fam_sp.index: v[fam_sp[f]] = True
    return v
all_og = list(og_fams)
pres = {o: vec(og_fams[o]) for o in all_og if og_fams[o]}
mc = np.array([tr.loc[s, "mc_state"] != "unicellular" and tr.loc[s, "reversal"] == 0 and tr.loc[s, "multicellular"] == 1 for s in sp])
uni = np.array([tr.loc[s, "mc_state"] == "unicellular" for s in sp])
print("species", len(sp), "multicellular", mc.sum(), "unicellular", uni.sum())
cx_og = {k: [o for o in v if o in pres] for k, v in members.items()}
# dedupe subunits with identical family sets (same family -> trivial J=1)
def dedupe(ogs):
    seen = {}; 
    for o in ogs: seen.setdefault(frozenset(og_fams[o]), o)
    return list(seen.values())
cx_u = {k: dedupe(v) for k, v in cx_og.items()}; cx_u = {k: v for k, v in cx_u.items() if len(v) >= 3}
cx_all = set(o for v in members.values() for o in v)
bg = [o for o in pres if o not in cx_all]; bg_all = list(pres)   # null pool = every mapped LECA OG (complex or not)
P = np.array([pres[o] for o in bg_all]); oidx = {o: i for i, o in enumerate(bg_all)}
def meanJ(B):  # B: k x n bool
    B = B.astype(np.float32); inter = B @ B.T; n = B.sum(1); union = n[:, None] + n[None, :] - inter
    J = np.where(union > 0, inter / np.maximum(union, 1), 0); k = len(B); return (J.sum() - np.trace(J)) / (k * (k - 1))
def analyse(mask, label):
    Pm = P[:, mask]; prev = Pm.sum(1); nsp = mask.sum()
    qs = np.unique(np.quantile(prev, np.linspace(0, 1, 21))); bins = np.clip(np.digitize(prev, qs[1:-1], right=True), 0, len(qs) - 2)
    binmem = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    out = []
    for k, ogs in cx_u.items():
        ii = np.array([oidx[o] for o in ogs]); obs = meanJ(Pm[ii]); null = np.empty(NPERM)
        for r in range(NPERM):
            pick = []
            for i in ii:
                pool = binmem[bins[i]]; pool = pool[~np.isin(pool, pick + list(ii))] if len(pool) > len(ii) + len(pick) + 1 else pool
                pick.append(int(rng.choice(pool)))
            null[r] = meanJ(Pm[pick])
        out.append(dict(cluster=k, complex=names[k], n_subunits=len(ogs), mean_prev=prev[ii].mean() / nsp, J=obs, null_mean=null.mean(), null_sd=null.std(),
                        z=(obs - null.mean()) / max(null.std(), 1e-9), p_perm=(1 + (null >= obs).sum()) / (NPERM + 1)))
    d = pd.DataFrame(out); d["set"] = label; return d
res = {}
for lab, mask in [("all", np.ones(len(sp), bool)), ("multicellular", mc), ("unicellular", uni)]:
    res[lab] = analyse(mask, lab); print(lab, mask.sum(), "done", flush=True)
a = res["all"].set_index("cluster"); m = res["multicellular"].set_index("cluster"); u = res["unicellular"].set_index("cluster")
c = a[["complex", "n_subunits", "J", "z", "p_perm"]].join(m[["J", "z", "p_perm"]].add_suffix("_mc")).join(u[["J", "z", "p_perm"]].add_suffix("_uni"))
c["dJ_mc_minus_uni"] = c.J_mc - c.J_uni; c["dz_mc_minus_uni"] = c.z_mc - c.z_uni
# subsample unicellular to size of mc species to check species-number sensitivity of J
ui = np.where(uni)[0]; sub = []
for r in range(20):
    mk = np.zeros(len(sp), bool); mk[rng.choice(ui, mc.sum(), replace=False)] = True
    Pm = P[:, mk]; sub.append([meanJ(Pm[[oidx[o] for o in ogs]]) for ogs in cx_u.values()])
c["J_uni_subsampled"] = np.mean(sub, axis=0); c["dJ_mc_minus_unisub"] = c.J_mc - c.J_uni_subsampled
c.to_csv(f"{OUT}/coretention_by_complex_{CUT}.tsv", sep="\t", float_format="%.4g")
s = []
s.append(f"complexes with >=3 distinct mapped subunit families: {len(c)}; species={len(sp)} (multicellular {mc.sum()}, unicellular {uni.sum()})")
for lab in ["all", "multicellular", "unicellular"]:
    d = res[lab]; s.append(f"{lab}: median J={d.J.median():.3f}, median null={d.null_mean.median():.3f}, median z={d.z.median():.2f}, perm p<0.05 for {(d.p_perm<0.05).sum()}/{len(d)} (expect {0.05*len(d):.1f}); z>2: {(d.z>2).sum()}")
w = wilcoxon(c.z_mc, c.z_uni); s.append(f"paired Wilcoxon z_mc vs z_uni: p={w.pvalue:.3g}; median dz={c.dz_mc_minus_uni.median():.2f}; complexes with dz>0: {(c.dz_mc_minus_uni>0).sum()}/{len(c)}")
w2 = wilcoxon(c.J_mc, c.J_uni_subsampled); s.append(f"J multicellular vs unicellular subsampled to {mc.sum()} species: median J_mc={c.J_mc.median():.3f} vs {c.J_uni_subsampled.median():.3f}; paired Wilcoxon p={w2.pvalue:.3g}")
low = c.sort_values("dz_mc_minus_uni").head(5); hi = c.sort_values("dz_mc_minus_uni").tail(5)
s.append("most MC-lower co-retention: " + "; ".join(f"{r.complex[:45]} (dz={r.dz_mc_minus_uni:.1f})" for r in low.itertuples()))
s.append("most MC-higher co-retention: " + "; ".join(f"{r.complex[:45]} (dz={r.dz_mc_minus_uni:.1f})" for r in hi.itertuples()))
open(f"{OUT}/coretention_summary_{CUT}.txt", "w").write("\n".join(s) + "\n"); print("\n".join(s))
