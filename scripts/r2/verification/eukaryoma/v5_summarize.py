"""Compare verifier per-branch results with the producer's; pooled agreement; >=80% mapped
restriction; depth analysis of the stems-vs-controls contrast."""
import numpy as np, pandas as pd
import statsmodels.api as sm
P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
O = f"{P}/results/r2/verification/eukaryoma"
v = pd.read_csv(f"{O}/v3_per_branch.tsv", sep="\t")
p = pd.concat([pd.read_csv(f"{P}/results/eukaryoma_xval/per_control_summary.tsv", sep="\t"), pd.read_csv(f"{P}/results/eukaryoma_xval/per_transition_summary.tsv", sep="\t")])
p = p[p.status == "exact"].drop_duplicates("pruned_stem")
bm = pd.read_csv(f"{P}/results/eukaryoma_xval/branch_mapping.tsv", sep="\t").drop_duplicates("branch_id")
p["child"] = p.pruned_stem.str.split(">").str[-1]
STEMS = ["Metazoa", "Dictyostelia", "Klebsormidiophyceae_Phragmoplastophyta", "Eurhodophytina",
         "Rhodymeniophycidae", "PX_clade"]
CTRL = list(bm[(bm.role == "control") & (bm.status == "exact")].bottom)
CTRL = [c for c in CTRL if c in set(v.child)]
pd.set_option("display.width", 250)
out = []

# ---- task 2: the 5 focal branches
FOC = ["Metazoa", "Klebsormidiophyceae_Phragmoplastophyta", "Viridiplantae", "Amoebozoa", "Fungi"]
pp = p.set_index("child"); vv = v.set_index("child")
rows = []
for c in FOC:
    for mine, theirs in [("A_gains", "A_gains_supported"), ("E_gains", "E_gains_supported"),
                         ("A_losses", "A_losses"), ("E_losses", "E_losses"), ("A_dups", "A_dups"), ("E_dups", "E_dups")]:
        a, b = vv.loc[c, mine], pp.loc[c, theirs]
        rows.append(dict(branch=c, quantity=mine, verifier=a, producer=b, rel_diff_pct=round(100 * (a / b - 1), 1),
                         ok=abs(a / b - 1) <= 0.05))
    for d in ("A2E_gain", "E2A_gain", "A2E_dup", "E2A_dup"):
        a, b = vv.loc[c, f"{d}_frac_same"], pp.loc[c, f"{d}_frac_same"]
        rows.append(dict(branch=c, quantity=f"{d}_frac_same", verifier=round(a, 4), producer=b,
                         diff_pp=round(100 * (a - b), 2), ok=abs(a - b) <= 0.02,
                         n_verifier=vv.loc[c, f"{d}_n_mapped"], n_producer=pp.loc[c, f"{d}_n"]))
foc = pd.DataFrame(rows)
foc.to_csv(f"{O}/v5_focal_branches.tsv", sep="\t", index=False)
print(foc.to_string())

# ---- pooled agreement (6 exact internal stems vs 26 exact controls), all and >=80% mapped
def pool(sel, d, tag=""):
    s = v[v.child.isin(sel)]
    return s[f"{d}{tag}_same"].sum(), s[f"{d}{tag}_n_mapped"].sum()
prow = []
for d in ("A2E_gain", "E2A_gain", "A2E_dup", "E2A_dup"):
    for grp, sel in (("stems6", STEMS), ("controls", CTRL)):
        r = dict(direction_event=d, group=grp, n_branches=len(sel))
        for tag in ("", "_ge80"):
            k, n = pool(sel, d, tag)
            r[f"same{tag}"] = k; r[f"n_mapped{tag}"] = n; r[f"frac_same{tag}"] = round(k / n, 4)
        prow.append(r)
pooled = pd.DataFrame(prow)
claims = {("A2E_gain", "stems6"): 0.290, ("A2E_gain", "controls"): 0.106, ("E2A_gain", "stems6"): 0.094,
          ("E2A_gain", "controls"): 0.142, ("A2E_dup", "stems6"): 0.038, ("A2E_dup", "controls"): 0.045,
          ("E2A_dup", "stems6"): 0.078, ("E2A_dup", "controls"): 0.050}
pooled["claimed"] = [claims[(a, b)] for a, b in zip(pooled.direction_event, pooled.group)]
pooled["diff_pp"] = (100 * (pooled.frac_same - pooled.claimed)).round(2)
pooled.to_csv(f"{O}/v5_pooled.tsv", sep="\t", index=False)
print(pooled.to_string())

# per-stem >=80% restriction
st = v[v.child.isin(STEMS + CTRL)][["child", "depth", "n_tips"] + [c for c in v if "frac_same" in c or c.endswith("n_mapped")]]
st.to_csv(f"{O}/v5_stems_controls_detail.tsv", sep="\t", index=False)
print(v[v.child.isin(STEMS)][["child", "A2E_gain_frac_same", "A2E_gain_ge80_frac_same", "A2E_gain_n_mapped",
                              "A2E_gain_ge80_n_mapped", "E2A_gain_frac_same", "E2A_gain_ge80_frac_same"]].round(3).to_string())

# ---- task 4: depth
v["role"] = np.where(v.child.isin(STEMS), "stem", np.where(v.child.isin(CTRL), "control", "other"))
print(v.groupby("role")[["depth", "n_tips"]].describe().round(1).to_string())
res = []
for d in ("A2E_gain", "E2A_gain", "A2E_dup", "E2A_dup"):
    x = v[v[f"{d}_n_mapped"] > 0].copy()
    x["logtips"] = np.log(x.n_tips)
    y = np.column_stack([x[f"{d}_same"], x[f"{d}_n_mapped"] - x[f"{d}_same"]])
    base = x.role != "stem"
    X = sm.add_constant(x[["depth", "logtips"]])
    g = sm.GLM(y[base.values], X[base.values], family=sm.families.Binomial()).fit()
    pr = g.predict(X)
    s = (x.role == "stem").values; c = (x.role == "control").values
    exp_stem = (pr[s] * x[f"{d}_n_mapped"][s]).sum() / x[f"{d}_n_mapped"][s].sum()
    exp_ctrl = (pr[c] * x[f"{d}_n_mapped"][c]).sum() / x[f"{d}_n_mapped"][c].sum()
    obs_s = x[f"{d}_same"][s].sum() / x[f"{d}_n_mapped"][s].sum()
    obs_c = x[f"{d}_same"][c].sum() / x[f"{d}_n_mapped"][c].sum()
    # model with stem indicator on stems + all non-stem internal branches
    x["stem"] = s.astype(float)
    g2 = sm.GLM(y, sm.add_constant(x[["depth", "logtips", "stem"]]), family=sm.families.Binomial()).fit()
    g3 = sm.GLM(y, sm.add_constant(x[["stem"]]), family=sm.families.Binomial()).fit()
    # depth-matched comparison: controls with depth within stems' range
    lo, hi = x.depth[s].min(), x.depth[s].max()
    cm = c & (x.depth >= lo).values & (x.depth <= hi).values
    obs_cm = x[f"{d}_same"][cm].sum() / x[f"{d}_n_mapped"][cm].sum() if cm.any() else np.nan
    res.append(dict(direction_event=d, obs_stems=round(obs_s, 4), obs_controls=round(obs_c, 4),
                    exp_stems_depth_size=round(exp_stem, 4), exp_controls_depth_size=round(exp_ctrl, 4),
                    frac_gap_explained=round((exp_stem - exp_ctrl) / (obs_s - obs_c), 2) if obs_s != obs_c else np.nan,
                    stem_logOR_raw=round(g3.params["stem"], 3), stem_logOR_adj=round(g2.params["stem"], 3),
                    depth_coef=round(g.params["depth"], 3), logtips_coef=round(g.params["logtips"], 3),
                    n_controls_depth_matched=int(cm.sum()), obs_controls_depth_matched=round(obs_cm, 4),
                    stem_depth_range=f"{lo}-{hi}", n_baseline=int(base.sum())))
dep = pd.DataFrame(res)
dep.to_csv(f"{O}/v5_depth.tsv", sep="\t", index=False)
print(dep.T.to_string())
print("stem depths", dict(zip(v.child[v.role == 'stem'], v.depth[v.role == 'stem'])))
print("control depths", sorted(v.depth[v.role == 'control']))
