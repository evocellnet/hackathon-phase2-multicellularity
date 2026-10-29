"""W7 step 2: duplication events x LECA complexes (stems vs controls)."""
import os, pickle, re, collections
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"; OUT = f"{ROOT}/results/r2/w7"; SCR = os.environ["SCR"]
CUT = os.environ.get("CUT", "cut_1202")
M = pickle.load(open(f"{SCR}/w7_map.pkl", "rb")); og_hogs = M["og_hogs"]; og_fams = M["og_fams"]; SP = M["SP"]
cl = pd.read_csv(os.environ["CLCSV"]); cl = cl.loc[:, ~cl.columns.str.startswith("Unnamed")]
t = pd.read_csv(f"{ROOT}/data/raw/marcotte/TableS1-LECA_OGs_annotated.csv", low_memory=False).set_index("ID")
cl = cl[cl[CUT] != 0]
MAXSZ = 40
allm = {k: list(v.ID) for k, v in cl.groupby(CUT)}
members = {k: v for k, v in allm.items() if len(v) <= MAXSZ}
print("clusters", len(allm), "kept (<=%d OGs)" % MAXSZ, len(members))
def cname(sub):
    c = collections.Counter()
    for col in ["cp_cmplx_name", "corum_cmplx_name", "go_cmplx_name"]:
        for s in sub[col].dropna():
            for n in re.split(r",\s(?=[A-Z0-9])|\|", s): c[n.strip()] += 1
    genes = [g.split(",")[0].strip() for g in sub.human_gene_names_primary.dropna()][:6]
    return (c.most_common(1)[0][0] if c and c.most_common(1)[0][1] >= 2 else "unnamed/mixed") + " [" + ",".join(genes) + "]"
names = {k: cname(cl[cl.ID.isin(v)]) for k, v in members.items()}
# complex -> HOG set (strict) and family set (coarse); and hog -> subunit OGs (for subunit listing)
cx_hogs = {k: set().union(*[og_hogs[o] for o in v]) for k, v in members.items()}
cx_fams = {k: set().union(*[og_fams[o] for o in v]) for k, v in members.items()}
cx_nmapped = {k: sum(1 for o in v if og_hogs[o]) for k, v in members.items()}
hog2og = collections.defaultdict(set)
for o, hs in og_hogs.items():
    if o in set(cl.ID):
        for h in hs: hog2og[h].add(o)
allcx_hogs = set(hog2og)
allcx_fams = set().union(*[og_fams[o] for o in cl.ID])
def _gn(o):
    for c in ["human_gene_names_primary", "arath_gene_names_primary"]:
        v = t.loc[o, c]
        if pd.notna(v): return str(v).split(",")[0].strip()
    v = t.loc[o, "human_protein_names"]
    return o + ("(" + str(v)[:25] + ")" if pd.notna(v) else "")
gname = {o: _gn(o) for o in cl.ID}
# testable HOGs: ancestral genes that contain >=1 gene of a Cox species in our set
genes = pd.read_parquet(f"{ROOT}/runs/analysis/euk/genes.parquet", columns=["species", "parent_idx"])
hogs = pd.read_parquet(f"{ROOT}/runs/analysis/euk/hogs.parquet", columns=["hog_id", "parent_idx", "family_id"])
hpar = hogs.parent_idx.to_numpy(); hid = hogs.hog_id.to_numpy()
cox_sp = set(SP.values())
gi = genes.parent_idx.to_numpy()[genes.species.isin(cox_sp).to_numpy()]
cur = np.unique(gi[gi >= 0]); mark = np.zeros(len(hogs), bool)
while len(cur):
    mark[cur] = True; cur = np.unique(hpar[cur][hpar[cur] >= 0]); cur = cur[~mark[cur]]
testable = set(hid[mark])
ev = pd.read_parquet(f"{ROOT}/runs/analysis/euk/events.parquet", columns=["branch_parent","branch_child","event","n_events","at_origin","hog_id","family_id","child_hog_ids","n_species"])
ev = ev[ev.event == "dup"].copy()
# retained: >=2 child copies each with >=2 extant species
nsp = dict(zip(hogs.hog_id, hogs.family_id))  # placeholder
hn = pd.read_parquet(f"{ROOT}/runs/analysis/euk/hogs.parquet", columns=["hog_id", "n_species"]); hn = dict(zip(hn.hog_id, hn.n_species))
def retained(s):
    k = 0
    for c in s.split(";"):
        c = c.split("|", 1)[-1]
        # child ids are "Level|id"; leaf gene copies count as 1 species
        if hn.get(s.split(";")[0].split("|")[0] + "|" + c, 1) >= 2: k += 1
    return k >= 2
ev["retained"] = ev.child_hog_ids.map(retained)
bv = pd.read_csv(f"{ROOT}/results/primary/euk/branch_values.tsv", sep="\t")
STEMS = {"Metazoa": ["Metazoa"], "Fungal hyphae": ["Blastocladiomycota_Amastigomycota", "Amastigomycota"], "Dictyostelia": ["Dictyostelia"],
         "Streptophyta": ["Klebsormidiophyceae_Phragmoplastophyta"], "Rhodophyta": ["Eurhodophytina"], "PX clade": ["PX_clade"],
         "Volvocine": ["TGV_clade", "Gonium_Volvocaceae"]}
ctrl = list(bv.branch[bv.in_U])
ctrl_broad = list(bv.branch[(~bv.is_terminal) & (bv.n_desc >= 2) & bv.U_exclusion.isin(["", "adjacent_to_stem", "class:other", "sister_tips<2"]) ])
ev["testable"] = ev.hog_id.isin(testable)
ev["cx_strict"] = ev.hog_id.isin(allcx_hogs)
ev["cx_family"] = ev.family_id.isin(allcx_fams)
def grp(branches, label):
    d = ev[ev.branch_child.isin(branches)].copy(); d["group"] = label; return d
groups = {k: grp(v, k) for k, v in STEMS.items()}
groups["Controls (U, 30 branches)"] = grp(ctrl, "ctrl"); groups["Controls broad"] = grp(ctrl_broad, "ctrlb")
rows = []
for k, d in groups.items():
    tt = d[d.testable]; r = d[d.retained]
    rows.append(dict(group=k, branches=d.branch_child.nunique(), dup_events=len(d), dup_events_weighted=int(d.n_events.sum()),
        testable=len(tt), testable_pct=100*len(tt)/max(len(d),1), strict_cx=int(tt.cx_strict.sum()),
        strict_cx_pct_of_testable=100*tt.cx_strict.mean() if len(tt) else np.nan,
        family_cx=int(d.cx_family.sum()), family_cx_pct_of_all=100*d.cx_family.mean() if len(d) else np.nan,
        retained_testable=int(tt.retained.sum()), retained_cx_pct=100*tt[tt.retained].cx_strict.mean() if tt.retained.sum() else np.nan))
frac = pd.DataFrame(rows); frac.to_csv(f"{OUT}/dup_fraction_by_group.tsv", sep="\t", index=False, float_format="%.3f")
# stem vs control fisher on strict, testable
C = groups["Controls (U, 30 branches)"]; C = C[C.testable]; Cb = groups["Controls broad"]; Cb = Cb[Cb.testable]
pw = []
for k in STEMS:
    d = groups[k]; d = d[d.testable]
    for cn, cc in [("U", C), ("broad", Cb)]:
        if len(d) < 5: continue
        tab = [[d.cx_strict.sum(), len(d)-d.cx_strict.sum()], [cc.cx_strict.sum(), len(cc)-cc.cx_strict.sum()]]
        o, p = fisher_exact(tab); pw.append(dict(stem=k, control=cn, stem_rate=d.cx_strict.mean(), ctrl_rate=cc.cx_strict.mean(), OR=o, p=p, n_stem=len(d), n_ctrl=len(cc)))
pw = pd.DataFrame(pw); pw.to_csv(f"{OUT}/dup_fraction_stem_vs_control_fisher.tsv", sep="\t", index=False, float_format="%.4g")
# per-complex enrichment
def enrich(stem_df, ctrl_df, label):
    out = []
    for k in members:
        if cx_nmapped[k] < 3: continue
        a = stem_df.hog_id.isin(cx_hogs[k]).sum(); b = len(stem_df) - a
        c = ctrl_df.hog_id.isin(cx_hogs[k]).sum(); d = len(ctrl_df) - c
        if a + c == 0: continue
        o, p = fisher_exact([[a, b], [c, d]])
        out.append(dict(cluster=k, complex=names[k], n_OG=len(members[k]), n_mapped=cx_nmapped[k], stem_dups=a, stem_n=len(stem_df), ctrl_dups=c, ctrl_n=len(ctrl_df),
                        stem_rate=a/len(stem_df), ctrl_rate=c/len(ctrl_df), OR=o, p=p))
    o = pd.DataFrame(out); o["q_BH"] = multipletests(o.p, method="fdr_bh")[1] if len(o) else []
    o = o.sort_values("p"); o.to_csv(f"{OUT}/complex_enrichment_{label}_{CUT}.tsv", sep="\t", index=False, float_format="%.4g"); return o
prim = pd.concat([groups[k][groups[k].testable] for k in STEMS if k != "Rhodophyta"])
e_pool = enrich(prim, C, "pooled_stems_vs_U")
e_met = enrich(groups["Metazoa"][groups["Metazoa"].testable], C, "Metazoa_vs_U")
e_met_b = enrich(groups["Metazoa"][groups["Metazoa"].testable], Cb, "Metazoa_vs_broad")
# Metazoa-stem duplicated subunits list
mz = groups["Metazoa"][groups["Metazoa"].testable]
L = []
for k in members:
    sub = collections.Counter(); nev = 0
    for h, ne, rt in zip(mz.hog_id, mz.n_events, mz.retained):
        if h in cx_hogs[k]:
            nev += 1
            for o in hog2og[h] & set(members[k]): sub[gname[o]] += 1
    if nev and cx_nmapped[k] >= 3: L.append(dict(cluster=k, complex=names[k], n_mapped=cx_nmapped[k], metazoa_dup_events=nev, ctrl_rate_per_event=C.hog_id.isin(cx_hogs[k]).mean(),
                          metazoa_rate=nev/len(mz), subunits_duplicated=";".join(f"{g}x{n}" for g, n in sub.most_common())))
L = pd.DataFrame(L).sort_values("metazoa_dup_events", ascending=False); L.to_csv(f"{OUT}/metazoa_stem_duplicated_complexes_{CUT}.tsv", sep="\t", index=False, float_format="%.4g")
pickle.dump(dict(names=names, members=members), open(f"{SCR}/w7_cx_{CUT}.pkl", "wb"))
pd.set_option("display.width", 250, "display.max_colwidth", 70)
print(frac.round(2).to_string()); print(pw.round(4).to_string())
print(e_pool.head(12)[["complex","n_mapped","stem_dups","stem_rate","ctrl_dups","ctrl_rate","OR","p","q_BH"]].to_string())
print(e_met.head(12)[["complex","n_mapped","stem_dups","stem_rate","ctrl_dups","ctrl_rate","OR","p","q_BH"]].to_string())
print(L.head(25)[["complex","n_mapped","metazoa_dup_events","metazoa_rate","ctrl_rate_per_event","subunits_duplicated"]].to_string())
print("n complexes", len(members), "with>=3 mapped", sum(v>=3 for v in cx_nmapped.values()), "sig pooled", (e_pool.q_BH<.05).sum(), "sig met", (e_met.q_BH<.05).sum())
