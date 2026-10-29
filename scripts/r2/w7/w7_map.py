"""W7 step 1: map Cox et al. LECA OGs -> our proteins -> round-1 Ancestra HOGs.
Writes results/r2/w7/og_map.parquet-like pickles (cache) and mapping-rate tables."""
import os, re, gzip, pickle, collections
import numpy as np, pandas as pd
ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = f"{ROOT}/results/r2/w7"; os.makedirs(OUT, exist_ok=True)
SP = {"human":"Homosap","caeel":"Caenele","drome":"Drosmel","nemve":"Nemavec","strpu":"Stropur","dicdi":"Dictdis",
      "yeast":"Sacccer","arath":"Arabtha","cerri":"Ceraric","selml":"Selamoe","chlre":"Chlarei","euggr":"Euglgra",
      "phatc":"Phaetri","tetts":"Tetrthe","tryb2":"Trypbru","plaf7":"Plasfal"}
t = pd.read_csv(f"{ROOT}/data/raw/marcotte/TableS1-LECA_OGs_annotated.csv", low_memory=False)
# complex clusters (xlsx exported to csv by hand in scratchpad)
cl = pd.read_csv(os.environ["CLCSV"]); cl = cl.loc[:, ~cl.columns.str.startswith("Unnamed")]
# idmaps
acc2gene = {}
rows = []
for col, code in SP.items():
    m = pd.read_csv(f"{ROOT}/data/proteomes/idmaps/{code}.tsv.gz", sep="\t", dtype=str)
    kept = dict(zip(m.original_id[m.status == "kept"], m.new_id[m.status == "kept"]))
    d = {}
    for o, n, s in zip(m.original_id, m.new_id, m.status):
        if s == "kept": d[o] = n
        else:
            mm = re.match(r"dropped_(?:identical_to|isoform_of):(.+)", s)
            if mm and mm.group(1) in kept: d[o] = kept[mm.group(1)]
    acc2gene[col] = d
genes = pd.read_parquet(f"{ROOT}/runs/analysis/euk/genes.parquet", columns=["gene","species","family_id","parent_idx"])
hogs = pd.read_parquet(f"{ROOT}/runs/analysis/euk/hogs.parquet", columns=["hog_id","parent_idx","family_id"])
hid = hogs.hog_id.to_numpy(); hpar = hogs.parent_idx.to_numpy()
gidx = pd.Series(np.arange(len(genes)), index=genes.gene.values)
gpar = genes.parent_idx.to_numpy(); gfam = genes.family_id.to_numpy()
def chain(gi):
    out = []; c = gpar[gi]
    while c >= 0:
        out.append(hid[c]); c = hpar[c]
    return out
# OG -> mapped genes
og_genes = {}; rate = collections.defaultdict(lambda: [0, 0, 0])  # species: n_acc, n_mapped_gene, n_with_hog
for og, *_ in t[["ID"]].itertuples(index=False): og_genes[og] = set()
for col in SP:
    s = t[["ID", col]].dropna()
    for og, v in zip(s.ID, s[col]):
        for a in v.split(","):
            a = a.strip().split("-")[0]
            if not a: continue
            rate[col][0] += 1
            g = acc2gene[col].get(a)
            if g is not None and g in gidx.index:
                rate[col][1] += 1; og_genes[og].add(g)
                if gpar[gidx[g]] >= 0: rate[col][2] += 1
r = pd.DataFrame(rate, index=["n_accessions", "n_mapped_to_our_protein", "n_with_ancestra_hog"]).T
r["pct_mapped"] = 100 * r.n_mapped_to_our_protein / r.n_accessions
r["pct_with_hog_of_mapped"] = 100 * r.n_with_ancestra_hog / r.n_mapped_to_our_protein
r.index = [f"{c} ({SP[c]})" for c in r.index]
r.to_csv(f"{OUT}/mapping_rates_by_species.tsv", sep="\t", float_format="%.2f")
# per-OG
rec = []; og_hogs = {}; og_fams = {}
for og, gs in og_genes.items():
    hs = set(); fs = set(); nh = 0
    for g in gs:
        gi = gidx[g]
        if gpar[gi] >= 0:
            nh += 1; hs.update(chain(gi))
        fs.add(gfam[gi])
    og_hogs[og] = hs; og_fams[og] = {f for f in fs if not str(f).startswith(tuple("ABCDEFGHIJKLMNOPQRSTUVWXYZ")) or str(f).startswith("HOG")}
    rec.append((og, len(gs), nh, len(hs)))
o = pd.DataFrame(rec, columns=["OG", "n_genes", "n_genes_with_hog", "n_hogs"]).set_index("OG")
o["any_cox_acc"] = t.set_index("ID")[list(SP)].notna().any(axis=1)
o["in_complex_map"] = o.index.isin(cl.ID)
o.to_csv(f"{OUT}/og_mapping.tsv", sep="\t")
S = []
for name, m in [("all LECA OGs (TableS1)", o), ("OGs in the CF-MS complex map", o[o.in_complex_map])]:
    S.append(f"{name}: n={len(m)}; with accession in our 16 shared species={m.any_cox_acc.sum()} ({100*m.any_cox_acc.mean():.1f}%); "
             f"mapped to >=1 of our proteins={(m.n_genes>0).sum()} ({100*(m.n_genes>0).mean():.1f}%); "
             f">=1 protein in an Ancestra HOG={(m.n_genes_with_hog>0).sum()} ({100*(m.n_genes_with_hog>0).mean():.1f}%)")
open(f"{OUT}/mapping_rates_summary.txt", "w").write("\n".join(S) + "\n")
print("\n".join(S)); print(r.round(1).to_string())
pickle.dump(dict(og_hogs=og_hogs, og_fams=og_fams, og_genes=og_genes, SP=SP), open(f"{os.environ['SCR']}/w7_map.pkl", "wb"))
