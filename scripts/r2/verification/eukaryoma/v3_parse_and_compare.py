"""Tasks 2-4: parse the independently filtered OrthoXMLs with mcell.hogs on the
guide tree (pruned by the parser to the 125 species), count gains/losses/dups,
and compute family-level agreement in both directions on every internal branch.
"""
import sys, json
import numpy as np, pandas as pd
sys.path.insert(0, "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity")
from mcell.hogs import parse_hogs

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = f"{P}/results/r2/verification/eukaryoma"
W = f"{OUT}/work"
GT = f"{P}/runs/euk_fast/stages/input/guide_tree.nwk"

R = {}
for m in ("ancestra", "eukaryoma"):
    R[m] = parse_hogs(f"{W}/{m}.filtered.orthoxml.gz", newick=GT)
    print(m, json.dumps({k: v for k, v in R[m].stats.items() if isinstance(v, (int, float))}), flush=True)
A, E = R["ancestra"], R["eukaryoma"]
T = A.tree
assert T.names == E.tree.names and T.parent == E.tree.parent

# tree comparison with the producer's pruned tree
sys.path.insert(0, f"{P}/external/ancestra-mvp-ai")
from ancestra.newick import parse_newick
prod = parse_newick(open(f"{P}/runs/eukaryoma_xval/tree_pruned.nwk").read())
def clades(root):
    d = {}
    for nd in root.postorder():
        d[nd.name] = frozenset(nd.leaf_names()) if nd.children else frozenset([nd.name])
    return d
mine = {T.names[i]: frozenset(T.names[l] for l in T.leaves if T.is_anc(i, l)) for i in range(len(T))}
pc = clades(prod)
tree_same = mine == pc
print("pruned tree identical to producer's (names+clades):", tree_same, "internal nodes", len(T.internal), flush=True)

# protein map
mp = pd.concat([pd.read_csv(f"{P}/runs/eukaryoma_xval/protein_map.tsv.gz", sep="\t"),
                pd.read_csv(f"{P}/runs/eukaryoma_xval/protein_map_rbh.tsv.gz", sep="\t")[["code", "our_id", "oma_id"]]])
a_idx = pd.Series(np.arange(len(A.genes)), index=A.genes.gene.astype(str))
e_idx = pd.Series(np.arange(len(E.genes)), index=E.genes.gene.astype(str))
mp = mp[mp.our_id.isin(a_idx.index) & mp.oma_id.isin(e_idx.index)]
a2e = np.full(len(A.genes), -1, np.int64); a2e[a_idx[mp.our_id].to_numpy()] = e_idx[mp.oma_id].to_numpy()
e2a = np.full(len(E.genes), -1, np.int64); e2a[e_idx[mp.oma_id].to_numpy()] = a_idx[mp.our_id].to_numpy()
print("A genes", len(A.genes), "mapped", (a2e >= 0).sum(), round((a2e >= 0).mean() * 100, 2),
      "| E genes", len(E.genes), "mapped", (e2a >= 0).sum(), flush=True)

def compare(X, Y, x2y, xat, yat, c):
    """events of method X on branch parent(c)->c checked in method Y."""
    hx_par = X.hogs.parent_idx.to_numpy(); hy_par = Y.hogs.parent_idx.to_numpy()
    nsp = X.hogs.n_species.to_numpy()
    g = np.nonzero(xat >= 0)[0]
    d = pd.DataFrame({"xh": xat[g], "yg": x2y[g]})
    d["mapped"] = d.yg >= 0
    d["yh"] = np.where(d.mapped, yat[np.where(d.mapped, d.yg, 0)], -2)
    # per X hog: n members, n mapped, majority Y hog (among members with Y hog at c)
    tot = d.groupby("xh").agg(n=("mapped", "size"), n_mapped=("mapped", "sum"))
    v = d[d.yh >= 0].groupby(["xh", "yh"]).size().rename("k").reset_index()
    v = v.sort_values(["xh", "k", "yh"], ascending=[True, False, True]).drop_duplicates("xh")
    tot["yh"] = v.set_index("xh").yh.reindex(tot.index).fillna(-1).astype(np.int64)
    tot["frac_mapped"] = tot.n_mapped / tot.n
    # gains: one event per family whose ancestral genes at c have no ancestor at p
    # (a top-level paralogGroup on p->c is one gain with k copies); vote over all members
    famx = X.hogs.family_id.to_numpy()
    dg = d[hx_par[d.xh.to_numpy()] < 0].copy(); dg["fam"] = famx[dg.xh.to_numpy()]
    gx = dg.groupby("fam").agg(n=("mapped", "size"), n_mapped=("mapped", "sum"))
    gsp = X.genes.species.cat.codes.to_numpy()[g][hx_par[d.xh.to_numpy()] < 0]
    gx["nsp"] = pd.Series(gsp, index=dg.index).groupby(dg.fam).nunique()
    gx = gx[gx.nsp >= 2].copy()
    vg = dg[dg.yh >= 0].groupby(["fam", "yh"]).size().rename("k").reset_index()
    vg = vg.sort_values(["fam", "k", "yh"], ascending=[True, False, True]).drop_duplicates("fam")
    gx["yh"] = vg.set_index("fam").yh.reindex(gx.index).fillna(-1).astype(np.int64)
    gx["frac_mapped"] = gx.n_mapped / gx.n
    gx["cls"] = np.where(gx.n_mapped == 0, "unmapped", np.where(gx.yh < 0, "shallower",
                np.where(hy_par[gx.yh.clip(lower=0)] < 0, "same", "deeper")))
    # duplications inside existing lineages: X hog at p with >=2 children at c
    ch = tot[hx_par[tot.index] >= 0].copy(); ch["par"] = hx_par[ch.index]
    rows = []
    for par, grp in ch.groupby("par"):
        if len(grp) < 2:
            continue
        S = set(grp.yh[grp.yh >= 0])
        nm = grp.n_mapped.sum()
        if nm == 0:
            cls = "unmapped"
        elif len(S) == 0:
            cls = "other"
        elif len(S) == 1:
            cls = "shallower"
        else:
            pars = [hy_par[s] for s in S]
            nn = [q for q in pars if q >= 0]
            if len(nn) != len(set(nn)):
                cls = "same"
            elif len(nn) >= 2:
                cls = "deeper"
            else:
                cls = "other"
        rows.append(dict(par=par, k=len(grp), cls=cls, frac_mapped=nm / grp.n.sum()))
    dx = pd.DataFrame(rows, columns=["par", "k", "cls", "frac_mapped"])
    return gx, dx

def summ(df, prefix, minfrac=0.0):
    df = df[df.frac_mapped >= minfrac] if minfrac > 0 else df
    n = len(df); un = (df.cls == "unmapped").sum(); nm = n - un
    out = {f"{prefix}_n": n, f"{prefix}_n_mapped": nm}
    for k in ("same", "deeper", "shallower", "other"):
        out[f"{prefix}_{k}"] = int((df.cls == k).sum())
    out[f"{prefix}_frac_same"] = (df.cls == "same").sum() / nm if nm else np.nan
    return out

rows = []
br = {m: R[m].branches.set_index("branch_child") for m in R}
for c in T.internal:
    if T.parent[c] < 0:
        continue
    cn = T.names[c]
    aat = A._gene_hog_at(cn); eat = E._gene_hog_at(cn)
    ga, da = compare(A, E, a2e, aat, eat, cn)
    ge, de = compare(E, A, e2a, eat, aat, cn)
    r = dict(parent=T.names[T.parent[c]], child=cn, depth=T.depth[c],
             n_tips=sum(1 for l in T.leaves if T.is_anc(c, l)),
             A_gains=len(ga), E_gains=len(ge),
             A_losses=int(br["ancestra"].loc[cn, "losses"]), E_losses=int(br["eukaryoma"].loc[cn, "losses"]),
             A_dups=int(br["ancestra"].loc[cn, "dups"]), E_dups=int(br["eukaryoma"].loc[cn, "dups"]),
             A_dup_events=len(da), E_dup_events=len(de))
    for mf, tag in ((0.0, ""), (0.8, "_ge80")):
        r.update(summ(ga, "A2E_gain" + tag, mf)); r.update(summ(ge, "E2A_gain" + tag, mf))
        r.update(summ(da, "A2E_dup" + tag, mf)); r.update(summ(de, "E2A_dup" + tag, mf))
    rows.append(r)
    print(cn, r["A_gains"], r["E_gains"], round(r["A2E_gain_frac_same"], 3) if r["A2E_gain_n_mapped"] else None, flush=True)
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/v3_per_branch.tsv", sep="\t", index=False)
for m in R:
    R[m].branches.to_csv(f"{W}/{m}_branches.tsv", sep="\t", index=False)
print("conservation violations", A.stats["conservation_violations"], E.stats["conservation_violations"])
