exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
import collections
H=cache.hogs; G=cache.genes; hits=cache.hits; ev=cache.events
tb=an.proteins_with(hits,["PF00907"])
def grp(s):
    for g in ["Metazoa","Choanoflagellata","Filasterea","Ichthyosporea","Holozoa","Holomycota"]:
        if s in ctx.clade(g): return g
    return s
gidx=pd.Series(np.arange(len(G)),index=G["gene"].to_numpy())
roots=collections.Counter(); rootsp=collections.defaultdict(collections.Counter)
for g in tb:
    p=G.parent_idx.iat[gidx[g]]
    r = "singleton" if p<0 else L.origin_row(p)
    roots[r]+=1; rootsp[r][grp(G.species.iat[gidx[g]])]+=1
fams=[]
for r,n in roots.most_common():
    if r=="singleton": print("singletons",dict(rootsp[r])); continue
    mg=L.member_genes([r])[r]
    print(f"root {H.hog_id.iat[r]} origin={L.origin_level(r)} n={len(mg)} Tbox={n} by group {dict(rootsp[r])}")
    fams.append(H.family_id.iat[r])
rows=an.pf_hogs(cache,["PF00907"]); s=set(rows.tolist())
print("\nPF-HOGs:",len(rows))
for h in rows:
    if L.level[h]==ctx.node("Metazoa") or L.level[h]==ctx.node("Choanozoa") or L.level[h]==ctx.node("Holozoa") or L.level[h]==ctx.node("Filozoa"):
        print("  ",H.hog_id.iat[h],"n",H.n_genes.iat[h],"origin",L.origin_level(h))
e=ev[ev.family_id.isin(fams)&(ev.branch_child.isin(["Metazoa","Choanozoa","Filozoa","Holozoa","Choanoflagellata"]))]
print(e[["branch_child","event","n_events","at_origin","hog_id","child_hog_ids"]].to_string())
