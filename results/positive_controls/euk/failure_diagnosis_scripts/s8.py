exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
import collections
H=cache.hogs; G=cache.genes; hits=cache.hits
gidx=pd.Series(np.arange(len(G)),index=G["gene"].to_numpy())
print("PF07657 (MNNL) in prereg hits:", (hits.pfam_acc=="PF07657").any())
mn=hits[hits.pfam_acc=="PF07657"].groupby("species").protein.nunique()
for g in ["Choanoflagellata","Filasterea","Ichthyosporea"]:
    print(" MNNL",g,{s:int(mn.get(s,0)) for s in sorted(ctx.clade(g))})
rx=pd.read_csv(f"{'/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad'}/relaxed_hits.tsv",sep="\t")
cho=ctx.clade("Choanoflagellata")
def grp(s):
    for g in ["Metazoa","Choanoflagellata","Filasterea","Ichthyosporea","Holozoa","Holomycota"]:
        if s in ctx.clade(g): return g
    return s
dsl=an.proteins_with(hits,["PF01414"])
def chain(h):
    out=[]
    while h>=0: out.append(h); h=L.parent[h]
    return out
for r_ in rx[(rx.pfam=="PF01414")&(rx.species.isin(cho|ctx.clade("Filasterea")|ctx.clade("Ichthyosporea")))].itertuples():
    p=G.parent_idx.iat[gidx[r_.protein]]
    if p<0: print(r_.species,r_.protein,"singleton"); continue
    ch=chain(p)
    # find the Choanozoa / Holozoa-level ancestor in chain
    desc=[]
    for x in ch:
        ms=L.member_species([x])[x]; mg=L.member_genes([x])[x]
        desc.append(f"{H.hog_id.iat[x]}[n={len(mg)},DSL_GA={len(mg&dsl)},Met={sum(grp(s)=='Metazoa' for s in ms)}]")
    print(r_.species, f"score={r_.seq_score}", " <- ".join(desc[:6]))
# Metazoan DSL PF-HOGs: topmost and their levels
rows=an.pf_hogs(cache,["PF01414"]); s=set(rows.tolist())
top=[h for h in rows if L.parent[h] not in s]
print("\nTop-most DSL PF-HOGs:")
for h in top:
    ms=collections.Counter(grp(x) for x in L.member_species([h])[h])
    print(" ",H.hog_id.iat[h],"origin",L.origin_level(h),dict(ms))
