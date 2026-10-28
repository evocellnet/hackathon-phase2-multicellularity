exec(open("load.py").read())
import collections
H=cache.hogs; G=cache.genes; ev=cache.events
print(H.origin.value_counts().head()); print(H.head(3).T)
pr=an.proteins_with(cache.hits,["PF01049"])
def chain(h):
    out=[]
    while h>=0:
        out.append(h); h=L.parent[h]
    return out
def frac(h,prots):
    mg=L.member_genes([h])[h]; return len(mg&prots),len(mg)
op=ctx.clade("Opisthokonta"); hol=ctx.clade("Holozoa")
for pf in ["PF01049","PF00362"]:
    prots=an.proteins_with(cache.hits,[pf])
    rows=an.pf_hogs(cache,[pf])
    # the topmost PF-HOGs (whose parent is not PF-HOG)
    s=set(rows.tolist())
    top=[h for h in rows if L.parent[h] not in s]
    print(f"\n=== {pf}: {len(rows)} PF-HOGs, {len(top)} top-most")
    for h in top:
        ch=chain(h)
        desc=" <- ".join(f"{H.hog_id.iat[x]}[{frac(x,prots)[0]}/{frac(x,prots)[1]}]" for x in ch)
        print("  ",desc)
    r=L.origin_row(rows[0]); ms=collections.Counter(L.member_species([r])[r])
    print("  root composition: Holozoa genes",sum(v for k,v in ms.items() if k in hol),
          "other Opisthokonta",sum(v for k,v in ms.items() if k in op-hol),
          "non-Opisthokonta",{k:v for k,v in ms.items() if k not in op})
    # events on the root lineage family at Metazoa stem / Choanozoa
    fam=H.family_id.iat[r]
    e=ev[(ev.family_id==fam)]
    print("  events in family by branch (top):"); print(e.groupby(["branch_child","event"])["n_events"].sum().sort_values(ascending=False).head(15).to_string())
