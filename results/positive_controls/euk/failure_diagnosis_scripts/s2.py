exec(open("load.py").read())
import collections
H=cache.hogs; G=cache.genes
print(H.columns.tolist()); print(G.columns.tolist()); print(cache.events.columns.tolist())
gidx=pd.Series(np.arange(len(G)),index=G["gene"].to_numpy())
def gene_lineage(g):
    """root HOG of the gene: (origin hog row, its level, deepest hog row)"""
    i=gidx[g]; p=G["parent_idx"].iat[i]
    if p<0: return None
    return p
def report(pf, extra=()):
    hits=cache.hits; prots=sorted(set(hits.loc[hits.pfam_acc==pf,"protein"]))
    print(f"\n######## {pf}: {len(prots)} hit proteins")
    pfrows=set(an.pf_hogs(cache,[pf]).tolist())
    roots=collections.defaultdict(list)
    for g in prots:
        p=gene_lineage(g)
        if p is None: roots["(singleton)"].append(g); continue
        roots[L.origin_row(p)].append(g)
    pr_all=an.proteins_with(hits,[pf])
    for r,gs in sorted(roots.items(), key=lambda kv:-len(kv[1])):
        if r=="(singleton)":
            print("  singletons (no HOG):", collections.Counter(G.loc[gidx[gs],'species']).most_common()); continue
        ms=L.member_species([r])[r]; mg=L.member_genes([r])[r]
        nhit=len(mg & pr_all)
        sc=collections.Counter(ms)
        hsp=collections.Counter(G.loc[gidx[gs],'species'])
        # PF-HOGs within this lineage
        lin_pf=[x for x in pfrows if L.origin_row(x)==r]
        print(f"  root {H.hog_id.iat[r]} origin={L.origin_level(r)} n_genes={len(mg)} n_species={len(sc)} hit-members={nhit} ({nhit/len(mg):.2f}); PF-HOG rows in lineage: {len(lin_pf)} at levels {sorted(set(T.names[L.level[x]] for x in lin_pf))[:12]}")
        print(f"     hit genes by species: {dict(hsp)}")
        print(f"     member species (top): {sc.most_common(25)}")
report("PF01049")
report("PF00362")
