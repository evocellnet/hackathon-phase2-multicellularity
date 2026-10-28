exec(open("load.py").read())
import collections
H=cache.hogs; ev=cache.events
met=ctx.clade("Metazoa"); cho=ctx.clade("Choanoflagellata")
def grp(s):
    for g in ["Metazoa","Choanoflagellata","Filasterea","Ichthyosporea","Holomycota"]:
        if s in ctx.clade(g): return g
    return "other:"+s
row=L.row_of
hits=cache.hits
def doms(genes):
    h=hits[hits.protein.isin(genes)]
    return h.pfam_acc.value_counts().head(8).to_dict()
pr=an.proteins_with(hits,["PF01049"]); pc_=an.proteins_with(hits,["PF00028"])
for hid in ["Choanozoa|HOG0058630.288","Choanozoa|HOG0058630.246","Choanozoa|HOG0058630.380","Choanozoa|HOG0064100.40","Choanozoa|HOG0012365.5","Filozoa|HOG0058630.244"]:
    r=row[hid]; mg=L.member_genes([r])[r]; ms=L.member_species([r])[r]
    c=collections.Counter(grp(s) for s in ms)
    print(hid, dict(c), "| PF01049:",len(mg&pr)," PF00028(cadherin rep):",len(mg&pc_))
    # children
    for ch in L.children(r):
        mgc=L.member_genes([ch])[ch]; msc=L.member_species([ch])[ch]
        print("    child",H.hog_id.iat[ch], dict(collections.Counter(grp(s) for s in msc)), "PF01049",len(mgc&pr))
    e=ev[(ev.hog_id==hid)|(ev.child_hog_ids.astype(str).str.contains(hid.split('|')[1]+"(?:;|$)",regex=True))]
print()
# events for the Metazoa-level PF-HOGs on the Metazoa stem
rows=an.pf_hogs(cache,["PF01049"]); ids=set(H.hog_id.to_numpy()[rows])
e=ev[ev.branch_child=="Metazoa"]
m=e[e.hog_id.isin(ids)|e.child_hog_ids.astype(str).map(lambda s:any(x in ids for x in s.split(";")))]
print(m[["branch_parent","branch_child","event","n_events","at_origin","hog_id","child_hog_ids"]].to_string())
print(ev[ev.event=="dup"].head(3).T)
