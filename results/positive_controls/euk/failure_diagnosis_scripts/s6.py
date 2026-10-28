exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
import collections
H=cache.hogs; G=cache.genes; hits=cache.hits
gidx=pd.Series(np.arange(len(G)),index=G["gene"].to_numpy())
pr=an.proteins_with(hits,["PF00362"])
def grp(s):
    for g in ["Metazoa","Choanoflagellata","Filasterea","Ichthyosporea","Holozoa","Holomycota"]:
        if s in ctx.clade(g): return g
    return s
def chain(h):
    out=[]
    while h>=0: out.append(h); h=L.parent[h]
    return out
for g in sorted(p for p in pr if p.startswith(("Thectra","Didycos","Pygsbif","Creofra","Sphaarc"))):
    p=G.parent_idx.iat[gidx[g]]
    print(g, "->", " <- ".join(f"{H.hog_id.iat[x]}" for x in chain(p)))
print()
for hid in ["Apusomonadida_Opisthokonta|HOG0060409.42","Opisthokonta|HOG0060409.43","Holozoa|HOG0060409.138","Holozoa|HOG0060409.252","Choanozoa|HOG0060409.253"]:
    r=L.row_of[hid]; mg=L.member_genes([r])[r]; ms=L.member_species([r])[r]
    print(hid, "n",len(mg),"hit",len(mg&pr), dict(collections.Counter(grp(s) for s in ms)))
    for ch in L.children(r):
        mgc=L.member_genes([ch])[ch]; msc=L.member_species([ch])[ch]
        if len(mgc)>2 or (mgc&pr): print("     child",H.hog_id.iat[ch],len(mgc),"hit",len(mgc&pr),dict(collections.Counter(grp(s) for s in msc)))
# what domains are on the non-integrin members of the Opisthokonta-level HOG
r=L.row_of["Opisthokonta|HOG0060409.43"]; mg=L.member_genes([r])[r]
non=mg-pr
print("non-integrin members:",len(non),"top Pfams:",hits[hits.protein.isin(non)].pfam_acc.value_counts().head(8).to_dict(), "no prereg hit:", len(non-set(hits.protein)))
print(hits[hits.protein.isin(non)].drop_duplicates("pfam_acc").set_index("pfam_acc").loc[hits[hits.protein.isin(non)].pfam_acc.value_counts().head(8).index,"pfam_name"].to_dict())
# Choanoflagellata crown / Holomycota state
for hid in ["Holozoa|HOG0060409.138","Holozoa|HOG0060409.252","Opisthokonta|HOG0060409.43","Apusomonadida_Opisthokonta|HOG0060409.42"]:
    r=L.row_of[hid]
    for node in ["Choanozoa","Choanoflagellata","Holomycota","Fungi","Opisthokonta","Holozoa"]:
        print(f"  {hid} present_at {node}: {L.present_at(r, ctx.node(node))}", end=";")
    print()
