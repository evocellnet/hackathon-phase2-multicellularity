exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
import collections
H=cache.hogs; hits=cache.hits
for hid in ["Holomycota|HOG0060409.44","Filozoa|HOG0060409.140"]:
    r=L.row_of[hid]; mg=sorted(L.member_genes([r])[r])
    print(hid, collections.Counter(L.member_species([r])[r]))
    if "Holomycota" in hid:
        for g in mg: print("   ",g, hits[hits.protein==g].pfam_acc.tolist())
    for node in ["Choanozoa","Choanoflagellata","Filozoa"]:
        print("   present_at",node, L.present_at(r, ctx.node(node)))
r=L.row_of["Opisthokonta|HOG0060409.43"]
for node in ["Holomycota","Fungi","Choanoflagellata"]:
    print("Opisthokonta|.43 present_at",node,L.present_at(r,ctx.node(node)))
print("Fungi stem: Holomycota -> Fungi; Fontalb in Holomycota not Fungi:", "Fontalb" in ctx.clade("Holomycota"), "Fontalb" in ctx.clade("Fungi"))
