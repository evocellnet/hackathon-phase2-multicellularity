exec(open("load.py").read())
# Opisthokonta + Apusozoa/near species list
op=ctx.clade("Opisthokonta") if "Opisthokonta" in T.index else set()
print("Apusomonadida_Opisthokonta" in T.index, sorted(ctx.clade("Apusomonadida_Opisthokonta")-op))
groups=["Metazoa","Choanoflagellata","Filasterea","Ichthyosporea","Holomycota","Fungi"]
for g in groups:
    if g in T.index: print(g, sorted(ctx.clade(g)))
print("Holozoa children / path:", path_to_root("Capsowc"))
print("Thectra path:", path_to_root("Thectra"))
dt=sp.set_index("code")["data_type_final"]
for pf in ["PF01049","PF00362","PF00907","PF01414","PF06816","PF07684","PF12499"]:
    hs=pc.hit_species(ctx,[pf])
    print("\n==",pf,"total species",len(hs))
    for g in groups:
        if g in T.index:
            cl=sorted(ctx.clade(g)); print(f"  {g}: "+", ".join(f"{s}:{int(hs.get(s,0))}" for s in cl))
    other=sorted(set(hs.index)-ctx.clade("Opisthokonta"))
    print("  non-opisthokont species with hits:", {s:int(hs[s]) for s in other})
