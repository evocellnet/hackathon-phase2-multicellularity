import pandas as pd
exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
d=pd.read_csv("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/relaxed_hits.tsv",sep="\t")
grp={}
for g in ["Metazoa","Choanoflagellata","Filasterea","Ichthyosporea","Holozoa","Holomycota"]:
    for s in ctx.clade(g): grp.setdefault(s,g)
d["group"]=d.species.map(lambda s:grp.get(s,"Apusozoa"))
for pf,g in d.groupby("pfam"):
    print("\n##",pf)
    nm=g.groupby(["group","species"]).agg(n_relaxed=("protein","size"),n_GA=("passes_GA","sum"),max_score=("seq_score","max"),min_E=("seq_evalue","min")).reset_index()
    print(nm[nm.group!="Metazoa"].to_string(index=False))
    print("Metazoa GA/relaxed:", ", ".join(f"{r.species}:{r.n_GA}/{r.n_relaxed}" for r in nm[nm.group=="Metazoa"].itertuples()))
