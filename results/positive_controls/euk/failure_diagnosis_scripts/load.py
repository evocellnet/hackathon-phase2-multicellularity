import sys, os
P="/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
sys.path.insert(0,P); os.chdir(P)
import pandas as pd, numpy as np
from mcell import analysis as an, controls as pc
p=an.run_paths("euk"); cache=an.load_cache("euk",p)
tr,traits,tree,bt=an.load_setup("euk",p)
sp=pd.read_csv(p["species"],sep="\t")
ctx=pc.Ctx(cache=cache,species=sp,transitions=tr,paths=p,bt=bt)
L=ctx.L; T=cache.tree
def path_to_root(s):
    i=T.index[s]; out=[]
    while i>=0: out.append(T.names[i]); i=T.parent[i]
    return out
