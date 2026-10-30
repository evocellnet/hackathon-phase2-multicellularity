"""Pick ~100 random sequence-novel families on internal event/control branches; write per-panel full AA/3Di FASTA + query lists."""
import pandas as pd, numpy as np, os, re
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
W=P+'/results/r2/verification/w5/work/spot'; os.makedirs(W,exist_ok=True)
F=pd.read_csv(P+'/results/r2/verification/w5/work/families_verifier.tsv',sep='\t')
Q=pd.read_csv(P+'/results/r2/verification/w5/work/queries_verifier.tsv.gz',sep='\t')
rng=np.random.default_rng(20261030)
F=F[(F.kind=='novel')&(~F.is_leaf)&((F.is_event)|(F.control))&(F.nhas3di>0)]
pick=pd.concat([g.sample(13,random_state=int(rng.integers(1e9))) for _,g in F.groupby('panel')])
pick.to_csv(W+'/picked.tsv',sep='\t',index=False)
src=pd.read_csv(P+'/data/structure3di/panel_sources.tsv',sep='\t',header=None,names=['sp','fa']).set_index('sp')['fa'].to_dict()
def readfa(fn):
    d={};k=None;buf=[]
    for l in open(fn):
        if l[0]=='>':
            if k: d[k]=''.join(buf)
            k=l[1:].split()[0];buf=[]
        else: buf.append(l.strip())
    if k: d[k]=''.join(buf)
    return d
def leaves_from_nwk(s):
    return re.findall(r'[(,]([A-Za-z0-9_]+):',s)
for R,g in pick.groupby('panel'):
    nwk=open(P+f'/results/r2/panels/{R}/tree_bl.nwk').read()
    sps=sorted(set(leaves_from_nwk(nwk)))
    with open(f'{W}/{R}.aa.fa','w') as fa, open(f'{W}/{R}.3di.fa','w') as f3, open(f'{W}/{R}.map.tsv','w') as mp:
        nmiss=0
        for sp in sps:
            fn3=f'{P}/data/structure3di/all/{sp}.3di.fa'
            if not os.path.exists(fn3) or sp not in src: print(R,'no 3di/aa for',sp); continue
            aa=readfa(src[sp]); di=readfa(fn3)
            for k,v in di.items():
                if k not in aa or len(aa[k])!=len(v): nmiss+=1; continue
                fa.write(f'>{k}\n{aa[k]}\n'); f3.write(f'>{k}\n{v.upper()}\n'); mp.write(f'{k}\t{sp}\n')
        print(R,len(sps),'species; len-mismatch/missing',nmiss)
    q=Q[(Q.panel==R)&Q.family_id.isin(g.family_id)&Q.has_3di]
    q=q.merge(g[['branch_child','family_id']],on=['branch_child','family_id'])
    q[['branch_child','family_id','gene','query']].to_csv(f'{W}/{R}.queries.tsv',sep='\t',index=False)
    open(f'{W}/{R}.qids','w').write('\n'.join(sorted(set(q.gene)))+'\n')
print(pick.groupby(['panel','is_event']).size())
