import pandas as pd, numpy as np, re, sys, os
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
E=P+'/results/r2/verification/w5/work/exh'
exec(open(P+'/results/r2/verification/w5/v2_spotcheck_compare.py').read().split('rows=[]')[0].split("pick=pd.read_csv")[0])
exec("def parse"+open(P+'/results/r2/verification/w5/v2_spotcheck_compare.py').read().split('def parse')[1].split('rows=[]')[0])
G=pd.read_csv(P+'/results/r2/verification/w5/work/glmm_input.tsv',sep='\t')
out=[]
for R in [f'R{k}' for k in range(1,9)]:
    if not os.path.exists(f'{E}/{R}.m8'): continue
    root,ch=parse(open(P+f'/results/r2/panels/{R}/tree_bl.nwk').read())
    q=pd.read_csv(f'{E}/{R}.queries.tsv',sep='\t')
    m=pd.read_csv(f'{E}/{R}.m8',sep='\t',header=None,names=['query','target','evalue','bits','qstart','qend','qlen','tstart','tend','tlen'])
    m=m[(m.evalue<1e-3)&((m.qend-m.qstart+1)/m.qlen>=0.5)]; m['tsp']=m.target.str.split('_').str[0]
    mq=m.merge(q.rename(columns={'gene':'query'}),on='query')
    resc=set()
    for b,g in mq.groupby('branch_child'):
        cl=lv(b,ch); gg=g[~g.tsp.isin(cl)]; resc|=set(zip([b]*len(gg),gg.family_id))
    fam=q.drop_duplicates(['branch_child','family_id'])[['branch_child','family_id','role']]
    fam['rescued']=[(b,f) in resc for b,f in zip(fam.branch_child,fam.family_id)]; fam['panel']=R
    out.append(fam)
D=pd.concat(out); D.to_csv(E+'/rescued.tsv',sep='\t',index=False)
print(D.groupby(['panel','role']).rescued.agg(['sum','size','mean']))
# recompute event vs control with exhaustive definition for panels done
G=G[G.panel.isin(D.panel.unique())].merge(D[['panel','branch_child','family_id','rescued']],on=['panel','branch_child','family_id'],how='left')
G['new_exh']=((G.new==1)&(G.rescued!=True)).astype(int)
for col in ['new','new_exh']:
    t=G.groupby('role')[col].agg(['sum','size','mean']); print(col); print(t)
