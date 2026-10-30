import pandas as pd, numpy as np, re
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
W=P+'/results/r2/verification/w5/work/spot'
pick=pd.read_csv(W+'/picked.tsv',sep='\t')
def parse(s):
    s=s.strip().rstrip(';'); i=0; ch={}; c=[0]
    def node():
        nonlocal i; kids=[]
        if s[i]=='(':
            i+=1
            while True:
                kids.append(node())
                if s[i]==',': i+=1; continue
                if s[i]==')': i+=1; break
        m=re.match(r"([^:,();]*)(:[^,();]*)?",s[i:]); n=m.group(1); i+=m.end()
        if not n: c[0]+=1; n=f'_n{c[0]}'
        ch[n]=kids; return n
    r=node(); return r,ch
def lv(n,ch): return {n} if not ch[n] else set().union(*[lv(k,ch) for k in ch[n]])
rows=[]; hitcmp=[]
for R,g in pick.groupby('panel'):
    root,ch=parse(open(P+f'/results/r2/panels/{R}/tree_bl.nwk').read())
    q=pd.read_csv(f'{W}/{R}.queries.tsv',sep='\t')
    m=pd.read_csv(f'{W}/{R}.m8',sep='\t',header=None,names=['query','target','evalue','bits','qstart','qend','qlen','tstart','tend','tlen'])
    m['tsp']=m.target.str.split('_').str[0]; m['qcov']=(m.qend-m.qstart+1)/m.qlen
    ph=pd.read_csv(P+f'/results/r2/w5/novelty/{R}/hits.tsv.gz',sep='\t')
    for _,f in g.iterrows():
        cl=lv(f.branch_child,ch); genes=set(q[(q.family_id==f.family_id)&(q.branch_child==f.branch_child)].gene)
        mm=m[m['query'].isin(genes)&~m.tsp.isin(cl)]
        ok=mm[(mm.evalue<1e-3)&(mm.qcov>=0.5)]
        pp=ph[(ph.branch_child==f.branch_child)&ph['query'].str.replace('Q|','',regex=False).isin(genes)]
        pok=pp[(pp.evalue<1e-3)&((pp.qend-pp.qstart+1)/pp.qlen>=0.5)]
        a=set(zip(ok['query'],ok.target)); b=set(zip(pok['query'].str.replace('Q|','',regex=False),pok.target))
        rows.append(dict(panel=R,branch=f.branch_child,family=f.family_id,role='event' if f.is_event else 'control',n_reps=len(genes),
            mine_old=len(ok)>0,prod_old=bool(f.old),n_hits_mine=len(a),n_hits_prod=len(b),jacc=len(a&b)/len(a|b) if a|b else 1.0))
D=pd.DataFrame(rows); D.to_csv(W+'/spotcheck_result.tsv',sep='\t',index=False)
print(pd.crosstab(D.mine_old,D.prod_old,rownames=['verifier_old'],colnames=['producer_old']))
print('agreement',(D.mine_old==D.prod_old).mean(),'n',len(D)); print('median jaccard of passing hit pairs',D.jacc.median(),'mean',D.jacc.mean())
print(D[D.mine_old!=D.prod_old].to_string())
