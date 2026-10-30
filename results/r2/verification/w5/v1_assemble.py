"""Independent verifier: assemble per-family table from producer's hit tables (own code)."""
import pandas as pd, numpy as np, re, json, gzip, os, sys
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
OUT=P+'/results/r2/verification/w5/work'
def parse_newick(s):
    s=s.strip().rstrip(';'); i=0; stack=[]; children={}; cnt=[0]
    def node():
        nonlocal i
        kids=[]
        if s[i]=='(':
            i+=1
            while True:
                kids.append(node())
                if s[i]==',': i+=1; continue
                if s[i]==')': i+=1; break
        m=re.match(r"([^:,();]*)(:[^,();]*)?",s[i:]); name=m.group(1); i+=m.end()
        if not name: cnt[0]+=1; name=f'_n{cnt[0]}'
        children[name]=kids; return name
    root=node(); return root,children
def leaves(n,ch,memo):
    if n in memo: return memo[n]
    r={n} if not ch[n] else set().union(*[leaves(c,ch,memo) for c in ch[n]])
    memo[n]=r; return r
# disorder lookup
def load_dis(sp):
    for d in ('euk','prok'):
        f=f'{P}/data/annotations/disorder/{d}/{sp}.disorder.tsv.gz'
        if os.path.exists(f): return pd.read_csv(f,sep='\t').set_index('protein')
    return None
src=pd.read_csv(P+'/data/structure3di/panel_sources.tsv',sep='\t',header=None,names=['sp','fa']).set_index('sp')['fa'].to_dict()
def lc_entropy(seq,w=12,thr=2.2):
    # fraction of residues covered by a 12-aa window with Shannon entropy < thr bits (SEG-like proxy)
    n=len(seq)
    if n<w: return 0.0
    mask=np.zeros(n,bool)
    from collections import Counter
    for k in range(n-w+1):
        c=Counter(seq[k:k+w]); p=np.array(list(c.values()))/w
        if -(p*np.log2(p)).sum()<thr: mask[k:k+w]=True
    return mask.mean()
allfam=[]; allq=[]
for R in [f'R{k}' for k in range(1,9)]:
    nd=P+f'/results/r2/w5/novelty/{R}'
    root,ch=parse_newick(open(P+f'/results/r2/panels/{R}/tree_bl.nwk').read()); memo={}
    allsp=leaves(root,ch,memo)
    bt=pd.read_csv(P+f'/results/r2/panels/{R}/branch_table.tsv',sep='\t').set_index('branch_child')
    s=pd.read_csv(nd+'/sample.tsv',sep='\t')
    h=pd.read_csv(nd+'/hits.tsv.gz',sep='\t')
    h['tsp']=h.target.str.split('_').str[0]
    h['qcov']=(h.qend-h.qstart+1)/h.qlen
    # target proteome sizes
    tg=pd.read_parquet(P+f'/results/r2/panels/{R}/work/ladder/db.targets.parquet',columns=['seqid','species','source','length'])
    panel_tg=tg[tg.source=='panel']; nprot=panel_tg.groupby('species').size()
    glen=tg.set_index('seqid')['length']
    gains=pd.read_parquet(P+f'/results/r2/panels/{R}/work/ladder/gains.parquet')
    gains=gains.set_index(['branch_child','family_id'])[['n_genes','n_species']]
    # check hits targets outside clade
    viol=0
    for b,hb in h.groupby('branch_child'):
        cl=leaves(b,ch,memo); viol+=hb.tsp.isin(cl).sum()
    qok=h[(h.evalue<1e-3)&(h.qcov>=0.5)]
    qold=set(zip(qok.branch_child,qok['query']))
    s['qold']=[(b,q) in qold for b,q in zip(s.branch_child,s['query'])]
    s['qhit_any']=s['query'].isin(set(h['query']))
    s['len']=s.gene.map(glen)
    # disorder
    dis=[];
    cache={}
    for sp,g in zip(s.species,s.gene):
        if sp not in cache: cache[sp]=load_dis(sp)
        d=cache[sp]; dis.append(d.loc[g,'frac_disordered'] if d is not None and g in d.index else np.nan)
    s['dis']=dis
    s['panel']=R
    allq.append(s)
    f=s.groupby(['branch_child','kind','family_id']).agg(n_reps=('gene','size'),old=('qold','any'),
        n_available=('n_available','first'),mlen=('len','mean'),dis=('dis','mean'),nhas3di=('has_3di','sum')).reset_index()
    f['panel']=R
    f=f.join(gains,on=['branch_child','family_id'])
    for col in ['depth','frac_tx','n_desc','is_leaf','is_event','control','events','bl','root_dist','inside_crown']:
        f[col]=f.branch_child.map(bt[col])
    f['n_out_sp']=[len(allsp-leaves(b,ch,memo)) for b in f.branch_child]
    f['n_out_prot']=[nprot[list(allsp-leaves(b,ch,memo))].sum() for b in f.branch_child]
    # compare with producer families.tsv
    pf=pd.read_csv(nd+'/families.tsv',sep='\t')
    m=f.merge(pf[['branch_child','kind','family_id','struct_old']],on=['branch_child','kind','family_id'],how='outer',indicator=True)
    agree=(m._merge=='both')&(m.old==m.struct_old)
    print(R,'fams',len(f),'producer',len(pf),'merge',m._merge.value_counts().to_dict(),'agree',agree.sum(),'target-in-clade violations',viol,
          'missing len',f.mlen.isna().sum(),'missing dis',f.dis.isna().sum(),'no3di reps',(s.has_3di==False).sum(),flush=True)
    allfam.append(f)
F=pd.concat(allfam); Q=pd.concat(allq)
F.to_csv(OUT+'/families_verifier.tsv',sep='\t',index=False); Q.to_csv(OUT+'/queries_verifier.tsv.gz',sep='\t',index=False)
