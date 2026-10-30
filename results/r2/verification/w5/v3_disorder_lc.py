"""metapredict for query reps lacking precomputed disorder; SEG-like low-complexity fraction for all query reps."""
import pandas as pd, numpy as np, os
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
W=P+'/results/r2/verification/w5/work'
Q=pd.read_csv(W+'/queries_verifier.tsv.gz',sep='\t')
src=pd.read_csv(P+'/data/structure3di/panel_sources.tsv',sep='\t',header=None,names=['sp','fa']).set_index('sp')['fa'].to_dict()
need=set(Q.gene)
seqs={}
for sp in Q.species.unique():
    fn=src.get(sp)
    if fn is None:
        for d in ('dense','euk','prok'):
            c=f'{P}/data/proteomes/{d}/{sp}.fa'
            if os.path.exists(c): fn=c
    k=None
    for l in open(fn):
        if l[0]=='>': k=l[1:].split()[0]; keep=k in need; 
            
        elif keep: seqs[k]=seqs.get(k,'')+l.strip()
print('seqs',len(seqs),'of',len(need),flush=True)
AA='ACDEFGHIKLMNPQRSTVWY'
idx={a:i for i,a in enumerate(AA)}
def lcfrac(s,w=12,thr=2.2):
    n=len(s)
    if n<w: return 0.0
    x=np.array([idx.get(c,20) for c in s]); oh=np.zeros((n+1,21)); oh[np.arange(1,n+1),x]=1; cs=oh.cumsum(0)
    cnt=cs[w:]-cs[:-w]; p=cnt/w
    with np.errstate(divide='ignore',invalid='ignore'): H=-np.nansum(np.where(p>0,p*np.log2(p),0),1)
    low=np.where(H<thr)[0]; mask=np.zeros(n,bool)
    for k in low: mask[k:k+w]=True
    return mask.mean()
g=sorted(need)
lc=pd.Series({k:lcfrac(seqs[k]) for k in g if k in seqs},name='lc')
print('lc done',flush=True)
miss=g
import torch; torch.set_num_threads(12)
import metapredict as meta
dd={}
import re
batch={k:re.sub('[^ACDEFGHIKLMNPQRSTVWY]','G',seqs[k].upper().rstrip('*')) for k in miss if k in seqs}
out=meta.predict_disorder(batch)
for k,v in out.items():
    sc=np.asarray(v[1] if isinstance(v,(tuple,list)) else v); dd[k]=(sc>=0.5).mean()
print('metapredict done',len(dd),flush=True)
pd.DataFrame({'gene':g}).assign(lc=lambda d:d.gene.map(lc),dis_new=lambda d:d.gene.map(dd),
    seq_len=lambda d:d.gene.map(lambda k:len(seqs.get(k,'')))).to_csv(W+'/query_lc_dis.tsv.gz',sep='\t',index=False)
