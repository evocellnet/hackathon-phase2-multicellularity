import pandas as pd, numpy as np
from scipy import stats
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
F=pd.read_csv(P+'/results/r2/verification/w5/work/families_verifier.tsv',sep='\t')
F=F[F.kind=='novel']
rows=[]
for R in sorted(F.panel.unique()):
    c=pd.read_parquet(P+f'/results/r2/panels/{R}/work/ladder/classified.parquet')
    nov=c[c.seq_novel==True]
    for b,g in F[F.panel==R].groupby('branch_child'):
        avail=nov[nov.branch_child==b].sort_values('family_id')
        ids=list(avail.family_id); pos={k:i for i,k in enumerate(ids)}
        samp=set(g.family_id)
        ranks=np.array([pos[k] for k in samp if k in pos])
        n=len(ids)
        # compare sampled vs unsampled family size (n_genes)
        av=avail.assign(s=avail.family_id.isin(samp))
        rows.append(dict(panel=R,branch=b,event=bool(g.is_event.iloc[0]) and not bool(g.is_leaf.iloc[0]),control=bool(g.control.iloc[0]),
            leaf=bool(g.is_leaf.iloc[0]),n_avail_ladder=n,n_avail_prod=int(g.n_available.iloc[0]),n_sampled=len(samp),n_sampled_found=len(ranks),
            capped=n>200,first200=(ranks.max()<len(samp)) if len(ranks) else None,
            mean_rank_frac=ranks.mean()/max(n-1,1) if len(ranks) else None,
            ngenes_samp=av[av.s].n_genes.mean(), ngenes_unsamp=av[~av.s].n_genes.mean() if (~av.s).any() else np.nan))
D=pd.DataFrame(rows)
D.to_csv(P+'/results/r2/verification/w5/work/sampling_check.tsv',sep='\t',index=False)
print('n_avail mismatches', (D.n_avail_ladder!=D.n_avail_prod).sum(), 'sampled-not-in-ladder', (D.n_sampled!=D.n_sampled_found).sum())
X=D[(D.event|D.control)]
print(X.groupby('event').agg(n=('branch','size'),capped=('capped','sum'),avail_med=('n_avail_ladder','median'),
    first200=('first200','sum'),mean_rank=('mean_rank_frac','mean')))
cap=X[X.capped]
print('capped branches: mean rank frac (0.5 expected if random)',cap.mean_rank_frac.describe().to_string())
print('capped: sampled vs unsampled n_genes', cap.ngenes_samp.mean(), cap.ngenes_unsamp.mean())
print(X[['panel','branch','event','n_avail_ladder','n_sampled','mean_rank_frac']].to_string())
