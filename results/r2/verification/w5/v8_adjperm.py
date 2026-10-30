"""Covariate-adjusted branch permutation: family logit (covariates + panel FE, no event term) -> branch mean (observed - expected) -> within-panel label permutation."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
W='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/verification/w5/work'
X=pd.read_csv(W+'/glmm_input.tsv',sep='\t')
X['llen']=np.log(X.len_m); X['lng']=np.log(X.n_genes); X['lprot']=np.log(X.n_out_prot)
res=[]
for name,form,sub in [('family covariates','new ~ C(panel) + llen + dis_m + lc_m + lng',X),
                      ('family+branch covariates','new ~ C(panel) + llen + dis_m + lc_m + lng + depth + lprot + frac_tx',X),
                      ('branch covariates only','new ~ C(panel) + depth + lprot + frac_tx',X)]:
    m=smf.logit(form,data=sub).fit(disp=0); X['res']=X.new-m.predict(X)
    B=X.groupby(['panel','branch_child']).agg(role=('role','first'),ekind=('ekind','first'),r=('res','mean')).reset_index()
    rng=np.random.default_rng(7)
    for kind in ['all','nested','origin']:
        Bk=B if kind=='all' else B[(B.role=='control')|(B.ekind==kind)]
        Bk=Bk.reset_index(drop=True); gs=[np.where(Bk.panel==p)[0] for p in Bk.panel.unique()]
        def st(roles):
            d=[]
            for gi in gs:
                r=Bk.r.values[gi]; ro=roles[gi]
                if (ro=='event').any() and (ro=='control').any(): d.append(r[ro=='event'].mean()-r[ro=='control'].mean())
            return np.mean(d)
        obs=st(Bk.role.values); s=np.empty(20000)
        for i in range(20000):
            ro=Bk.role.values.copy()
            for gi in gs: ro[gi]=rng.permutation(ro[gi])
            s[i]=st(ro)
        line=f'{name:28s} {kind:7s}: mean residual diff (event-control, panel-averaged) {obs:+.3f}; one-sided p={(1+(s>=obs).sum())/20001:.4f}'
        print(line); res.append(line)
open(W+'/adjperm_results.txt','w').write('\n'.join(res)+'\n')
