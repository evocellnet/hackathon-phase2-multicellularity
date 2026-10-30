"""Re-derive event vs control with exhaustive-search definition of 'structurally new'."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'; W=P+'/results/r2/verification/w5/work'
G=pd.read_csv(W+'/glmm_input.tsv',sep='\t'); D=pd.read_csv(W+'/exh/rescued.tsv',sep='\t')
G=G.merge(D[['panel','branch_child','family_id','rescued']],on=['panel','branch_child','family_id'],how='left')
assert G[G.new==1].rescued.notna().all()
G['new_exh']=((G.new==1)&(G.rescued!=True)).astype(int)
G.to_csv(W+'/glmm_input_exh.tsv',sep='\t',index=False)
out=[]
def p(s): print(s); out.append(s)
p('Rescue of producer-"new" families by exhaustive search, by panel and role (fraction rescued):')
p(G[G.new==1].groupby(['panel','role']).rescued.mean().unstack().round(3).to_string())
r=G[G.new==1].groupby('role').rescued.agg(['sum','size','mean']); p(r.to_string())
def mh(G,col):
    Rs=Ss=PR=PSQR=QS=0
    for R,g in G.groupby('panel'):
        e=g[g.role=='event'][col]; c=g[g.role=='control'][col]
        if len(e)==0: continue
        a,b,cc,d=e.sum(),len(e)-e.sum(),c.sum(),len(c)-c.sum(); n=a+b+cc+d
        R_=a*d/n; S_=b*cc/n; Rs+=R_; Ss+=S_; Pi=(a+d)/n; Qi=(b+cc)/n; PR+=Pi*R_; PSQR+=Pi*S_+Qi*R_; QS+=Qi*S_
    OR=Rs/Ss; se=np.sqrt(PR/(2*Rs**2)+PSQR/(2*Rs*Ss)+QS/(2*Ss**2)); return OR,OR*np.exp(-1.96*se),OR*np.exp(1.96*se)
for col in ['new','new_exh']:
    e=G[G.role=='event'][col]; c=G[G.role=='control'][col]
    p(f'{col}: event {e.sum()}/{len(e)}={e.mean():.3f} control {c.sum()}/{len(c)}={c.mean():.3f}; MH OR %.2f [%.2f, %.2f]'%mh(G,col))
    for k in ['nested','origin']:
        Gk=G[(G.role=='control')|(G.ekind==k)]; p(f'   {k}: MH OR %.2f [%.2f, %.2f]'%mh(Gk,col))
    p('   per panel OR: '+'; '.join(f'{R} %.2f'%mh(g,col)[0] for R,g in G.groupby('panel')))
    # branch permutation equal weight
    B=G.groupby(['panel','branch_child']).agg(role=('role','first'),ekind=('ekind','first'),f=(col,'mean')).reset_index()
    rng=np.random.default_rng(3)
    for kind in ['all','nested','origin']:
        Bk=(B if kind=='all' else B[(B.role=='control')|(B.ekind==kind)]).reset_index(drop=True)
        gs=[np.where(Bk.panel==q)[0] for q in Bk.panel.unique()]
        def st(ro):
            d=[Bk.f.values[gi][ro[gi]=='event'].mean()-Bk.f.values[gi][ro[gi]=='control'].mean() for gi in gs if (ro[gi]=='event').any()]
            return np.mean(d)
        obs=st(Bk.role.values); s=np.array([st(np.concatenate([rng.permutation(Bk.role.values[gi]) for gi in gs])[np.argsort(np.concatenate(gs))]) for _ in range(20000)])
        p(f'   branch-perm {kind}: diff {obs:+.3f} one-sided p={(1+(s>=obs).sum())/20001:.4f}')
open(W+'/exh_stats.txt','w').write('\n'.join(out)+'\n')
