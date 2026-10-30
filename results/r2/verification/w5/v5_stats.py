"""Item 1/4: per-branch fractions, MH OR by panel, branch-level permutation, per-panel ORs, leave-one-panel-out."""
import pandas as pd, numpy as np
from scipy import stats
P='/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
W=P+'/results/r2/verification/w5/work'
F=pd.read_csv(W+'/families_verifier.tsv',sep='\t'); F=F[(F.kind=='novel')&(F.nhas3di>0)].copy()
ev=pd.concat([pd.read_csv(P+f'/results/r2/panels/R{k}/event_table.tsv',sep='\t') for k in range(1,9)]).set_index('event_id')['kind']
def ekind(s):
    if not isinstance(s,str): return None
    ks=[ev.get(e.strip()) for e in s.replace(';',',').split(',') if e.strip()]
    for k in ('nested','origin_complex','origin'):
        if k in ks: return 'nested' if k=='nested' else 'origin'
    return ks[0] if ks else None
F['ekind']=F.events.map(ekind)
F['role']=np.where(F.control,'control',np.where((~F.is_leaf)&F.ekind.isin(['origin','nested']),'event',np.where((~F.is_leaf)&F.ekind.notna(),'other_internal_event','other')))
F['new']=(~F.old).astype(int)
F.to_csv(W+'/families_roles.tsv',sep='\t',index=False)
B=F[F.role.isin(['event','control'])].groupby(['panel','branch_child','role','ekind'],dropna=False).agg(n=('new','size'),n_new=('new','sum'),
   depth=('depth','first'),frac_tx=('frac_tx','first'),n_out_sp=('n_out_sp','first'),n_out_prot=('n_out_prot','first'),n_desc=('n_desc','first')).reset_index()
B['frac_new']=B.n_new/B.n
B.to_csv(W+'/branch_level.tsv',sep='\t',index=False)
print(B.groupby(['panel','role']).size().unstack())
def mh(tabs):
    # tabs: list of (a,b,c,d) a=event new, b=event old, c=ctrl new, d=ctrl old ; RGB variance
    num=den=0; P_=Q_=R_=S_=0; Rs=Ss=0; PR=PS_QR=QS=0
    for a,b,c,d in tabs:
        n=a+b+c+d; R=a*d/n; S=b*c/n; Rs+=R; Ss+=S
        Pi=(a+d)/n; Qi=(b+c)/n
        PR+=Pi*R; PS_QR+=Pi*S+Qi*R; QS+=Qi*S
    OR=Rs/Ss; v=PR/(2*Rs**2)+PS_QR/(2*Rs*Ss)+QS/(2*Ss**2); se=np.sqrt(v)
    return OR,OR*np.exp(-1.96*se),OR*np.exp(1.96*se)
def tabs_of(B):
    t=[]
    for R,g in B.groupby('panel'):
        e=g[g.role=='event']; c=g[g.role=='control']
        if len(e)==0 or len(c)==0: continue
        t.append((e.n_new.sum(),e.n.sum()-e.n_new.sum(),c.n_new.sum(),c.n.sum()-c.n_new.sum()))
    return t
out=[]
def line(s): print(s); out.append(s)
e=B[B.role=='event']; c=B[B.role=='control']
line(f'pooled: events {e.n_new.sum()}/{e.n.sum()}={e.n_new.sum()/e.n.sum():.4f} ({len(e)} branches); controls {c.n_new.sum()}/{c.n.sum()}={c.n_new.sum()/c.n.sum():.4f} ({len(c)} branches)')
OR,lo,hi=mh(tabs_of(B)); line(f'MH OR {OR:.3f} [{lo:.3f},{hi:.3f}]')
line('per panel:')
for R,g in B.groupby('panel'):
    t=tabs_of(g)[0]; a,b,c_,d=t
    o=(a*d)/(b*c_); se=np.sqrt(1/a+1/b+1/c_+1/d)
    line(f'  {R} ev {a}/{a+b}={a/(a+b):.3f} ctl {c_}/{c_+d}={c_/(c_+d):.3f} OR {o:.2f} [{o*np.exp(-1.96*se):.2f},{o*np.exp(1.96*se):.2f}] n_ev_br={len(g[g.role=="event"])} n_ctl_br={len(g[g.role=="control"])}')
line('leave-one-panel-out MH OR:')
for R in sorted(B.panel.unique()):
    OR_,lo_,hi_=mh(tabs_of(B[B.panel!=R])); line(f'  without {R}: {OR_:.3f} [{lo_:.3f},{hi_:.3f}]')
# branch-level permutation within panel
rng=np.random.default_rng(1)
def stat_mh(Bp): return np.log(mh(tabs_of(Bp))[0])
def stat_br(Bp):  # mean over panels of (mean event branch frac - mean ctrl branch frac)
    d=[]
    for R,g in Bp.groupby('panel'):
        d.append(g[g.role=='event'].frac_new.mean()-g[g.role=='control'].frac_new.mean() if (g.role=='event').any() else np.nan)
    return np.nanmean(d)
obs1,obs2=stat_mh(B),stat_br(B)
NP=20000; s1=np.empty(NP); s2=np.empty(NP)
Bp=B.copy()
groups=[np.where(B.panel==R)[0] for R in sorted(B.panel.unique())]
roles=B.role.values
for i in range(NP):
    r=roles.copy()
    for gi in groups: r[gi]=rng.permutation(r[gi])
    Bp['role']=r; s1[i]=stat_mh(Bp); s2[i]=stat_br(Bp)
line(f'branch permutation (within panel, {NP}x): MH logOR stat obs {obs1:.3f} p_one-sided={(1+(s1>=obs1).sum())/(NP+1):.4f} two-sided={(1+(np.abs(s1)>=abs(obs1)).sum())/(NP+1):.4f}')
line(f'   equal-branch-weight mean diff obs {obs2:.4f} p_one-sided={(1+(s2>=obs2).sum())/(NP+1):.4f} two-sided={(1+(np.abs(s2)>=abs(obs2)).sum())/(NP+1):.4f}')
# by kind
for k in ['origin','nested']:
    Bk=B[(B.role=='control')|(B.ekind==k)]
    OR_,lo_,hi_=mh(tabs_of(Bk)); ee=Bk[Bk.role=='event']
    line(f'kind {k}: {len(ee)} events, {ee.n_new.sum()}/{ee.n.sum()}={ee.n_new.sum()/ee.n.sum():.3f}; MH OR {OR_:.3f} [{lo_:.3f},{hi_:.3f}]')
    # permutation for kind
    Bk=Bk.reset_index(drop=True); obs=stat_br(Bk); gs=[np.where(Bk.panel==R)[0] for R in Bk.panel.unique()]
    rr=Bk.role.values; s=np.empty(5000); Bq=Bk.copy()
    for i in range(5000):
        r=rr.copy()
        for gi in gs: r[gi]=rng.permutation(r[gi])
        Bq['role']=r; s[i]=stat_br(Bq)
    line(f'   branch-perm equal-weight diff {obs:.3f}, one-sided p={(1+(s>=obs).sum())/5001:.4f}')
line('event branches:')
line(B[B.role=='event'][['panel','branch_child','ekind','n','n_new','frac_new','depth','n_out_sp']].to_string(index=False))
open(W+'/stats_item1_4.txt','w').write('\n'.join(out)+'\n')
