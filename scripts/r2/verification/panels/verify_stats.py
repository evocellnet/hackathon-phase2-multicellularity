"""Independent recomputation of R2-P1 statistics from per-branch tables (no reuse of scripts/r2/panels or mcell.branchstats)."""
import sys, json, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats, optimize
P = '/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
OUT = P + '/results/r2/verification/panels'
panels = ['R1','R2','R3','R4','R5','R6','R7','R8']
B = []; E = []
for R in panels:
    b = pd.read_csv(f'{P}/results/r2/panels/{R}/branch_table.tsv', sep='\t'); b['panel'] = R; B.append(b)
    e = pd.read_csv(f'{P}/results/r2/panels/{R}/event_table.tsv', sep='\t'); E.append(e)
B = pd.concat(B, ignore_index=True); E = pd.concat(E, ignore_index=True)
def tf(x): return x.astype(str) == 'True'
B['control'] = tf(B['control'])
covs = ['log_bl','log_n','frac_tx','mean_busco','log_dout']
B['log_bl'] = np.log(B['bl']); B['log_n'] = np.log(B['n_desc']); B['log_dout'] = np.log(B['dist_outgroup'])
C = B[B.control].copy().reset_index(drop=True)
print('controls by panel', C.groupby('panel').size().to_dict(), 'total', len(C))
mu_c = C[covs].mean(); sd_c = C[covs].std(ddof=0)
def design(df):
    X = pd.DataFrame({'const': 1.0}, index=df.index)
    for R in panels[1:]: X['panel_'+R] = (df['panel'] == R).astype(float)
    for c in covs:
        X[c] = (df[c] - mu_c[c]) / (sd_c[c] if sd_c[c] > 0 else 1.0)
    return X
from scipy.special import gammaln
def nb_nll(th, X, y):
    b = th[:-1]; a = np.exp(th[-1]); mu = np.exp(X @ b); r = 1/a
    return -(gammaln(y+r) - gammaln(r) - gammaln(y+1) + r*np.log(r/(r+mu)) + y*np.log(mu/(r+mu))).sum()
def num_hess(f, x, h=1e-4):
    n = len(x); H = np.zeros((n,n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ej = np.zeros(n); ei[i] = h; ej[j] = h
            H[i,j] = H[j,i] = (f(x+ei+ej) - f(x+ei-ej) - f(x-ei+ej) + f(x-ei-ej)) / (4*h*h)
    return H
def fit(df):
    X = design(df); cols = list(X.columns); X = X.values.astype(float); y = df['seq_novel'].astype(float).values
    p0 = sm.GLM(y, X, family=sm.families.Poisson()).fit().params
    th0 = np.r_[p0, np.log(0.5)]
    r = optimize.minimize(nb_nll, th0, args=(X, y), method='BFGS', options={'gtol':1e-8, 'maxiter':5000})
    r = optimize.minimize(nb_nll, r.x, args=(X, y), method='Nelder-Mead', options={'xatol':1e-9,'fatol':1e-10,'maxiter':200000,'maxfev':200000}) if not r.success else r
    H = num_hess(lambda t: nb_nll(t, X, y), r.x)
    cov = np.linalg.inv(H); k = X.shape[1]
    return r.x[:k], cov[:k,:k], float(np.exp(r.x[-1])), cols
beta, V, alpha, names = fit(C)
print('alpha', alpha); print(pd.Series(beta, index=names).round(4).to_dict())
# LOO control residuals
C['s'] = np.nan
for i in range(len(C)):
    b_, V_, a_, _ = fit(C.drop(i))
    x = design(C.iloc[[i]]).values[0]; yhat = np.exp(x @ b_)
    C.loc[i, 's'] = np.log((C.loc[i,'seq_novel'] + 0.5) / (yhat + 0.5))
# events: branch row with branch_child == M
def depth_class(d): return 'A' if d <= 3 else ('B' if d <= 6 else 'C')
def matched(R, d):
    pc = C[C.panel == R]; cls = depth_class(d)
    order = {'A': ['A','B','C'], 'B': ['B','A','C'], 'C': ['C','B','A']}[cls]
    sel = []
    for c in order:
        sel += list(pc.index[pc['depth'].map(depth_class) == c])
        if len(sel) >= 5: break
    return pc.loc[sel]
rows = []
for _, ev in E.iterrows():
    br = B[(B.panel == ev.panel) & (B.branch_child == ev.M)]
    assert len(br) == 1, ev.event_id
    br = br.iloc[0]
    if pd.isna(br['seq_novel']): continue
    x = design(br.to_frame().T.astype({c: float for c in covs})).values[0].astype(float)
    yhat = float(np.exp(x @ beta)); s = np.log((br.seq_novel + 0.5) / (yhat + 0.5))
    var = 1/(yhat+0.5) + alpha + x @ V @ x
    mc = matched(ev.panel, br.depth)
    p = (1 + (mc.s >= s).sum()) / (len(mc) + 1)
    rows.append(dict(panel=ev.panel, event_id=ev.event_id, kind=ev.kind, frm=ev.from_state, to=ev.to_state, M=ev.M,
                     tip=bool(br.is_leaf), depth=int(br.depth), seq_novel=br.seq_novel, yhat=yhat, log_rr=s, var=var,
                     rr=np.exp(s), lo=np.exp(s-1.96*np.sqrt(var)), hi=np.exp(s+1.96*np.sqrt(var)), n_matched=len(mc), p_perm=p))
R_ = pd.DataFrame(rows)
def reml(y, v):
    y = np.asarray(y); v = np.asarray(v); k = len(y)
    def nll(t2):
        w = 1/(v+t2); mu = (w*y).sum()/w.sum()
        return 0.5*(np.log(v+t2).sum() + np.log(w.sum()) + (w*(y-mu)**2).sum())
    r = optimize.minimize_scalar(nll, bounds=(0, 100), method='bounded', options={'xatol':1e-10})
    t2 = r.x if nll(r.x) < nll(0) else 0.0
    w = 1/(v+t2); mu = (w*y).sum()/w.sum()
    q = (w*(y-mu)**2).sum()/(k-1); se = np.sqrt(max(1,q)/w.sum()); t = stats.t.ppf(0.975, k-1)
    wf = 1/v; muf = (wf*y).sum()/wf.sum(); Q = (wf*(y-muf)**2).sum(); I2 = max(0,(Q-(k-1))/Q)
    out = dict(k=k, tau2=t2, mu=mu, q=q, RR=np.exp(mu), lo=np.exp(mu-t*se), hi=np.exp(mu+t*se), I2=I2)
    if k >= 3:
        tp = stats.t.ppf(0.975, k-2); h = tp*np.sqrt(t2+se**2); out.update(PI_lo=np.exp(mu-h), PI_hi=np.exp(mu+h))
    return out
p1 = ['R4_astrephomene','R4_volvox_carteri','R4_volvox_africanus','R6_fucophycidae']
S = R_.set_index('event_id').loc[p1]
Z = stats.norm.ppf(1 - S.p_perm).sum()/np.sqrt(4); pst = 1 - stats.norm.cdf(Z)
meta = reml(S.log_rr, S['var'])
sc = R_[(R_.frm=='simple')&(R_.to=='complex')&(~R_.tip)&(R_.panel.isin(['R2','R3','R4','R5','R6','R8']))&(R_.kind=='nested')]
meta_int = reml(sc.log_rr, sc['var'])
R_.to_csv(OUT+'/verify_event_residuals.tsv', sep='\t', index=False)
C[['panel','branch_child','depth','seq_novel','s']].to_csv(OUT+'/verify_control_loo_residuals.tsv', sep='\t', index=False)
res = dict(alpha=alpha, n_controls=len(C), p1=S[['rr','lo','hi','p_perm','n_matched']].to_dict('index'),
           stouffer_Z=Z, stouffer_p=pst, p1_meta=meta, internal_simple_complex=dict(meta_int, events=list(sc.event_id)))
json.dump(res, open(OUT+'/verify_stats.json','w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
