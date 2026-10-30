"""Dependence-aware sensitivity for confirmatory test R2-P1 (MUST FIX 3, round-2 final review).

The four confirmatory simple->complex events are not independent: three (Astrephomene, Volvox carteri,
V. africanus/reticuliferus) sit in panel R4 and are compared with the same six R4 control branches.
Here each panel contributes ONE effect:
  * R4 effect = mean residual (log rate ratio) of its 3 events. Permutation p: enumerate all C(6,3) = 20 ways of
    drawing 3 control branches without replacement from the same matched R4 control set; p = (1 + #{control mean >=
    observed mean}) / (N + 1), the same "+1" convention as the per-event test.
  * R6 effect = Fucophycidae residual and its per-event permutation p (5 matched controls).
  * Stouffer over the 2 panels; pooled RR with one effect per panel (REML tau^2, modified Hartung-Knapp, q*=max(1,q)).
Input: spec-literal (no branch-length floor) residuals from the independent verifier
  results/r2/verification/panels/verify_event_residuals.tsv, verify_control_loo_residuals.tsv
and, for comparison only, the same files under clip1e-4/ (producer's 1e-4 branch-length floor).
This is a sensitivity analysis; the confirmatory result (spec-literal Stouffer p 0.0504) is not changed.
"""
import itertools, json, math
import numpy as np, pandas as pd
from scipy import stats, optimize

P = '/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity'
OUT = P + '/results/r2/panels_sensitivity_onepanel'
R4_EV = ['R4_astrephomene', 'R4_volvox_carteri', 'R4_volvox_africanus']
R6_EV = 'R6_fucophycidae'


def depth_class(d): return 'A' if d <= 3 else ('B' if d <= 6 else 'C')


def matched(C, R, d):  # identical rule to scripts/r2/verification/panels/verify_stats.py
    pc = C[C.panel == R]; cls = depth_class(d)
    order = {'A': ['A', 'B', 'C'], 'B': ['B', 'A', 'C'], 'C': ['C', 'B', 'A']}[cls]
    sel = []
    for c in order:
        sel += list(pc.index[pc['depth'].map(depth_class) == c])
        if len(sel) >= 5: break
    return pc.loc[sel]


def meta_mhk(y, v):  # REML + modified Hartung-Knapp, as in verify_stats.reml
    y = np.asarray(y, float); v = np.asarray(v, float); k = len(y)
    def nll(t2):
        w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
        return 0.5 * (np.log(v + t2).sum() + np.log(w.sum()) + (w * (y - mu) ** 2).sum())
    r = optimize.minimize_scalar(nll, bounds=(0, 100), method='bounded', options={'xatol': 1e-10})
    t2 = r.x if nll(r.x) < nll(0) else 0.0
    w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
    q = (w * (y - mu) ** 2).sum() / (k - 1); se = math.sqrt(max(1, q) / w.sum()); t = stats.t.ppf(0.975, k - 1)
    return dict(k=k, tau2=float(t2), RR=math.exp(mu), lo=math.exp(mu - t * se), hi=math.exp(mu + t * se))


def run(sub):
    d = P + '/results/r2/verification/panels' + sub
    E = pd.read_csv(d + '/verify_event_residuals.tsv', sep='\t').set_index('event_id')
    C = pd.read_csv(d + '/verify_control_loo_residuals.tsv', sep='\t')
    r4 = E.loc[R4_EV]
    sets = [set(matched(C, 'R4', int(r4.loc[e, 'depth'])).index) for e in R4_EV]
    assert sets[0] == sets[1] == sets[2], 'R4 events do not share one matched control set'
    ctrl = C.loc[sorted(sets[0])]
    obs = float(r4.log_rr.mean())
    combs = [float(np.mean(c)) for c in itertools.combinations(ctrl.s.values, 3)]
    n_ge = int(sum(m >= obs - 1e-12 for m in combs))
    p_r4 = (1 + n_ge) / (len(combs) + 1)
    p_r4_exact = n_ge / len(combs)  # without the +1 (reported for transparency; 0 is not usable in Stouffer)
    e6 = E.loc[R6_EV]; p_r6 = float(e6.p_perm)
    z = [stats.norm.isf(p_r4), stats.norm.isf(p_r6)]
    Z = sum(z) / math.sqrt(2); p_st = float(stats.norm.sf(Z))
    # one effect per panel: R4 = mean of 3 log RRs; variance under two dependence extremes
    v = r4['var'].values
    var_indep = float(v.sum() / 9)                  # event errors independent
    var_corr = float(np.mean(np.sqrt(v)) ** 2)      # event errors perfectly correlated (conservative)
    w = 1 / v; iv_mu = float((w * r4.log_rr).sum() / w.sum()); iv_var = float(1 / w.sum())  # producer's convention
    pooled = {
        'mean_effect_var_perfectly_correlated': meta_mhk([obs, e6.log_rr], [var_corr, e6['var']]),
        'mean_effect_var_independent': meta_mhk([obs, e6.log_rr], [var_indep, e6['var']]),
        'iv_weighted_effect (pooled.py secondary convention)': meta_mhk([iv_mu, e6.log_rr], [iv_var, e6['var']]),
    }
    return dict(
        R4_events={e: dict(log_rr=float(r4.loc[e, 'log_rr']), rr=float(r4.loc[e, 'rr']), p_perm=float(r4.loc[e, 'p_perm'])) for e in R4_EV},
        R4_controls=dict(zip(ctrl.branch_child, ctrl.s.round(4))),
        R4_panel_mean_log_rr=obs, R4_panel_rr=math.exp(obs),
        R4_n_combinations=len(combs), R4_max_control_mean=max(combs), R4_n_combinations_ge_observed=n_ge,
        R4_p_perm=p_r4, R4_p_perm_without_plus1=p_r4_exact,
        R6_log_rr=float(e6.log_rr), R6_rr=float(e6.rr), R6_p_perm=p_r6, R6_n_matched=int(e6.n_matched),
        stouffer_2panel_Z=Z, stouffer_2panel_p=p_st,
        holm_step2_adjusted=min(1.0, 2 * p_st),
        pooled_one_effect_per_panel=pooled)


res = {'spec_literal_no_floor (primary)': run(''), 'producer_floor_1e-4 (comparison)': run('/clip1e-4')}
json.dump(res, open(OUT + '/onepanel_sensitivity.json', 'w'), indent=1)
L = []
for k, r in res.items():
    L.append(f'== {k} ==')
    for e, x in r['R4_events'].items():
        L.append(f'  {e:22s} RR {x["rr"]:6.2f}  per-event p {x["p_perm"]:.3f}')
    L.append(f'  R4 controls (LOO residuals): {r["R4_controls"]}')
    L.append(f'  R4 panel effect: mean log RR {r["R4_panel_mean_log_rr"]:.3f} (RR {r["R4_panel_rr"]:.2f}); '
             f'{r["R4_n_combinations_ge_observed"]} of {r["R4_n_combinations"]} control triples have mean >= observed '
             f'(max control-triple mean {r["R4_max_control_mean"]:.3f}); p = {r["R4_p_perm"]:.4f} (= 1/21)')
    L.append(f'  R6 Fucophycidae: RR {r["R6_rr"]:.2f}, p = {r["R6_p_perm"]:.4f} [{r["R6_n_matched"]} controls]')
    L.append(f'  Stouffer over 2 panels: Z = {r["stouffer_2panel_Z"]:.3f}, one-sided p = {r["stouffer_2panel_p"]:.4f}; '
             f'x2 (Holm second step) = {r["holm_step2_adjusted"]:.3f}')
    for nm, m in r['pooled_one_effect_per_panel'].items():
        L.append(f'  pooled RR, one effect per panel [{nm}]: {m["RR"]:.2f} (95% CI {m["lo"]:.2f}-{m["hi"]:.1f}), tau2 {m["tau2"]:.2f}, k=2')
    L.append('')
hdr = __doc__.strip().splitlines()
txt = '\n'.join(['ONE-EFFECT-PER-PANEL SENSITIVITY FOR CONFIRMATORY TEST R2-P1', 'Script: scripts/r2/panels_sensitivity_onepanel.py', ''] + hdr[2:] + [''] + L)
open(OUT + '/summary_raw.txt', 'w').write(txt)
print(txt)
