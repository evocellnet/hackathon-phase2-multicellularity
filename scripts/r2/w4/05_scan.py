"""Task C: per-family scans, multicellular clade vs sister clade, for each primary / secondary eukaryotic origin.

Model per family f and event e: count_sf ~ NB(mu_sf, alpha_f), log mu_sf = c_{f,group(s)} + log(proteome size_s)
(two-group GLM with offset; the two group log-rates are separable). alpha_f by Cox-Reid adjusted profile likelihood
on a 41-point log grid (1e-4 .. 1e2). Test: LRT of one shared rate vs two rates at alpha_f -> chi2(1); BH within event.
Families: Pfam families present in >= 2 species of the event, transposon mask (data/r2_sets/transposon_mask.tsv) removed.
Log fold change for ranking / ashr: log rate ratio with 0.5 pseudocount on group totals, SE from NB information.
Recurrence: number of primary origins (7) with q < 0.1 in the same direction; secondary events counted separately.
Null: per-event permutation of mc/sister labels among the event's species (group sizes kept), same permutation for all
families in an event, events permuted independently; B = 1000, seed 11.
Outputs: results/r2/w4/scan/event_<id>.tsv, recurrence_observed.tsv, recurrence_null_summary.tsv, perm_recurrence.npz
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(__file__))
from common import P, OUT, out, load_transitions, SEED  # noqa: E402

ALPHAS = np.logspace(-4, 2, 41)
NEWTON = 8
B = 1000
QCUT = 0.1


def nb_ll(y, mu, a):
    r = 1.0 / a
    with np.errstate(divide="ignore", invalid="ignore"):
        t = gammaln(y + r) - gammaln(r) - gammaln(y + 1) + y * np.log(a * mu / (1 + a * mu)) - r * np.log1p(a * mu)
    return np.where(mu > 0, t, np.where(y == 0, 0.0, -np.inf))


def group_fit(Y, N, a):
    """Y: F x m counts of one group, N: m sizes. Returns c (F,), mu (F x m). c=-inf if all zero."""
    tot = Y.sum(1)
    pos = tot > 0
    c = np.full(Y.shape[0], -np.inf)
    c[pos] = np.log(tot[pos] / N.sum())
    cp = c[pos]
    Yp = Y[pos]
    for _ in range(NEWTON):
        mu = np.exp(cp)[:, None] * N[None, :]
        w = mu / (1 + a * mu)
        score = ((Yp - mu) / (1 + a * mu)).sum(1)
        info = w.sum(1)
        cp = cp + score / info
    c[pos] = cp
    mu = np.where(np.isfinite(c)[:, None], np.exp(np.where(np.isfinite(c), c, 0))[:, None] * N[None, :], 0.0)
    return c, mu


def fit_event(Y, N, g, want_lfc=False):
    """Y: F x n, N: n, g: n bool (True = multicellular clade). Returns dict of arrays."""
    Y1, Y0, N1, N0 = Y[:, g], Y[:, ~g], N[g], N[~g]
    best = np.full(Y.shape[0], -np.inf)
    a_hat = np.full(Y.shape[0], np.nan)
    for a in ALPHAS:
        c1, mu1 = group_fit(Y1, N1, a)
        c0, mu0 = group_fit(Y0, N0, a)
        ll = nb_ll(Y1, mu1, a).sum(1) + nb_ll(Y0, mu0, a).sum(1)
        i1 = (mu1 / (1 + a * mu1)).sum(1)
        i0 = (mu0 / (1 + a * mu0)).sum(1)
        cr = -0.5 * (np.where(i1 > 0, np.log(np.where(i1 > 0, i1, 1)), 0) + np.where(i0 > 0, np.log(np.where(i0 > 0, i0, 1)), 0))
        apl = ll + cr
        upd = apl > best
        best[upd] = apl[upd]
        a_hat[upd] = a
    # full and null fits at a_hat (per family: group families by alpha value)
    llf = np.zeros(Y.shape[0]); lln = np.zeros(Y.shape[0])
    c1v = np.zeros(Y.shape[0]); c0v = np.zeros(Y.shape[0])
    for a in np.unique(a_hat):
        idx = np.where(a_hat == a)[0]
        c1, mu1 = group_fit(Y1[idx], N1, a)
        c0, mu0 = group_fit(Y0[idx], N0, a)
        llf[idx] = nb_ll(Y1[idx], mu1, a).sum(1) + nb_ll(Y0[idx], mu0, a).sum(1)
        cn, mun = group_fit(Y[idx], N, a)
        lln[idx] = nb_ll(Y[idx], mun, a).sum(1)
        c1v[idx] = c1; c0v[idx] = c0
    lr = np.clip(2 * (llf - lln), 0, None)
    p = chi2.sf(lr, 1)
    tot1, tot0 = Y1.sum(1), Y0.sum(1)
    sign = np.where(tot1 == 0, -1, np.where(tot0 == 0, 1, np.sign(c1v - c0v)))
    res = dict(alpha=a_hat, lr=lr, p=p, sign=sign)
    if want_lfc:
        r1 = (tot1 + 0.5) / N1.sum(); r0 = (tot0 + 0.5) / N0.sum()
        mu1 = r1[:, None] * N1[None, :]; mu0 = r0[:, None] * N0[None, :]
        i1 = (mu1 / (1 + a_hat[:, None] * mu1)).sum(1); i0 = (mu0 / (1 + a_hat[:, None] * mu0)).sum(1)
        res.update(lfc=np.log(r1) - np.log(r0), se=np.sqrt(1 / i1 + 1 / i0), lfc_mle=c1v - c0v,
                   mean_mc=Y1.mean(1), mean_sis=Y0.mean(1), prev_mc=(Y1 > 0).mean(1), prev_sis=(Y0 > 0).mean(1))
    return res


def bh(p):
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def sig_vectors(Y, N, g):
    r = fit_event(Y, N, g)
    q = bh(r["p"])
    return (q < QCUT) & (r["sign"] > 0), (q < QCUT) & (r["sign"] < 0)


def perm_worker(args):
    Y, N, perms = args
    n1 = perms.shape[1]
    E = np.zeros((len(perms), Y.shape[0]), dtype=bool); C = np.zeros_like(E)
    for k, idx in enumerate(perms):
        g = np.zeros(len(N), dtype=bool); g[idx] = True
        E[k], C[k] = sig_vectors(Y, N, g)
    return E, C


def main():
    counts = pd.read_csv(f"{OUT}/matrix/counts.tsv.gz", sep="\t", index_col=0)
    sp = pd.read_csv(f"{OUT}/matrix/species.tsv", sep="\t").set_index("code")
    mask = set(pd.read_csv(f"{P}/data/r2_sets/transposon_mask.tsv", sep="\t").pfam)
    tr = load_transitions()
    ev = tr[(tr.primary_eligible == 1) | (tr.analysis_set == "secondary")].copy()
    ev["tier"] = np.where(ev.primary_eligible == 1, "primary", "secondary")
    fams_all = [f for f in counts.columns if f not in mask]
    fidx = {f: i for i, f in enumerate(fams_all)}
    F = len(fams_all)
    print(f"{F} families after mask ({len(mask)} masked); events:", ev.transition_id.tolist(), flush=True)
    f_obs = out("scan/recurrence_observed.tsv"); f_null = out("scan/recurrence_null_summary.tsv"); f_npz = out("scan/perm_recurrence.npz")
    f_events = out("scan/events_used.tsv")
    obs_E = {}; obs_C = {}; perm_E = {}; perm_C = {}; tested = {}
    rows_ev = []
    for ei, (_, e) in enumerate(ev.iterrows()):
        codes = [c for c in e.mc_list + e.sis_list if c in counts.index]
        g = np.array([c in e.mc_list for c in codes])
        sub = counts.loc[codes, fams_all].to_numpy(dtype=float).T
        keep = (sub > 0).sum(1) >= 2
        Y = sub[keep]; N = sp.loc[codes, "n_proteins"].to_numpy(dtype=float)
        fam_e = [f for f, k in zip(fams_all, keep) if k]
        r = fit_event(Y, N, g, want_lfc=True)
        q = bh(r["p"])
        df = pd.DataFrame(dict(pfam=fam_e, n_mc=g.sum(), n_sister=(~g).sum(), mean_count_mc=r["mean_mc"], mean_count_sister=r["mean_sis"],
                               prev_mc=r["prev_mc"], prev_sister=r["prev_sis"], lfc=r["lfc"], se=r["se"], lfc_mle=r["lfc_mle"],
                               alpha=r["alpha"], lr=r["lr"], p=r["p"], q_bh=q, direction=r["sign"]))
        df.to_csv(out(f"scan/event_{e.transition_id}.tsv"), sep="\t", index=False)
        ii = np.array([fidx[f] for f in fam_e])
        tested[e.transition_id] = ii
        oe = np.zeros(F, bool); oc = np.zeros(F, bool)
        oe[ii] = (q < QCUT) & (r["sign"] > 0); oc[ii] = (q < QCUT) & (r["sign"] < 0)
        obs_E[e.transition_id], obs_C[e.transition_id] = oe, oc
        tr_t = transcriptome = int(sum(sp.loc[c, "data_type"] == "transcriptome" for c in codes))
        rows_ev.append(dict(event=e.transition_id, tier=e.tier, label=e.label, n_mc=int(g.sum()), n_sister=int((~g).sum()),
                            n_transcriptome=tr_t, mc_species=",".join(c for c in codes if c in e.mc_list),
                            sister_species=",".join(c for c in codes if c not in e.mc_list), n_families_tested=len(fam_e),
                            n_sig_expanded=int(oe.sum()), n_sig_contracted=int(oc.sum())))
        print(rows_ev[-1], flush=True)
        # permutations
        rng = np.random.default_rng([SEED, ei])
        n = len(codes); n1 = int(g.sum())
        perms = np.array([rng.permutation(n)[:n1] for _ in range(B)])
        chunks = np.array_split(perms, 8)
        with ProcessPoolExecutor(8) as ex:
            parts = list(ex.map(perm_worker, [(Y, N, ch) for ch in chunks]))
        PE = np.concatenate([p[0] for p in parts]); PC = np.concatenate([p[1] for p in parts])
        fullE = np.zeros((B, F), bool); fullC = np.zeros((B, F), bool)
        fullE[:, ii] = PE; fullC[:, ii] = PC
        perm_E[e.transition_id], perm_C[e.transition_id] = fullE, fullC
        print(e.transition_id, "perm mean sig exp/con", PE.sum(1).mean(), PC.sum(1).mean(), flush=True)
    pd.DataFrame(rows_ev).to_csv(f_events, sep="\t", index=False)
    prim = ev.transition_id[ev.tier == "primary"].tolist(); sec = ev.transition_id[ev.tier == "secondary"].tolist()
    rec = {}
    for nm, ids in (("primary", prim), ("secondary", sec)):
        for d, O, Pm in (("expanded", obs_E, perm_E), ("contracted", obs_C, perm_C)):
            rec[(nm, d, "obs")] = sum(O[i].astype(int) for i in ids)
            rec[(nm, d, "perm")] = sum(Pm[i].astype(np.int8) for i in ids)
    # per-family table
    fam_df = pd.DataFrame(dict(pfam=fams_all))
    for nm in ("primary", "secondary"):
        for d in ("expanded", "contracted"):
            o = rec[(nm, d, "obs")]; pm = rec[(nm, d, "perm")]
            fam_df[f"rec_{nm}_{d}"] = o
            fam_df[f"perm_p_{nm}_{d}"] = ((pm >= o[None, :]).sum(0) + 1) / (B + 1)
    for i in prim + sec:
        fam_df[f"tested_{i}"] = 0
        fam_df.loc[tested[i], f"tested_{i}"] = 1
        fam_df[f"sig_{i}"] = obs_E[i].astype(int) - obs_C[i].astype(int)
    fam_df["n_primary_tested"] = fam_df[[f"tested_{i}" for i in prim]].sum(1)
    fam_df.to_csv(f_obs, sep="\t", index=False)
    # global null: number of families with recurrence >= k
    rows = []
    for nm, kmax in (("primary", len(prim)), ("secondary", len(sec))):
        for d in ("expanded", "contracted"):
            o = rec[(nm, d, "obs")]; pm = rec[(nm, d, "perm")]
            for k in range(1, kmax + 1):
                ob = int((o >= k).sum()); nl = (pm >= k).sum(1)
                rows.append(dict(tier=nm, direction=d, k=k, observed_n_families=ob, null_mean=nl.mean(), null_q025=np.quantile(nl, 0.025),
                                 null_q975=np.quantile(nl, 0.975), null_q95=np.quantile(nl, 0.95), p_perm=((nl >= ob).sum() + 1) / (B + 1)))
    pd.DataFrame(rows).to_csv(f_null, sep="\t", index=False)
    np.savez_compressed(f_npz, pfam=np.array(fams_all), **{f"{nm}_{d}": rec[(nm, d, "perm")] for nm in ("primary", "secondary") for d in ("expanded", "contracted")})
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
