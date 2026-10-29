#!/usr/bin/env python
"""POST HOC robustness checks of three positive EXPLORATORY results (euk). Not pre-registered.

    python scripts/robustness/exploratory_posthoc.py --part {loss,arch,coopt,all} [--no-registry]

1. loss  : convergent Pfam-domain LOSS (pfam_wide.py, flag 'lost', G_k = n domains lost in >= k transitions).
           Checks: reproduction; leave-one-transition-out (LOTO); species filters (genome-only, BUSCO C >= 80,
           both + not low_quality); dependence of per-contrast loss counts on n_M / n_S / BUSCO across the 30
           control pseudo-transitions; nulls that hold the per-transition loss counts fixed (uniform and
           control-propensity-weighted draws within each transition's eligible domains); matched-control nulls
           with controls matched on (log n_M, log n_S) or on log n_lost; functional enrichment of k >= 2 domains.
2. arch  : H11 new multi-domain architectures on stems vs controls (architectures.py).
           Checks: reproduction; LOTO; per-transition u; GLM + mean BUSCO of descendants; genome-only and
           genome + BUSCO >= 80 annotation sets (origins recomputed, n_desc/n_sister recounted); LOTO per variant.
3. coopt : co-option index C_i (cooption.py), versions pfam_expanded and toolkit_present.
           Checks: C vs n_M across controls; M of T_metazoa / T_streptophyta subsampled to 3, 4, 5 species
           (100 subsamples each, seed 11; dropped M species removed from the data), u vs the matched controls'
           full-data C and vs controls whose M is subsampled to the same size (20 subsamples per control).

Outputs: results/robustness/exploratory/*.tsv; registry rows tier 'exploratory', notes start 'posthoc'.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import warnings

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from mcell import branchstats as bs  # noqa: E402
from mcell import explore as ex  # noqa: E402

warnings.simplefilter("ignore")
OUT = os.path.join(ROOT, "results", "robustness", "exploratory")
TAB = os.path.join(ROOT, "runs", "analysis", "euk", "pfam_full")
SEED = 11
T0 = time.time()


def log(m):
    print(f"[posthoc {time.time() - T0:7.1f}s] {m}", file=sys.stderr, flush=True)


def save(df, name):
    df.to_csv(os.path.join(OUT, name), sep="\t", index=False)
    log(f"wrote {name}")


def species_info():
    sf = pd.read_csv(os.path.join(ROOT, "data", "species_final.tsv"), sep="\t")
    sf = sf.set_index("code")
    return sf[["data_type_final", "busco_C", "low_quality"]]


def species_filters(sf, allsp):
    g = set(sf.index[sf["data_type_final"] == "genome"]) & set(allsp)
    b = set(sf.index[sf["busco_C"] >= 80]) & set(allsp)
    hq = set(sf.index[sf["low_quality"] == 0]) & set(allsp)
    return {"all": set(allsp), "genome_only": g, "busco80": b, "genome_busco80_hq": g & b & hq}


def reg(row, args):
    if args.no_registry:
        return
    r = dict(row)
    r["tier"] = "exploratory"
    r["holm_p"] = "NA"
    r["notes"] = "posthoc; label=exploratory_posthoc; " + r.get("notes", "")
    bs.append_registry(os.path.join(ROOT, "results", "analysis_registry.tsv"), r)


# =====================================================================================================
# 1. convergent loss
# =====================================================================================================

def ok_contrast(con, keep, min_M=2, min_S=1):
    r = con.restrict(keep)
    return len(r.M) >= min(min_M, len(con.M)) and len(r.S) >= min_S


def lost_flags(X, cons):
    tabs = {c.cid: ex.domain_contrast(X, c) for c in cons}
    F = pd.DataFrame({k: t.set_index("pfam_acc")["status"] == "lost" for k, t in tabs.items()}).fillna(False)
    E = pd.DataFrame({k: t.set_index("pfam_acc")["frac_S"] >= 0.5 for k, t in tabs.items()}).fillna(False)
    return F.astype(bool), E.astype(bool), tabs


def conv(Ft, Fc, sets, n_perm=10000):
    sets = {t: [x for x in sets[t] if x in Fc.columns] for t in Ft.columns}
    per, glob, info = ex.convergence_test(Ft, Fc, sets, n_perm=n_perm, seed=SEED)
    out = {}
    for _, g in glob.iterrows():
        out[f"G{g['k_min']}_obs"] = int(g["n_domains_obs"])
        out[f"G{g['k_min']}_{g['scheme']}_null_mean"] = round(float(g["null_mean"]), 2)
        out[f"G{g['k_min']}_{g['scheme']}_p"] = float(g["p_perm"])
    out["matched_scheme"] = info.get("matched", "")
    return out, per


def fixed_margin_null(Ft, E, w, n_draws=2000, seed=SEED, chunk=250):
    """Per transition i draw |L_i| domains without replacement from its eligible set E_i, with weights w
    (Efraimidis-Spirakis keys; w = 1 -> uniform). Returns null arrays of G2, G3 and pairwise overlaps."""
    rng = np.random.default_rng(seed)
    dom = Ft.index
    tids = list(Ft.columns)
    idx = {t: np.flatnonzero(E.reindex(dom).fillna(False)[t].to_numpy()) for t in tids}
    n = {t: int(Ft[t].sum()) for t in tids}
    G2, G3 = np.zeros(n_draws, int), np.zeros(n_draws, int)
    pairs = [(a, b) for i, a in enumerate(tids) for b in tids[i + 1:]]
    PO = np.zeros((n_draws, len(pairs)), int)
    wv = w.reindex(dom).fillna(w.median()).to_numpy(float)
    for s in range(0, n_draws, chunk):
        m = min(chunk, n_draws - s)
        K = np.zeros((m, len(dom)), np.int16)
        Ls = {}
        for t in tids:
            e = idx[t]
            if n[t] == 0 or len(e) == 0:
                Ls[t] = np.zeros((m, len(dom)), bool)
                continue
            keys = np.log(rng.random((m, len(e)))) / wv[e][None, :]
            top = np.argpartition(-keys, n[t] - 1, axis=1)[:, :n[t]]
            L = np.zeros((m, len(dom)), bool)
            np.put_along_axis(L, e[top], True, axis=1)
            K += L
            Ls[t] = L
        G2[s:s + m] = (K >= 2).sum(1)
        G3[s:s + m] = (K >= 3).sum(1)
        for j, (a, b) in enumerate(pairs):
            PO[s:s + m, j] = (Ls[a] & Ls[b]).sum(1)
    return G2, G3, pairs, PO


def nearest_sets(feat_t, feat_c, k=10):
    Z = pd.concat([feat_t, feat_c])
    mu, sd = feat_c.mean(), feat_c.std(ddof=1).replace(0, 1)
    zt, zc = (feat_t - mu) / sd, (feat_c - mu) / sd
    out = {}
    for t in zt.index:
        d = np.sqrt(((zc - zt.loc[t]) ** 2).sum(1))
        out[t] = list(d.sort_values().index[:k])
    return out


CATS = {
    "cilium_flagellum": r"cili|flagell|axonem|dynein|intraflagellar|\bIFT|radial.spoke|bardet|BBS|tektin|centriol|basal.body|"
                        r"TCTEX|nexin|DM10|hydin|CFAP|transition.zone|tubby|kinesin|MORN|PACRG|ODA|IQ.calmodulin",
    "DNA_repair_recomb": r"repair|\bRAD\d|Rad\d|MutS|MutL|mismatch|glycosylase|photolyase|excision|Fanconi|BRCT|Ku70|Ku80|"
                         r"recombinase|helicase|XPG|XPC|UvrD|meio|synaptonemal|DMC1|Spo11|HORMA",
    "metabolism_enzyme": r"synthase|synthetase|dehydrogenase|reductase|oxidase|oxygenase|isomerase|mutase|aldolase|lyase|"
                         r"carboxylase|decarboxylase|transaminase|aminotransferase|hydratase|dehydratase|epimerase|"
                         r"phosphatase|hydrolase|esterase|lipase|transferase|kinase.*(sugar|carbohydrate|glycer|pyruv)|"
                         r"biosynth|catabol|metabol|NAD|FAD|CoA",
    "transporter_channel": r"transport|permease|channel|carrier|pump|porin|ABC|MFS|symporter|antiporter|exchanger",
    "TE_viral": r"transpos|integrase|retrotrans|reverse.transcr|\bRVT|Tnp|gag|viral|virus|phage|DDE",
    "photosynth_plastid": r"photosyst|chlorophyll|thylakoid|rubisco|ribulose|light.harvest|plastid|chloroplast|PsbP|"
                          r"phycobil|carotenoid",
    "DUF_UPF": r"^DUF|^UPF|Domain of unknown function|Uncharacteri[sz]ed",
}


def enrichment(hits, bg, ann_sets: dict):
    """Fisher (greater) of set membership among hits vs background; ann_sets: name -> set of accs."""
    hits, bg = set(hits), set(bg)
    rows = []
    for name, S in ann_sets.items():
        a = len(hits & S)
        b = len(hits) - a
        c = len((bg - hits) & S)
        d = len(bg - hits) - c
        if a + c < 3:
            continue
        orr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(dict(category=name, n_hits_in_cat=a, n_hits=len(hits), n_bg_in_cat=a + c, n_bg=len(bg),
                         frac_hits=a / max(1, len(hits)), frac_bg=(a + c) / max(1, len(bg)), odds_ratio=orr, p=p))
    R = pd.DataFrame(rows)
    if len(R):
        R["q_bh"] = ex.bh(R["p"].to_numpy())
        R = R.sort_values("p")
    return R


def part_loss(ctx, X, sf, args):
    log("part 1: convergent loss")
    allsp = list(X.index)
    filt = species_filters(sf, allsp)
    sets = {t: ctx.match.sets[t] for t in ctx.tids}
    # ---- baseline (all species)
    Ft, Et, tabs = lost_flags(X, ctx.trans)
    Fc, Ec, ctabs = lost_flags(X, ctx.ctrls)
    Fc.columns = [c.replace("ctrl:", "") for c in Fc.columns]
    Ec.columns = Fc.columns
    rows = []
    base, per = conv(Ft, Fc, sets)
    rows.append(dict(analysis="baseline_reproduction", dropped="", species_set="all", n_transitions=Ft.shape[1],
                     n_controls=Fc.shape[1], n_lost_per_transition=json.dumps({k: int(v) for k, v in Ft.sum().items()}),
                     **base))
    log(f"baseline {base}")
    # ---- LOTO
    for t in ctx.tids:
        keep = [x for x in Ft.columns if x != t]
        o, _ = conv(Ft[keep], Fc, sets)
        rows.append(dict(analysis="LOTO", dropped=t, species_set="all", n_transitions=len(keep), n_controls=Fc.shape[1],
                         n_lost_per_transition="", **o))
    # drop both reduced-proteome transitions
    keep = [x for x in Ft.columns if x not in ("T_dictyostelia", "T_rhodophyta")]
    o, _ = conv(Ft[keep], Fc, sets)
    rows.append(dict(analysis="LOTO2", dropped="T_dictyostelia;T_rhodophyta", species_set="all", n_transitions=len(keep),
                     n_controls=Fc.shape[1], n_lost_per_transition="", **o))
    # ---- species filters
    cov_rows = []
    for fname, keep_sp in filt.items():
        if fname == "all":
            continue
        Xs = X.loc[[s for s in allsp if s in keep_sp]]
        tr = [c for c in ctx.trans if ok_contrast(c, Xs.index)]
        cr = [c for c in ctx.ctrls if ok_contrast(c, Xs.index)]
        for c in ctx.trans + ctx.ctrls:
            r = c.restrict(Xs.index)
            cov_rows.append(dict(species_set=fname, contrast=c.cid, n_M=len(r.M), n_S=len(r.S), n_out=len(r.outside),
                                 evaluable=c in tr or c in cr))
        ft, _, _ = lost_flags(Xs, tr)
        fc, _, _ = lost_flags(Xs, cr)
        fc.columns = [c.replace("ctrl:", "") for c in fc.columns]
        o, _ = conv(ft, fc, sets)
        rows.append(dict(analysis="species_filter", dropped=",".join(c.cid for c in ctx.trans if c not in tr),
                         species_set=fname, n_transitions=ft.shape[1], n_controls=fc.shape[1],
                         n_lost_per_transition=json.dumps({k: int(v) for k, v in ft.sum().items()}), **o))
        for t in [c.cid for c in tr if c.cid in ("T_dictyostelia", "T_rhodophyta")]:
            kk = [x for x in ft.columns if x != t]
            o, _ = conv(ft[kk], fc, sets)
            rows.append(dict(analysis="species_filter+LOTO", dropped=t, species_set=fname, n_transitions=len(kk),
                             n_controls=fc.shape[1], n_lost_per_transition="", **o))
        log(f"filter {fname}: {rows[-1]}")
    save(pd.DataFrame(cov_rows), "loss_species_filter_coverage.tsv")

    # ---- dependence of loss counts on n_M, n_S, completeness (controls), transitions vs prediction
    def feats(c):
        r = c.restrict(X.index)
        bM = sf.reindex(r.M)["busco_C"].astype(float)
        bS = sf.reindex(r.S)["busco_C"].astype(float)
        txM = (sf.reindex(r.M)["data_type_final"] == "transcriptome").mean()
        return dict(n_M=len(r.M), n_S=len(r.S), busco_M=bM.mean(), min_busco_M=bM.min(), busco_S=bS.mean(),
                    frac_tx_M=txM)
    F = []
    for c in ctx.trans:
        F.append(dict(contrast=c.cid, kind="transition", n_lost=int(Ft[c.cid].sum()), n_eligible=int(Et[c.cid].sum()),
                      **feats(c)))
    for c in ctx.ctrls:
        b = c.branch
        F.append(dict(contrast=b, kind="control", n_lost=int(Fc[b].sum()), n_eligible=int(Ec[b].sum()), **feats(c)))
    F = pd.DataFrame(F)
    F["frac_lost_of_eligible"] = F["n_lost"] / F["n_eligible"]
    Cc = F[F["kind"] == "control"]
    sp_rows = []
    for v in ("n_M", "n_S", "busco_M", "min_busco_M", "busco_S", "frac_tx_M"):
        for y in ("n_lost", "frac_lost_of_eligible"):
            rho, p = stats.spearmanr(Cc[v], Cc[y])
            sp_rows.append(dict(y=y, x=v, n=len(Cc), spearman_rho=rho, p=p))
    import statsmodels.api as sm
    design = ["log_n_M", "log_n_S", "busco_M", "frac_tx_M"]
    F["log_n_M"] = np.log(F["n_M"])
    F["log_n_S"] = np.log(F["n_S"])
    Xd = sm.add_constant(F.loc[F["kind"] == "control", design])
    fit = sm.OLS(np.log(F.loc[F["kind"] == "control", "n_lost"] + 1), Xd).fit()
    F["pred_log_n_lost"] = fit.predict(sm.add_constant(F[design], has_constant="add"))
    F["resid_log_n_lost"] = np.log(F["n_lost"] + 1) - F["pred_log_n_lost"]
    F["pred_n_lost"] = np.exp(F["pred_log_n_lost"]) - 1
    cres = F.loc[F["kind"] == "control", "resid_log_n_lost"].to_numpy()
    F["resid_percentile_vs_controls"] = [((cres < r).sum() + 0.5 * (cres == r).sum() + 0.5) / (len(cres) + 1)
                                         for r in F["resid_log_n_lost"]]
    for k_, v_ in fit.params.items():
        sp_rows.append(dict(y="OLS log(n_lost+1) on controls", x=k_, n=int(fit.nobs), spearman_rho=v_,
                            p=fit.pvalues[k_]))
    sp_rows.append(dict(y="OLS log(n_lost+1) on controls", x="R2", n=int(fit.nobs), spearman_rho=fit.rsquared, p=np.nan))
    save(F, "loss_counts_per_contrast.tsv")
    save(pd.DataFrame(sp_rows), "loss_count_dependence.tsv")

    # ---- matched-on-size / matched-on-loss-count null
    ft_ = F[F["kind"] == "transition"].set_index("contrast")
    fc_ = F[F["kind"] == "control"].set_index("contrast")
    for lab, cols in (("matched_nM_nS", ["log_n_M", "log_n_S"]), ("matched_nM_nS_busco", ["log_n_M", "log_n_S", "busco_M"]),
                      ("matched_n_lost", ["log_n_lost"])):
        if lab == "matched_n_lost":
            ft_["log_n_lost"] = np.log(ft_["n_lost"] + 1)
            fc_["log_n_lost"] = np.log(fc_["n_lost"] + 1)
        ns = nearest_sets(ft_[cols], fc_[cols], k=10)
        o, _ = conv(Ft, Fc, ns)
        rows.append(dict(analysis=f"null_{lab}(k=10)", dropped="", species_set="all", n_transitions=Ft.shape[1],
                         n_controls=Fc.shape[1], n_lost_per_transition=json.dumps(ns), **o))
        log(f"{lab}: {o}")

    # ---- fixed-margin nulls
    elig_c = Ec.sum(1)
    lost_c = Fc.reindex(Ec.index).fillna(False).sum(1)
    w = ((lost_c + 0.5) / (elig_c + 1)).reindex(Ft.index.union(Ec.index))
    w = w.fillna(float(w.median()))
    fm_rows, pair_rows = [], []
    G2o, G3o = int((Ft.sum(1) >= 2).sum()), int((Ft.sum(1) >= 3).sum())
    for wlab, ww in (("uniform", pd.Series(1.0, index=w.index)), ("control_propensity", w)):
        for drop in [None] + list(ctx.tids) + ["T_dictyostelia+T_rhodophyta"]:
            if drop is None:
                cols = list(Ft.columns)
            elif drop == "T_dictyostelia+T_rhodophyta":
                cols = [c for c in Ft.columns if c not in ("T_dictyostelia", "T_rhodophyta")]
            else:
                cols = [c for c in Ft.columns if c != drop]
            ft = Ft[cols]
            g2o, g3o = int((ft.sum(1) >= 2).sum()), int((ft.sum(1) >= 3).sum())
            G2, G3, pairs, PO = fixed_margin_null(ft, Et[cols], ww, n_draws=args.n_fm)
            fm_rows.append(dict(weights=wlab, dropped=drop or "", n_transitions=len(cols), G2_obs=g2o,
                                G2_null_mean=G2.mean(), G2_null_q95=np.quantile(G2, 0.95),
                                G2_p=(1 + (G2 >= g2o).sum()) / (len(G2) + 1), G3_obs=g3o, G3_null_mean=G3.mean(),
                                G3_null_q95=np.quantile(G3, 0.95), G3_p=(1 + (G3 >= g3o).sum()) / (len(G3) + 1),
                                obs_over_null_G2=g2o / max(G2.mean(), 1e-9), n_draws=len(G2)))
            if drop is None:
                for j, (a_, b_) in enumerate(pairs):
                    o_ = int((Ft[a_] & Ft[b_]).sum())
                    pair_rows.append(dict(weights=wlab, t1=a_, t2=b_, overlap_obs=o_, overlap_null_mean=PO[:, j].mean(),
                                          p=(1 + (PO[:, j] >= o_).sum()) / (len(PO) + 1)))
            log(f"fixed-margin {wlab} drop={drop}: {fm_rows[-1]['G2_obs']} vs {fm_rows[-1]['G2_null_mean']:.1f} "
                f"p={fm_rows[-1]['G2_p']:.3g}")
    save(pd.DataFrame(fm_rows), "loss_fixed_margin_null.tsv")
    save(pd.DataFrame(pair_rows), "loss_pairwise_overlap.tsv")
    R = pd.DataFrame(rows)
    save(R, "loss_convergence_checks.tsv")

    # ---- functional characterisation
    names_tab = pd.read_csv(os.path.join(TAB, "names.tsv"), sep="\t")
    clans = pd.read_csv(os.path.join(ROOT, "data", "raw", "pfam", "Pfam-A.clans.tsv.gz"), sep="\t", header=None,
                        names=["pfam_acc", "clan", "clan_name", "pfam_name", "desc"])
    clans = clans.set_index("pfam_acc")
    nm = clans["pfam_name"].to_dict()
    nm.update({a: n for a, n in zip(names_tab["pfam_acc"], names_tab["pfam_name"]) if a not in nm})
    desc = clans["desc"].fillna("").to_dict()
    k = Ft.sum(1)
    hits = set(k.index[k >= 2])
    ke = Et.sum(1)
    bg1 = set(ke.index[ke >= 2])              # eligible to be lost in >= 2 transitions
    bg2 = set(ke.index[ke >= 1])              # present (>= 50 % of sisters) at >= 1 P
    allacc = sorted(bg2 | hits)
    ann = {}
    for c, pat in CATS.items():
        rx = re.compile(pat, re.I)
        ann[f"kw:{c}"] = {a for a in allacc if rx.search(f"{nm.get(a, '')} | {desc.get(a, '')}")}
    h12 = pd.read_csv(os.path.join(ROOT, "data", "h12_sets.tsv"), sep="\t")
    for mname, g in h12.groupby("module"):
        ann[f"h12:{mname}"] = set(g["pfam_acc"].str.split(".").str[0])
    ann["h12:cilia_all"] = set(h12.loc[h12["module"] != "phagotrophy_exploratory", "pfam_acc"].str.split(".").str[0])
    tk = pd.read_csv(os.path.join(ROOT, "data", "toolkit_pfam.tsv"), sep="\t")
    for cat, g in tk.groupby("category"):
        ann[f"toolkit:{cat}"] = set(g["pfam_acc"].str.split(".").str[0])
    cl = clans["clan"].dropna()
    cl_name = clans.dropna(subset=["clan"]).drop_duplicates("clan").set_index("clan")["clan_name"].to_dict()
    for c_, g in cl.groupby(cl):
        ann[f"clan:{c_}:{cl_name.get(c_, '')}"] = set(g.index)
    ann["no_clan"] = set(allacc) - set(cl.index)
    E1 = enrichment(hits, bg1, ann)
    E1.insert(0, "background", "eligible_in>=2_transitions")
    E2 = enrichment(hits, bg2, ann)
    E2.insert(0, "background", "present_at>=1_P")
    save(pd.concat([E1, E2]), "loss_functional_enrichment.tsv")
    # per-domain table of convergently lost domains
    D = pd.DataFrame({"pfam_acc": sorted(hits)})
    D["pfam_name"] = D["pfam_acc"].map(nm)
    D["desc"] = D["pfam_acc"].map(desc)
    D["clan"] = D["pfam_acc"].map(clans["clan"].to_dict())
    D["k_lost"] = D["pfam_acc"].map(k)
    D["transitions"] = [";".join(Ft.columns[Ft.loc[a].to_numpy()]) for a in D["pfam_acc"]]
    D["n_controls_lost"] = D["pfam_acc"].map(lost_c).fillna(0).astype(int)
    D["n_controls_eligible"] = D["pfam_acc"].map(elig_c).fillna(0).astype(int)
    for c in CATS:
        D[f"kw_{c}"] = D["pfam_acc"].isin(ann[f"kw:{c}"])
    save(D.sort_values(["k_lost", "pfam_acc"], ascending=[False, True]), "loss_convergent_domains.tsv")
    # which transitions contribute to k>=2
    contrib = pd.DataFrame({t: [int((Ft.loc[list(hits), t]).sum())] for t in Ft.columns})
    save(contrib.T.reset_index().rename(columns={"index": "transition", 0: "n_of_k2_domains_lost_here"}),
         "loss_k2_contribution.tsv")
    # propensity: are convergently lost domains generally loss-prone among controls?
    lp = pd.DataFrame(dict(group=["k>=2 lost", "background eligible>=2"],
                           median_control_loss_rate=[float(w.reindex(list(hits)).median()),
                                                     float(w.reindex(list(bg1)).median())],
                           mean_control_loss_rate=[float(w.reindex(list(hits)).mean()),
                                                   float(w.reindex(list(bg1)).mean())]))
    mw = stats.mannwhitneyu(w.reindex(list(hits)).dropna(), w.reindex(list(bg1 - hits)).dropna(), alternative="greater")
    lp["mwu_p_hits_gt_bg"] = mw.pvalue
    save(lp, "loss_control_propensity.tsv")

    # ---- convergence beyond the cilium: drop domains with a ciliary phylogenetic profile
    cil_markers = [a for a in ann["h12:cilia_all"] if a in X.columns]
    ciliated = (X[cil_markers] > 0).mean(1) >= 0.5
    P_ = (X > 0)
    a_ = P_[ciliated.to_numpy()].sum(0)
    b_ = P_[~ciliated.to_numpy()].sum(0)
    n1, n0 = int(ciliated.sum()), int((~ciliated).sum())
    p1, p0 = a_ / n1, b_ / n0
    pall = (a_ + b_) / (n1 + n0)
    phi = (p1 - p0) * np.sqrt(n1 * n0) / ((n1 + n0) * np.sqrt(np.clip(pall * (1 - pall), 1e-12, None)))
    cil_dom = set(phi.index[(phi >= 0.4)]) | ann["kw:cilium_flagellum"] | ann["h12:cilia_all"]
    save(pd.DataFrame(dict(pfam_acc=phi.index, pfam_name=phi.index.map(nm), phi_ciliated=phi.to_numpy(),
                           frac_in_ciliated=p1.to_numpy(), frac_in_nonciliated=p0.to_numpy(),
                           excluded=phi.index.isin(list(cil_dom)))).sort_values("phi_ciliated", ascending=False),
         "loss_ciliary_profile.tsv")
    log(f"ciliated species {n1}/{n1 + n0}; ciliary-profile domains excluded: {len(cil_dom)}; "
        f"of k>=2 hits: {len(hits & cil_dom)}")
    nc_rows = []
    keepd = [d for d in Ft.index.union(Fc.index) if d not in cil_dom]
    Ftn, Fcn = Ft.reindex(keepd).fillna(False).astype(bool), Fc.reindex(keepd).fillna(False).astype(bool)
    for drop in [None, "T_dictyostelia", "T_fungi_hyphae"]:
        cols = [c for c in Ftn.columns if c != drop]
        o, _ = conv(Ftn[cols], Fcn, sets)
        G2, G3, _, _ = fixed_margin_null(Ftn[cols], Et.reindex(keepd).fillna(False)[cols], w, n_draws=args.n_fm)
        G2u, _, _, _ = fixed_margin_null(Ftn[cols], Et.reindex(keepd).fillna(False)[cols],
                                         pd.Series(1.0, index=w.index), n_draws=args.n_fm)
        nc_rows.append(dict(dropped=drop or "", n_domains_excluded=len(cil_dom), n_transitions=len(cols),
                            n_lost_per_transition=json.dumps({k: int(v) for k, v in Ftn[cols].sum().items()}), **o,
                            fm_uniform_null_mean=G2u.mean(), fm_uniform_p=(1 + (G2u >= o["G2_obs"]).sum()) / (len(G2u) + 1),
                            fm_propensity_null_mean=G2.mean(),
                            fm_propensity_p=(1 + (G2 >= o["G2_obs"]).sum()) / (len(G2) + 1)))
        log(f"non-ciliary drop={drop}: {nc_rows[-1]}")
    NC = pd.DataFrame(nc_rows)
    save(NC, "loss_convergence_nonciliary.tsv")
    for _, r in NC.iterrows():
        reg(dict(hypothesis="Q2 Pfam-wide convergence (lost) robustness",
                 statistic="G2 after excluding ciliary-profile domains (phi>=0.4 with ciliated species, keyword or H12 "
                           "cilia set)" + (f"; LOTO -{r['dropped']}" if r["dropped"] else ""),
                 n=int(r["n_transitions"]), k=30, estimate=float(r["G2_obs"]), p=float(r["G2_permute_p"]), seed=SEED,
                 n_draws=10000, scheme="permute",
                 notes=f"p_matched={r['G2_matched_p']:.3g}; fixed-margin uniform p={r['fm_uniform_p']:.3g}; "
                       f"propensity p={r['fm_propensity_p']:.3g}; excluded={r['n_domains_excluded']}"), args)

    fmr = pd.DataFrame(fm_rows)
    for _, r in R.iterrows():
        if r["analysis"] in ("LOTO", "LOTO2", "species_filter", "null_matched_nM_nS(k=10)", "null_matched_n_lost(k=10)"):
            reg(dict(hypothesis="Q2 Pfam-wide convergence (lost) robustness",
                     statistic=f"G2 = n domains lost in >= 2 transitions; {r['analysis']}",
                     n=int(r["n_transitions"]), k=int(r["n_controls"]), estimate=float(r["G2_obs"]),
                     p=float(r["G2_permute_p"] if not r["analysis"].startswith("null_") else r["G2_matched_p"]),
                     seed=SEED, n_draws=10000, scheme="permute" if not r["analysis"].startswith("null_") else "matched",
                     notes=f"dropped={r['dropped']}; species_set={r['species_set']}; "
                           f"null_mean_permute={r.get('G2_permute_null_mean')}; p_matched={r.get('G2_matched_p')}"), args)
    for _, r in fmr[fmr["dropped"].isin(["", "T_dictyostelia", "T_dictyostelia+T_rhodophyta"])].iterrows():
        reg(dict(hypothesis="Q2 Pfam-wide convergence (lost) robustness",
                 statistic=f"G2 vs fixed-margin null ({r['weights']} weights within eligible domains)",
                 n=int(r["n_transitions"]), estimate=float(r["G2_obs"]), p=float(r["G2_p"]), seed=SEED,
                 n_draws=int(r["n_draws"]), scheme="fixed_margin",
                 notes=f"dropped={r['dropped']}; null_mean={r['G2_null_mean']:.1f}"), args)


# =====================================================================================================
# 2. architectures
# =====================================================================================================

def arch_branch_counts(ctx, arch, annotated):
    org = ex.architecture_origins(arch, ctx.tree, annotated)
    return ex.branch_architecture_counts(org, ctx.bt, ctx.tree, annotated), org


def bt_for_species(ctx, keep, sf):
    """copy of bt with n_desc / n_sister recounted on ``keep`` and mean BUSCO of kept descendants."""
    bt = ctx.bt.copy()
    tips = ex.clade_tips(ctx.tree)
    keep = set(keep)
    nd, ns, mb = [], [], []
    for b in bt.index:
        edges = bs.path_edges(bt, b) if b in bt.index else [b]
        first = edges[0]
        u = bt.at[first, "parent"]
        d = tips[b] & keep
        s = (tips[u] - tips[first]) & keep if u in tips else set()
        nd.append(len(d))
        ns.append(len(s))
        mb.append(float(sf.reindex(sorted(d))["busco_C"].mean()) if d else np.nan)
    bt["n_desc_kept"] = nd
    bt["n_sister_kept"] = ns
    bt["log_n_desc_kept"] = np.log(np.maximum(bt["n_desc_kept"], 0.5))
    bt["log_n_sister_kept"] = np.log(np.maximum(bt["n_sister_kept"], 0.5))
    bt["mean_busco_kept"] = mb
    return bt


def h11_test(y, bt, stems, sets, pool, covs, label, n_draws):
    res = bs.count_residual_statistic(y, bt, pool=pool, covariates=covs)
    s = res["s"]
    out = []
    for vname, vals in (("residual", s), ("raw", y)):
        r = bs.relabel_test(vals, sets, stems=stems, alternative="greater", n_draws=n_draws, seed=SEED)
        out.append(dict(variant=label, version=vname, covariates=",".join(res.attrs.get("covariates", covs))
                        if vname == "residual" else "", glm_method=res.attrs.get("method") if vname == "residual" else "",
                        n=r.n, T_obs=r.T_obs, p=r.p, u=json.dumps({k: round(float(v), 3) for k, v in r.u.items()})))
    return out, s


def part_arch(ctx, sf, args):
    log("part 2: architectures")
    arch = pd.read_parquet(os.path.join(TAB, "arch.parquet"))
    arch["species"] = arch["species"].astype(str)
    allsp = sorted(arch["species"].unique())
    filt = species_filters(sf, allsp)
    stems = ctx.stems
    sets = {t: list(ctx.match.sets[t]) for t in ctx.tids}
    rows, per = [], []
    for fname in ("all", "genome_only", "genome_busco80_hq", "busco80"):
        keep = sorted(filt[fname])
        br, _ = arch_branch_counts(ctx, arch, keep)
        bt = bt_for_species(ctx, keep, sf)
        for stat in ("new_arch_multi", "new_arch_all"):
            y = br[stat].astype(float)
            # a branch with < 2 kept descendant species cannot carry a 'new' architecture (needs >= 2 species):
            # structural zero -> excluded (NaN) rather than counted
            y[bt["n_desc_kept"] < 2] = np.nan
            variants = [("base_GLM", list(bs.RESID_COVARIATES), bs.RESID_COVARIATES)]
            variants.append(("GLM+mean_busco_desc", list(bs.RESID_COVARIATES) + ["mean_busco_kept"], None))
            if fname != "all":
                variants.append(("GLM_recounted_n", ["log_bl", "log_n_desc_kept", "log_n_sister_kept", "mean_busco_kept"],
                                 None))
            for vlab, covs, _ in variants:
                res, s = h11_test(y, bt, stems, sets, ctx.pool, covs, vlab, args.n_draws)
                for r in res:
                    rows.append(dict(species_set=fname, statistic=stat, dropped="", **r))
                if stat == "new_arch_multi":
                    for t, b in stems.items():
                        cn = s.reindex(sets[t]).dropna().to_numpy()
                        rn = y.reindex(sets[t]).dropna().to_numpy()
                        per.append(dict(species_set=fname, variant=vlab, transition=t, stem=b, y=y[b],
                                        n_controls_eval=int(len(rn)),
                                        median_y_controls=float(np.median(rn)),
                                        n_desc_kept=int(bt.at[b, "n_desc_kept"]), s_resid=s[b],
                                        u_resid=bs.u_score(s[b], cn) if np.isfinite(s[b]) else np.nan,
                                        u_raw=bs.u_score(y[b], rn) if np.isfinite(y[b]) else np.nan))
                # LOTO (residual + raw), every variant
                for drop in ctx.tids:
                    st2 = {t: b for t, b in stems.items() if t != drop}
                    s2 = {t: v for t, v in sets.items() if t != drop}
                    for vname, vals in (("residual", s), ("raw", y)):
                        r = bs.relabel_test(vals, s2, stems=st2, alternative="greater", n_draws=args.n_draws, seed=SEED)
                        rows.append(dict(species_set=fname, statistic=stat, dropped=drop, variant=vlab, version=vname,
                                         covariates="", glm_method="", n=r.n, T_obs=r.T_obs, p=r.p, u=""))
            log(f"arch {fname} {stat} done")
    R = pd.DataFrame(rows)
    save(R, "arch_checks.tsv")
    save(pd.DataFrame(per), "arch_per_transition.tsv")
    for _, r in R[(R["statistic"] == "new_arch_multi")].iterrows():
        if r["dropped"] not in ("", "T_metazoa"):
            continue
        reg(dict(hypothesis="H11 new domain architectures robustness",
                 statistic=f"new_arch_multi {r['version']}; {r['variant']}; species={r['species_set']}"
                           + (f"; LOTO -{r['dropped']}" if r["dropped"] else ""),
                 n=int(r["n"]), k=20, estimate=float(r["T_obs"]), p=float(r["p"]), seed=SEED, n_draws=args.n_draws,
                 scheme="permute", notes=f"estimate=T_obs (mean u); one-sided greater; covariates={r['covariates']}"),
            args)


# =====================================================================================================
# 3. co-option
# =====================================================================================================

def part_coopt(ctx, X, args):
    log("part 3: co-option")
    tk = pd.read_csv(os.path.join(ROOT, "data", "toolkit_pfam.tsv"), sep="\t")
    tk = set(tk.loc[tk["role"].isin(["member", "conditional_member"]), "pfam_acc"].str.split(".").str[0])
    versions = {"pfam_expanded": (None, "expanded"), "toolkit_present": (tk, "present")}
    rng = np.random.default_rng(SEED)
    full, ctrlC = [], {}
    for v, (fams, rule) in versions.items():
        cc = {}
        for con in ctx.trans + ctx.ctrls:
            r = ex.coption_index(X, con, ctx.dtype, families=fams, rule=rule)
            full.append(dict(version=v, kind=con.kind, contrast=con.cid, branch=con.branch, n_M=r["n_M"], n_S=r["n_S"],
                             n_out=r["n_out"], n_F=r["n_F"], n_P=r["n_P"], C=r["C"]))
            if con.kind == "control":
                cc[con.branch] = r["C"]
        ctrlC[v] = pd.Series(cc)
    FU = pd.DataFrame(full)
    save(FU, "coopt_full_data.tsv")
    dep = []
    for v in versions:
        c = FU[(FU["version"] == v) & (FU["kind"] == "control")].dropna(subset=["C"])
        for x in ("n_M", "n_out", "n_F"):
            rho, p = stats.spearmanr(c[x], c["C"])
            dep.append(dict(version=v, x=x, n_controls=len(c), spearman_rho=rho, p=p))
    save(pd.DataFrame(dep), "coopt_C_vs_size_controls.tsv")

    def sub_C(con, m, fams, rule, rng_):
        Ms = list(con.restrict(X.index).M)
        if len(Ms) <= m:
            r = ex.coption_index(X, con, ctx.dtype, families=fams, rule=rule)
            return r["C"], r["n_F"]
        pick = set(rng_.choice(Ms, size=m, replace=False))
        Xs = X.drop(index=[s for s in Ms if s not in pick])
        r = ex.coption_index(Xs, con, ctx.dtype, families=fams, rule=rule)
        return r["C"], r["n_F"]

    sizes = (3, 4, 5)
    # controls subsampled to the same M size (n_ctrl_sub draws, median C per control)
    ctrl_sub = {}
    for v, (fams, rule) in versions.items():
        for m in sizes:
            vals = {}
            for con in ctx.ctrls:
                cs = [sub_C(con, m, fams, rule, rng)[0] for _ in range(args.n_ctrl_sub if len(con.M) > m else 1)]
                vals[con.branch] = float(np.nanmedian(cs)) if np.isfinite(cs).any() else np.nan
            ctrl_sub[(v, m)] = pd.Series(vals)
            log(f"controls subsampled {v} m={m}")
    rows, sumr = [], []
    target = ["T_metazoa", "T_streptophyta"]
    for v, (fams, rule) in versions.items():
        fullT = FU[(FU["version"] == v) & (FU["kind"] == "transition")].set_index("contrast")
        for m in sizes:
            for t in target:
                con = next(c for c in ctx.trans if c.cid == t)
                N = ctx.match.sets[t]
                cfull = ctrlC[v].reindex(N).dropna().to_numpy()
                csub = ctrl_sub[(v, m)].reindex(N).dropna().to_numpy()
                for i in range(args.n_sub):
                    C, nF = sub_C(con, m, fams, rule, rng)
                    rows.append(dict(version=v, transition=t, m=m, rep=i, n_F=nF, C=C,
                                     u_vs_full_controls=bs.u_score(C, cfull) if np.isfinite(C) else np.nan,
                                     u_vs_size_matched_controls=bs.u_score(C, csub) if np.isfinite(C) else np.nan))
                d = pd.DataFrame([r for r in rows if r["version"] == v and r["transition"] == t and r["m"] == m])
                sumr.append(dict(version=v, transition=t, m=m, n_sub=len(d), C_full=float(fullT.at[t, "C"]),
                                 C_median=d["C"].median(), C_q025=d["C"].quantile(0.025), C_q975=d["C"].quantile(0.975),
                                 median_n_F=d["n_F"].median(),
                                 median_C_controls_full=float(np.median(cfull)),
                                 median_C_controls_size_matched=float(np.median(csub)),
                                 frac_C_below_ctrl_median_full=float((d["C"] < np.median(cfull)).mean()),
                                 frac_C_below_ctrl_median_sizematched=float((d["C"] < np.median(csub)).mean()),
                                 median_u_full=d["u_vs_full_controls"].median(),
                                 median_u_sizematched=d["u_vs_size_matched_controls"].median(),
                                 frac_u_min_sizematched=float((d["u_vs_size_matched_controls"] <= 0.5 / (len(csub) + 1) + 1e-9).mean()),
                                 frac_u_le_0_1_sizematched=float((d["u_vs_size_matched_controls"] <= 0.1).mean())))
                log(f"coopt {v} {t} m={m}: {sumr[-1]['C_median']:.3f} vs ctrl {sumr[-1]['median_C_controls_size_matched']:.3f}")
    save(pd.DataFrame(rows), "coopt_subsamples.tsv")
    S = pd.DataFrame(sumr)
    # recompute T = mean u across the 7 transitions with M subsampled (others full / unchanged; their M <= ... )
    trows = []
    for v, (fams, rule) in versions.items():
        fullT = FU[(FU["version"] == v) & (FU["kind"] == "transition")].set_index("contrast")
        for m in sizes:
            us = []
            for t in ctx.tids:
                N = ctx.match.sets[t]
                csub = ctrl_sub[(v, m)].reindex(N).dropna().to_numpy()
                if t in target:
                    d = S[(S["version"] == v) & (S["transition"] == t) & (S["m"] == m)]
                    us.append(float(d["median_u_sizematched"].iloc[0]))
                else:
                    con = next(c for c in ctx.trans if c.cid == t)
                    cs = [sub_C(con, m, fams, rule, rng)[0] for _ in range(args.n_ctrl_sub if len(con.M) > m else 1)]
                    us.append(bs.u_score(float(np.nanmedian(cs)), csub))
            trows.append(dict(version=v, m=m, T_mean_u=float(np.mean(us)),
                              u=json.dumps({t: round(x, 3) for t, x in zip(ctx.tids, us)})))
    save(S, "coopt_subsample_summary.tsv")
    save(pd.DataFrame(trows), "coopt_T_size_matched.tsv")
    for _, r in S.iterrows():
        reg(dict(hypothesis=f"Q3 co-option index ({r['version']}) robustness",
                 statistic=f"C_i with M subsampled to {r['m']} species ({r['n_sub']} subsamples); vs controls subsampled "
                           f"to the same size; estimate = median C_sub - median C(N_i, size-matched)",
                 n=int(r["n_sub"]), k=20, estimate=float(r["C_median"] - r["median_C_controls_size_matched"]),
                 ci_low=float(r["C_q025"]), ci_high=float(r["C_q975"]), seed=SEED, n_draws=int(r["n_sub"]), scheme="subsample",
                 notes=f"transition={r['transition']}; descriptive, no p; median_u_sizematched={r['median_u_sizematched']:.3g}; "
                       f"C_full={r['C_full']:.3g}; ci = 2.5/97.5% of subsample C"), args)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all", choices=["loss", "arch", "coopt", "all"])
    ap.add_argument("--n-draws", type=int, default=10000)
    ap.add_argument("--n-fm", type=int, default=2000)
    ap.add_argument("--n-sub", type=int, default=100)
    ap.add_argument("--n-ctrl-sub", type=int, default=20)
    ap.add_argument("--no-registry", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    ctx = ex.load_context("euk", "default", need_cache=False)
    sf = species_info()
    man = json.load(open(os.path.join(TAB, "manifest.json")))
    ann = man["annotated"]
    X = ex.count_matrix(pd.read_parquet(os.path.join(TAB, "counts.parquet")), ann)
    if a.part in ("loss", "all"):
        part_loss(ctx, X, sf, a)
    if a.part in ("arch", "all"):
        part_arch(ctx, sf, a)
    if a.part in ("coopt", "all"):
        part_coopt(ctx, X, a)
    log("done")


if __name__ == "__main__":
    sys.exit(main())
