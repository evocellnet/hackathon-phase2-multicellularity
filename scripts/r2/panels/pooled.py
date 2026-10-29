#!/usr/bin/env python
"""Pooled residual model, R2-P1 (confirmatory) and secondary meta-analyses over the finished panels
(task steps 6-7; METHODS_DEFINITIONS.txt items 5-9).

    pooled.py --tag interim_5panels      # -> results/r2/panels/pooled/<tag>/ (refuses to overwrite)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, PANELS, OUT, log, guard, write_done, paths  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from mcell import branchstats as bs  # noqa: E402

COVS = ["log_bl", "log_n", "frac_tx", "mean_busco", "log_dout"]
TIP_COVS = ["log_bl", "frac_tx", "mean_busco", "log_dout"]
P1_MAP = {"R6_fucophycidae": "Fucophycidae"}
MIN_MATCH = 5


def depth_class(d):
    return "A" if d <= 3 else ("B" if d <= 6 else "C")


# ----------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------

def load(panels):
    B, E = [], []
    for R in panels:
        b = pd.read_csv(os.path.join(OUT, R, "branch_table.tsv"), sep="\t")
        tf = os.path.join(OUT, R, "tip_gain_ladder.tsv")
        if os.path.exists(tf):   # add-on tip ladder (DEVIATIONS D-P6)
            t = pd.read_csv(tf, sep="\t").set_index("branch_child")
            m = b["branch_child"].isin(t.index) & ~b["laddered"].astype(bool)
            for c in ("seq_novel", "seq_novel_tx_sens", "seq_novel_genomes_only", "split_repaired", "gains_after_repair"):
                b.loc[m, c] = t.loc[b.loc[m, "branch_child"], c].to_numpy()
            b.loc[m, "laddered"] = True
            b["tip_control"] = (b["is_leaf"].astype(bool) & ~b["is_event"].astype(bool) & (b["state_child"] == "uni")
                                & b["laddered"].astype(bool))
        e = pd.read_csv(os.path.join(OUT, R, "event_table.tsv"), sep="\t")
        B.append(b)
        E.append(e)
    B = pd.concat(B, ignore_index=True)
    E = pd.concat(E, ignore_index=True)
    B["key"] = B["panel"] + ":" + B["branch_child"]
    B = B.set_index("key")
    B["log_bl"] = np.log(B["bl"].clip(lower=1e-4))
    B["log_n"] = np.log(B["n_desc"])
    B["log_dout"] = np.log(B["dist_outgroup"].clip(lower=1e-4))
    B["depth_class"] = B["depth"].map(depth_class)
    B["y"] = B["seq_novel"]
    return B, E


def design(B, rows, covs, ref_rows, panels):
    """Intercept + panel dummies (first panel = reference) + covariates standardised on ref_rows."""
    covs = [c for c in covs if B.loc[ref_rows, c].astype(float).std(ddof=1) > 0]   # drop covariates constant in the pool
    Z = B.loc[rows, covs].astype(float)
    mu, sd = B.loc[ref_rows, covs].mean(), B.loc[ref_rows, covs].std(ddof=1)
    Z = (Z - mu) / sd
    cols = [np.ones(len(rows))]
    names = ["intercept"]
    for P in panels[1:]:
        cols.append((B.loc[rows, "panel"] == P).astype(float).to_numpy())
        names.append(f"panel_{P}")
    for c in covs:
        cols.append(Z[c].to_numpy())
        names.append(c)
    return pd.DataFrame(np.column_stack(cols), index=rows, columns=names)


def cov_beta(X, fit):
    mu = np.exp(np.clip(X @ fit["beta"], -30, 30))
    if fit["method"] == "nb":
        W = mu / (1 + fit["alpha"] * mu)
        return np.linalg.pinv((X.T * W) @ X)
    return fit["phi"] * np.linalg.pinv((X.T * mu) @ X)


def fit_model(B, ctrl, covs, panels, all_rows, loo=True):
    """Pooled NB GLM on ctrl rows; predictions (LOO for ctrl rows) and per-row variance of log RR."""
    panels = [P for P in panels if (B.loc[ctrl, "panel"] == P).any()]
    D = design(B, all_rows, covs, ctrl, panels)
    Xc, yc = D.loc[ctrl].to_numpy(), B.loc[ctrl, "y"].astype(float).to_numpy()
    fit = bs.fit_count_glm(Xc, yc, "nb")
    V = cov_beta(Xc, fit)
    X = D.to_numpy()
    eta = X @ fit["beta"]
    yhat = pd.Series(np.exp(np.clip(eta, -30, 30)), index=all_rows)
    vpred = pd.Series(np.einsum("ij,jk,ik->i", X, V, X), index=all_rows)
    how = pd.Series("full", index=all_rows)
    fits = []
    if loo:
        for i, b in enumerate(ctrl):
            keep = [k for k in range(len(ctrl)) if k != i]
            f = bs.fit_count_glm(Xc[keep], yc[keep], "nb", start=fit)
            xb = D.loc[b].to_numpy()
            yhat[b] = float(np.exp(np.clip(xb @ f["beta"], -30, 30)))
            how[b] = "loo"
            fits.append(f["method"])
    in_panel = B.loc[all_rows, "panel"].isin(panels)
    yhat[~in_panel] = np.nan
    y = B.loc[all_rows, "y"].astype(float)
    s = np.log((y + 0.5) / (yhat + 0.5))
    samp = (1 / (yhat + 0.5) + fit["alpha"]) if fit["method"] == "nb" else fit["phi"] / (yhat + 0.5)
    var = samp + vpred
    out = pd.DataFrame({"y": y, "yhat": yhat, "s": s, "var": var, "fit": how})
    coef = pd.DataFrame({"term": D.columns, "beta": fit["beta"], "se": np.sqrt(np.diag(V))})
    info = dict(covariates_used=[c for c in D.columns if not c.startswith(("intercept", "panel_"))], method=fit["method"], alpha=fit["alpha"], phi=fit["phi"], n_controls=len(ctrl),
                panels=panels, fallback=fit["fallback_reason"], n_loo_fallback=int(sum(m != "nb" for m in fits)),
                covariates=covs)
    return out, coef, info


# ----------------------------------------------------------------------------
# meta-analysis
# ----------------------------------------------------------------------------

def meta_mhk(y, v, alpha=0.05):
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    if k == 0:
        return {}
    if k == 1:
        return dict(k=1, mu=float(y[0]), se=float(math.sqrt(v[0])), ci_low=float(y[0] - 1.96 * math.sqrt(v[0])),
                    ci_high=float(y[0] + 1.96 * math.sqrt(v[0])), tau2=np.nan, q=np.nan, pi_low=np.nan, pi_high=np.nan,
                    p_one_sided=float(stats.norm.sf(y[0] / math.sqrt(v[0]))))
    re = bs.random_effects(y, v, method="REML", ci_method="z")
    tau2 = re["tau2"]
    w = 1 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    se = math.sqrt(max(1.0, q) / np.sum(w))
    t = stats.t.ppf(1 - alpha / 2, k - 1)
    pi = (np.nan, np.nan)
    if k >= 3:
        tp = stats.t.ppf(1 - alpha / 2, k - 2)
        h = tp * math.sqrt(tau2 + se ** 2)
        pi = (mu - h, mu + h)
    return dict(k=k, mu=mu, se=se, ci_low=mu - t * se, ci_high=mu + t * se, tau2=tau2, q=q, I2=re["I2"],
                pi_low=pi[0], pi_high=pi[1], p_one_sided=float(stats.t.sf(mu / se, k - 1)))


def rr(d):
    """add rate-ratio columns to a meta dict"""
    if not d:
        return d
    for a, b in (("mu", "RR"), ("ci_low", "RR_low"), ("ci_high", "RR_high"), ("pi_low", "PI_RR_low"), ("pi_high", "PI_RR_high")):
        d[b] = float(np.exp(d[a])) if d.get(a) is not None and np.isfinite(d[a]) else np.nan
    return d


# ----------------------------------------------------------------------------
# permutation p (P1)
# ----------------------------------------------------------------------------

def matched_controls(B, ctrl_rows, panel, depth):
    c = B.loc[ctrl_rows]
    c = c[c["panel"] == panel]
    ec = depth_class(depth)
    order = {"A": [["A"], ["B"], ["C"]], "B": [["B"], ["A", "C"]], "C": [["C"], ["B"], ["A"]]}[ec]
    use, widened = [], False
    for k, cls in enumerate(order):
        use += list(c.index[c["depth_class"].isin(cls)])
        if len(use) >= MIN_MATCH:
            break
        if k + 1 < len(order):
            widened = True
    return use, ec, widened


def perm_p(s_event, s_ctrl):
    s_ctrl = np.asarray(s_ctrl, float)
    s_ctrl = s_ctrl[np.isfinite(s_ctrl)]
    return (1 + int(np.sum(s_ctrl >= s_event))) / (len(s_ctrl) + 1), len(s_ctrl)


def smd(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 1 or len(b) < 2:
        return np.nan
    va = np.var(a, ddof=1) if len(a) > 1 else 0.0
    sd = math.sqrt((va + np.var(b, ddof=1)) / 2) if len(a) > 1 else np.std(b, ddof=1)
    return float((np.mean(a) - np.mean(b)) / sd) if sd > 0 else np.nan


# ----------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", default=None, help="testing only: write here instead of pooled/<tag>")
    ap.add_argument("--panels", default=None, help="comma list; default = all panels with event_table.tsv")
    a = ap.parse_args(argv)
    panels = a.panels.split(",") if a.panels else [R for R in PANELS if os.path.exists(os.path.join(OUT, R, "event_table.tsv"))]
    od = a.outdir or os.path.join(OUT, "pooled", a.tag)
    if os.path.exists(od):
        raise SystemExit(f"refusing to overwrite {od}")
    os.makedirs(od)
    B, E = load(panels)
    st = {R: pd.read_csv(paths(R)["species"], sep="\t").set_index("code") for R in panels}
    ctrl = list(B.index[B["control"].astype(bool)])
    internal_or_laddered = list(B.index[B["laddered"].astype(bool)])
    log(f"panels {panels}: {len(ctrl)} control branches; per panel {B.loc[ctrl].groupby('panel').size().to_dict()}")
    # ---- primary pooled model
    res, coef, info = fit_model(B, ctrl, COVS, panels, internal_or_laddered)
    coef.to_csv(os.path.join(od, "glm_coefficients.tsv"), sep="\t", index=False)
    br = B.loc[internal_or_laddered, ["panel", "branch_child", "branch_parent", "is_leaf", "depth", "depth_class", "n_desc", "bl",
                                     "frac_tx", "mean_busco", "dist_outgroup", "control", "tip_control", "inside_crown",
                                     "events", "hog_gains", "supported_gains", "seq_novel", "dups", "losses"]].join(res)
    br.to_csv(os.path.join(od, "branch_residuals.tsv"), sep="\t")
    sd_ctrl = br.loc[ctrl].groupby("panel")["s"].std(ddof=1)
    # ---- tip-control sensitivity model (panels with tip controls)
    tipc = list(B.index[B["tip_control"].astype(bool)])
    tip_panels = sorted(B.loc[tipc, "panel"].unique()) if tipc else []
    tres = None
    if len(tipc) >= 6:
        rows_t = [k for k in internal_or_laddered if B.at[k, "is_leaf"] and B.at[k, "panel"] in tip_panels]
        tres, tcoef, tinfo = fit_model(B, tipc, TIP_COVS, tip_panels, rows_t)
        tcoef.to_csv(os.path.join(od, "glm_coefficients_tipcontrols.tsv"), sep="\t", index=False)
        info["tip_model"] = tinfo
    # ---- events
    p1 = pd.read_csv(os.path.join(ROOT, "data", "r2_prereg", "p1_events.tsv"), sep="\t")
    p1_role = dict(zip(p1["event_id"], p1["role"]))
    erows = []
    for e in E.itertuples():
        k = f"{e.panel}:{e.M}"
        if k not in br.index:
            continue
        r = br.loc[k]
        sp = st[e.panel]
        mspecies = [x for x in str(pd.read_csv(paths(e.panel)["transitions"], sep="\t").set_index("event_id").at[e.event_id, "M_species"]).split(",")]
        msp = sp.reindex(mspecies)
        msp = msp[msp["mc_state"] == "complex"] if (msp["mc_state"] == "complex").any() else msp
        Dv = msp["D"].astype(float)
        Tv = msp["T"].astype(float)
        mc, ec, widened = matched_controls(B, ctrl, e.panel, int(r["depth"]))
        pp, nm = perm_p(r["s"], br.loc[mc, "s"]) if np.isfinite(r["s"]) else (np.nan, 0)
        row = dict(panel=e.panel, event_id=e.event_id, kind=e.kind, from_state=e.from_state, to_state=e.to_state,
                   confidence=e.confidence, M=e.M, is_tip=bool(r["is_leaf"]), depth=int(r["depth"]), depth_class=ec,
                   n_desc=int(r["n_desc"]), bl=r["bl"], hog_gains=r["hog_gains"], supported_gains=r["supported_gains"],
                   seq_novel=r["y"], yhat=r["yhat"], log_rr=r["s"], var=r["var"],
                   rr=np.exp(r["s"]), rr_low=np.exp(r["s"] - 1.96 * np.sqrt(r["var"])),
                   rr_high=np.exp(r["s"] + 1.96 * np.sqrt(r["var"])), dups=r["dups"], losses=r["losses"],
                   n_matched=nm, matched_widened=widened, p_perm=pp,
                   smd_depth_vs_matched=smd([r["depth"]], B.loc[mc, "depth"]),
                   ctrl_resid_sd_panel=sd_ctrl.get(e.panel, np.nan),
                   D_only=bool(Tv.mean() < 0.5), T_frac=float(Tv.mean()),
                   p1_role=p1_role.get(e.event_id, ""))
        if tres is not None and r["is_leaf"] and k in tres.index:
            tc = [x for x in tipc if B.at[x, "panel"] == e.panel]
            row.update(tip_log_rr=tres.at[k, "s"], tip_var=tres.at[k, "var"],
                       tip_p_perm=perm_p(tres.at[k, "s"], tres.loc[tc, "s"])[0], tip_n_controls=len(tc))
        erows.append(row)
    ER = pd.DataFrame(erows)
    ER.to_csv(os.path.join(od, "event_residuals.tsv"), sep="\t", index=False)
    # ---- P1
    out = dict(panels=panels, model=info, n_controls_by_panel=B.loc[ctrl].groupby("panel").size().to_dict(),
               depth_smd_all_events_vs_controls=smd(ER["depth"], B.loc[ctrl, "depth"]))
    conf = ER[ER["p1_role"] == "confirmatory"].copy()
    missing = sorted(set(p1.loc[p1["role"] == "confirmatory", "event_id"]) - set(conf["event_id"]))
    out["p1_mapping"] = {"R6_fucophycidae": "node Fucophycidae of the R6 tree (= MRCA of sampled Fucophycidae; "
                                            "checked in p1_mapping.txt)"}
    out["p1_missing_events"] = missing
    if len(conf):
        zs = stats.norm.isf(conf["p_perm"].to_numpy())
        Z = float(zs.sum() / math.sqrt(len(zs)))
        out["p1"] = dict(k=len(conf), events=conf["event_id"].tolist(), p_perm=conf["p_perm"].round(4).tolist(),
                         n_matched=conf["n_matched"].tolist(), stouffer_Z=Z, p_stouffer_one_sided=float(stats.norm.sf(Z)),
                         complete=not missing)
        out["p1_meta"] = rr(meta_mhk(conf["log_rr"], conf["var"]))
        sens = ER[ER["event_id"].isin(conf["event_id"].tolist() + ["R4_pleodorina"])]
        out["p1_meta_with_pleodorina"] = rr(meta_mhk(sens["log_rr"], sens["var"]))
        if "tip_log_rr" in ER.columns:
            ct = conf.copy()
            m = ct["tip_log_rr"].notna()
            ct.loc[m, "p_tipsens"] = ct.loc[m, "tip_p_perm"]
            ct.loc[~m, "p_tipsens"] = ct.loc[~m, "p_perm"]
            zt = stats.norm.isf(ct["p_tipsens"].to_numpy())
            out["p1_tip_sensitivity"] = dict(p=ct["p_tipsens"].round(4).tolist(),
                                             stouffer_Z=float(zt.sum() / math.sqrt(len(zt))),
                                             p_stouffer_one_sided=float(stats.norm.sf(zt.sum() / math.sqrt(len(zt)))))
            lv = np.where(m, ct["tip_log_rr"], ct["log_rr"])
            vv = np.where(m, ct["tip_var"], ct["var"])
            out["p1_tip_sensitivity_meta"] = rr(meta_mhk(lv, vv))
    # ---- secondary
    nested = ER[(ER["kind"] == "nested") & (ER["from_state"] == "simple") & (ER["to_state"] == "complex")
                & ER["panel"].str.split("_").str[0].isin(["R2", "R3", "R4", "R5", "R6", "R8"]) & ER["log_rr"].notna()]
    out["secondary_all_nested"] = rr(meta_mhk(nested["log_rr"], nested["var"]))
    out["secondary_all_nested"]["events"] = nested["event_id"].tolist()
    pp = []
    for P, g in nested.groupby("panel"):
        w = 1 / g["var"]
        pp.append((float((w * g["log_rr"]).sum() / w.sum()), float(1 / w.sum())))
    if pp:
        out["secondary_one_per_panel"] = rr(meta_mhk([x for x, _ in pp], [v for _, v in pp]))
    tonly = nested[~nested["D_only"]]
    out["secondary_T_only"] = rr(meta_mhk(tonly["log_rr"], tonly["var"]))
    out["secondary_T_only"]["events"] = tonly["event_id"].tolist()
    internal_nested = nested[~nested["is_tip"]]
    out["secondary_nested_internal_only"] = rr(meta_mhk(internal_nested["log_rr"], internal_nested["var"]))
    out["secondary_nested_internal_only"]["events"] = internal_nested["event_id"].tolist()
    # permutation-p Stouffer over all nested (descriptive)
    zs = stats.norm.isf(nested["p_perm"].to_numpy())
    out["secondary_all_nested_stouffer"] = dict(Z=float(zs.sum() / math.sqrt(len(zs))) if len(zs) else np.nan,
                                                p=float(stats.norm.sf(zs.sum() / math.sqrt(len(zs)))) if len(zs) else np.nan)
    orig = ER[ER["kind"].isin(["origin", "origin_complex"]) & ER["log_rr"].notna()]
    out["descriptive_origins"] = rr(meta_mhk(orig["log_rr"], orig["var"]))
    out["descriptive_origins"]["events"] = orig["event_id"].tolist()
    with open(os.path.join(od, "p1_and_meta.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda x: None if (isinstance(x, float) and not np.isfinite(x)) else str(x))
    write_done(od, [os.path.join(od, f) for f in os.listdir(od)], dict(panels=panels))
    log(json.dumps({k: out.get(k) for k in ("p1", "p1_meta")}, default=str))


if __name__ == "__main__":
    main()
