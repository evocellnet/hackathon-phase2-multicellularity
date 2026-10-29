#!/usr/bin/env python
"""W3 pooled: recurrent Pfam domains in sequence-novel and duplicated families across simple->complex events,
against a permutation null in which each event is replaced by a matched control branch of its own panel
(internal events: depth-class-matched internal controls, as P1; tip events: unicellular-state tip controls of the
same panel are not annotated, so tip events are listed descriptively and excluded from the calibrated test).

    w3_pooled.py --tag <pooled tag>      # reads pooled/<tag>/branch_residuals.tsv; writes pooled/<tag>/w3/
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, SEEDS, log, paths, write_done  # noqa: E402
from pooled import matched_controls, depth_class  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

N_PER_SEED = 200
N_RARE = 50


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    a = ap.parse_args(argv)
    pdir = os.path.join(OUT, "pooled", a.tag)
    od = os.path.join(pdir, "w3")
    if os.path.exists(od):
        raise SystemExit(f"refusing to overwrite {od}")
    br = pd.read_csv(os.path.join(pdir, "branch_residuals.tsv"), sep="\t", index_col=0)
    ER = pd.read_csv(os.path.join(pdir, "event_residuals.tsv"), sep="\t")
    panels = sorted(br["panel"].unique())
    W = []
    for R in panels:
        f = os.path.join(OUT, R, "w3_families.tsv.gz")
        if os.path.exists(f):
            W.append(pd.read_csv(f, sep="\t", keep_default_na=False))
    if not W:
        raise SystemExit("no W3 family tables")
    os.makedirs(od)
    W = pd.concat(W, ignore_index=True)
    W["key"] = W["panel"] + ":" + W["branch_child"]
    have = set(W["panel"])
    dom = {}   # key -> {"novel": set, "dup": set}
    for (k, kind), g in W[W["pfam"] != ""].groupby(["key", "kind"]):
        s = set()
        for x in g["pfam"]:
            s |= set(x.split(";"))
        dom.setdefault(k, {})[kind] = s

    def dset(k, which):
        d = dom.get(k, {})
        if which == "both":
            return d.get("seq_novel", set()) | d.get("dup", set())
        return d.get(which, set())

    sc = ER[(ER["from_state"] == "simple") & (ER["to_state"] == "complex") & (ER["kind"] == "nested")
            & ER["panel"].isin(have)].copy()
    sc["key"] = sc["panel"] + ":" + sc["M"]
    ctrl = list(br.index[br["control"].astype(bool)])
    tipc = list(br.index[br["tip_control"].astype(bool)]) if "tip_control" in br.columns else []
    pools, calibrated = {}, []
    for e in sc.itertuples():
        if e.is_tip:   # non-event tips carry no W3 annotation (w3_panel.py scope): tip events descriptive only
            continue
        pool, _, _ = matched_controls(br, ctrl, e.panel, int(e.depth))
        if len(pool) >= 3:
            pools[e.key] = pool
            calibrated.append(e.key)
    names = {}
    for acc, nm in zip(W.loc[W["pfam"] != "", "pfam"], W.loc[W["pfam"] != "", "pfam_names"]):
        names.update(zip(acc.split(";"), nm.split(";")))
    # ---- per-event summary (all simple->complex events, descriptive)
    rows = []
    for e in sc.itertuples():
        g = W[W["key"] == e.key]
        nv, dp = g[g["kind"] == "seq_novel"], g[g["kind"] == "dup"]
        topd = pd.Series([x for s in g.loc[g["pfam"] != "", "pfam"] for x in s.split(";")]).value_counts().head(8)
        cats = pd.Series([x for s in g.loc[g["categories"] != "", "categories"] for x in s.split(";")]).value_counts()
        rows.append(dict(panel=e.panel, event_id=e.event_id, M=e.M, is_tip=e.is_tip, calibrated=e.key in calibrated,
                         n_null_pool=len(pools.get(e.key, [])), seq_novel_families=len(nv),
                         seq_novel_with_pfam=int((nv["pfam"] != "").sum()), dup_units=len(dp),
                         dup_events=int(dp["n_dup_events"].astype(int).sum()), dup_with_pfam=int((dp["pfam"] != "").sum()),
                         top_domains="; ".join(f"{names.get(k, k)}({k}):{v}" for k, v in topd.items()),
                         categories="; ".join(f"{k}:{v}" for k, v in cats.items())))
    pd.DataFrame(rows).to_csv(os.path.join(od, "w3_event_summary.tsv"), sep="\t", index=False)
    # ---- recurrence vs permutation null (calibrated events only)
    out_rows, summ = [], []
    for which in ("both", "seq_novel", "dup"):
        obs = pd.Series([d for k in calibrated for d in dset(k, which)]).value_counts()
        draws_k = []
        stat_null = {2: [], 3: []}
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            for _ in range(N_PER_SEED):
                cnt = {}
                used = set()
                for k in calibrated:
                    pool = [x for x in pools[k] if x not in used] or pools[k]
                    c = pool[int(rng.integers(len(pool)))]
                    used.add(c)
                    for d in dset(c, which):
                        cnt[d] = cnt.get(d, 0) + 1
                s = pd.Series(cnt, dtype=float)
                draws_k.append(s)
                for t in (2, 3):
                    stat_null[t].append(int((s >= t).sum()))
        nd = len(draws_k)
        cand = obs[obs >= 2]
        ge = {d: 0 for d in cand.index}
        for s in draws_k:
            for d in cand.index:
                if s.get(d, 0) >= cand[d]:
                    ge[d] += 1
        for t in (2, 3):
            o = int((obs >= t).sum())
            arr = np.array(stat_null[t])
            summ.append(dict(units=which, threshold=t, n_events=len(calibrated), observed_recurrent_domains=o,
                             null_mean=float(arr.mean()), null_q95=float(np.quantile(arr, 0.95)),
                             p_perm=(1 + int((arr >= o).sum())) / (nd + 1), n_draws=nd))
        if len(cand):
            pv = np.array([(1 + ge[d]) / (nd + 1) for d in cand.index])
            q = bh(pv)
            for (d, kobs), p_, q_ in zip(cand.items(), pv, q):
                evs = [k for k in calibrated if d in dset(k, which)]
                out_rows.append(dict(units=which, pfam=d, name=names.get(d, d), k_events=int(kobs),
                                     events=";".join(ER.set_index(ER["panel"] + ":" + ER["M"]).loc[evs, "event_id"]),
                                     p_perm=p_, q_bh=q_))
    # ---- rarefied recurrence (size-adjusted): every branch contributes at most N_RARE annotated units, drawn at
    # random in each draw, for events and controls alike; paired draws -> p = P(R_null >= R_obs)
    units = {}
    for (k, kind), g in W[W["pfam"] != ""].groupby(["key", "kind"]):
        units.setdefault(k, {})[kind] = [set(x.split(";")) for x in g["pfam"]]

    def rare(k, which, rng):
        u = units.get(k, {})
        lst = (u.get("seq_novel", []) + u.get("dup", [])) if which == "both" else u.get(which, [])
        if len(lst) > N_RARE:
            lst = [lst[i] for i in rng.choice(len(lst), N_RARE, replace=False)]
        out = set()
        for x in lst:
            out |= x
        return out
    rrows = []
    for which in ("both", "seq_novel", "dup"):
        ro, rn = {2: [], 3: []}, {2: [], 3: []}
        for seed in SEEDS:
            rng = np.random.default_rng(seed + 1000)
            for _ in range(N_PER_SEED):
                co, cn, used = {}, {}, set()
                for k in calibrated:
                    for d in rare(k, which, rng):
                        co[d] = co.get(d, 0) + 1
                    pool = [x for x in pools[k] if x not in used] or pools[k]
                    c = pool[int(rng.integers(len(pool)))]
                    used.add(c)
                    for d in rare(c, which, rng):
                        cn[d] = cn.get(d, 0) + 1
                for t in (2, 3):
                    ro[t].append(sum(v >= t for v in co.values()))
                    rn[t].append(sum(v >= t for v in cn.values()))
        for t in (2, 3):
            a_, b_ = np.array(ro[t]), np.array(rn[t])
            rrows.append(dict(units=which, threshold=t, n_rare=N_RARE, n_events=len(calibrated),
                              observed_mean=float(a_.mean()), null_mean=float(b_.mean()),
                              p_paired=(1 + int((b_ >= a_).sum())) / (len(a_) + 1)))
    pd.DataFrame(rrows).to_csv(os.path.join(od, "w3_recurrence_rarefied.tsv"), sep="\t", index=False)
    pd.DataFrame(summ).to_csv(os.path.join(od, "w3_recurrence_summary.tsv"), sep="\t", index=False)
    rec = pd.DataFrame(out_rows).sort_values(["units", "k_events", "p_perm"], ascending=[True, False, True])
    rec.to_csv(os.path.join(od, "w3_recurrent_domains.tsv"), sep="\t", index=False)
    # matrix for the heatmap: domains recurrent in 'both' (top 40) x events; 1 novel, 2 dup, 3 both
    top = rec[rec["units"] == "both"].head(40)
    M = pd.DataFrame(0, index=top["pfam"], columns=[k for k in sc["key"]])
    for k in M.columns:
        for d in M.index:
            M.at[d, k] = (1 if d in dset(k, "seq_novel") else 0) + (2 if d in dset(k, "dup") else 0)
    M.columns = [ER.set_index(ER["panel"] + ":" + ER["M"]).at[k, "event_id"] for k in M.columns]
    M.insert(0, "name", [names.get(d, d) for d in M.index])
    M.to_csv(os.path.join(od, "w3_heatmap_matrix.tsv"), sep="\t")
    with open(os.path.join(od, "README.txt"), "w") as fh:
        fh.write(f"W3 recurrence over {len(calibrated)} calibrated simple->complex events ({', '.join(calibrated)}); "
                 f"null: {len(SEEDS)} seeds x {N_PER_SEED} draws, each event replaced by a matched control of its panel "
                 f"(distinct within a draw where possible). Units: seq_novel = sequence-novel gained families, dup = "
                 f"non-origin duplications (parent-level HOG), both = union. Domains after transposon mask.\n")
    write_done(od, [os.path.join(od, f) for f in os.listdir(od)])
    log(pd.DataFrame(summ).to_string())


if __name__ == "__main__":
    main()
