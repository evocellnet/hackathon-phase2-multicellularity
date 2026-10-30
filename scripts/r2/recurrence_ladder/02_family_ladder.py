"""D17 recurrence ladder, family level (round-2 W4 Pfam scan; transposon-masked; BH q < 0.1 within event).
Null: within each event, the event's significant families (per direction) are re-assigned at random among the families
tested in that event (family shuffle preserving per-event hit numbers); events independent. B = 10,000, seed 11.
Per-family recurrence p: exact Poisson-binomial with per-event hit probability k_e / n_e (tested events only).
Outputs: family_ladder.tsv, family_recurrence_all.tsv.gz, recurrent_families_k2.tsv, set_families_ladder.tsv
"""
import gzip
import os
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import P, R, R1_ORIGINS, R1_NESTED, LABEL, out, poibin_pmf, bh  # noqa: E402

B, SEED, QCUT = 10000, 11, 0.1
SCAN = f"{R}/r2/w4/scan"
EVENTS = R1_ORIGINS + list(R1_NESTED)

# ---- per-event significance matrix ----
fam = None
tested, sig = {}, {}
tags = {}
for e in EVENTS:
    d = pd.read_csv(f"{SCAN}/event_{e}.tsv", sep="\t")
    s = np.where(d.q_bh < QCUT, np.sign(d.direction), 0).astype(int)
    tested[e] = pd.Series(1, index=d.pfam)
    sig[e] = pd.Series(s, index=d.pfam)
allf = sorted(set().union(*[set(v.index) for v in tested.values()]))
T = pd.DataFrame({e: tested[e].reindex(allf).fillna(0).astype(int) for e in EVENTS}, index=allf)
S = pd.DataFrame({e: sig[e].reindex(allf).fillna(0).astype(int) for e in EVENTS}, index=allf)
# cross-check against the published recurrence table
obs = pd.read_csv(f"{SCAN}/recurrence_observed.tsv", sep="\t").set_index("pfam")
chk = {e: int((obs[f"sig_{e}"].reindex(allf).fillna(0).astype(int) != S[e]).sum()) for e in EVENTS}
assert sum(chk.values()) == 0, chk
F = len(allf)
print(f"{F} families; per-event mismatches vs recurrence_observed.tsv: {chk}", flush=True)

UNIT = {e: e for e in R1_ORIGINS}
UNIT.update(R1_NESTED)
LADDERS = {"7 origins": R1_ORIGINS, "3 nested": list(R1_NESTED),
           "7 origins collapsed (origin OR nested)": EVENTS}

rng = np.random.default_rng(SEED)
fam_rows = pd.DataFrame(index=allf)
tier_rows = []
null_store = {}
for dname, sgn in (("expanded", 1), ("contracted", -1)):
    H = (S == sgn).astype(int)
    # null hits per event: B x F boolean per event, generated in chunks
    null_hits = {}
    for e in EVENTS:
        idx = np.where(T[e].to_numpy() == 1)[0]
        k = int(H[e].sum())
        M = np.zeros((B, F), dtype=bool)
        if k > 0:
            for b0 in range(0, B, 500):
                nb = min(500, B - b0)
                r = rng.random((nb, len(idx)))
                pick = np.argpartition(r, k - 1, axis=1)[:, :k]
                rows_ = np.repeat(np.arange(b0, b0 + nb), k)
                M[rows_, idx[pick.ravel()]] = True
        null_hits[e] = M
    for lname, evs in LADDERS.items():
        units = sorted(set(UNIT[e] for e in evs)) if "collapsed" in lname else evs
        # observed recurrence per family
        if "collapsed" in lname:
            ob = np.zeros(F, int)
            nl = np.zeros((B, F), np.int8)
            pis = np.zeros((F, len(units)))
            for j, u in enumerate(units):
                es = [e for e in evs if UNIT[e] == u]
                ob += H[es].max(axis=1).to_numpy()
                hit = np.zeros((B, F), bool)
                miss = np.ones(F)
                for e in es:
                    hit |= null_hits[e]
                    miss *= 1 - T[e].to_numpy() * H[e].sum() / max(T[e].sum(), 1)
                nl += hit
                pis[:, j] = 1 - miss
        else:
            ob = H[evs].sum(axis=1).to_numpy()
            nl = np.zeros((B, F), np.int8)
            for e in evs:
                nl += null_hits[e]
            pis = np.column_stack([T[e].to_numpy() * H[e].sum() / max(T[e].sum(), 1) for e in evs])
        # exact per-family p
        pfam_p = np.array([poibin_pmf(pis[i])[ob[i]:].sum() if ob[i] > 0 else 1.0 for i in range(F)])
        pfam_exp = pis.sum(1)
        tag = f"{lname}|{dname}"
        fam_rows[f"r[{tag}]"] = ob
        fam_rows[f"p[{tag}]"] = pfam_p
        fam_rows[f"q[{tag}]"] = bh(pfam_p)
        fam_rows[f"E[{tag}]"] = pfam_exp
        for k in range(1, len(units) + 1):
            n_null = (nl >= k).sum(1)
            o = int((ob >= k).sum())
            # how many families at >= k reach BH q < 0.1
            tier_rows.append(dict(ladder=lname, direction=dname, tier_k=k, n_units=len(units), observed_n_families=o,
                                  null_mean=n_null.mean(), null_q025=np.quantile(n_null, .025),
                                  null_q975=np.quantile(n_null, .975), p_perm=((n_null >= o).sum() + 1) / (B + 1),
                                  obs_over_expected=o / n_null.mean() if n_null.mean() > 0 else np.nan,
                                  excess=o - n_null.mean(),
                                  n_families_q01=int(((ob >= k) & (fam_rows[f"q[{tag}]"].to_numpy() < 0.1)).sum()),
                                  B=B, seed=SEED, source="recomputed from results/r2/w4/scan/event_T_*.tsv"))
        null_store[tag] = (nl, ob)
    del null_hits

TR = pd.DataFrame(tier_rows)
# Holm across the 2 directions x 7 tiers of the primary ladder is not meaningful (nested tiers); report BH over tiers as info
pd.set_option("display.width", 220)
print(TR.drop(columns=["source", "B", "seed"]).to_string(index=False))

# ---- annotations ----
cl = pd.read_csv(f"{P}/data/raw/pfam/Pfam-A.clans.tsv.gz", sep="\t", header=None,
                 names=["pfam", "clan", "clan_name", "name", "description"]).set_index("pfam")
go = {}
with open(f"{P}/data/raw/go/pfam2go") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        pf = line.split()[0].replace("Pfam:", "")
        term = line.split(" > ")[1].strip().replace("GO:", "", 1)
        go.setdefault(pf, []).append(term)
sets = {}
for fn, nm in (("elm_motif_binding.tsv", "motif-binding"), ("cell_death.tsv", "cell-death"),
               ("tf_celltype.tsv", "cell-type TF")):
    s = pd.read_csv(f"{P}/data/r2_sets/{fn}", sep="\t")
    if "set" in s.columns:
        s = s[s.set == "primary"]
    for pf in s.pfam:
        sets.setdefault(pf, []).append(nm)
tk = pd.read_csv(f"{P}/data/toolkit_pfam.tsv", sep="\t")
tkcol = [c for c in tk.columns if c.lower() in ("pfam", "pfam_acc", "accession")][0]
catcol = [c for c in tk.columns if c != tkcol][0]
toolkit = tk.groupby(tkcol)[catcol].apply(lambda x: ";".join(sorted(set(map(str, x))))).to_dict()

fam_rows["name"] = cl.name.reindex(allf).to_numpy()
fam_rows["description"] = cl.description.reindex(allf).to_numpy()
fam_rows["go_terms"] = [" | ".join(go.get(f, [])[:6]) for f in allf]
fam_rows["frozen_sets"] = [";".join(sets.get(f, [])) for f in allf]
fam_rows["round1_toolkit"] = [toolkit.get(f, "") for f in allf]
fam_rows["n_events_tested"] = T.sum(1).to_numpy()
for dname, sgn in (("expanded", 1), ("contracted", -1)):
    fam_rows[f"origins_{dname}"] = [";".join(LABEL[e] for e in R1_ORIGINS if S.at[f, e] == sgn) for f in allf]
    fam_rows[f"nested_{dname}"] = [";".join(LABEL[e] for e in R1_NESTED if S.at[f, e] == sgn) for f in allf]
fam_rows.index.name = "pfam"
fam_rows.reset_index().to_csv(out("family_recurrence_all.tsv.gz"), sep="\t", index=False, compression="gzip")

# ---- recurrent families (>= 2 of the 7 origins), per direction ----
rec = []
for dname in ("expanded", "contracted"):
    tag = f"7 origins|{dname}"
    sub = fam_rows[fam_rows[f"r[{tag}]"] >= 2]
    for pf, x in sub.iterrows():
        rec.append(dict(pfam=pf, direction=dname, n_origins=int(x[f"r[{tag}]"]),
                        expected_origins=round(x[f"E[{tag}]"], 3), recurrence_p=x[f"p[{tag}]"],
                        recurrence_q_bh=x[f"q[{tag}]"], exceeds_chance_q01=int(x[f"q[{tag}]"] < 0.1),
                        n_origins_tested=int(T.loc[pf, R1_ORIGINS].sum()), origins=x[f"origins_{dname}"],
                        nested_same_direction=x[f"nested_{dname}"],
                        opposite_direction_origins=x["origins_" + ("contracted" if dname == "expanded" else "expanded")],
                        name=x["name"], description=x.description, frozen_sets=x.frozen_sets,
                        round1_toolkit=x.round1_toolkit, go_terms=x.go_terms))
REC = pd.DataFrame(rec).sort_values(["direction", "n_origins", "recurrence_p"], ascending=[False, False, True])
REC.to_csv(out("recurrent_families_k2.tsv"), sep="\t", index=False)
print(REC.groupby(["direction", "n_origins"]).agg(n=("pfam", "size"), n_q01=("exceeds_chance_q01", "sum")))

# ---- frozen sets: set families significant in >= 2 origins, vs the same null ----
srows = []
for setname in ("motif-binding", "cell-death", "cell-type TF"):
    mem = np.array([setname in sets.get(f, []) for f in allf])
    for dname in ("expanded", "contracted"):
        for lname in ("7 origins", "7 origins collapsed (origin OR nested)"):
            nl, ob = null_store[f"{lname}|{dname}"]
            for k in (1, 2, 3):
                o = int(((ob >= k) & mem).sum())
                n_null = (nl[:, mem] >= k).sum(1)
                fams = [f"{fam_rows.at[f, 'name']}({int(ob[i])}: {fam_rows.at[f, 'origins_' + dname]}"
                        + (f" + {fam_rows.at[f, 'nested_' + dname]}" if "collapsed" in lname and
                           fam_rows.at[f, 'nested_' + dname] else "") + ")"
                        for i, f in enumerate(allf) if mem[i] and ob[i] >= k] if k >= 2 else []
                srows.append(dict(set=setname, n_set_families_scanned=int(mem.sum()), direction=dname, ladder=lname,
                                  tier_k=k, observed=o, null_mean=n_null.mean(), p_perm=((n_null >= o).sum() + 1) / (B + 1),
                                  families=";".join(fams)))
SR = pd.DataFrame(srows)
SR.to_csv(out("set_families_ladder.tsv"), sep="\t", index=False)
print(SR[SR.ladder == "7 origins"].drop(columns=["families"]).to_string(index=False))

# ---- quoted existing ladders (not recomputed) ----
q = []
cg = pd.read_csv(f"{R}/exploratory/euk/pfam_wide/convergence_global.tsv", sep="\t")
for _, x in cg[cg.scheme == "permute"].iterrows():
    q.append(dict(ladder=f"round-1 Pfam-wide convergence ({x.flag}; Mann-Whitney flags; branch-label permutation)",
                  direction=x.flag, tier_k=x.k_min, n_units=7, observed_n_families=x.n_domains_obs, null_mean=x.null_mean,
                  null_q025=np.nan, null_q975=x.null_q95, p_perm=x.p_perm,
                  source="results/exploratory/euk/pfam_wide/convergence_global.tsv (null_q975 column = null q95)"))
w3 = pd.read_csv(f"{R}/r2/panels/pooled/final_8panels/w3/w3_recurrence_summary.tsv", sep="\t")
for _, x in w3.iterrows():
    q.append(dict(ladder=f"panel simple->complex events, domains gained/duplicated ({x.units}; event-label null)",
                  direction="gained_or_duplicated" if x.units == "both" else "gained", tier_k=x.threshold,
                  n_units=x.n_events, observed_n_families=x.observed_recurrent_domains, null_mean=x.null_mean,
                  null_q025=np.nan, null_q975=x.null_q95, p_perm=x.p_perm,
                  source="results/r2/panels/pooled/final_8panels/w3/w3_recurrence_summary.tsv (null_q975 = q95)"))
w3r = pd.read_csv(f"{R}/r2/panels/pooled/final_8panels/w3/w3_recurrence_rarefied.tsv", sep="\t")
for _, x in w3r.iterrows():
    q.append(dict(ladder=f"panel simple->complex events, RAREFIED ({x.units}; paired null)", direction=x.units,
                  tier_k=x.threshold, n_units=x.n_events, observed_n_families=x.observed_mean, null_mean=x.null_mean,
                  null_q025=np.nan, null_q975=np.nan, p_perm=x.p_paired,
                  source="results/r2/panels/pooled/final_8panels/w3/w3_recurrence_rarefied.tsv"))
FL = pd.concat([TR, pd.DataFrame(q)], ignore_index=True)
FL.to_csv(out("family_ladder.tsv"), sep="\t", index=False)
