#!/usr/bin/env python
"""W4 cell-death follow-up, step 2 (descriptive, exploratory): copy-number tables per species and per lineage,
which families carry the multicellular-minus-sister excess in each lineage, and the domain architectures of the
proteins that carry each cell-death domain.

Lineages = the 10 eukaryotic events used by the W4 Pfam-wide scans (results/r2/w4/scan/events_used.tsv: 7 primary
origins + 3 nested simple->complex events), each multicellular (MC) clade vs its sister clade from
data/transitions.tsv. All data types are used, as in the scans (transcriptome species flagged in the tables).

Outputs (results/r2/w4_celldeath/):
  copy_number_species.tsv        187 species x 22 families (+ state, data type, proteome size, primary-148 flag)
  copy_number_by_state.tsv       mean copies per family in unicellular / simple / complex species (148 primary set)
  lineage_family_summary.tsv     event x family: mean copies MC vs sister, per 10k proteins, W4 scan lfc/q
  lineage_drivers.tsv            event x family: share of the total MC-sister cell-death excess (raw and size-scaled)
  architectures_top.tsv          event x group x family: 5 most frequent simplified architectures
  architectures_top_global.tsv   family: 10 most frequent architectures over all 187 proteomes, with lineage spread
  stand_classes.tsv              event x group: NB-ARC/NACHT (STAND) proteins by N-terminal and C-terminal class
  stand_excess.tsv               event: MC-sister excess of STAND proteins split by architecture class
"""
import os
import re
import sys

import numpy as np
import pandas as pd

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = os.path.join(P, "results/r2/w4_celldeath")
CD = pd.read_csv(os.path.join(P, "data/r2_sets/cell_death.tsv"), sep="\t")
SHORT = {"PF00656": "Caspase", "PF31181": "Metacaspase", "PF01650": "Legumain (VPE)", "PF00452": "Bcl-2",
         "PF02180": "BH4", "PF00653": "BIR", "PF00619": "CARD", "PF00531": "Death", "PF01335": "DED",
         "PF00931": "NB-ARC", "PF05729": "NACHT", "PF06985": "HET", "PF14479": "HeLo", "PF04598": "Gasdermin",
         "PF26164": "Gasdermin (bact./fungal)", "PF01823": "MACPF", "PF01582": "TIR", "PF13676": "TIR_2",
         "PF01027": "Bax inhibitor-1", "PF06943": "zf-LSD1", "PF14721": "AIF_C", "PF02452": "MazF"}
STAND = {"PF00931", "PF05729"}


def out(name):
    path = os.path.join(OUT, name)
    if os.path.exists(path):
        sys.exit(f"REFUSE: {path} exists")
    return path


# ---------------- architecture simplification ----------------
def dom_class(name):
    base = re.sub(r"\(x\d+\)$", "", name)
    if base.startswith("LRR") or base in ("LRR_RPS2",):
        return "LRR"
    if base.startswith("TPR") or base in ("Sel1", "TPR_MalT"):
        return "TPR"
    if "WD40" in base or base.startswith("prop_"):
        return "WD40"
    if base.startswith("Ank"):
        return "ANK"
    if base.startswith("HEAT") or base.startswith("ARM"):
        return "HEAT/ARM"
    if base.startswith("WHD_") or base in ("APAF1_C", "NLRC4_HD2", "NACHT_sigma", "WHD_APAF1"):
        return None   # winged-helix/helical STAND sub-domains: part of the NB-ARC/NACHT module, dropped
    if base.startswith("AAA") and base != "AAA":
        return "AAA-like"
    return base


def simplify(arch):
    """Ordered domain classes, repeats collapsed to one token (LRR, TPR, WD40, ANK)."""
    toks = [dom_class(t) for t in arch.split(" | ")]
    toks = [t for t in toks if t]
    res = []
    for t in toks:
        if not res or res[-1] != t:
            res.append(t)
    return "–".join(res)


NT_CLASS = {"TIR": "TIR", "TIR_2": "TIR", "Rx_N": "CC (Rx_N)", "RPW8": "RPW8", "HET": "HET", "HeLo": "HeLo",
            "HeLo_like": "HeLo", "CARD": "CARD/Death/DED/PYD", "Death": "CARD/Death/DED/PYD",
            "DED": "CARD/Death/DED/PYD", "PYRIN": "CARD/Death/DED/PYD", "Patatin": "Patatin",
            "PNP_UDP_1": "PNP-UDP phosphorylase", "BIR": "BIR", "Pkinase": "kinase", "PK_Tyr_Ser-Thr": "kinase",
            "NPHP3_N": "NPHP3_N", "Goodbye": "Goodbye", "SesB": "SesB"}
REPEATS = ("LRR", "TPR", "WD40", "ANK", "HEAT/ARM")


def stand_class(arch):
    toks = [dom_class(t) for t in arch.split(" | ")]
    toks = [t for t in toks if t]
    idx = [i for i, t in enumerate(toks) if t in ("NB-ARC", "NACHT")]
    core = "NB-ARC" if "NB-ARC" in toks else "NACHT"
    i0, i1 = idx[0], idx[-1]
    nterm = [NT_CLASS.get(t, "other") for t in toks[:i0]]
    nt = next((n for n in nterm if n != "other"), "other" if nterm else "none")
    cts = [t for t in toks[i1 + 1:] if t in REPEATS]
    ct = "+".join(sorted(set(cts))) if cts else ("other" if toks[i1 + 1:] else "none")
    nlr_like = bool(cts)
    return core, nt, ct, nlr_like


def main():
    tr = pd.read_csv(os.path.join(P, "data/species_traits.tsv"), sep="\t")
    tr = tr[tr.kingdom == "euk"].set_index("code")
    sp = pd.read_csv(os.path.join(P, "results/r2/w4/matrix/species.tsv"), sep="\t").set_index("code")
    # primary set as in W4 / verification: genome-based proteomes with BUSCO C >= 50% (148)
    prim = set(tr[(tr.data_type == "genome") & (tr.busco_C >= 50)].index)
    prot = pd.read_csv(os.path.join(OUT, "celldeath_proteins.tsv.gz"), sep="\t")
    prot["simple"] = prot.arch.map(simplify)
    ev = pd.read_csv(os.path.join(P, "results/r2/w4/scan/events_used.tsv"), sep="\t")
    ex = prot.assign(fam=prot.cd_accs.str.split(";")).explode("fam")

    # ---------- 1. species x family copy numbers ----------
    cn = ex.groupby(["species", "fam"]).size().unstack(fill_value=0).reindex(
        index=tr.index, columns=CD.pfam, fill_value=0)
    cn.columns = [f"{a}|{SHORT[a]}" for a in cn.columns]
    meta = pd.DataFrame({"species": tr.species, "mc_state": tr.mc_state, "multicellular": tr.multicellular,
                         "data_type": tr.data_type, "busco_C": tr.busco_C, "n_proteins": sp.n_proteins,
                         "primary148": [c in prim for c in tr.index]})
    tab = meta.join(cn)
    tab["total_celldeath_proteins"] = ex.groupby("species").protein.nunique().reindex(tr.index).fillna(0).astype(int)
    tab.to_csv(out("copy_number_species.tsv"), sep="\t", index_label="code")
    assert len(prim) == 148, len(prim)

    # by state (148 primary): unicellular includes reversals and grey zone coded 0 (round-1 binary coding)
    t148 = tab[tab.primary148].copy()
    t148["state3"] = np.where(t148.multicellular == 0, "unicellular/coded0", t148.mc_state)
    famcols = list(cn.columns)
    rows = []
    for st, g in t148.groupby("state3"):
        for f in famcols + ["total_celldeath_proteins"]:
            rows.append(dict(state=st, n_species=len(g), family=f, mean_copies=g[f].mean(),
                             median_copies=g[f].median(), prevalence=(g[f] > 0).mean(),
                             mean_per_10k_proteins=(g[f] / g.n_proteins * 1e4).mean()))
    pd.DataFrame(rows).to_csv(out("copy_number_by_state.tsv"), sep="\t", index=False, float_format="%.4g")

    # ---------- 2. lineage (event) summaries ----------
    ashr = {}
    for e in ev.event:
        a = pd.read_csv(os.path.join(P, f"results/r2/w4/scan/ashr_{e}.tsv"), sep="\t").set_index("pfam")
        s = pd.read_csv(os.path.join(P, f"results/r2/w4/scan/event_{e}.tsv"), sep="\t").set_index("pfam")
        ashr[e] = s.join(a)
    rows, drv = [], []
    for r in ev.itertuples():
        mc, sis = r.mc_species.split(","), r.sister_species.split(",")
        A, B = tab.loc[mc], tab.loc[sis]
        for acc, f in zip(CD.pfam, famcols):
            sc = ashr[r.event].loc[acc] if acc in ashr[r.event].index else None
            ma, mb = A[f].mean(), B[f].mean()
            na, nb = (A[f] / A.n_proteins * 1e4).mean(), (B[f] / B.n_proteins * 1e4).mean()
            rows.append(dict(event=r.event, label=r.label, tier=r.tier, n_mc=len(mc), n_sister=len(sis),
                             pfam=acc, family=SHORT[acc], mean_mc=ma, mean_sister=mb, diff=ma - mb,
                             log2_ratio=np.log2((ma + 0.5) / (mb + 0.5)), per10k_mc=na, per10k_sister=nb,
                             log2_ratio_per10k=np.log2((na + 0.05) / (nb + 0.05)),
                             prev_mc=(A[f] > 0).mean(), prev_sister=(B[f] > 0).mean(),
                             scan_lfc_ash=None if sc is None else sc.lfc_ash,
                             scan_q=None if sc is None else sc.q_bh,
                             scan_sig=None if sc is None else bool(sc.q_bh < 0.1)))
        d = pd.DataFrame([x for x in rows if x["event"] == r.event])
        pos = d["diff"].clip(lower=0).sum()
        size_ratio = A.n_proteins.mean() / B.n_proteins.mean()
        for x in d.itertuples():
            # size-scaled excess: MC mean minus sister mean scaled up to MC proteome size
            exp_mc = x.per10k_sister * A.n_proteins.mean() / 1e4
            drv.append(dict(event=r.event, label=r.label, pfam=x.pfam, family=x.family, diff=x.diff,
                            share_of_positive_excess=x.diff / pos if pos > 0 else np.nan,
                            total_mc=A.total_celldeath_proteins.mean(), total_sister=B.total_celldeath_proteins.mean(),
                            proteome_mc=A.n_proteins.mean(), proteome_sister=B.n_proteins.mean(),
                            proteome_ratio=size_ratio, size_scaled_excess=x.mean_mc - exp_mc))
    lfs = pd.DataFrame(rows)
    lfs.to_csv(out("lineage_family_summary.tsv"), sep="\t", index=False, float_format="%.4g")
    dr = pd.DataFrame(drv)
    dr["share_of_positive_size_scaled_excess"] = dr.groupby("event").size_scaled_excess.transform(
        lambda s: s / s.clip(lower=0).sum())
    dr.sort_values(["event", "diff"], ascending=[True, False]).to_csv(
        out("lineage_drivers.tsv"), sep="\t", index=False, float_format="%.4g")

    # ---------- 3. architectures ----------
    grp_rows, st_rows, se_rows = [], [], []
    stand = prot[prot.cd_accs.str.contains("PF00931|PF05729")].copy()
    sc = stand.arch.map(stand_class)
    stand["core"], stand["nterm"], stand["cterm"], stand["nlr_like"] = zip(*sc)
    stand["stand_class"] = stand.nterm + " / " + stand.core + " / " + stand.cterm
    stand.to_csv(os.path.join(OUT, "stand_proteins.tsv.gz"), sep="\t", index=False, compression="gzip")
    for r in ev.itertuples():
        for grp, codes in (("MC", r.mc_species.split(",")), ("sister", r.sister_species.split(","))):
            sub = ex[ex.species.isin(codes)]
            for acc in CD.pfam:
                s = sub[sub.fam == acc]
                if s.empty:
                    continue
                vc = s.simple.value_counts()
                for rank, (a, n) in enumerate(vc.head(5).items(), 1):
                    grp_rows.append(dict(event=r.event, group=grp, n_species=len(codes), pfam=acc,
                                         family=SHORT[acc], n_proteins_family=len(s), rank=rank, architecture=a,
                                         n_proteins=n, frac=n / len(s), per_species=n / len(codes),
                                         n_species_with=s[s.simple == a].species.nunique()))
            ss = stand[stand.species.isin(codes)]
            for (core, nt, ct, nl), g in ss.groupby(["core", "nterm", "cterm", "nlr_like"]):
                st_rows.append(dict(event=r.event, group=grp, n_species=len(codes), core=core, nterm=nt, cterm=ct,
                                    nlr_like=nl, n_proteins=len(g), per_species=len(g) / len(codes),
                                    n_species_with=g.species.nunique()))
        mc, sis = r.mc_species.split(","), r.sister_species.split(",")
        sm, ssis = stand[stand.species.isin(mc)], stand[stand.species.isin(sis)]
        cls = lambda d: np.where(d.nlr_like, "repeat-containing (NLR-like)",  # noqa: E731
                                 np.where(d.nterm != "none", "effector domain, no repeats", "STAND only / other"))
        a = pd.Series(cls(sm)).value_counts() / len(mc)
        b = pd.Series(cls(ssis)).value_counts() / len(sis)
        for k in sorted(set(a.index) | set(b.index)):
            se_rows.append(dict(event=r.event, label=r.label, stand_type=k, per_species_mc=a.get(k, 0),
                                per_species_sister=b.get(k, 0), diff=a.get(k, 0) - b.get(k, 0)))
        se_rows.append(dict(event=r.event, label=r.label, stand_type="TOTAL", per_species_mc=len(sm) / len(mc),
                            per_species_sister=len(ssis) / len(sis), diff=len(sm) / len(mc) - len(ssis) / len(sis)))
    pd.DataFrame(grp_rows).to_csv(out("architectures_top.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(st_rows).sort_values(["event", "group", "n_proteins"], ascending=[True, True, False]).to_csv(
        out("stand_classes.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(se_rows).to_csv(out("stand_excess.tsv"), sep="\t", index=False, float_format="%.4g")

    # global top architectures with species / state spread
    g_rows = []
    ex2 = ex.join(tab[["mc_state", "multicellular"]], on="species")
    for acc in CD.pfam:
        s = ex2[ex2.fam == acc]
        for a, n in s.simple.value_counts().head(10).items():
            t = s[s.simple == a]
            g_rows.append(dict(pfam=acc, family=SHORT[acc], n_proteins_family=len(s), architecture=a, n_proteins=n,
                               frac=n / len(s), n_species=t.species.nunique(),
                               n_species_mc=t[t.multicellular == 1].species.nunique(),
                               n_species_coded0=t[t.multicellular == 0].species.nunique(),
                               top_species=";".join(t.species.value_counts().head(3).index)))
    pd.DataFrame(g_rows).to_csv(out("architectures_top_global.tsv"), sep="\t", index=False, float_format="%.4g")
    print("done")


if __name__ == "__main__":
    main()
