"""Re-score round-1 positive controls (mcell.controls) on a dense panel (task step 5).

The round-1 rules are used unchanged; node names absent from a panel tree are mapped through
ALIASES (recorded in the output) or the control is reported 'not_applicable'. Pfam hits = the
pre-registered subset (round-1 files for round-1 species, pyhmmer GA on the same subset HMM file
for new panel species). Ancestra-event controls (gain/loss/duplication placement) are flagged.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from common import ROOT, paths, guard, log, kingdom, base
from mcell import analysis as an
from mcell import controls as pc

#: panel -> controls that refer to clades sampled in the panel
RELEVANT = {"R1": ["C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08"],
            "R2": ["C13", "C14"], "R3": ["C09", "C10"], "R4": ["C11", "C12"], "R5": ["C09"], "R6": ["C09"],
            "R7": ["C16"], "R8": ["P1"]}
ANCESTRA_EVENT = {"C01", "C02", "C03", "C06", "C07", "C10", "C13"}
#: round-1 node name -> panel node name (only where the round-1 name does not exist in the panel)
ALIASES = {"R1": {"Fungi": "Holomycota", "Apusomonadida_Opisthokonta": "Eukaryota"},
           "R2": {}, "R3": {}, "R4": {}, "R5": {}, "R6": {}, "R7": {}, "R8": {}}


class PanelCtx(pc.Ctx):
    alias: dict = {}

    def node(self, name):
        return super().node(self.alias.get(name, name))

    def clade(self, name):
        return self.L.clade(self.alias.get(name, name))


def c09_panel(ctx):
    """C09 restricted to the transition pair(s) present in the panel."""
    pairs = [("Streptophyta", "Klebsormidiophyceae_Phragmoplastophyta"), ("Phaeophyceae_Schizocladia", "Phaeophyceae"),
             ("Rhodophytina", "Eurhodophytina")]
    T = ctx.cache.tree
    pairs = [(P, M) for P, M in pairs if P in T.index and M in T.index]
    lq = ctx.low_quality()
    res, parts = {}, []
    for P, M in pairs:
        Mset = ctx.clade(M)
        Sset = ctx.clade(P) - Mset
        if not (Mset - lq) or not (Sset - lq):
            parts.append(f"{P}: not_evaluable")
            continue
        rows = pc.co_option_hogs(ctx, ["PF05920"], P, [Mset])
        res[P] = bool(rows)
        parts.append(f"{P}->{M}: {len(rows)} PF-HOG(PF05920) present at P with a member in {M}")
    if not res:
        return pc.Score("C09", "not_evaluable", "; ".join(parts) or "no C09 pair in panel")
    return pc.Score("C09", pc.yn(all(res.values())), "; ".join(parts), res)


def p1_panel(ctx, st):
    """P1 with the panel's own state coding: filamentous = simple or complex, heterocystous = Nostocales clade,
    unicellular = mc_state unicellular."""
    fil = list(st.loc[st["mc_state"] != "unicellular", "code"])
    het = sorted(ctx.clade("Nostocales"))
    uni = list(st.loc[st["mc_state"] == "unicellular", "code"])
    parts, passed = [], 0
    for pf, grp, gname in (("PF18460", fil, "filamentous"), ("PF28069", het, "heterocystous"), ("PF24301", fil, "filamentous")):
        hs = pc.hit_species(ctx, [pf])
        ft = np.mean([hs.get(s, 0) > 0 for s in grp])
        fu = np.mean([hs.get(s, 0) > 0 for s in uni])
        ok = bool(ft >= 0.8 and fu <= 0.2)
        passed += ok
        parts.append(f"{pf}: {gname} {ft:.2f} (n={len(grp)}), unicellular {fu:.2f} (n={len(uni)}) -> {'pass' if ok else 'fail'}")
    return pc.Score("P1", pc.yn(passed >= 2), f"{passed}/3 Pfams pass; " + "; ".join(parts))


def cmd_controls(a):
    R = a.R
    p = paths(R)
    dst = os.path.join(p["out"], "positive_controls.tsv")
    guard(dst)
    cp = dict(orthoxml=p["orthoxml"], tree=p["tree"], cache=p["cache"],
              toolkit=os.path.join(ROOT, "data", "toolkit_pfam.tsv"), h12_sets=os.path.join(ROOT, "data", "h12_sets.tsv"),
              pfam=os.path.join(p["work"], "pfam_prereg"), transitions=os.path.join(ROOT, "data", "transitions.tsv"))
    if not os.path.exists(os.path.join(p["cache"], "DONE.json")):
        an.build_cache(R, paths=cp, levels="all", log=log)
    cache = an.load_cache(R, paths=cp)
    st = pd.read_csv(p["species"], sep="\t")
    sp = st.rename(columns={})[["code", "data_type", "busco_C", "mc_state"]].copy()
    sp["low_quality"] = (sp["busco_C"] < 50).astype(int)
    cls = pd.read_parquet(os.path.join(p["ladder"], "classified.parquet"))
    pq = pd.read_parquet(os.path.join(p["ladder"], "per_query.parquet"))
    PanelCtx.alias = ALIASES[base(R)]
    ctx = PanelCtx(cache=cache, species=sp, transitions=pd.read_csv(p["transitions"], sep="\t"), paths=cp,
                   ladder_gains=cls, ladder_per_query=pq)
    ref = pd.read_csv(os.path.join(ROOT, "results", "positive_controls", kingdom(R), "controls.tsv"), sep="\t")
    ref = ref.set_index("control")
    rows = []
    for cid in RELEVANT[base(R)]:
        try:
            if cid == "C09":
                s = c09_panel(ctx)
            elif cid == "P1":
                s = p1_panel(ctx, st)
            else:
                s = pc.EUK_RULES[cid](ctx)
        except KeyError as e:
            s = pc.Score(cid, "not_applicable", f"missing node/species in panel: {e}")
        rows.append(dict(panel=R, control=cid, ancestra_event_control=cid in ANCESTRA_EVENT,
                         round1_recovered=ref.loc[cid, "recovered"] if cid in ref.index else "",
                         panel_recovered=s.recovered, observed=s.observed,
                         components=";".join(f"{k}={v}" for k, v in s.components.items()), notes=s.notes,
                         aliases=";".join(f"{k}->{v}" for k, v in ALIASES[base(R)].items())))
    out = pd.DataFrame(rows)
    out.to_csv(dst, sep="\t", index=False)
    log(f"wrote {dst}: " + ", ".join(f"{r.control}:{r.panel_recovered}" for r in out.itertuples()))
