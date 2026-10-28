"""Positive-control scoring (data/positive_controls.tsv, pass_rule column; node mapping
data/positive_controls_nodes.tsv incl. pi_resolution; terms of data/prereg_definitions.tsv).

Every rule returns a :class:`Score` with ``recovered`` in
{'yes', 'no', 'not_evaluable', 'undetermined'}:

* ``not_evaluable`` -- the pass rule itself says so (e.g. 'loss on stem(C)' when C
  has < 2 genome-based species); counts as not recovered in the >= 11/16 criterion.
* ``undetermined`` -- an input the rule needs is not available yet (gain ladder,
  anchor protein, ...); re-run the scorer when it exists.

Primitive terms (prereg_definitions.tsv) are implemented on :class:`mcell.analysis.Cache`:
``PF-HOG(X)`` = :func:`mcell.analysis.pf_hogs`; ``present_at``, ``origin``,
``loss on stem`` = :class:`mcell.analysis.Lineage`; ``hit(X)`` = a protein with a
gathering-threshold hit in the pre-registered Pfam subset.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd

from . import analysis as an
from . import branchstats as bs

LEAK = 2


@dataclass
class Score:
    control: str
    recovered: str
    observed: str
    components: dict = field(default_factory=dict)
    notes: str = ""


@dataclass
class Ctx:
    cache: an.Cache
    species: pd.DataFrame              # species_final(_prok).tsv
    transitions: pd.DataFrame
    paths: dict
    ladder_gains: pd.DataFrame | None = None
    ladder_per_query: pd.DataFrame | None = None
    bt: pd.DataFrame | None = None     # branch table (C15)

    @property
    def L(self) -> an.Lineage:
        return self.cache.lineage

    def node(self, name: str) -> int:
        if name not in self.cache.tree.index:
            raise KeyError(f"node {name!r} not in the species tree")
        return self.cache.tree.index[name]

    def clade(self, name: str) -> set:
        return self.L.clade(name)

    def genome_species(self, sps: Iterable[str]) -> list:
        dt = self.species.set_index("code")["data_type_final" if "data_type_final" in self.species.columns
                                             else "data_type"].astype(str)
        return [s for s in sps if dt.get(s, "") == "genome"]

    def low_quality(self) -> set:
        if "low_quality" not in self.species.columns:
            return set()
        return set(self.species.loc[self.species["low_quality"] == 1, "code"])


# =============================================================================
# primitives
# =============================================================================

def hit_species(ctx: Ctx, pfams: Iterable[str]) -> pd.Series:
    """species -> number of proteins with hit(X) (any listed Pfam)."""
    h = ctx.cache.hits
    pf = {an._acc(x) for x in pfams}
    return h[h["pfam_acc"].isin(pf)].groupby("species")["protein"].nunique()


def leak(ctx: Ctx, pfams, crown: str) -> list:
    """Species outside(crown) with hit(X)."""
    hs = hit_species(ctx, pfams)
    inside = ctx.clade(crown)
    return sorted(s for s in hs.index if s not in inside and hs[s] > 0)


def gains_on_stem(ctx: Ctx, pfams, crown: str) -> list:
    """Rows of PF-HOG(X) whose origin(H) = crown(C) (gain on stem(C))."""
    L = ctx.L
    rows = an.pf_hogs(ctx.cache, pfams)
    X = ctx.node(crown)
    return [int(h) for h in rows if L.level[L.origin_row(h)] == X]


def ladder_status(ctx: Ctx, rows: Sequence[int]) -> str:
    """Sequence-novel status (ladder tier 3) of the families of these HOGs, if available."""
    if ctx.ladder_gains is None or not len(rows):
        return "gain ladder not available" if ctx.ladder_gains is None else ""
    fam = set(ctx.cache.hogs["family_id"].to_numpy()[list(rows)])
    g = ctx.ladder_gains[ctx.ladder_gains["family_id"].isin(fam)]
    if g.empty:
        return "families not in gain ladder"
    return f"ladder: {int(g['seq_novel'].sum())}/{len(g)} gained families sequence-novel"


def co_option_hogs(ctx: Ctx, pfams, at: str, groups: Sequence[set], origin_at_or_above: str | None = None) -> list:
    """PF-HOG(X) rows H with present_at(at), >= 1 member in every species set of
    ``groups`` and (optionally) origin(H) = ``origin_at_or_above`` or an ancestor."""
    L, T = ctx.L, ctx.cache.tree
    X = ctx.node(at)
    rows = an.pf_hogs(ctx.cache, pfams)
    if not len(rows):
        return []
    # necessary condition: the HOG's level is an ancestor of >= 1 species of every group
    end = np.asarray(T.end)
    for g in groups:
        gi = np.array(sorted(T.index[s] for s in g if s in T.index))
        if not len(gi):
            return []
        lv = L.level[np.asarray(rows, dtype=np.int64)]
        lo = np.searchsorted(gi, lv, side="left")
        hi = np.searchsorted(gi, end[lv], side="right")
        rows = [h for h, a, b in zip(rows, lo, hi) if b > a]
        if not rows:
            return []
    if origin_at_or_above is not None:
        O = ctx.node(origin_at_or_above)
        rows = [h for h in rows if T.is_anc(L.level[L.origin_row(h)], O)]
    rows = [h for h in rows if L.present_at(h, X)]
    if not rows:
        return []
    ms = L.member_species(rows)
    return [h for h in rows if all(set(ms[h]) & g for g in groups)]


def describe(ctx: Ctx, rows: Sequence[int], k: int = 3) -> str:
    """'n_families families: id (origin X), ...' -- HOG rows of one lineage at several levels share a family."""
    L = ctx.L
    ids = ctx.cache.hogs["hog_id"].to_numpy()
    fam = ctx.cache.hogs["family_id"].to_numpy()
    rows = list(rows)
    nf = len({fam[h] for h in rows})
    return f"{nf} famil{'y' if nf == 1 else 'ies'}: " + ", ".join(
        f"{ids[h]} (origin {L.origin_level(h)})" for h in rows[:k]) + (f" +{len(rows) - k} more HOGs" if len(rows) > k else "")


def yn(x: bool) -> str:
    return "yes" if x else "no"


# =============================================================================
# eukaryote controls
# =============================================================================

def rule_gain(ctx: Ctx, cid: str, pfam: str, crown: str, companion: Sequence[str] = ()) -> Score:
    g = gains_on_stem(ctx, [pfam], crown)
    lk = leak(ctx, [pfam], crown)
    a, b = bool(g), len(lk) <= LEAK
    obs = (f"(a) PF-HOG({pfam}) with origin on stem({crown}): {len(g)}"
           + (f" [{describe(ctx, g)}]" if g else "")
           + f"; (b) species outside {crown} with hit: {len(lk)}" + (f" ({', '.join(lk[:8])})" if lk else "")
           + f"; {ladder_status(ctx, g)}")
    for c in companion:
        gc = gains_on_stem(ctx, [c], crown)
        obs += f"; companion {c} (not required): {len(gc)} gain(s) on stem, {len(leak(ctx, [c], crown))} outside species"
    notes = f"outside species with hits (possible contamination/annotation FP): {', '.join(lk)}" if lk and b else ""
    return Score(cid, yn(a and b), obs, dict(a=a, b=b), notes)


def c04(ctx: Ctx) -> Score:
    filo = ctx.clade("Choanoflagellata") | ctx.clade("Filasterea")
    rows = co_option_hogs(ctx, ["PF00028"], "Choanozoa", [ctx.clade("Metazoa"), filo], origin_at_or_above="Choanozoa")
    return Score("C04", yn(bool(rows)),
                 f"PF-HOG(PF00028) present at Choanozoa with origin at/above Choanozoa and members in Metazoa and "
                 f"Choanoflagellata/Filasterea: {len(rows)}" + (f" [{describe(ctx, rows)}]" if rows else ""))


def c05(ctx: Ctx) -> Score:
    arch = an.architecture_proteins(ctx.cache.hits, "PF00017", "PF07714")
    grp = {"Metazoa": ctx.clade("Metazoa"), "Choanoflagellata": ctx.clade("Choanoflagellata"),
           "Filasterea": ctx.clade("Filasterea")}
    cnt = {k: int(arch["species"].isin(v).sum()) for k, v in grp.items()}
    a = all(v >= 1 for v in cnt.values())
    filo = grp["Choanoflagellata"] | grp["Filasterea"]
    rows = co_option_hogs(ctx, ["PF00102"], "Choanozoa", [grp["Metazoa"], filo])
    b = bool(rows)
    return Score("C05", yn(a and b),
                 f"(a) SH2+PK_Tyr_Ser-Thr proteins: " + ", ".join(f"{k} {v}" for k, v in cnt.items())
                 + f"; (b) PF-HOG(PF00102) present at Choanozoa with Metazoa + non-metazoan Filozoa members: {len(rows)}",
                 dict(a=a, b=b), "TK catalytic-loop motif check not implemented (secondary, descriptive only)")


def _loss_eval(ctx: Ctx, crown: str):
    n_gen = len(ctx.genome_species(ctx.clade(crown)))
    return n_gen >= 2, n_gen


def c06(ctx: Ctx) -> Score:
    L, T = ctx.L, ctx.cache.tree
    rows = an.pf_hogs(ctx.cache, ["PF00362"])
    O = ctx.node("Apusomonadida_Opisthokonta")
    rows = [h for h in rows if T.is_anc(L.level[L.origin_row(h)], O) and T.is_anc(L.level[h], ctx.node("Thectra"))]
    ms = L.member_species(rows) if rows else {}
    cand = [h for h in rows if "Thectra" in ms[h]]
    okF, nF = _loss_eval(ctx, "Fungi")
    okC, nC = _loss_eval(ctx, "Choanoflagellata")
    F, Cn = ctx.node("Fungi"), ctx.node("Choanoflagellata")
    lostF = [h for h in cand if L.lost_on_stem(h, F)]
    lostC = [h for h in cand if L.lost_on_stem(h, Cn)]
    both = [h for h in cand if h in set(lostF) and h in set(lostC)]
    lostHolo = [h for h in cand if L.lost_on_stem(h, ctx.node("Holomycota"))]
    obs = (f"(a) PF-HOG(PF00362) with origin at/above Apusomonadida_Opisthokonta and a Thecamonas member: {len(cand)}"
           + (f" [{describe(ctx, cand)}]" if cand else "")
           + f"; (b) lost on stem(Fungi): {len(lostF)}; (c) lost on stem(Choanoflagellata): {len(lostC)}; "
           f"same HOG (a)-(c): {len(both)}; reported: lost on stem(Holomycota): {len(lostHolo)}")
    if not (okF and okC):
        return Score("C06", "not_evaluable", obs, notes=f"genome species: Fungi {nF}, Choanoflagellata {nC} (<2)")
    notes = ""
    if cand and not lostF and lostHolo:
        notes = "integrin loss maps to the Holomycota stem (audit caveat in positive_controls_nodes.tsv)"
    return Score("C06", yn(bool(both)), obs, dict(a=bool(cand), b=bool(lostF), c=bool(lostC)), notes)


def c07(ctx: Ctx) -> Score:
    L = ctx.L
    nonmet = ctx.clade("Ichthyosporea") | ctx.clade("Filasterea")
    cand = co_option_hogs(ctx, ["PF00907"], "Holozoa", [nonmet, ctx.clade("Metazoa")])
    okC, nC = _loss_eval(ctx, "Choanoflagellata")
    Cn = ctx.node("Choanoflagellata")
    lost = [h for h in cand if L.lost_on_stem(h, Cn)]
    pf = set(ctx.cache.hogs["hog_id"].to_numpy()[an.pf_hogs(ctx.cache, ["PF00907"])])
    ev = ctx.cache.events
    d = ev[(ev["event"] == "dup") & (ev["branch_child"] == "Metazoa")]
    dm = d[d["hog_id"].isin(pf) | d["child_hog_ids"].map(lambda s: any(x in pf for x in str(s).split(";")))]
    c = len(dm) > 0
    fungi_members = sum(bool(set(m) & ctx.clade("Fungi")) for m in L.member_species(cand).values()) if cand else 0
    obs = (f"(a) PF-HOG(PF00907) present at Holozoa with Ichthyosporea/Filasterea + Metazoa members: {len(cand)}"
           + (f" [{describe(ctx, cand)}]" if cand else "")
           + f"; (b) of these lost on stem(Choanoflagellata): {len(lost)}; (c) duplication events on stem(Metazoa) in "
           f"PF-HOG(PF00907): {len(dm)} ({int(dm['at_origin'].sum()) if len(dm) else 0} at origin); reported: "
           f"candidates with fungal members {fungi_members}")
    if not okC:
        return Score("C07", "not_evaluable", obs, notes=f"Choanoflagellata genome species {nC} (<2)")
    return Score("C07", yn(bool(lost) and c), obs, dict(a=bool(cand), b=bool(lost), c=c))


def c08(ctx: Ctx) -> Score:
    g = [ctx.clade("Metazoa"), ctx.clade("Choanoflagellata")]
    a = co_option_hogs(ctx, ["PF01414"], "Choanozoa", g)
    b = co_option_hogs(ctx, ["PF06816", "PF07684"], "Choanozoa", g)
    return Score("C08", yn(bool(a) and bool(b)),
                 f"(a) DSL PF-HOGs present at Choanozoa with Metazoa + Choanoflagellata members: {len(a)}; "
                 f"(b) NOD/NODP PF-HOGs: {len(b)}", dict(a=bool(a), b=bool(b)), "each half reported separately")


def c09(ctx: Ctx) -> Score:
    pairs = [("Streptophyta", "Klebsormidiophyceae_Phragmoplastophyta"), ("Phaeophyceae_Schizocladia", "Phaeophyceae"),
             ("Rhodophytina", "Eurhodophytina")]
    lq = ctx.low_quality()
    res, parts, ne = {}, [], []
    for P, M in pairs:
        Mset = ctx.clade(M)
        Sset = ctx.clade(P) - Mset
        if not (Mset - lq) or not (Sset - lq):
            ne.append(P)
            parts.append(f"{P}: not_evaluable (no adequate-quality species on one side)")
            continue
        rows = co_option_hogs(ctx, ["PF05920"], P, [Mset])
        res[P] = bool(rows)
        parts.append(f"{P}->{M}: {len(rows)} PF-HOG(PF05920) present at P with a member in {M}")
    ok = bool(res) and all(res.values())
    notes = ""
    if ne:
        notes = f"scored on the remaining transitions; omitted: {', '.join(ne)}"
    lqn = {P: sorted((ctx.clade(P)) & lq) for P, _ in pairs}
    notes += ("; " if notes else "") + "low-quality species: " + "; ".join(f"{k}: {','.join(v) or '-'}" for k, v in lqn.items())
    return Score("C09", yn(ok), "; ".join(parts), res, notes)


def c10(ctx: Ctx) -> Score:
    crown = "Anydrophyta"
    g = gains_on_stem(ctx, ["PF03514"], crown)
    lk = leak(ctx, ["PF03514"], crown)
    a, b = bool(g), len(lk) <= LEAK
    obs = (f"(a) PF-HOG(PF03514) with origin on stem({crown}): {len(g)}" + (f" [{describe(ctx, g)}]" if g else "")
           + f"; (b) outside species with hit: {len(lk)}" + (f" ({', '.join(lk[:8])})" if lk else ""))
    if not (a and b):
        return Score("C10", "no", obs + "; (c) not needed", dict(a=a, b=b))
    fam = set(ctx.cache.hogs["family_id"].to_numpy()[g])
    pq = ctx.ladder_per_query
    if pq is None and ctx.ladder_gains is None:
        return Score("C10", "undetermined", obs + "; (c) gain ladder DIAMOND hits not available",
                     dict(a=a, b=b), "re-run when results/gainladder/euk/ exists")
    if pq is not None:
        q = pq[pq["family_id"].isin(fam)]
        doms = q["best_out_domain"].dropna().astype(str)
    else:
        q = ctx.ladder_gains[ctx.ladder_gains["family_id"].isin(fam)]
        doms = q["best_out_domain"].dropna().astype(str)
    if len(q) == 0:
        return Score("C10", "undetermined", obs + "; (c) gained families not found in the gain ladder", dict(a=a, b=b))
    if len(doms) == 0:
        return Score("C10", "yes", obs + "; (c) not_evaluable: no outside hit at e<=1e-5", dict(a=a, b=b, c="not_evaluable"),
                     "scored on (a) and (b) with flag (pass rule)")
    c = bool((doms == "prok").any())
    return Score("C10", yn(c), obs + f"; (c) best outside hit domains of representatives: {doms.value_counts().to_dict()}",
                 dict(a=a, b=b, c=c))


def c11(ctx: Ctx) -> Score:
    P = "Chlamydomonas_TGV"
    rows = co_option_hogs(ctx, ["PF01858", "PF01857"], P, [{"Chlarei"}, {"Volvcar"}], origin_at_or_above=P)
    tgv = [h for h in rows if ctx.L.present_at(h, ctx.node("TGV_clade"))]
    return Score("C11", yn(bool(rows)),
                 f"PF-HOG(PF01858|PF01857) with Chlarei + Volvcar members, present at {P}, origin at/above {P}: "
                 f"{len(rows)}" + (f" [{describe(ctx, rows)}]" if rows else "") + f"; reported: also present at TGV_clade: {len(tgv)}",
                 notes="P = Chlamydomonas_TGV (D7; single volvocine primary transition; pi_resolution)")


def c12(ctx: Ctx) -> Score:
    hs = hit_species(ctx, ["PF12499"])
    v, c = int(hs.get("Volvcar", 0)), int(hs.get("Chlarei", 0))
    return Score("C12", yn(v >= 2 * c and v - c >= 5), f"proteins with hit(PF12499): Volvcar {v}, Chlarei {c}")


def c13(ctx: Ctx, hex1: str = "Neurcra_P87252") -> Score:
    L, T = ctx.L, ctx.cache.tree
    G = ctx.cache.genes
    gi = np.flatnonzero(G["gene"].to_numpy() == hex1)
    if not len(gi):
        return Score("C13", "undetermined", f"anchor {hex1} (UniProt P87252, HEX1_NEUCR) not in the run",
                     notes="fallback A. nidulans HexA needs a reciprocal-best-hit search (not run)")
    hits = ctx.cache.hits
    nc = hits[(hits["species"] == "Neurcra") & (hits["pfam_acc"] == "PF01287")]
    eif = sorted(set(nc["protein"]) - {hex1})
    p = int(G["parent_idx"].to_numpy()[gi[0]])
    if p < 0:
        return Score("C13", "no", f"Hex1 ({hex1}) is a species-specific gene (no HOG)", notes=f"eIF5A candidates: {eif}")
    hhex = L.origin_row(p)
    mem = L.member_genes([hhex])[hhex]
    sps = L.member_species([hhex])[hhex]
    has_eif = bool(set(eif) & mem)
    pez, taph = ctx.clade("Pezizomycotina"), ctx.clade("Taphrinomycotina")
    usp = set(sps)
    out_sp = sorted(usp - pez)
    frac = len(usp & pez) / max(1, len(usp))
    org = L.origin_level(hhex)
    obs = (f"H_hex = {ctx.cache.hogs['hog_id'].iat[hhex]} (origin {org}, {len(mem)} genes, {len(usp)} species, "
           f"{frac:.2f} Pezizomycotina species, outside {len(out_sp)}); N. crassa PF01287 proteins other than Hex1: "
           f"{eif}; H_hex contains eIF5A candidate: {has_eif}")
    if not has_eif:
        ok = org == "Pezizomycotina" and len(out_sp) <= LEAK
        return Score("C13", yn(ok), obs + " -> branch (i)", dict(branch="i"),
                     "'>=90% Pezizomycotina species' read with the leak tolerance '(<=2 outside species)'")
    sacch = ctx.node("saccharomyceta")
    at_p = L.at_level(hhex, sacch)
    ids = set(ctx.cache.hogs["hog_id"].to_numpy()[at_p]) if at_p else set()
    ev = ctx.cache.events
    dups = ev[(ev["event"] == "dup") & (ev["branch_child"] == "Pezizomycotina") & ev["hog_id"].isin(ids)]
    cnt = pd.Series(sps).value_counts()
    med_p = float(np.median([cnt.get(s, 0) for s in sorted(pez)]))
    med_t = float(np.median([cnt.get(s, 0) for s in sorted(taph)]))
    ok = len(dups) > 0 and med_p >= 2 and med_t <= 1
    return Score("C13", yn(ok), obs + f" -> branch (ii): duplications on stem(Pezizomycotina) {int(dups['n_events'].sum()) if len(dups) else 0}; "
                 f"median copies Pezizomycotina {med_p}, Taphrinomycotina {med_t}", dict(branch="ii"),
                 "Hex1 identified by its UniProt accession in the proteome ids (P87252) instead of a reciprocal best hit; "
                 "eIF5A = the other N. crassa PF01287 protein(s)")


def c14(ctx: Ctx) -> Score:
    hs = hit_species(ctx, ["PF01185"])
    ag = sorted(ctx.clade("Agaricales"))
    non = sorted(ctx.clade("Dikarya") - set(ag))
    sac = sorted(ctx.clade("Saccharomycotina"))
    ma = float(np.median([hs.get(s, 0) for s in ag]))
    mn = float(np.median([hs.get(s, 0) for s in non]))
    sac_hit = [s for s in sac if hs.get(s, 0) > 0]
    a = ma >= max(3, 2 * mn)
    b = len(sac_hit) <= 1
    return Score("C14", yn(a and b), f"(a) median hit(PF01185) per species: Agaricales {ma} ({', '.join(f'{s}:{int(hs.get(s, 0))}' for s in ag)}), "
                 f"other Dikarya {mn} (n={len(non)}); (b) Saccharomycotina species with hits: {sac_hit or 'none'}",
                 dict(a=a, b=b), "Agaricomycetes = node Agaricales (pi_resolution)")


def c15(ctx: Ctx) -> Score:
    if ctx.bt is None:
        return Score("C15", "undetermined", "branch table not available")
    bt = ctx.bt
    stem = bs.stem_branches(bt).get("T_fungi_hyphae")
    if stem is None:
        return Score("C15", "undetermined", "T_fungi_hyphae stem not found")
    wgd = an.wgd_flags(ctx.cache, bt, ctx.paths.get("qc_flags"))
    excl = set(wgd.index[wgd["wgd_flag"]]) - set(bs.path_edges(bt, stem))
    parts, calls, ne = [], {}, []
    for cls in ("kinase_ptyr", "adhesion_ecm", "signalling"):
        comp = bs.compositional_statistic(ctx.cache.events, ctx.cache.hogs, ctx.cache.class_hogs(cls), bt,
                                          min_d=10, exclude=excl)
        if not np.isfinite(comp.at[stem, "s"]):
            ne.append(cls)
            parts.append(f"{cls}: {comp.at[stem, 'reason']}")
            continue
        pool = [b for b in bt.index[bt["in_U"]] if np.isfinite(comp.at[b, "s"])]
        m = bs.match_controls(bt, {"T_fungi_hyphae": stem}, pool=pool, extra=comp[["log_d"]])
        e = bs.effect_table(comp["s"], m, scale="difference").loc["T_fungi_hyphae"]
        oe = float(np.exp(e["est"]))
        exp_call = e["percentile"] >= 0.95 and oe >= 1.5
        calls[cls] = exp_call
        parts.append(f"{cls}: d={comp.at[stem, 'd']}, a={comp.at[stem, 'a']}, percentile {e['percentile']:.3f} vs k={m.k}, "
                     f"O/E (odds ratio) {oe:.2f} -> {'EXPANSION' if exp_call else 'no expansion'}")
    obs = f"stem {bs.path_edges(bt, stem)}; " + "; ".join(parts)
    if any(calls.values()):
        return Score("C15", "no", obs, calls)
    if ne:
        return Score("C15", "not_evaluable", obs, calls, "statistic undefined (d<10) for: " + ", ".join(ne))
    return Score("C15", "yes", obs, calls, "O/E = exp(adjusted-logit difference vs matched controls) (operational)")


def c16(ctx: Ctx) -> Score:
    arch = an.architecture_proteins(ctx.cache.hits, "PF00109", "PF00698")
    n = arch.groupby("species").size()
    dic = sorted(ctx.clade("Dictyostelia"))
    oth = sorted(ctx.clade("Amoebozoa") - set(dic))
    og = ctx.genome_species(oth)
    comp = og if len(og) >= 2 else oth
    md = float(np.median([n.get(s, 0) for s in dic]))
    mo = float(np.median([n.get(s, 0) for s in comp]))
    return Score("C16", yn(md >= 10 and md >= 2 * mo),
                 f"PKS (PF00109+PF00698) proteins per species: Dictyostelia median {md} "
                 f"({', '.join(f'{s}:{int(n.get(s, 0))}' for s in dic)}); non-Dictyostelia Amoebozoa genomes median {mo} "
                 f"({', '.join(f'{s}:{int(n.get(s, 0))}' for s in comp)})")


# =============================================================================
# prokaryote controls
# =============================================================================

def _prok_groups(ctx: Ctx):
    s = ctx.species
    cy = s[s["lineage_group"].astype(str).str.startswith("Bacteria-Cyanobacteria-")]
    roles = cy["role"].astype(str).str.split(";")
    fil = cy[roles.map(lambda r: "multicellular" in r)]
    het = fil[fil["transition_id"].astype(str).str.contains("T_cyano_heterocyst")]
    uni = cy[roles.map(lambda r: bool({"closest_unicellular_sister", "outgroup"} & set(r)))
             & (cy["multicellularity_type"].astype(str) == "none")]
    return list(fil["code"]), list(het["code"]), list(uni["code"])


def p1(ctx: Ctx) -> Score:
    fil, het, uni = _prok_groups(ctx)
    parts, passed = [], 0
    for pf, grp, gname in (("PF18460", fil, "filamentous"), ("PF28069", het, "heterocystous"), ("PF24301", fil, "filamentous")):
        hs = hit_species(ctx, [pf])
        ft = np.mean([hs.get(s, 0) > 0 for s in grp]) if grp else np.nan
        fu = np.mean([hs.get(s, 0) > 0 for s in uni]) if uni else np.nan
        ok = bool(ft >= 0.8 and fu <= 0.2)
        passed += ok
        parts.append(f"{pf}: {gname} {ft:.2f} (n={len(grp)}), unicellular {fu:.2f} (n={len(uni)}) -> {'pass' if ok else 'fail'}")
    return Score("P1", yn(passed >= 2), f"{passed}/3 Pfams pass; " + "; ".join(parts),
                 notes=f"unicellular group: {', '.join(uni)}")


def p2(ctx: Ctx) -> Score:
    s = ctx.species
    sp = s[s["transition_id"].astype(str).str.contains("T_streptomyces")
           & (s["lineage_group"].astype(str) == "Bacteria-Actinomycetota-Streptomycetales")]["code"].tolist()
    if len(sp) < 3:
        return Score("P2", "not_evaluable", f"only {len(sp)} species ({sp})")
    parts, oks = [], []
    for pf in ("PF25848", "PF21179"):
        hs = hit_species(ctx, [pf])
        f = np.mean([hs.get(x, 0) > 0 for x in sp])
        oks.append(f >= 0.8)
        parts.append(f"{pf}: {f:.2f} of {len(sp)} species ({', '.join(f'{x}:{int(hs.get(x, 0))}' for x in sp)})")
    return Score("P2", yn(all(oks)), "; ".join(parts))


EUK_RULES: dict[str, Callable[[Ctx], Score]] = {
    "C01": lambda c: rule_gain(c, "C01", "PF00110", "Metazoa"),
    "C02": lambda c: rule_gain(c, "C02", "PF00019", "Metazoa", companion=["PF00688"]),
    "C03": lambda c: rule_gain(c, "C03", "PF01049", "Metazoa"),
    "C04": c04, "C05": c05, "C06": c06, "C07": c07, "C08": c08, "C09": c09, "C10": c10, "C11": c11,
    "C12": c12, "C13": c13, "C14": c14, "C15": c15, "C16": c16,
}
PROK_RULES: dict[str, Callable[[Ctx], Score]] = {"P1": p1, "P2": p2}


def score_all(ctx: Ctx, controls: pd.DataFrame, rules: dict) -> pd.DataFrame:
    rows = []
    for _, r in controls.iterrows():
        cid = r["control_id"]
        fn = rules.get(cid)
        try:
            s = fn(ctx) if fn else Score(cid, "undetermined", "no rule implemented")
        except KeyError as e:
            s = Score(cid, "undetermined", f"missing node/species: {e}")
        rows.append(dict(control=cid, expected=f"{r['expected_event']}: {r['expected_node_description']}",
                         confidence=r["confidence"], recovered=s.recovered, observed=s.observed,
                         components=";".join(f"{k}={v}" for k, v in s.components.items()), notes=s.notes))
    return pd.DataFrame(rows)
