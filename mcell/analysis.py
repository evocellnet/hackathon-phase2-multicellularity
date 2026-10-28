"""Shared layer of the analysis scripts (scripts/analysis/*.py).

1. **Cache** -- :func:`build_cache` parses an Ancestra run once (:mod:`mcell.hogs`)
   and writes parquet tables to ``runs/analysis/<run>/``: branch summary, HOGs,
   genes, per-(branch, HOG) events, HOG membership and per-species copy number
   at chosen levels (integer rows into hogs/genes), the pre-registered Pfam hits, gene- and HOG-level class
   annotation and per-branch class composition of gains / duplications / losses.
   :func:`load_cache` reads it back (:class:`Cache`).
2. **Pfam classes** -- :class:`ClassSpec` (from data/toolkit_pfam.tsv and
   data/h12_sets.tsv), :func:`protein_flags`, :func:`hog_class_table`.
   A HOG is in class C if >= 50 % of its member genes carry a class domain
   (DESIGN §4.6). Classes: ``tf_dbd`` (toolkit tf_dbd; C2H2 is its own category
   ``c2h2`` and is NOT part of tf_dbd), ``adhesion_ecm`` (members, plus the
   conditional LRR members LRR_1/LRR_8 only on proteins without an NB-ARC, NACHT
   or TIR hit), ``signalling``, ``adh_sig`` (protein carries an adhesion_ecm or a
   signalling domain; then the 50 % rule; H5), ``kinase_ptyr``, ``pcd``,
   ``cell_cycle``; flagellar modules ``mod_<module>`` from h12_sets.tsv
   (role == member; Dynein_heavy and Kinesin are excluded there).
3. **Lineage queries** -- :class:`Lineage` (origin of a HOG, presence at a node,
   members) for the positive-control rules (data/prereg_definitions.tsv).
4. **Branch statistics** -- thin wrappers that feed :mod:`mcell.branchstats`:
   :func:`h1h2_counts` (sequence-novel gains from the gain ladder),
   :func:`wgd_flags`, :func:`h12_module_hogs`, :func:`run_relabel_suite`.

Choices where the specification is silent are listed in results/README_run_primary.md.
"""
from __future__ import annotations

import glob
import gzip
import hashlib
import json
import math
import os
import time
import warnings
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import branchstats as bs
from .hogs import HogResult, SpeciesTree, load_species_tree, parse_hogs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#: class columns used for compositions (in this order)
CLASS_COLS = ("tf_dbd", "adhesion_ecm", "signalling", "adh_sig", "kinase_ptyr", "pcd", "cell_cycle", "c2h2")
#: the four flagellar modules of H12 (DESIGN §4.6)
H12_MODULES = ("IFT-B", "BBSome", "axoneme", "transition_zone")
MIN_FRAC = 0.5


def _mod_col(m: str) -> str:
    return "mod_" + m.replace("-", "_")


# =============================================================================
# paths
# =============================================================================

def run_paths(run: str, root: str = ROOT) -> dict:
    """Standard locations for run ``run`` ('euk' or 'prok')."""
    rd = os.path.join(root, "runs", f"{run}_fast")
    ox = [p for p in (os.path.join(rd, "hogs.orthoxml.gz"), os.path.join(rd, "hogs.orthoxml")) if os.path.exists(p)]
    return dict(
        run=run, run_dir=rd, orthoxml=ox[0] if ox else os.path.join(rd, "hogs.orthoxml"),
        tree=os.path.join(rd, "stages", "input", "guide_tree.nwk"),
        guide=os.path.join(root, "data", f"guide_tree_{run}.nwk"),
        bl=os.path.join(root, "data", f"species_tree_bl_{run}.nwk"),
        cache=os.path.join(root, "runs", "analysis", run),
        pfam=os.path.join(root, "data", "annotations", "pfam_prereg", run),
        pfam_full=os.path.join(root, "data", "annotations", "pfam", run),
        ladder=os.path.join(root, "results", "gainladder", run, "branch_gain_ladder.tsv"),
        ladder_gains=os.path.join(root, "results", "gainladder", run, "gains_classified.tsv.gz"),
        ladder_per_query=os.path.join(root, "runs", "gainladder", run, "per_query.parquet"),
        ladder_reps=os.path.join(root, "runs", "gainladder", run, "reps.parquet"),
        transitions=os.path.join(root, "data", "transitions.tsv"),
        traits=os.path.join(root, "data", "species_traits.tsv"),
        species=os.path.join(root, "data", "species_final.tsv" if run == "euk" else "species_final_prok.tsv"),
        toolkit=os.path.join(root, "data", "toolkit_pfam.tsv"),
        h12_sets=os.path.join(root, "data", "h12_sets.tsv"),
        h12_dir=os.path.join(root, "data", "h12_expected_direction_frozen.tsv"),
        qc_flags=os.path.join(root, "results", "qc", "flags.tsv"),
        registry=os.path.join(root, "results", "analysis_registry.tsv"),
    )


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# =============================================================================
# Pfam classes
# =============================================================================

def _acc(x: str) -> str:
    return str(x).split(".")[0]


@dataclass
class ClassSpec:
    """Pre-registered Pfam class definitions.

    ``members[c]``: Pfam accessions that put a protein in class c;
    ``conditional[c]``: accessions that count only if the protein carries none of
    ``exclusion`` (LRR rule of toolkit_pfam.tsv); ``modules[m][component]``:
    accessions of each component of a flagellar module (h12_sets.tsv; the three
    BBS2 domain models form one component 'BBS2', every other Pfam is its own
    component)."""
    members: dict
    conditional: dict
    exclusion: set
    modules: dict

    def all_pfams(self) -> set:
        out = set(self.exclusion)
        for d in (self.members, self.conditional):
            for v in d.values():
                out |= set(v)
        for comps in self.modules.values():
            for v in comps.values():
                out |= set(v)
        return out


def load_class_spec(toolkit: str, h12_sets: str) -> ClassSpec:
    tk = pd.read_csv(toolkit, sep="\t")
    tk["pfam_acc"] = tk["pfam_acc"].map(_acc)
    cat_map = {"c2h2_separate": "c2h2"}
    members, conditional = {}, {}
    for _, r in tk.iterrows():
        c = cat_map.get(r["category"], r["category"])
        if r["role"] == "member":
            members.setdefault(c, set()).add(r["pfam_acc"])
        elif r["role"] == "conditional_member":
            conditional.setdefault(c, set()).add(r["pfam_acc"])
    exclusion = set(tk.loc[tk["role"] == "exclusion_marker", "pfam_acc"])
    h = pd.read_csv(h12_sets, sep="\t")
    h = h[h["role"] == "member"]
    modules = {}
    for _, r in h.iterrows():
        comp = "BBS2" if str(r["pfam_name"]).startswith("BBS2_") else _acc(r["pfam_acc"])
        modules.setdefault(r["module"], {}).setdefault(comp, set()).add(_acc(r["pfam_acc"]))
    return ClassSpec(members=members, conditional=conditional, exclusion=exclusion, modules=modules)


def read_pfam_dir(d: str, species: Iterable[str] | None = None, pfams: set | None = None) -> pd.DataFrame:
    """Hits (protein, species, pfam_acc, env_from, dom_score) from ``<d>/<code>.pfam.tsv.gz``.

    ``species`` restricts the files read; missing species are listed in ``.attrs['missing']``."""
    files = {os.path.basename(f).split(".pfam")[0]: f for f in glob.glob(os.path.join(d, "*.pfam.tsv.gz"))}
    want = sorted(files) if species is None else list(species)
    rows, missing = [], []
    for sp in want:
        f = files.get(sp)
        if f is None:
            missing.append(sp)
            continue
        t = pd.read_csv(f, sep="\t", usecols=lambda c: c in ("protein", "pfam_acc", "env_from", "dom_score"))
        if t.empty:
            continue
        t["pfam_acc"] = t["pfam_acc"].map(_acc)
        if pfams is not None:
            t = t[t["pfam_acc"].isin(pfams)]
        t["species"] = sp
        rows.append(t)
    out = (pd.concat(rows, ignore_index=True) if rows else
           pd.DataFrame(columns=["protein", "pfam_acc", "env_from", "dom_score", "species"]))
    out.attrs["missing"] = missing
    out.attrs["n_files"] = len(want) - len(missing)
    return out


def protein_flags(hits: pd.DataFrame, spec: ClassSpec) -> pd.DataFrame:
    """Bool table indexed by protein: class columns (:data:`CLASS_COLS`), one
    ``mod_<module>`` column per module (any member domain) and one
    ``comp|<module>|<component>`` column per module component."""
    by = hits.groupby("pfam_acc")["protein"].agg(lambda x: set(x)).to_dict()

    def has(pf):
        s = set()
        for p in pf:
            s |= by.get(p, set())
        return s

    prots = pd.Index(sorted(set(hits["protein"])), name="protein")
    excl = has(spec.exclusion)
    sets = {}
    for c in ("tf_dbd", "signalling", "kinase_ptyr", "pcd", "cell_cycle", "c2h2"):
        sets[c] = has(spec.members.get(c, ()))
    sets["adhesion_ecm"] = has(spec.members.get("adhesion_ecm", ())) | (has(spec.conditional.get("adhesion_ecm", ())) - excl)
    sets["adh_sig"] = sets["adhesion_ecm"] | sets["signalling"]
    out = pd.DataFrame(index=prots)
    for c in CLASS_COLS:
        out[c] = prots.isin(sets[c])
    for m, comps in spec.modules.items():
        allp = set()
        for k, pf in comps.items():
            s = has(pf)
            allp |= s
            out[f"comp|{m}|{k}"] = prots.isin(s)
        out[_mod_col(m)] = prots.isin(allp)
    return out


def _depth_of_hogs(hogs: pd.DataFrame, tree: SpeciesTree) -> np.ndarray:
    return np.asarray(tree.depth, dtype=np.int32)[hogs["level_idx"].to_numpy()]


def accumulate(hogs: pd.DataFrame, genes: pd.DataFrame, w: np.ndarray, tree: SpeciesTree) -> np.ndarray:
    """Per-HOG sum of a gene weight over all member genes (any level).

    ``w`` has one value per row of ``genes``. Genes add to their parent HOG
    (``genes.parent_idx``); HOG sums are then pushed to parent HOGs from the
    deepest level upwards (a parent HOG is always at a shallower level)."""
    n = len(hogs)
    w = np.asarray(w, dtype=np.float64)
    gp = genes["parent_idx"].to_numpy()
    ok = gp >= 0
    acc = np.bincount(gp[ok], weights=w[ok], minlength=n)
    par = hogs["parent_idx"].to_numpy()
    dep = _depth_of_hogs(hogs, tree)
    for d in sorted(np.unique(dep), reverse=True):
        rows = np.flatnonzero((dep == d) & (par >= 0))
        if len(rows):
            acc += np.bincount(par[rows], weights=acc[rows], minlength=n)
    return acc


def hog_class_table(hogs: pd.DataFrame, genes: pd.DataFrame, gflags: pd.DataFrame, tree: SpeciesTree,
                    min_frac: float = MIN_FRAC) -> pd.DataFrame:
    """Per HOG (rows of ``hogs``): n_genes_check, n_<class> member counts,
    ``in_<class>`` (>= ``min_frac`` of members), ``mod_<module>`` (module HOG)
    and ``components`` (';'-joined '<module>|<component>' carried by >= min_frac
    of the members of module HOGs)."""
    out = pd.DataFrame({"hog_id": hogs["hog_id"].to_numpy(), "level": hogs["level"].astype(str).to_numpy()})
    ng = accumulate(hogs, genes, np.ones(len(genes)), tree)
    out["n_genes"] = ng.astype(np.int64)
    if not np.array_equal(out["n_genes"].to_numpy(), hogs["n_genes"].to_numpy()):
        raise AssertionError("accumulated member counts differ from hogs.n_genes (bug)")
    for c in CLASS_COLS:
        a = accumulate(hogs, genes, gflags[c].to_numpy(), tree)
        out[f"n_{c}"] = a.astype(np.int32)
        out[f"in_{c}"] = a >= min_frac * ng
    comp_cols = [c for c in gflags.columns if c.startswith("comp|")]
    comps: dict = {}
    for m in sorted({c.split("|")[1] for c in comp_cols}):
        mc = _mod_col(m)
        a = accumulate(hogs, genes, gflags[mc].to_numpy(), tree)
        inm = a >= min_frac * ng
        out[mc] = inm
        idx = np.flatnonzero(inm)
        if not len(idx):
            continue
        for c in [c for c in comp_cols if c.split("|")[1] == m]:
            a = accumulate(hogs, genes, gflags[c].to_numpy(), tree)
            for i in idx[a[idx] >= min_frac * ng[idx]]:
                comps.setdefault(int(i), []).append(c[5:])
    col = np.full(len(hogs), "", dtype=object)
    for i, v in comps.items():
        col[i] = ";".join(v)
    out["components"] = col
    return out


# =============================================================================
# cache
# =============================================================================

def _auto_levels(tree: SpeciesTree, transitions: pd.DataFrame, kingdom: str, extra: Iterable[str] = ()) -> list:
    want = set(extra)
    t = transitions[transitions["kingdom"] == kingdom]
    for _, r in t.iterrows():
        for col in ("M", "P"):
            v = r.get(col)
            if isinstance(v, str):
                want |= {x.strip() for x in v.split(",")}
        sp = r.get("stem_path")
        if isinstance(sp, str):
            want |= {x.strip() for x in sp.split(">")}
    return [tree.names[i] for i in tree.internal if tree.names[i] in want]


CONTROL_LEVELS = ("Metazoa", "Choanozoa", "Filozoa", "Holozoa", "Apusomonadida_Opisthokonta", "Opisthokonta",
                  "Fungi", "Holomycota", "Choanoflagellata", "Anydrophyta", "Coleochaetophyceae_Anydrophyta",
                  "Chlamydomonas_TGV", "TGV_clade", "Gonium_Volvocaceae", "Pezizomycotina", "saccharomyceta",
                  "Agaricales", "Eumycota", "Blastocladiomycota_Amastigomycota", "Amastigomycota", "Streptophyta",
                  "Phaeophyceae_Schizocladia", "Rhodophytina", "Eurhodophytina", "Streptomycetaceae", "Nostocales")


def build_cache(run: str, paths: dict | None = None, levels: str = "auto", log=print) -> dict:
    """Parse the run once and write the cache (see module doc). Returns file paths."""
    p = paths or run_paths(run)
    out = p["cache"]
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    res = parse_hogs(p["orthoxml"], newick=p["tree"] if os.path.exists(p["tree"]) else None,
                     progress=lambda n: log(f"{n:,} families"))
    bad = res.check_conservation()
    if not bad.empty:
        raise SystemExit(f"conservation law violated on {len(bad)} branches (parser bug)")
    log(f"parsed {res.stats['n_genes']:,} genes, {len(res.hogs):,} HOGs in {time.time() - t0:.0f}s")
    files = res.write(out, fmt="parquet")
    T = res.tree
    # ---- Pfam (pre-registered subset) and classes
    spec = load_class_spec(p["toolkit"], p["h12_sets"])
    species = [T.names[s] for s in T.leaves]
    hits = read_pfam_dir(p["pfam"], species)
    hits = hits[hits["protein"].isin(set(res.genes["gene"]))].reset_index(drop=True)
    hits.to_parquet(os.path.join(out, "pfam_hits.parquet"), index=False)
    files["pfam_hits"] = os.path.join(out, "pfam_hits.parquet")
    pf = protein_flags(hits[hits["pfam_acc"].isin(spec.all_pfams())], spec)
    gflags = pf.reindex(res.genes["gene"].to_numpy(), fill_value=False).astype(bool).reset_index(drop=True)
    gcl = pd.concat([res.genes[["gene", "species"]].reset_index(drop=True),
                     gflags[[c for c in gflags.columns if not c.startswith("comp|")]]], axis=1)
    gcl.to_parquet(os.path.join(out, "gene_classes.parquet"), index=False)
    hc = hog_class_table(res.hogs, res.genes, gflags, T)
    hc.to_parquet(os.path.join(out, "hog_classes.parquet"), index=False)
    files.update(gene_classes=os.path.join(out, "gene_classes.parquet"),
                 hog_classes=os.path.join(out, "hog_classes.parquet"))
    log(f"classes annotated ({len(hits):,} Pfam hits; {len(hits.attrs.get('missing', []))} species without Pfam file)")
    # ---- per-branch class composition
    bce = branch_class_events(res.events, hc)
    bce.to_csv(os.path.join(out, "branch_class_events.tsv"), sep="\t", index=False)
    files["branch_class_events"] = os.path.join(out, "branch_class_events.tsv")
    # ---- membership / copy number at chosen levels
    tr = pd.read_csv(p["transitions"], sep="\t")
    if levels == "all":
        lv = [T.names[i] for i in T.internal]
    elif levels == "auto":
        lv = _auto_levels(T, tr, run, CONTROL_LEVELS)
    else:
        lv = [x for x in levels.split(",") if x]
    # members as integer rows (hog_row -> hogs.parquet row, gene_row -> genes.parquet row): compact for euk
    lvl_codes, hr, gr = [], [], []
    for i, x in enumerate(lv):
        idx = res._gene_hog_at(x)
        k = np.flatnonzero(idx >= 0)
        lvl_codes.append(np.full(len(k), i, dtype=np.int32))
        hr.append(idx[k].astype(np.int64))
        gr.append(k.astype(np.int64))
    cat = pd.CategoricalDtype(lv)
    mem = pd.DataFrame({"level": pd.Categorical.from_codes(np.concatenate(lvl_codes) if lv else np.array([], np.int32),
                                                           dtype=cat),
                        "hog_row": np.concatenate(hr) if lv else np.array([], np.int64),
                        "gene_row": np.concatenate(gr) if lv else np.array([], np.int64)})
    mem.to_parquet(os.path.join(out, "members.parquet"), index=False)
    sp_codes = res.genes["species"].cat.codes.to_numpy()
    cn = (pd.DataFrame({"level": mem["level"], "hog_row": mem["hog_row"],
                        "species": pd.Categorical.from_codes(sp_codes[mem["gene_row"].to_numpy()],
                                                             dtype=res.genes["species"].dtype)})
          .groupby(["level", "hog_row", "species"], observed=True).size().rename("n_copies").reset_index())
    cn["hog_id"] = res.hogs["hog_id"].to_numpy()[cn["hog_row"].to_numpy()]
    cn.to_parquet(os.path.join(out, "copy_number.parquet"), index=False)
    del mem, cn
    files.update(members=os.path.join(out, "members.parquet"), copy_number=os.path.join(out, "copy_number.parquet"))
    info = dict(run=run, orthoxml=p["orthoxml"], orthoxml_size=os.path.getsize(p["orthoxml"]), tree=p["tree"],
                tree_names=list(T.names), species=species, levels=lv, stats=res.stats,
                pfam_dir=p["pfam"], pfam_missing_species=hits.attrs.get("missing", []),
                n_pfam_hits=int(len(hits)), wall_s=round(time.time() - t0, 1),
                finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    with open(os.path.join(out, "run_info.json"), "w") as fh:
        json.dump(info, fh, indent=1, default=str)
    with open(os.path.join(out, "DONE.json"), "w") as fh:
        json.dump(dict(status="done", finished=info["finished"], files=files), fh, indent=1)
    return files


def branch_class_events(events: pd.DataFrame, hc: pd.DataFrame) -> pd.DataFrame:
    """Per branch (branch_child): number of HOGs gained / duplicated (non-origin) /
    lost, overall and per class. Gains use the gained HOG, duplications and
    losses the HOG at the parent node. Gains on leaf branches have no HOG row
    (species-specific genes) and count only in the totals."""
    cls = hc.set_index("hog_id")[[f"in_{c}" for c in CLASS_COLS]]
    ev = events[(events["event"] != "dup") | ~events["at_origin"].astype(bool)]
    out = []
    for e, name in (("gain", "gain"), ("dup", "dup"), ("loss", "loss")):
        x = ev[ev["event"] == e][["branch_child", "hog_id"]].drop_duplicates()
        x = x.join(cls, on="hog_id")
        g = x.groupby("branch_child")
        t = pd.DataFrame({f"{name}_hogs": g.size()})
        for c in CLASS_COLS:
            t[f"{name}_{c}"] = g[f"in_{c}"].sum().astype(int)
        out.append(t)
    df = pd.concat(out, axis=1).fillna(0).astype(int)
    return df.reset_index().rename(columns={"index": "branch_child"})


@dataclass
class Cache:
    run: str
    dir: str
    tree: SpeciesTree
    branches: pd.DataFrame
    hogs: pd.DataFrame
    genes: pd.DataFrame
    events: pd.DataFrame
    hog_classes: pd.DataFrame
    gene_classes: pd.DataFrame
    hits: pd.DataFrame
    branch_class_events: pd.DataFrame
    info: dict
    _lineage: object = None

    @property
    def lineage(self) -> "Lineage":
        if self._lineage is None:
            self._lineage = Lineage(self)
        return self._lineage

    def class_hogs(self, col: str) -> set:
        """hog ids in class ``col`` ('tf_dbd', ...) or module (``mod_...``)."""
        c = col if col.startswith("mod_") else f"in_{col}"
        return set(self.hog_classes.loc[self.hog_classes[c].to_numpy(), "hog_id"])

    def as_hogresult(self) -> HogResult:
        return HogResult(tree=self.tree, branches=self.branches, hogs=self.hogs, genes=self.genes,
                         events=self.events, stats=self.info.get("stats", {}))


def load_cache(run: str, paths: dict | None = None) -> Cache:
    p = paths or run_paths(run)
    d = p["cache"]
    if not os.path.exists(os.path.join(d, "DONE.json")):
        raise SystemExit(f"no analysis cache in {d}; run scripts/analysis/build_branch_events.py --run {run}")
    info = json.load(open(os.path.join(d, "run_info.json")))
    tree = load_species_tree(p["tree"], info["species"]) if os.path.exists(p["tree"]) else None
    if tree is None or list(tree.names) != list(info["tree_names"]):
        raise SystemExit("species tree differs from the one used to build the cache")
    rd = lambda n: pd.read_parquet(os.path.join(d, n))  # noqa: E731
    return Cache(run=run, dir=d, tree=tree,
                 branches=pd.read_csv(os.path.join(d, "branch_summary.tsv"), sep="\t", keep_default_na=False,
                                      na_values=[]),
                 hogs=rd("hogs.parquet"), genes=rd("genes.parquet"), events=rd("events.parquet"),
                 hog_classes=rd("hog_classes.parquet"), gene_classes=rd("gene_classes.parquet"),
                 hits=rd("pfam_hits.parquet"),
                 branch_class_events=pd.read_csv(os.path.join(d, "branch_class_events.tsv"), sep="\t"),
                 info=info)


# =============================================================================
# lineage queries (positive controls; data/prereg_definitions.tsv)
# =============================================================================

class Lineage:
    """HOG lineage navigation on a :class:`Cache` (rows of ``cache.hogs``)."""

    def __init__(self, c: Cache):
        self.c = c
        self.T = c.tree
        self.level = c.hogs["level_idx"].to_numpy()
        self.parent = c.hogs["parent_idx"].to_numpy()
        n = len(self.level)
        order = np.argsort(self.parent, kind="stable")
        ps = self.parent[order]
        self._ch = order
        self._lo = np.searchsorted(ps, np.arange(n), side="left")
        self._hi = np.searchsorted(ps, np.arange(n), side="right")
        self.row_of = pd.Series(np.arange(n), index=c.hogs["hog_id"].to_numpy())
        self.gene_parent = c.genes["parent_idx"].to_numpy()
        self.gene_species = c.genes["species"].astype(str).to_numpy()
        self.gene_name = c.genes["gene"].to_numpy()

    def node(self, name: str) -> int:
        return self.T.index[name]

    def children(self, h: int) -> np.ndarray:
        return self._ch[self._lo[h]:self._hi[h]]

    def origin_row(self, h: int) -> int:
        while self.parent[h] >= 0:
            h = self.parent[h]
        return int(h)

    def origin_level(self, h: int) -> str:
        """origin(H): node at the lower end of the branch carrying the gain of H's lineage."""
        return self.T.names[self.level[self.origin_row(h)]]

    def at_level(self, h: int, X: int) -> list:
        """HOG rows of H's lineage at node X (internal): its ancestor at X if X is
        above H's level, its descendants at X if X is below, [h] if equal."""
        L = self.level[h]
        T = self.T
        if L == X:
            return [int(h)]
        if T.is_anc(X, L):
            x = h
            while x >= 0 and self.level[x] != X:
                x = self.parent[x]
            return [int(x)] if x >= 0 else []
        if T.is_anc(L, X):
            out, stack = [], [h]
            while stack:
                y = stack.pop()
                for z in self.children(y):
                    lz = self.level[z]
                    if lz == X:
                        out.append(int(z))
                    elif T.is_anc(lz, X):
                        stack.append(z)
            return out
        return []

    def present_at(self, h: int, X: int) -> bool:
        """present_at(N) for H's lineage. For a leaf X: a member gene in that species."""
        if self.T.is_leaf[X]:
            return self.T.names[X] in set(self.member_species([h])[h])
        return len(self.at_level(h, X)) > 0

    def lost_on_stem(self, h: int, crown: int) -> bool:
        """'loss on stem(C)' for H: lineage present at parent(crown), absent at crown."""
        par = self.T.parent[crown]
        if par < 0:
            return False
        return self.present_at(h, par) and not self.present_at(h, crown)

    def members(self, rows: Sequence[int]) -> pd.DataFrame:
        """Long table (hog_row, gene_row) of the member genes of the given HOG rows."""
        rows = np.unique(np.asarray(list(rows), dtype=np.int64))
        mark = np.zeros(len(self.level), dtype=bool)
        mark[rows] = True
        cur = self.gene_parent
        gi = np.flatnonzero(cur >= 0)
        cur = cur[gi]
        oh, og = [], []
        while len(cur):
            hit = mark[cur]
            if hit.any():
                oh.append(cur[hit])
                og.append(gi[hit])
            nxt = self.parent[cur]
            ok = nxt >= 0
            cur, gi = nxt[ok], gi[ok]
        if not oh:
            return pd.DataFrame({"hog_row": np.array([], dtype=np.int64), "gene_row": np.array([], dtype=np.int64)})
        return pd.DataFrame({"hog_row": np.concatenate(oh), "gene_row": np.concatenate(og)})

    def member_species(self, rows: Sequence[int]) -> dict:
        m = self.members(rows)
        m["species"] = self.gene_species[m["gene_row"].to_numpy()]
        d = m.groupby("hog_row")["species"].agg(list).to_dict()
        return {int(r): d.get(int(r), []) for r in rows}

    def member_genes(self, rows: Sequence[int]) -> dict:
        m = self.members(rows)
        m["gene"] = self.gene_name[m["gene_row"].to_numpy()]
        d = m.groupby("hog_row")["gene"].agg(set).to_dict()
        return {int(r): d.get(int(r), set()) for r in rows}

    def clade(self, name: str) -> set:
        i = self.T.index[name]
        return {self.T.names[s] for s in self.T.leaves if self.T.is_anc(i, s)}


def proteins_with(hits: pd.DataFrame, pfams: Iterable[str]) -> set:
    pf = {_acc(x) for x in pfams}
    return set(hits.loc[hits["pfam_acc"].isin(pf), "protein"])


def pf_hogs(c: Cache, pfams: Iterable[str], min_frac: float = MIN_FRAC, min_genes: int = 2) -> np.ndarray:
    """PF-HOG(X): HOG rows with >= ``min_genes`` members of which >= 50 % have a
    hit to any Pfam in ``pfams`` (data/prereg_definitions.tsv)."""
    prots = proteins_with(c.hits, pfams)
    w = c.genes["gene"].isin(prots).to_numpy()
    if not w.any():
        return np.array([], dtype=np.int64)
    a = accumulate(c.hogs, c.genes, w, c.tree)
    ng = c.hogs["n_genes"].to_numpy()
    return np.flatnonzero((ng >= min_genes) & (a >= min_frac * ng))


def architecture_proteins(hits: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Proteins with hit(A) and hit(B), A starting (envelope start) before B.
    Returns (protein, species)."""
    a, b = _acc(a), _acc(b)
    ha = hits[hits["pfam_acc"] == a].groupby("protein").agg(a_min=("env_from", "min"), species=("species", "first"))
    hb = hits[hits["pfam_acc"] == b].groupby("protein")["env_from"].max().rename("b_max")
    j = ha.join(hb, how="inner")
    j = j[j["a_min"] < j["b_max"]]
    return j.reset_index()[["protein", "species"]]


# =============================================================================
# branch statistics inputs
# =============================================================================

def load_setup(run: str, paths: dict | None = None):
    """(transitions of this kingdom, traits, ete3 guide tree, branch table with
    lengths from data/species_tree_bl_<run>.nwk) -- as in scripts/phylo/h9_rates.py."""
    p = paths or run_paths(run)
    tr = pd.read_csv(p["transitions"], sep="\t")
    traits = pd.read_csv(p["traits"], sep="\t")
    tree = bs.load_tree(p["guide"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bl = bs.branch_lengths_from_tree(tree, p["bl"] if os.path.exists(p["bl"]) else None)
        bt = bs.build_branch_table(tree, traits, tr[tr["kingdom"] == run], branch_lengths=bl)
    return tr[tr["kingdom"] == run].copy(), traits, tree, bt


def h1h2_counts(ladder_path: str, bt: pd.DataFrame, column: str = "seq_novel") -> pd.Series:
    """Sequence-novel gains per branch from results/gainladder/<run>/branch_gain_ladder.tsv,
    summed over path-stem edges (D7). Branches absent from the ladder (root
    children are present; leaf branches are not) get NaN."""
    lad = pd.read_csv(ladder_path, sep="\t")
    y = lad.set_index("branch_child")[column].astype(float).reindex(bt.index)
    return bs.path_sum(y, bt)


def wgd_flags(c: Cache, bt: pd.DataFrame, qc_path: str | None, ratio: float = 3.0) -> pd.DataFrame:
    """WGD flags for the H4/H5 exclusion (DESIGN §4.6).

    Rule (fixed before any eukaryote result was seen): an internal branch is
    flagged if its duplication rate d_b / n_parent (d_b = HOGs at the parent
    with >= 1 non-origin duplication on the branch; n_parent = HOGs at the
    parent) is >= ``ratio`` x the median rate over all internal, non-root-child
    branches; or if it is named in results/qc/flags.tsv (column ``code``,
    ``branch`` or ``node``) with an issue mentioning duplication or WGD.
    Terminal branches are reported but play no role (not in U, never stems)."""
    ev = c.events
    d = (ev[(ev["event"] == "dup") & ~ev["at_origin"].astype(bool)].groupby("branch_child")["hog_id"].nunique())
    br = c.branches.set_index("branch_child")
    df = pd.DataFrame(index=bt.index)
    df["n_parent"] = br["n_parent"].reindex(df.index)
    df["d"] = d.reindex(df.index).fillna(0).astype(int)
    df["dup_rate"] = df["d"] / df["n_parent"].where(df["n_parent"] > 0)
    ref = ~bt["is_terminal"] & ~bt["root_child"]
    med = float(df.loc[ref, "dup_rate"].median())
    df["median_rate_internal"] = med
    df["ratio_to_median"] = df["dup_rate"] / med if med > 0 else np.nan
    df["flag_dup_excess"] = (df["ratio_to_median"] >= ratio) & ~bt["is_terminal"]
    df["flag_qc"] = False
    df["qc_issue"] = ""
    if qc_path and os.path.exists(qc_path):
        q = pd.read_csv(qc_path, sep="\t")
        key = next((k for k in ("branch", "node", "code", "branch_child") if k in q.columns), None)
        if key is not None:
            txt = q.astype(str).agg(" ".join, axis=1).str.lower()
            sel = q[txt.str.contains("wgd|duplicat|polyploid")]
            for _, r in sel.iterrows():
                b = str(r[key])
                if b in df.index:
                    df.at[b, "flag_qc"] = True
                    df.at[b, "qc_issue"] = str(r.get("issue", ""))[:120]
    df["wgd_flag"] = df["flag_dup_excess"] | df["flag_qc"]
    return df


def h12_module_hogs(c: Cache, modules: Sequence[str] = H12_MODULES, min_present: float = 0.5) -> tuple:
    """(module_hogs, coverage) for :func:`mcell.branchstats.module_loss_statistic`.

    module_hogs: (level, hog_id, module) for module HOGs at every internal
    node where the module is 'present', i.e. at least ``min_present`` of the
    module's components (data/h12_sets.tsv members; the three BBS2 domain
    models count as one component) are carried by >= 50 % of the members of
    some module HOG at that node. Elsewhere the module is skipped (§4.6:
    "modules with < 50 % of HOGs present at u skipped")."""
    spec = load_class_spec(run_paths(c.run)["toolkit"], run_paths(c.run)["h12_sets"])
    hc = c.hog_classes
    rows, cov = [], []
    for m in modules:
        col = _mod_col(m)
        if col not in hc.columns:
            continue
        ncomp = len(spec.modules.get(m, {}))
        sub = hc.loc[hc[col].to_numpy(), ["hog_id", "level", "components"]]
        for lvl, g in sub.groupby("level"):
            comps = set()
            for s in g["components"]:
                comps |= {x.split("|", 1)[1] for x in s.split(";") if x.startswith(m + "|")}
            ok = ncomp > 0 and len(comps) >= min_present * ncomp
            cov.append(dict(level=lvl, module=m, n_hogs=len(g), n_components=len(comps),
                            n_components_total=ncomp, present=ok))
            if ok:
                rows.append(pd.DataFrame({"level": lvl, "hog_id": g["hog_id"].to_numpy(), "module": m}))
    mh = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["level", "hog_id", "module"])
    return mh, pd.DataFrame(cov)


def supergroups(transitions: pd.DataFrame, species: pd.DataFrame) -> dict:
    """transition_id -> eukaryote supergroup of its multicellular species (majority of
    the first lineage_group token: Holozoa/Holomycota -> Opisthokonta;
    Stramenopiles/Alveolata/Rhizaria -> SAR; Archaeplastida, Amoebozoa as is)."""
    mp = {"Holozoa": "Opisthokonta", "Holomycota": "Opisthokonta", "Stramenopiles": "SAR",
          "Alveolata": "SAR", "Rhizaria": "SAR"}
    lg = species.set_index("code")["lineage_group"].astype(str).str.split("-").str[0]
    out = {}
    for _, r in transitions.iterrows():
        sps = [s for s in str(r.get("mc_species", "")).split(",") if s in lg.index]
        if not sps:
            out[r["transition_id"]] = "NA"
            continue
        g = lg.loc[sps].map(lambda x: mp.get(x, x)).value_counts()
        out[r["transition_id"]] = g.index[0]
    return out


# =============================================================================
# one relabelling suite (primary statistic + sensitivities + effects + LOO)
# =============================================================================

@dataclass
class Suite:
    name: str
    stems: dict
    match: object
    relabel: object
    relabel_spec: object | None
    effects: pd.DataFrame
    pooled: dict
    pooled_hksj: dict
    loo: pd.DataFrame
    excluded_stems: dict = field(default_factory=dict)


def run_relabel_suite(name: str, values: pd.Series, bt: pd.DataFrame, stems: Mapping[str, str],
                      alternative: str, effect_values: pd.Series | None = None, scale: str = "difference",
                      extra=None, signs: Mapping[str, float] | None = None, groups: Mapping[str, str] | None = None,
                      n_draws: int = bs.N_DRAWS, seed: int = bs.SEED, spec_sensitivity: bool = True,
                      pool: Sequence[str] | None = None) -> Suite:
    """Match, relabel ('permute', primary; 'spec' sensitivity), effect table,
    REML pooling (z and HKSJ), leave-one-transition-out and (with ``groups``)
    leave-one-supergroup-out.

    Stems with a missing value are excluded (``excluded_stems``); the pool is U
    restricted to branches with a finite value (and finite ``extra``)."""
    values = values.astype(float)
    excl = {t: "value_missing" for t, b in stems.items() if not np.isfinite(values.get(b, np.nan))}
    st = {t: b for t, b in stems.items() if t not in excl}
    if pool is None:
        pool = list(bt.index[bt["in_U"]])
    pool = [b for b in pool if np.isfinite(values.get(b, np.nan))]
    if extra is not None:
        ex = pd.DataFrame(extra) if isinstance(extra, pd.Series) else extra
        good = ex.reindex(bt.index).notna().all(axis=1)
        pool = [b for b in pool if good.get(b, False)]
        for t, b in list(st.items()):
            if not good.get(b, False):
                excl[t] = "extra_covariate_missing"
                st.pop(t)
    if not st:
        raise ValueError(f"{name}: no stem with a value")
    m = bs.match_controls(bt, st, pool=pool, extra=extra)
    sg = {t: signs[t] for t in st} if signs else None
    r = bs.relabel_test(values, m, alternative=alternative, n_draws=n_draws, seed=seed, scheme="permute", signs=sg)
    rs = None
    if spec_sensitivity and len(st) > 1:
        try:
            rs = bs.relabel_test(values, m, alternative=alternative, n_draws=n_draws, seed=seed, scheme="spec",
                                 signs=sg)
        except Exception as e:  # report, do not abort
            warnings.warn(f"{name}: spec scheme failed: {e}")
    ev = values if effect_values is None else effect_values.astype(float)
    eff = bs.effect_table(ev, m, scale=scale, signs=sg)
    eff["u"] = r.u.reindex(eff.index)
    ok = eff.dropna(subset=["est"]) if "est" in eff.columns else eff.iloc[0:0]
    pooled = bs.random_effects(ok["est"], ok["var"]) if len(ok) else {}
    pooled_h = bs.random_effects(ok["est"], ok["var"], ci_method="hksj") if len(ok) > 1 else {}

    def fn(keep):
        out = {}
        k2 = [t for t in keep if t in ok.index]
        if k2:
            re = bs.random_effects(ok.loc[k2, "est"], ok.loc[k2, "var"])
            out.update(mu=re["mu"], ci_low=re["ci_low"], ci_high=re["ci_high"], tau2=re["tau2"])
        rr = bs.relabel_test(values, m, alternative=alternative, n_draws=n_draws, seed=seed,
                             signs=sg, transitions=keep)
        out.update(T_obs=rr.T_obs, p=rr.p)
        return out

    tids = list(st)
    loo = bs.leave_one_out(tids, fn)
    loo.insert(0, "kind", "transition")
    if groups:
        los = bs.leave_one_out(tids, fn, groups={t: groups.get(t, "NA") for t in tids})
        los.insert(0, "kind", "supergroup")
        loo = pd.concat([loo, los[los["left_out"] != "none"]], ignore_index=True)
    return Suite(name=name, stems=st, match=m, relabel=r, relabel_spec=rs, effects=eff, pooled=pooled,
                 pooled_hksj=pooled_h, loo=loo, excluded_stems=excl)


def matched_sets_table(m) -> pd.DataFrame:
    rows = []
    for t, C in m.sets.items():
        for rank, (cc, d) in enumerate(zip(C, m.distances[t])):
            rows.append(dict(transition_id=t, stem=m.stems[t], rank=rank + 1, control=cc, distance=float(d)))
    return pd.DataFrame(rows)


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and not math.isfinite(p)):
        return "NA"
    return f"{p:.4g}"
