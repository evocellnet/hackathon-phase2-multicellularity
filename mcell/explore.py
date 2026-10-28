"""Shared layer of the EXPLORATORY analyses (DESIGN §1 Q2-Q5, §3.3, §4.3-4.6; scripts/analysis/
{pfam_wide,cooption,architectures,protein_traits,reversals}.py). Nothing here is part of the primary family.

1. **Contrasts** -- :class:`Contrast`: a transition (M clade vs sister species S, P, path edges) or a
   pseudo-transition built from a control branch b = (u -> v) of the pool U (M = tips(v), S = tips(u) - tips(v)).
   :func:`transition_contrasts`, :func:`control_contrasts`, :func:`control_pool` (U, or a relaxed pool when
   |U| is too small, e.g. on the prokaryote tree), :func:`select_transitions`.
2. **Pfam tables** -- :func:`build_pfam_tables` reads full Pfam annotations species by species
   (data/annotations/pfam/<run>/<code>.pfam.tsv.gz; species without a file are listed, never treated as zeros)
   into per-species domain counts, per-protein domain sets and ordered architectures.
3. **Pfam-wide contrasts** -- :func:`domain_contrast` (expansion: one-sided Mann-Whitney on per-species counts,
   exact by enumeration when small, BH within the contrast, or log2 fold >= 1 with presence in >= 50 % of M;
   orthology-free gains and losses), :func:`convergence_test` (k transitions vs pseudo-transitions, permutation).
4. **HOG view** -- :func:`hog_domain_events` (HOGs gained / duplicated / lost on a stem that carry a domain in
   >= 50 % of their members).
5. **Co-option index** -- :func:`coption_index` (DESIGN §4.6, descriptive).
6. **Architectures** -- :func:`architecture_origins`, :func:`branch_architecture_counts` (H11).
7. **Paired traits** -- :func:`gene_hog_at`, :func:`single_copy_hogs`, :func:`paired_test`,
   :func:`two_sample_test`, Hodges-Lehmann estimators (H9 companion, H10).
8. **Reversals** -- :func:`reversal_test`.
"""
from __future__ import annotations

import glob
import itertools
import json
import math
import os
import warnings
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from . import branchstats as bs

SEED = 11
N_PERM = 10_000
TOOLKIT_CLASSES = ("tf_dbd", "adhesion_ecm", "signalling", "kinase_ptyr", "pcd", "cell_cycle", "c2h2")


# =============================================================================
# small statistics helpers
# =============================================================================

def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(r, 1.0)
    q[ok] = out
    return q


def hl_paired(d, alpha: float = 0.05, max_pairs: int = 20_000_000, seed: int = SEED) -> tuple:
    """One-sample (paired) Hodges-Lehmann estimate = median of the Walsh averages, with the
    distribution-free CI from the signed-rank critical value (normal approximation).
    For n(n+1)/2 > ``max_pairs`` the Walsh averages are sampled (seeded) and the CI ranks
    are applied as quantiles. Returns (est, lo, hi)."""
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    n = len(d)
    if n == 0:
        return (np.nan, np.nan, np.nan)
    N = n * (n + 1) // 2
    z = stats.norm.ppf(1 - alpha / 2)
    K = int(math.floor(n * (n + 1) / 4 - z * math.sqrt(n * (n + 1) * (2 * n + 1) / 24)))
    if N <= max_pairs:
        i, j = np.triu_indices(n)
        w = np.sort((d[i] + d[j]) / 2)
        est = float(np.median(w))
        if K < 1:
            return (est, np.nan, np.nan)
        return (est, float(w[K - 1]), float(w[N - K]))
    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, max_pairs)
    j = rng.integers(0, n, max_pairs)
    w = (d[i] + d[j]) / 2
    est = float(np.median(w))
    if K < 1:
        return (est, np.nan, np.nan)
    return (est, float(np.quantile(w, (K - 1) / N)), float(np.quantile(w, (N - K) / N)))


def hl_two_sample(x, y, alpha: float = 0.05) -> tuple:
    """Two-sample Hodges-Lehmann shift (median of x_i - y_j) with the Mann-Whitney CI
    (normal approximation). Returns (est, lo, hi); CI NaN when too few observations."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    m, n = len(x), len(y)
    if m == 0 or n == 0:
        return (np.nan, np.nan, np.nan)
    D = np.sort((x[:, None] - y[None, :]).ravel())
    est = float(np.median(D))
    z = stats.norm.ppf(1 - alpha / 2)
    K = int(math.floor(m * n / 2 - z * math.sqrt(m * n * (m + n + 1) / 12)))
    if K < 1:
        return (est, np.nan, np.nan)
    return (est, float(D[K - 1]), float(D[m * n - K]))


def paired_test(d, alpha: float = 0.05) -> dict:
    """Wilcoxon signed-rank (two-sided, zeros dropped) + paired Hodges-Lehmann on differences d."""
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    out = dict(n=int(len(d)), median_diff=float(np.median(d)) if len(d) else np.nan,
               frac_positive=float(np.mean(d > 0)) if len(d) else np.nan)
    est, lo, hi = hl_paired(d, alpha)
    out.update(hl=est, hl_lo=lo, hl_hi=hi)
    if len(d) >= 2 and np.any(d != 0):
        out["p"] = float(stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided").pvalue)
    else:
        out["p"] = np.nan
    return out


def two_sample_test(x, y, alpha: float = 0.05) -> dict:
    """Mann-Whitney (two-sided) + two-sample Hodges-Lehmann for x (M side) vs y (S side)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    out = dict(n_M=int(len(x)), n_S=int(len(y)),
               median_M=float(np.median(x)) if len(x) else np.nan,
               median_S=float(np.median(y)) if len(y) else np.nan)
    est, lo, hi = hl_two_sample(x, y, alpha)
    out.update(hl=est, hl_lo=lo, hl_hi=hi)
    if len(x) >= 1 and len(y) >= 1 and len(x) + len(y) >= 3:
        out["p"] = float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
    else:
        out["p"] = np.nan
    return out


def mwu_greater(XM: np.ndarray, XS: np.ndarray, max_exact: int = 2000) -> np.ndarray:
    """One-sided Mann-Whitney p (M > S) for every column of XM (nM x D) vs XS (nS x D).

    Exact conditional permutation distribution (enumeration of all C(nM+nS, nM) labellings,
    ties handled through mid-ranks) when the number of labellings is <= ``max_exact``,
    else the normal approximation with tie correction and continuity correction.
    Columns that are constant across all species get p = NaN (untestable)."""
    XM = np.asarray(XM, dtype=float)
    XS = np.asarray(XS, dtype=float)
    nM, nS = XM.shape[0], XS.shape[0]
    D = XM.shape[1]
    p = np.full(D, np.nan)
    if nM == 0 or nS == 0 or D == 0:
        return p
    X = np.vstack([XM, XS])
    var = X.max(0) != X.min(0)
    if not var.any():
        return p
    Xv = X[:, var]
    R = stats.rankdata(Xv, axis=0)
    W_obs = R[:nM].sum(0)
    N = nM + nS
    ncomb = math.comb(N, nM)
    if ncomb <= max_exact:
        C = np.zeros((ncomb, N))
        for r, c in enumerate(itertools.combinations(range(N), nM)):
            C[r, list(c)] = 1.0
        W = C @ R
        pv = (W >= W_obs[None, :] - 1e-9).mean(0)
    else:
        mu = nM * (N + 1) / 2
        # tie correction per column
        Sx = np.sort(Xv, axis=0)
        brk = np.vstack([np.ones((1, Sx.shape[1]), bool), Sx[1:] != Sx[:-1]])
        gid = np.cumsum(brk, axis=0) - 1 + np.arange(Sx.shape[1])[None, :] * N
        cnt = np.bincount(gid.ravel(), minlength=N * Sx.shape[1]).astype(float)
        tie = (cnt ** 3 - cnt).reshape(Sx.shape[1], N).sum(1)
        sd = np.sqrt(nM * nS / 12 * ((N + 1) - tie / (N * (N - 1))))
        zz = (W_obs - mu - 0.5) / np.where(sd > 0, sd, np.nan)
        pv = stats.norm.sf(zz)
    p[var] = pv
    return p


def fisher_or(a: int, b: int, c: int, d: int) -> dict:
    """2x2 table [[a, b], [c, d]] (rows: group of interest / background; cols: event / no event).
    Fisher exact p (two-sided and greater), sample OR, Haldane-corrected log OR with Woolf CI."""
    t = np.array([[a, b], [c, d]])
    p2 = float(stats.fisher_exact(t, alternative="two-sided")[1])
    pg = float(stats.fisher_exact(t, alternative="greater")[1])
    ah, bh_, ch, dh = (x + 0.5 for x in (a, b, c, d))
    lor = math.log(ah * dh / (bh_ * ch))
    se = math.sqrt(1 / ah + 1 / bh_ + 1 / ch + 1 / dh)
    return dict(odds_ratio=float(np.exp(lor)), or_lo=float(np.exp(lor - 1.96 * se)),
                or_hi=float(np.exp(lor + 1.96 * se)), p_two_sided=p2, p_greater=pg,
                rate_group=a / (a + b) if a + b else np.nan, rate_background=c / (c + d) if c + d else np.nan)


# =============================================================================
# 1. contrasts
# =============================================================================

def clade_tips(tree) -> dict:
    """node name -> frozenset of tip names (ete3 tree)."""
    tips = {}
    for n in tree.traverse("postorder"):
        tips[n.name] = frozenset([n.name]) if n.is_leaf() else frozenset().union(*(tips[c.name] for c in n.children))
    return tips


@dataclass
class Contrast:
    cid: str                    # transition_id, or 'ctrl:<branch>'
    kind: str                   # 'transition' | 'control'
    branch: str                 # stem branch id (= M node) or control branch id
    P: str
    M: list                     # species of the new-state clade (tips of M)
    S: list                     # sister species
    outside: list               # every species not below M
    outside_P: list             # every species not below P
    edges: list = field(default_factory=list)   # stem edges (child names), D7 path stems

    def restrict(self, species: Iterable[str]) -> "Contrast":
        keep = set(species)
        return Contrast(self.cid, self.kind, self.branch, self.P, [s for s in self.M if s in keep],
                        [s for s in self.S if s in keep], [s for s in self.outside if s in keep],
                        [s for s in self.outside_P if s in keep], list(self.edges))


def _split(x, sep=","):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    return [y.strip() for y in str(x).replace(";", sep).split(sep) if y.strip()]


def select_transitions(tr: pd.DataFrame, run: str, which: str = "default") -> list:
    """transition ids of this kingdom. 'default': analysis_set 'primary' if the kingdom has any
    (eukaryotes: the 7 primary origins), else 'secondary' (prokaryotes: the 4 descriptive
    transitions of PREREG); 'primary' / 'secondary' / 'exploratory' / 'all' (every non-reversal,
    non-grey transition with an internal M) or a comma list of ids."""
    t = tr[tr["kingdom"] == run]
    if which == "default":
        which = "primary" if (t["analysis_set"] == "primary").any() else "secondary"
    if which in ("primary", "secondary", "exploratory"):
        return list(t.loc[t["analysis_set"] == which, "transition_id"])
    if which == "all":
        ok = ~t["to_state"].isin(["grey", "reversal"]) & ~t["analysis_set"].isin(["dropped_no_data"])
        return list(t.loc[ok, "transition_id"])
    return [x for x in which.split(",") if x]


def transition_contrasts(tr: pd.DataFrame, tree, bt: pd.DataFrame, tids: Sequence[str],
                         require_internal: bool = True) -> list:
    """One :class:`Contrast` per transition id. M species = tips below M (union for a list of M tips);
    S = ``sister_species`` (else tips(P) - M); P = first node of ``stem_path`` (else the table's P,
    else parent(M)); edges = D7 path edges from the branch table."""
    tips = clade_tips(tree)
    nodes = {n.name: n for n in tree.traverse()}
    allsp = sorted(tips[tree.name])
    out = []
    rows = tr.set_index("transition_id")
    for tid in tids:
        r = rows.loc[tid]
        Ms = [m for m in _split(r.get("M")) if m in nodes]
        if not Ms:
            warnings.warn(f"{tid}: M not in tree; skipped")
            continue
        if require_internal and (len(Ms) != 1 or nodes[Ms[0]].is_leaf()):
            warnings.warn(f"{tid}: M is not a single internal node; skipped")
            continue
        Msp = sorted(frozenset().union(*(tips[m] for m in Ms)))
        sp = r.get("stem_path")
        if isinstance(sp, str) and ">" in sp:
            P = sp.split(">")[0].strip()
        elif isinstance(r.get("P"), str) and r.get("P") in nodes:
            P = r.get("P")
        else:
            P = nodes[Ms[0]].up.name
        S = [s for s in _split(r.get("sister_species")) if s in tips[tree.name] and s not in Msp]
        if not S:
            S = sorted(tips[P] - set(Msp))
        b = Ms[0]
        own = _own_path(sp, b, nodes)
        if own is not None:
            edges = own
        else:
            edges = bs.path_edges(bt, b) if (bt is not None and b in bt.index) else [b]
        out.append(Contrast(tid, "transition", b, P, Msp, sorted(S), sorted(set(allsp) - set(Msp)),
                            sorted(set(allsp) - tips[P]), list(edges)))
    return out


def _own_path(sp, m, nodes):
    """Edges (child names) of this transition's own stem_path 'P>X>M' if it is a valid chain ending at M."""
    if not isinstance(sp, str) or ">" not in sp:
        return None
    path = [x.strip() for x in sp.split(">")]
    ok = (path[-1] == m and all(x in nodes for x in path)
          and all(nodes[path[i + 1]].up is nodes[path[i]] for i in range(len(path) - 1)))
    return path[1:] if ok else None


def patch_path_stems(bt: pd.DataFrame, contrasts: Sequence[Contrast]) -> pd.DataFrame:
    """Copy of ``bt`` in which every selected transition's stem row carries that transition's own path
    (stem_path, summed branch_length, log_bl). Needed because build_branch_table keeps the first path among
    transitions sharing an M (T_volvocales_gv shadows the primary T_volvocine_tgv path)."""
    bt = bt.copy()
    for c in contrasts:
        if c.kind != "transition" or c.branch not in bt.index:
            continue
        if ";".join(c.edges) != (bt.at[c.branch, "stem_path"] or c.branch):
            bt.at[c.branch, "stem_path"] = ";".join(c.edges)
            el = bt.loc[c.edges, "edge_branch_length"]
            bt.at[c.branch, "branch_length"] = float(el.sum(min_count=len(c.edges)))
            bl = bt.at[c.branch, "branch_length"]
            bt.at[c.branch, "log_bl"] = np.log(bl) if bl > 0 else np.nan
            bt.attrs.setdefault("patched_paths", {})[c.cid] = ";".join(c.edges)
    return bt


def control_contrasts(tree, bt: pd.DataFrame, pool: Sequence[str]) -> list:
    """Pseudo-transitions from control branches b = (u -> v): M = tips(v), S = tips(u) - tips(v), P = u."""
    tips = clade_tips(tree)
    allsp = set(tips[tree.name])
    out = []
    for b in pool:
        u = bt.at[b, "parent"]
        M = sorted(tips[b])
        S = sorted(tips[u] - tips[b])
        out.append(Contrast(f"ctrl:{b}", "control", b, u, M, S, sorted(allsp - tips[b]),
                            sorted(allsp - tips[u]), [b]))
    return out


def control_pool(bt: pd.DataFrame, min_size: int = 10) -> tuple:
    """(pool, label). U (DESIGN §4.6) if |U| >= ``min_size``; else the relaxed pool: internal
    within-unicellular branches, not root children, >= 2 tips below and >= 2 sister tips
    (adjacency to stems allowed). The relaxed pool is only for dry runs / descriptive use."""
    U = list(bt.index[bt["in_U"]])
    if len(U) >= min_size:
        return U, "U"
    rel = bt[(bt["branch_class"] == "within_unicellular") & ~bt["is_terminal"] & ~bt["root_child"]
             & (bt["n_desc"] >= 2) & (bt["n_sister"] >= 2)]
    return list(rel.index), f"relaxed(|U|={len(U)}<{min_size})"


def match_sets(bt: pd.DataFrame, contrasts: Sequence[Contrast], pool: Sequence[str], k: int = bs.K_MAX):
    """Matched control sets N_i (mcell.branchstats.match_controls on the §4.6 covariates)."""
    stems = {c.cid: c.branch for c in contrasts}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = bs.match_controls(bt, stems, pool=list(pool), k=k)
    return m


def control_self_sets(bt: pd.DataFrame, pool: Sequence[str], k: int = bs.K_MAX) -> dict:
    """For every control c, its own matched set in pool \\ {c} (used when a flag needs a reference set)."""
    out = {}
    for c in pool:
        rest = [x for x in pool if x != c]
        if len(rest) < 2:
            out[c] = rest
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = bs.match_controls(bt, {c: c}, pool=rest, k=k)
        out[c] = m.sets[c]
    return out


def evaluable(c: Contrast, annotated: Iterable[str], min_out_cov: float = 0.8) -> bool:
    """A contrast is evaluable on a partial annotation if >= min(2, |M|) M species and >= 1 sister species are
    annotated and >= ``min_out_cov`` of the species outside M are annotated (gains/presence need them)."""
    r = c.restrict(annotated)
    return (len(r.M) >= min(2, len(c.M)) and len(r.S) >= 1
            and len(r.outside) >= min_out_cov * max(1, len(c.outside)))


def species_data_types(traits: pd.DataFrame) -> dict:
    """code -> 'genome' | 'transcriptome' (data_type of data/species_traits.tsv)."""
    if traits is None or "data_type" not in traits.columns:
        return {}
    return {str(k): str(v).lower() for k, v in zip(traits["code"], traits["data_type"]) if isinstance(v, str)}


def present_rule(counts: pd.DataFrame, species: Sequence[str], dtype: Mapping[str, str]) -> pd.Series:
    """Presence per column of ``counts`` (species x families) among ``species``: >= 1 genome species
    or >= 2 transcriptome species with count > 0 (DESIGN §4.6; unknown data type counts as genome)."""
    sp = [s for s in species if s in counts.index]
    if not sp:
        return pd.Series(False, index=counts.columns)
    sub = counts.loc[sp] > 0
    tx = np.array([dtype.get(s, "genome") == "transcriptome" for s in sp])
    g = sub.loc[~tx].sum(0) if (~tx).any() else pd.Series(0, index=counts.columns)
    t = sub.loc[tx].sum(0) if tx.any() else pd.Series(0, index=counts.columns)
    return (g >= 1) | (t >= 2)


# =============================================================================
# 2. Pfam tables
# =============================================================================

def pfam_files(d: str) -> dict:
    """code -> finished per-species Pfam file (``.pfam.tsv.gz``; ``.tmp`` files are in progress)."""
    return {os.path.basename(f).split(".pfam")[0]: f for f in glob.glob(os.path.join(d, "*.pfam.tsv.gz"))}


def _resolve_overlaps(d: pd.DataFrame, min_frac: float = 0.5) -> pd.DataFrame:
    from .gainladder import _resolve_overlaps as ro
    return ro(d, min_frac)


def architectures_from_hits(h: pd.DataFrame) -> pd.DataFrame:
    """(protein, arch, n_domains) for one species' hits (columns protein, pfam_acc, ali_from, ali_to,
    dom_score). Same rule as mcell.gainladder.pfam_architectures: overlapping hits resolved by score,
    accessions ordered by ali_from, tandem repeats collapsed (A-A-B -> A-B; A-B-A kept)."""
    if h.empty:
        return pd.DataFrame(columns=["protein", "arch", "n_domains"])
    d = h[["protein", "pfam_acc", "ali_from", "ali_to", "dom_score"]].copy()
    d = _resolve_overlaps(d).sort_values(["protein", "ali_from"], kind="mergesort")
    nd = d.groupby("protein")["pfam_acc"].nunique().rename("n_domains")
    prev = d.groupby("protein")["pfam_acc"].shift()
    d = d[d["pfam_acc"] != prev]
    a = d.groupby("protein", sort=False)["pfam_acc"].agg("-".join).rename("arch").to_frame()
    return a.join(nd).reset_index()


def build_pfam_tables(pfam_dir: str, species: Sequence[str], out_dir: str, force: bool = False, log=print) -> dict:
    """Per-species pass over the full Pfam annotation; writes to ``out_dir``:
    counts.parquet (species, pfam_acc, n_proteins, n_hits), protein_domains.parquet (protein, species,
    pfam_acc; unique pairs), arch.parquet (protein, species, arch, n_domains), names.tsv (pfam_acc,
    pfam_name) and manifest.json (annotated / missing species). Rebuilt only when the set of finished
    per-species files changed (or ``force``)."""
    os.makedirs(out_dir, exist_ok=True)
    files = pfam_files(pfam_dir)
    ann = [s for s in species if s in files]
    missing = [s for s in species if s not in files]
    sig = {s: [os.path.getsize(files[s]), int(os.path.getmtime(files[s]))] for s in ann}
    man_p = os.path.join(out_dir, "manifest.json")
    if not force and os.path.exists(man_p):
        man = json.load(open(man_p))
        if man.get("signature") == sig and all(os.path.exists(os.path.join(out_dir, f)) for f in
                                               ("counts.parquet", "protein_domains.parquet", "arch.parquet")):
            return man
    cnt, pdom, arch, names = [], [], [], {}
    for i, s in enumerate(ann):
        h = pd.read_csv(files[s], sep="\t", usecols=["protein", "pfam_acc", "pfam_name", "ali_from", "ali_to",
                                                     "dom_score"])
        if h.empty:
            continue
        h["pfam_acc"] = h["pfam_acc"].astype(str).str.split(".").str[0]
        names.update(dict(zip(h["pfam_acc"], h["pfam_name"])))
        u = h[["protein", "pfam_acc"]].drop_duplicates()
        c = u.groupby("pfam_acc").size().rename("n_proteins").to_frame()
        c["n_hits"] = h.groupby("pfam_acc").size()
        c = c.reset_index()
        c.insert(0, "species", s)
        cnt.append(c)
        u = u.assign(species=s)
        pdom.append(u)
        a = architectures_from_hits(h)
        a["species"] = s
        arch.append(a)
        if (i + 1) % 25 == 0:
            log(f"pfam tables: {i + 1}/{len(ann)} species")
    C = pd.concat(cnt, ignore_index=True) if cnt else pd.DataFrame(columns=["species", "pfam_acc", "n_proteins",
                                                                            "n_hits"])
    PD = pd.concat(pdom, ignore_index=True) if pdom else pd.DataFrame(columns=["protein", "pfam_acc", "species"])
    A = pd.concat(arch, ignore_index=True) if arch else pd.DataFrame(columns=["protein", "arch", "n_domains",
                                                                              "species"])
    for df in (C, PD, A):
        for col in ("species", "pfam_acc"):
            if col in df.columns:
                df[col] = df[col].astype("category")
    C.to_parquet(os.path.join(out_dir, "counts.parquet"), index=False)
    PD.to_parquet(os.path.join(out_dir, "protein_domains.parquet"), index=False)
    A.to_parquet(os.path.join(out_dir, "arch.parquet"), index=False)
    pd.DataFrame({"pfam_acc": list(names), "pfam_name": list(names.values())}).to_csv(
        os.path.join(out_dir, "names.tsv"), sep="\t", index=False)
    man = dict(pfam_dir=pfam_dir, annotated=ann, missing=missing, signature=sig, n_hits_rows=int(len(PD)))
    with open(man_p, "w") as fh:
        json.dump(man, fh, indent=1)
    return man


def count_matrix(counts: pd.DataFrame, species: Sequence[str], value: str = "n_proteins") -> pd.DataFrame:
    """species x pfam_acc matrix (annotated species only; zeros where a domain is absent)."""
    c = counts[counts["species"].astype(str).isin(set(species))]
    X = c.pivot_table(index="species", columns="pfam_acc", values=value, aggfunc="sum", fill_value=0,
                      observed=True)
    X.index = X.index.astype(str)
    X.columns = X.columns.astype(str)
    return X.reindex(index=list(species), fill_value=0).astype(np.int32)


# =============================================================================
# 3. Pfam-wide contrasts
# =============================================================================

STATUS_CODE = {"gained": 2, "expanded": 1, "": 0, "contracted": -1, "lost": -2}


def domain_contrast(X: pd.DataFrame, c: Contrast, min_per_side: int = 2, fold: float = 1.0,
                    min_frac_M: float = 0.5, min_frac_S_lost: float = 0.5, q_max: float = 0.05,
                    max_exact: int = 2000) -> pd.DataFrame:
    """Per Pfam domain, M species vs sister species S (only species present in X, i.e. annotated).

    Columns: n_M, n_S, n_out (annotated species of each set), mean_M, mean_S, frac_M, frac_S,
    n_out_present, log2fc = log2((mean_M + 1)/(mean_S + 1)), p_mwu (one-sided M > S; NaN with
    < ``min_per_side`` species on a side), q (BH within the contrast, over tested domains), and status:
      gained     present in >= 50 % of M, absent from every annotated species outside M;
      lost       absent from every M species, present in >= 50 % of S;
      expanded   not gained, and (q <= q_max with log2fc > 0) or (log2fc >= 1 and present in >= 50 % of M);
      contracted not lost, log2fc <= -1 and present in >= 50 % of S.
    Domains absent from M, S and outside are omitted."""
    c = c.restrict(X.index)
    XM = X.loc[c.M].to_numpy(float)
    XS = X.loc[c.S].to_numpy(float)
    XO = X.loc[c.outside].to_numpy(float) if c.outside else np.zeros((0, X.shape[1]))
    mM = XM.mean(0) if len(c.M) else np.full(X.shape[1], np.nan)
    mS = XS.mean(0) if len(c.S) else np.full(X.shape[1], np.nan)
    fM = (XM > 0).mean(0) if len(c.M) else np.full(X.shape[1], np.nan)
    fS = (XS > 0).mean(0) if len(c.S) else np.full(X.shape[1], np.nan)
    nout = (XO > 0).sum(0)
    anyM = (XM > 0).any(0) if len(c.M) else np.zeros(X.shape[1], bool)
    keep = anyM | ((XS > 0).any(0) if len(c.S) else False) | (nout > 0)
    df = pd.DataFrame({"pfam_acc": X.columns, "n_M": len(c.M), "n_S": len(c.S), "n_out": len(c.outside),
                       "mean_M": mM, "mean_S": mS, "frac_M": fM, "frac_S": fS, "n_out_present": nout})
    df["log2fc"] = np.log2((df["mean_M"] + 1) / (df["mean_S"] + 1))
    df["p_mwu"] = np.nan
    if len(c.M) >= min_per_side and len(c.S) >= min_per_side:
        tested = keep & (anyM | ((XS > 0).any(0)))
        idx = np.flatnonzero(tested)
        if len(idx):
            df.loc[idx, "p_mwu"] = mwu_greater(XM[:, idx], XS[:, idx], max_exact=max_exact)
    df["q"] = bh(df["p_mwu"].to_numpy())
    gained = (df["frac_M"] >= min_frac_M) & (df["n_out_present"] == 0) & (len(c.outside) > 0)
    lost = ~anyM & (df["frac_S"] >= min_frac_S_lost) & (len(c.M) > 0)
    exp_t = (df["q"] <= q_max) & (df["log2fc"] > 0)
    exp_f = (df["log2fc"] >= fold) & (df["frac_M"] >= min_frac_M)
    expanded = ~gained & (exp_t | exp_f)
    contracted = ~lost & (df["log2fc"] <= -fold) & (df["frac_S"] >= min_frac_S_lost)
    st = np.full(len(df), "", dtype=object)
    st[contracted.to_numpy()] = "contracted"
    st[expanded.to_numpy()] = "expanded"
    st[lost.to_numpy()] = "lost"
    st[gained.to_numpy()] = "gained"
    df["status"] = st
    df["expanded_by_test"] = exp_t & ~gained
    df["expanded_by_fold"] = exp_f & ~gained
    df.insert(0, "contrast", c.cid)
    return df[keep].reset_index(drop=True)


def status_matrix(tables: Mapping[str, pd.DataFrame], domains: Sequence[str] | None = None,
                  col: str = "status") -> pd.DataFrame:
    """domain x contrast matrix of status codes (gained 2, expanded 1, contracted -1, lost -2, 0)."""
    cols = {}
    for cid, t in tables.items():
        s = t.set_index("pfam_acc")[col]
        cols[cid] = s.map(STATUS_CODE).fillna(0).astype(int) if col == "status" else s
    M = pd.DataFrame(cols)
    if domains is not None:
        M = M.reindex(domains)
    return M.fillna(0)


def _draws(rng, n, n_cols, n_perm, sets_idx=None, max_tries=200):
    """(draws, scheme): n column indices per draw. ``sets_idx`` (list of index arrays, one per transition) ->
    'matched' (one per set, distinct, rejection); else n distinct of ``n_cols`` ('permute'/'pool')."""
    if sets_idx is not None:
        got, tries = [], 0
        while sum(len(g) for g in got) < n_perm and tries < max_tries:
            J = np.stack([s[rng.integers(0, len(s), size=n_perm)] for s in sets_idx], axis=1)
            if n > 1:
                J = J[~(np.diff(np.sort(J, 1), axis=1) == 0).any(1)]
            got.append(J)
            tries += 1
        if sum(len(g) for g in got) >= n_perm:
            return np.concatenate(got)[:n_perm], "matched"
        return None, "matched_infeasible"
    if n_cols >= n:
        return np.argsort(rng.random((n_perm, n_cols)), axis=1)[:, :n], "distinct"
    return rng.integers(0, n_cols, size=(n_perm, n)), "with_replacement"


def _null_counts(F, draws, k_obs, kmins, chunk):
    ge = np.zeros(F.shape[0], dtype=np.int64)
    ksum = np.zeros(F.shape[0], dtype=np.float64)
    G = {k: [] for k in kmins}
    for a in range(0, len(draws), chunk):
        K = F[:, draws[a:a + chunk]].sum(-1)            # D x m
        ge += (K >= k_obs[:, None]).sum(1)
        ksum += K.sum(1)
        for k in kmins:
            G[k].append((K >= k).sum(0))
    return ge, ksum / len(draws), {k: np.concatenate(v) for k, v in G.items()}


def convergence_test(obs: pd.DataFrame, ctrl: pd.DataFrame, sets: Mapping[str, Sequence[str]] | None = None,
                     n_perm: int = N_PERM, seed: int = SEED, kmins: Sequence[int] = (2, 3),
                     chunk: int = 500) -> tuple:
    """Convergence of a per-domain binary flag across independent transitions.

    obs: domain x transition (bool); ctrl: domain x control pseudo-transition (bool, the same flag computed on
    child clades of the control pool). Statistics: per domain k = number of flagged transitions; global G_k =
    number of domains flagged in >= k transitions.
    Null 'permute' (primary, as DESIGN v3.1 D6): the n transitions and all controls are exchangeable; each draw
    takes n distinct columns of [transitions | controls] at random and recomputes k and G_k (per domain this is
    the exact hypergeometric conditional test, here by simulation). Sensitivity 'matched' (``sets``): one control
    per transition from its matched set N_i, distinct across transitions (controls only; anti-conservative for
    domains rarely flagged among the controls -- reported, not used for q).
    p = (1 + #{null >= obs}) / (n_perm + 1); per-domain p only for k_obs >= 1, BH over those (q).
    Returns (per_domain, global, info)."""
    tids = list(obs.columns)
    n = len(tids)
    dom = obs.index.union(ctrl.index)
    O = obs.reindex(dom).fillna(False).astype(bool)
    Cdf = ctrl.reindex(dom).fillna(False).astype(bool)
    ccols = list(Cdf.columns)
    cidx = {c: i for i, c in enumerate(ccols)}
    k_obs = O.sum(1).to_numpy()
    rng = np.random.default_rng(seed)
    Fall = np.hstack([O.to_numpy(np.int8), Cdf.to_numpy(np.int8)])
    dr, sch = _draws(rng, n, Fall.shape[1], n_perm)
    ge, kmean, G = _null_counts(Fall, dr, k_obs, kmins, chunk)
    per = pd.DataFrame({"pfam_acc": dom, "k_obs": k_obs, "k_null_mean": kmean,
                        "n_controls_flagged": Cdf.sum(1).to_numpy(),
                        "transitions": [";".join(t for t, v in zip(tids, row) if v) for row in O.to_numpy()]})
    per["p_perm"] = np.where(k_obs >= 1, (1 + ge) / (n_perm + 1), np.nan)
    per["q"] = bh(per["p_perm"].to_numpy())
    rows = []
    for k in kmins:
        go = int((k_obs >= k).sum())
        rows.append(dict(k_min=k, scheme="permute", n_domains_obs=go, null_mean=float(G[k].mean()),
                         null_q95=float(np.quantile(G[k], 0.95)),
                         p_perm=float((1 + (G[k] >= go).sum()) / (n_perm + 1))))
    info = dict(scheme="permute", draw=sch, n_transitions=n, n_controls=len(ccols), n_perm=n_perm, seed=seed)
    # sensitivity: matched controls only
    if sets is not None and all(len([x for x in sets.get(t, []) if x in cidx]) >= 1 for t in tids):
        S = [np.array([cidx[x] for x in sets[t] if x in cidx]) for t in tids]
        dm, schm = _draws(np.random.default_rng(seed), n, len(ccols), n_perm, S)
        if dm is not None:
            gm, _, Gm = _null_counts(Cdf.to_numpy(np.int8), dm, k_obs, kmins, chunk)
            per["p_matched"] = np.where(k_obs >= 1, (1 + gm) / (n_perm + 1), np.nan)
            for k in kmins:
                go = int((k_obs >= k).sum())
                rows.append(dict(k_min=k, scheme="matched", n_domains_obs=go, null_mean=float(Gm[k].mean()),
                                 null_q95=float(np.quantile(Gm[k], 0.95)),
                                 p_perm=float((1 + (Gm[k] >= go).sum()) / (n_perm + 1))))
        info["matched"] = schm
    return per.sort_values(["k_obs", "p_perm"], ascending=[False, True]).reset_index(drop=True), \
        pd.DataFrame(rows), info


# =============================================================================
# 4. HOG view
# =============================================================================

def members_of(hog_parent: np.ndarray, gene_parent: np.ndarray, rows: Iterable[int]) -> pd.DataFrame:
    """(hog_row, gene_row) for the member genes of the given HOG rows (any level)."""
    rows = np.unique(np.asarray(list(rows), dtype=np.int64))
    mark = np.zeros(len(hog_parent), dtype=bool)
    if len(rows):
        mark[rows] = True
    cur = gene_parent.copy()
    gi = np.flatnonzero(cur >= 0)
    cur = cur[gi]
    oh, og = [], []
    while len(cur):
        hit = mark[cur]
        if hit.any():
            oh.append(cur[hit])
            og.append(gi[hit])
        nxt = hog_parent[cur]
        ok = nxt >= 0
        cur, gi = nxt[ok], gi[ok]
    if not oh:
        return pd.DataFrame({"hog_row": np.array([], np.int64), "gene_row": np.array([], np.int64)})
    return pd.DataFrame({"hog_row": np.concatenate(oh), "gene_row": np.concatenate(og)})


def gene_domain_codes(genes: pd.DataFrame, protein_domains: pd.DataFrame) -> tuple:
    """(gene_row array, domain code array, domain names) for the genes of the HOG tables."""
    gi = pd.Series(np.arange(len(genes)), index=genes["gene"].to_numpy())
    pdm = protein_domains[protein_domains["protein"].isin(gi.index)]
    cat = pd.Categorical(pdm["pfam_acc"].astype(str))
    return gi.loc[pdm["protein"].to_numpy()].to_numpy(), cat.codes.astype(np.int64), np.asarray(cat.categories)


def hog_domains(rows: np.ndarray, hogs: pd.DataFrame, genes: pd.DataFrame, gdom: tuple,
                min_frac: float = 0.5) -> pd.DataFrame:
    """(hog_row, pfam_acc) for HOG rows whose members carry the domain in >= ``min_frac`` of genes."""
    if len(rows) == 0:
        return pd.DataFrame(columns=["hog_row", "pfam_acc"])
    m = members_of(hogs["parent_idx"].to_numpy(), genes["parent_idx"].to_numpy(), rows)
    g_rows, d_codes, names = gdom
    gd = pd.DataFrame({"gene_row": g_rows, "d": d_codes})
    j = m.merge(gd, on="gene_row")
    if j.empty:
        return pd.DataFrame(columns=["hog_row", "pfam_acc"])
    cnt = j.groupby(["hog_row", "d"]).size().rename("n").reset_index()
    ng = hogs["n_genes"].to_numpy()[cnt["hog_row"].to_numpy()]
    cnt = cnt[cnt["n"] >= min_frac * ng]
    return pd.DataFrame({"hog_row": cnt["hog_row"].to_numpy(), "pfam_acc": names[cnt["d"].to_numpy()]})


def hog_domain_events(events: pd.DataFrame, hogs: pd.DataFrame, genes: pd.DataFrame, gdom: tuple,
                      edges: Sequence[str], row_of: pd.Series | None = None) -> pd.DataFrame:
    """Per Pfam domain: number of HOGs gained / duplicated (non-origin duplications; the HOG at the
    parent node) / lost on the given stem edges that carry the domain in >= 50 % of their members."""
    if row_of is None:
        row_of = pd.Series(np.arange(len(hogs)), index=hogs["hog_id"].to_numpy())
    ev = events[events["branch_child"].isin(set(edges))]
    sel = {"gain": ev[ev["event"] == "gain"],
           "dup": ev[(ev["event"] == "dup") & ~ev["at_origin"].astype(bool)],
           "loss": ev[ev["event"] == "loss"]}
    rows = {k: np.unique(row_of.reindex(v["hog_id"].unique()).dropna().to_numpy().astype(np.int64))
            for k, v in sel.items()}
    allr = np.unique(np.concatenate([r for r in rows.values()])) if rows else np.array([], np.int64)
    hd = hog_domains(allr, hogs, genes, gdom)
    out = []
    for k, r in rows.items():
        x = hd[hd["hog_row"].isin(set(r.tolist()))].groupby("pfam_acc").size().rename(k)
        out.append(x)
    df = pd.concat(out, axis=1).fillna(0).astype(int) if out else pd.DataFrame()
    for k in ("gain", "dup", "loss"):
        if k not in df.columns:
            df[k] = 0
    df.index.name = "pfam_acc"
    return df[["gain", "dup", "loss"]]


# =============================================================================
# 5. co-option index
# =============================================================================

def coption_index(X: pd.DataFrame, c: Contrast, dtype: Mapping[str, str], families: Iterable[str] | None = None,
                  rule: str = "expanded", fold: float = 1.0, min_frac_M: float = 0.5,
                  extra_present: Iterable[str] = ()) -> dict:
    """C_i of DESIGN §4.6 for one contrast on a species x family count matrix X.

    F = families with log2((mean_M+1)/(mean_S+1)) >= ``fold`` present in >= 50 % of M species
    (rule 'expanded'), or present in >= 50 % of M species (rule 'present'; toolkit version),
    optionally restricted to ``families``. Present at P: >= 1 genome or >= 2 transcriptome species
    with the family in S U outgroups (= every species outside M); present at P's parent: the same
    among species outside P. ``extra_present`` = families counted as present at P regardless (HOG
    version: split-repaired gains). Returns n_F, n_P, C (= n_P / n_F), n_Pparent, C_Pparent, F."""
    c = c.restrict(X.index)
    cols = list(X.columns) if families is None else [f for f in X.columns if f in set(families)]
    Y = X[cols]
    out = dict(contrast=c.cid, n_M=len(c.M), n_S=len(c.S), n_out=len(c.outside))
    if not c.M or (rule == "expanded" and not c.S) or not cols:
        out.update(n_F=0, n_P=0, C=np.nan, n_Pparent=0, C_Pparent=np.nan, F=[])
        return out
    mM = Y.loc[c.M].mean(0)
    fM = (Y.loc[c.M] > 0).mean(0)
    if rule == "expanded":
        mS = Y.loc[c.S].mean(0)
        sel = (np.log2((mM + 1) / (mS + 1)) >= fold) & (fM >= min_frac_M)
    elif rule == "present":
        sel = fM >= min_frac_M
    else:
        raise ValueError(rule)
    F = list(sel.index[sel.to_numpy()])
    pres = present_rule(Y[F], c.outside, dtype) if F else pd.Series(dtype=bool)
    presP = present_rule(Y[F], c.outside_P, dtype) if F else pd.Series(dtype=bool)
    ex = set(extra_present)
    if ex and F:
        pres = pres | pd.Series([f in ex for f in F], index=F)
    nF = len(F)
    out.update(n_F=nF, n_P=int(pres.sum()) if nF else 0, C=float(pres.mean()) if nF else np.nan,
               n_Pparent=int(presP.sum()) if nF else 0, C_Pparent=float(presP.mean()) if nF else np.nan, F=F)
    return out


# =============================================================================
# 6. architectures (H11)
# =============================================================================

def node_masks(tree) -> tuple:
    """(tip index dict, list of (node name, depth, bitmask)) sorted deepest first."""
    tips = clade_tips(tree)
    leaves = sorted(tips[tree.name])
    ti = {s: i for i, s in enumerate(leaves)}
    depth = {}
    for n in tree.traverse("preorder"):
        depth[n.name] = 0 if n.up is None else depth[n.up.name] + 1
    nodes = []
    for n in tree.traverse():
        m = 0
        for s in tips[n.name]:
            m |= 1 << ti[s]
        nodes.append((n.name, depth[n.name], m))
    nodes.sort(key=lambda x: -x[1])
    return ti, nodes


def architecture_origins(arch: pd.DataFrame, tree, annotated: Sequence[str], min_species: int = 2,
                         domain_species: Mapping[str, frozenset] | None = None) -> pd.DataFrame:
    """Per architecture (string over annotated species): n_species, n_domains (distinct accessions),
    origin = MRCA of its species (the architecture is 'new on the branch into origin': present in >= 2
    descendant species, absent from every annotated species outside), multi (>= 2 distinct domains),
    recombined (multi and every component domain occurs in some annotated species outside the origin
    clade = a new combination of old domains)."""
    a = arch[arch["species"].astype(str).isin(set(annotated))]
    sp = a.groupby("arch", observed=True)["species"].agg(lambda x: frozenset(map(str, x)))
    nd = a.groupby("arch", observed=True)["n_domains"].max()
    ti, nodes = node_masks(tree)
    tips = clade_tips(tree)
    ann = set(annotated)
    cache = {}
    origin = []
    for s in sp.to_numpy():
        if len(s) < min_species:
            origin.append("")
            continue
        if s not in cache:
            m = 0
            for x in s:
                m |= 1 << ti[x]
            cache[s] = next(nm for nm, _, mk in nodes if (mk & m) == m)
        origin.append(cache[s])
    out = pd.DataFrame({"arch": sp.index.astype(str), "n_species": [len(s) for s in sp.to_numpy()],
                        "n_domains": nd.reindex(sp.index).to_numpy(), "origin": origin})
    out["multi"] = out["n_domains"] >= 2
    if domain_species is None:
        ds = {}
        for ar, s in zip(sp.index.astype(str), sp.to_numpy()):
            for d in set(ar.split("-")):
                ds.setdefault(d, set()).update(s)
        domain_species = {k: frozenset(v) for k, v in ds.items()}
    rec = []
    for ar, o, mu in zip(out["arch"], out["origin"], out["multi"]):
        if not o or not mu:
            rec.append(False)
            continue
        inside = tips[o]
        rec.append(all(len(set(domain_species.get(d, frozenset())) & (ann - inside)) > 0 for d in set(ar.split("-"))))
    out["recombined"] = rec
    return out


def branch_architecture_counts(origins: pd.DataFrame, bt: pd.DataFrame, tree, annotated: Sequence[str]) -> pd.DataFrame:
    """Per branch (child name): new architectures (all / multi-domain / recombined), path-summed for D7
    stems, plus Pfam coverage (annotated fraction of descendant and of outside species)."""
    o = origins[origins["origin"] != ""]
    df = pd.DataFrame(index=bt.index)
    df["new_arch_all"] = o.groupby("origin").size().reindex(bt.index).fillna(0)
    df["new_arch_multi"] = o[o["multi"]].groupby("origin").size().reindex(bt.index).fillna(0)
    df["new_arch_recombined"] = o[o["recombined"]].groupby("origin").size().reindex(bt.index).fillna(0)
    for col in ("new_arch_all", "new_arch_multi", "new_arch_recombined"):
        df[col] = bs.path_sum(df[col].astype(float), bt)
    tips = clade_tips(tree)
    allsp = tips[tree.name]
    ann = set(annotated)
    df["cov_desc"] = [len(tips[b] & ann) / len(tips[b]) for b in df.index]
    df["cov_outside"] = [len((allsp - tips[b]) & ann) / max(1, len(allsp - tips[b])) for b in df.index]
    return df


# =============================================================================
# 7. HOG membership at a level, paired traits
# =============================================================================

def gene_hog_at(hogs: pd.DataFrame, genes: pd.DataFrame, level_idx: int, max_depth: int = 200) -> np.ndarray:
    """For every gene, the hogs row of its ancestral gene at ``level_idx`` (-1 if none)."""
    hl = hogs["level_idx"].to_numpy()
    hp = hogs["parent_idx"].to_numpy()
    cur = genes["parent_idx"].to_numpy().copy()
    out = np.full(len(cur), -1, dtype=np.int64)
    for _ in range(max_depth):
        ok = cur >= 0
        if not ok.any():
            break
        safe = np.where(ok, cur, 0)
        hit = ok & (hl[safe] == level_idx)
        out[hit] = cur[hit]
        cur = np.where(ok & ~hit, hp[safe], -1)
    return out


def hog_at_level(hogs: pd.DataFrame, rows: np.ndarray, level_idx: int, max_depth: int = 200) -> np.ndarray:
    """For HOG rows (any level at or below ``level_idx``), their ancestor row at ``level_idx`` (-1 if none)."""
    hl = hogs["level_idx"].to_numpy()
    hp = hogs["parent_idx"].to_numpy()
    cur = np.asarray(rows, dtype=np.int64).copy()
    out = np.full(len(cur), -1, dtype=np.int64)
    for _ in range(max_depth):
        ok = cur >= 0
        if not ok.any():
            break
        safe = np.where(ok, cur, 0)
        hit = ok & (hl[safe] == level_idx)
        out[hit] = cur[hit]
        cur = np.where(ok & ~hit, hp[safe], -1)
    return out


def single_copy_hogs(memb: pd.DataFrame, M: Sequence[str], S: Sequence[str], min_frac: float = 0.5,
                     single_copy: bool = True) -> pd.DataFrame:
    """HOGs (``memb``: hog_row, gene, species at one level) present in >= ``min_frac`` (and >= 1) of the M
    species and of the S species, and (``single_copy``) with exactly one gene in every M U S species that
    has it. Returns the member rows of those HOGs restricted to M U S, with column side ('M'/'S')."""
    side = {s: "M" for s in M}
    side.update({s: "S" for s in S})
    m = memb[memb["species"].isin(side)].copy()
    m["side"] = m["species"].map(side)
    per = m.groupby(["hog_row", "species"], observed=True).size().rename("n").reset_index()
    per["side"] = per["species"].map(side)
    g = per.groupby("hog_row")
    nM = per[per["side"] == "M"].groupby("hog_row").size()
    nS = per[per["side"] == "S"].groupby("hog_row").size()
    mx = g["n"].max()
    ok = pd.DataFrame({"nM": nM, "nS": nS, "mx": mx}).fillna(0)
    sel = (ok["nM"] >= max(1, math.ceil(min_frac * len(M)))) & (ok["nS"] >= max(1, math.ceil(min_frac * len(S))))
    if single_copy:
        sel &= ok["mx"] == 1
    keep = set(ok.index[sel])
    return m[m["hog_row"].isin(keep)]


def per_hog_side_medians(m: pd.DataFrame, trait: pd.Series) -> pd.DataFrame:
    """hog_row x (M, S): median trait over the member genes of each side (genes without a value dropped)."""
    x = m.assign(v=trait.reindex(m["gene"]).to_numpy()).dropna(subset=["v"])
    return x.pivot_table(index="hog_row", columns="side", values="v", aggfunc="median").reindex(columns=["M", "S"])


# =============================================================================
# 8. reversals
# =============================================================================

def descendant_mark(hogs: pd.DataFrame, genes: pd.DataFrame, tree_names: Sequence[str], is_leaf: Sequence[bool],
                    node: str, level_idx: int) -> np.ndarray:
    """Bool per hogs row at ``level_idx``: has >= 1 descendant ancestral gene at ``node`` (internal node:
    a HOG row at that level; leaf: a gene of that species)."""
    names = list(tree_names)
    ni = names.index(node)
    n = len(hogs)
    out = np.zeros(n, dtype=bool)
    if is_leaf[ni]:
        g = genes[genes["species"].astype(str) == node]
        start = g["parent_idx"].to_numpy()
        start = start[start >= 0]
    else:
        start = np.flatnonzero(hogs["level_idx"].to_numpy() == ni)
    anc = hog_at_level(hogs, start, level_idx)
    anc = anc[anc >= 0]
    out[anc] = True
    return out


def reversal_test(hogs: pd.DataFrame, genes: pd.DataFrame, events: pd.DataFrame, tree_names: Sequence[str],
                  is_leaf: Sequence[bool], prev_edges: Sequence[str], rev_parent: str, rev_child: str) -> dict:
    """Are HOGs gained / duplicated on the preceding transition stem (``prev_edges``) lost preferentially
    on the reversal branch rev_parent -> rev_child, relative to the other HOGs present at rev_parent?

    A HOG (ancestral gene) at rev_parent is 'stem-gained' if an ancestor of it (or itself) is a HOG
    gained on a stem edge, 'stem-duplicated' if an ancestor is the parent-node HOG of a non-origin
    duplication on a stem edge. It is lost on the reversal branch if it has no descendant at rev_child.
    Returns counts and Fisher tests for gained, duplicated and gained-or-duplicated vs the rest."""
    names = list(tree_names)
    li = names.index(rev_parent)
    row_of = pd.Series(np.arange(len(hogs)), index=hogs["hog_id"].to_numpy())
    ev = events[events["branch_child"].isin(set(prev_edges))]
    G = row_of.reindex(ev.loc[ev["event"] == "gain", "hog_id"].unique()).dropna().astype(np.int64).to_numpy()
    D = row_of.reindex(ev.loc[(ev["event"] == "dup") & ~ev["at_origin"].astype(bool), "hog_id"].unique()) \
        .dropna().astype(np.int64).to_numpy()
    rows = np.flatnonzero(hogs["level_idx"].to_numpy() == li)
    hp = hogs["parent_idx"].to_numpy()
    mg = np.zeros(len(hogs), bool)
    mg[G] = True
    md = np.zeros(len(hogs), bool)
    md[D] = True
    fg = np.zeros(len(rows), bool)
    fd = np.zeros(len(rows), bool)
    cur = rows.copy()
    for _ in range(500):
        ok = cur >= 0
        if not ok.any():
            break
        safe = np.where(ok, cur, 0)
        fg |= ok & mg[safe]
        fd |= ok & md[safe]
        cur = np.where(ok, hp[safe], -1)
    present = descendant_mark(hogs, genes, names, is_leaf, rev_child, li)[rows]
    lost = ~present
    out = dict(rev_parent=rev_parent, rev_child=rev_child, n_hogs_parent=int(len(rows)),
               n_stem_gained=int(fg.sum()), n_stem_dup=int(fd.sum()), n_lost=int(lost.sum()))
    for name, grp in (("gained", fg), ("dup", fd), ("gained_or_dup", fg | fd)):
        a = int((grp & lost).sum())
        b = int((grp & ~lost).sum())
        c = int((~grp & lost).sum())
        d = int((~grp & ~lost).sum())
        f = fisher_or(a, b, c, d) if (a + b) > 0 and (c + d) > 0 else dict(
            odds_ratio=np.nan, or_lo=np.nan, or_hi=np.nan, p_two_sided=np.nan, p_greater=np.nan,
            rate_group=np.nan, rate_background=np.nan)
        out[name] = dict(a_lost=a, b_kept=b, c_bg_lost=c, d_bg_kept=d, **f)
    return out


# =============================================================================
# context for the scripts
# =============================================================================

@dataclass
class Context:
    run: str
    paths: dict
    tr: pd.DataFrame
    traits: pd.DataFrame
    tree: object                 # ete3 guide tree
    bt: pd.DataFrame
    tids: list
    trans: list                  # transition Contrasts
    pool: list
    pool_label: str
    ctrls: list                  # control Contrasts
    match: object                # MatchResult (transitions vs pool)
    dtype: dict
    cache: object = None

    @property
    def stems(self) -> dict:
        return {c.cid: c.branch for c in self.trans}

    @property
    def species(self) -> list:
        return sorted(clade_tips(self.tree)[self.tree.name])


def load_context(run: str, which: str = "default", need_cache: bool = True, pool_min: int = 10) -> Context:
    from . import analysis as an
    p = an.run_paths(run)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tr, traits, tree, bt = an.load_setup(run, p)
        tids = select_transitions(pd.read_csv(p["transitions"], sep="\t"), run, which)
        trans = transition_contrasts(tr, tree, bt, tids)
    bt = patch_path_stems(bt, trans)
    pool, label = control_pool(bt, pool_min)
    ctrls = control_contrasts(tree, bt, pool)
    m = match_sets(bt, trans, pool)
    cache = an.load_cache(run, p) if need_cache else None
    return Context(run=run, paths=p, tr=tr, traits=traits, tree=tree, bt=bt, tids=[c.cid for c in trans],
                   trans=trans, pool=pool, pool_label=label, ctrls=ctrls, match=m,
                   dtype=species_data_types(traits), cache=cache)


def out_dir(run: str, name: str, root: str | None = None) -> str:
    from .analysis import ROOT
    d = os.path.join(root or ROOT, "results", "exploratory", run, name)
    os.makedirs(d, exist_ok=True)
    return d


def write_summary(d: str, summary: dict) -> None:
    with open(os.path.join(d, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


# =============================================================================
# misc
# =============================================================================

def prereg_ok(root: str) -> bool:
    import sys
    sys.path.insert(0, os.path.join(root, "scripts", "prereg"))
    try:
        import check_prereg
        return bool(check_prereg.verify(strict=False))
    except Exception as e:  # pragma: no cover
        warnings.warn(f"PREREG check failed to run: {e}")
        return False


def registry_row(hypothesis: str, statistic: str, n, estimate=np.nan, ci=(np.nan, np.nan), p=np.nan,
                 k="", notes: str = "", seed="", n_draws="", scheme="", inputs="") -> dict:
    return dict(hypothesis=hypothesis, statistic=statistic, n=n, k=k, estimate=estimate, ci_low=ci[0],
                ci_high=ci[1], p=p, p_mcse="", holm_p="NA", tier="exploratory", input_md5s=inputs,
                seed=seed, n_draws=n_draws, scheme=scheme, notes="label=exploratory; " + notes)
