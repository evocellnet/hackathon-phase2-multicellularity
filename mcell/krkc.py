"""Radical vs conservative amino-acid replacements (Kr/Kc proxy) on fixed gene trees.

Within-protein constraint proxy for the H9 wording rule (DESIGN §4.6: "relaxed constraint" only if a
within-protein ratio supports it). Amino-acid version of Hanada et al. 2007 / Zhang 2000:

* classifications (Zhang 2000): **charge** {R,H,K} (+) | {D,E} (−) | rest (neutral); **polarity**
  polar {R,N,D,C,Q,E,G,H,K,S,T,Y} | non-polar {A,I,L,M,F,P,W,V}. A replacement is *radical* if it
  changes class, *conservative* otherwise;
* ancestral states: Fitch parsimony on the rooted (pruned guide) gene tree, generalised to polytomies
  (bottom-up set = states shared by the maximum number of children; for two children this is the
  usual intersection-else-union); top-down: a node keeps its parent's state when that state is in its
  set, otherwise a state is drawn uniformly from its set (seeded). Gaps/X/ambiguous = missing
  (wildcard leaf: it never causes a change);
* per edge, counts of radical (n_R) and conservative (n_C) replacements.

For a unit (P, M, S) and a gene: M side = the gene-tree edges from P down to M (stem, incl. D7 paths)
plus every edge below M; S side = the edges from P to its S children plus everything below them.
Statistic per transition: log Kr/Kc ratio
    LR = log( (ΣR_M + ½)/(ΣC_M + ½) ) − log( (ΣR_S + ½)/(ΣC_S + ½) )
pooled over genes, with a gene-bootstrap CI. LR > 0 = proportionally more radical change on the M
side (consistent with relaxed constraint); the same sequences and sites are compared on both sides,
so site opportunities largely cancel. Caveats: parsimony undercounts multiple hits on long branches,
and the radical fraction rises with divergence, so a faster side can show a higher fraction even
without relaxed constraint -- compare with the H9 controls, which share the bias.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
CHARGE = {**{a: 1 for a in "RHK"}, **{a: -1 for a in "DE"}}
POLAR = set("RNDCQEGHKSTY")
SCHEMES = ("charge", "polarity")


def class_vector(scheme: str) -> np.ndarray:
    if scheme == "charge":
        return np.array([CHARGE.get(a, 0) for a in AA])
    if scheme == "polarity":
        return np.array([1 if a in POLAR else 0 for a in AA])
    raise ValueError(scheme)


def radical_matrix(scheme: str) -> np.ndarray:
    c = class_vector(scheme)
    return c[:, None] != c[None, :]


def encode(seq: str) -> np.ndarray:
    """(L, 20) bool state sets; non-standard characters = all states (missing)."""
    out = np.ones((len(seq), 20), dtype=bool)
    for i, ch in enumerate(seq.upper()):
        j = AA_IDX.get(ch)
        if j is not None:
            out[i] = False
            out[i, j] = True
    return out


def fitch_states(gt, seqs: Mapping[str, str], seed: int = 11) -> dict:
    """Most-parsimonious state (int 0..19, per site) at every node of a rates.GeneTree."""
    tips = sorted(gt.tips())
    L = len(seqs[tips[0]])
    sets = {}
    order = []
    stack = [gt.root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(gt.children.get(n, []))
    for n in reversed(order):
        if gt.is_leaf(n):
            sets[n] = encode(seqs[n]) if n in seqs else np.ones((L, 20), dtype=bool)
        else:
            cnt = sum(sets[c].astype(np.int16) for c in gt.children[n])
            sets[n] = cnt == cnt.max(axis=1, keepdims=True)
    rng = np.random.default_rng(seed)
    state = {}
    for n in order:
        S = sets[n]
        pick = np.argmax(S * rng.random(S.shape), axis=1)
        if n == gt.root:
            state[n] = pick
        else:
            ps = state[gt.parent[n]]
            keep = S[np.arange(L), ps]
            state[n] = np.where(keep, ps, pick)
    return state


def edge_counts(gt, seqs: Mapping[str, str], schemes: Sequence[str] = SCHEMES, seed: int = 11) -> pd.DataFrame:
    """Per edge (child node): n_changes and radical/conservative counts per scheme."""
    st = fitch_states(gt, seqs, seed)
    R = {s: radical_matrix(s) for s in schemes}
    rows = []
    for c, p in gt.parent.items():
        a, b = st[p], st[c]
        ch = a != b
        row = dict(node=c, parent=p, n_changes=int(ch.sum()))
        for s in schemes:
            rad = R[s][a[ch], b[ch]]
            row[f"R_{s}"] = int(rad.sum())
            row[f"C_{s}"] = int((~rad).sum())
        rows.append(row)
    return pd.DataFrame(rows).set_index("node")


def side_nodes(gt, P: str, M: str, S: Sequence[str], guide_idx) -> tuple[list, list] | None:
    """Gene-tree edges (child names) on the M and S sides of P, or None if the unit is not identifiable
    (same inclusion rules as mcell.rates.h9_gene_stat)."""
    from mcell import rates
    d = rates.h9_gene_stat(gt, P, M, S, guide_idx)
    if not d["included"]:
        return None
    tipsS = frozenset().union(*(guide_idx[s][1] for s in S))

    def below(n):
        out, stack = [], [n]
        while stack:
            x = stack.pop()
            out.append(x)
            stack.extend(gt.children.get(x, []))
        return out

    m_side, n = below(M), gt.parent[M]
    while n != P:
        m_side.append(n)
        n = gt.parent[n]
    s_side = []
    for c in gt.children[P]:
        if gt.tips(c) <= tipsS:
            s_side.extend(below(c))
    return m_side, s_side


def stem_nodes(gt, P: str, M: str, S: Sequence[str], guide_idx) -> tuple[list, list] | None:
    """Only the stem edges: P ⇝ M path on the M side, P → S-children edges on the S side."""
    sn = side_nodes(gt, P, M, S, guide_idx)
    if sn is None:
        return None
    m_path, n = [M], gt.parent[M]
    while n != P:
        m_path.append(n)
        n = gt.parent[n]
    tipsS = frozenset().union(*(guide_idx[s][1] for s in S))
    return m_path, [c for c in gt.children[P] if gt.tips(c) <= tipsS]


def log_ratio(RM, CM, RS, CS):
    return float(np.log((RM + 0.5) / (CM + 0.5)) - np.log((RS + 0.5) / (CS + 0.5)))


def unit_table(counts: Mapping[str, pd.DataFrame], gene_trees, guide_idx, units, scope: str = "clade",
               schemes: Sequence[str] = SCHEMES) -> pd.DataFrame:
    """Long table unit × gene with R/C sums per side and scheme. ``scope`` 'clade' or 'stem'."""
    fn = side_nodes if scope == "clade" else stem_nodes
    rows = []
    for uid, P, M, S in units:
        for g, ec in counts.items():
            sn = fn(gene_trees[g], P, M, S, guide_idx)
            if sn is None:
                continue
            ms, ss = sn
            row = dict(unit=uid, gene=g)
            for s in schemes:
                row[f"RM_{s}"] = int(ec.loc[ms, f"R_{s}"].sum())
                row[f"CM_{s}"] = int(ec.loc[ms, f"C_{s}"].sum())
                row[f"RS_{s}"] = int(ec.loc[ss, f"R_{s}"].sum())
                row[f"CS_{s}"] = int(ec.loc[ss, f"C_{s}"].sum())
            rows.append(row)
    return pd.DataFrame(rows)


def summarise(tab: pd.DataFrame, schemes: Sequence[str] = SCHEMES, n_boot: int = 2000, seed: int = 11) -> pd.DataFrame:
    """Per unit and scheme: pooled LR, gene-bootstrap 95 % CI, radical fractions per side."""
    rng = np.random.default_rng(seed)
    rows = []
    for uid, d in tab.groupby("unit", sort=False):
        n = len(d)
        bi = rng.integers(0, n, size=(n_boot, n))
        row = dict(unit=uid, n_genes=n)
        for s in schemes:
            a = d[[f"RM_{s}", f"CM_{s}", f"RS_{s}", f"CS_{s}"]].to_numpy(float)
            tot = a.sum(0)
            row[f"LR_{s}"] = log_ratio(*tot)
            bt = a[bi].sum(1)
            lr = np.log((bt[:, 0] + .5) / (bt[:, 1] + .5)) - np.log((bt[:, 2] + .5) / (bt[:, 3] + .5))
            row[f"ci_lo_{s}"], row[f"ci_hi_{s}"] = np.quantile(lr, [0.025, 0.975])
            row[f"fracR_M_{s}"] = tot[0] / max(tot[0] + tot[1], 1)
            row[f"fracR_S_{s}"] = tot[2] / max(tot[2] + tot[3], 1)
            row[f"n_changes_M_{s}"] = int(tot[0] + tot[1])
            row[f"n_changes_S_{s}"] = int(tot[2] + tot[3])
        rows.append(row)
    return pd.DataFrame(rows).set_index("unit")
