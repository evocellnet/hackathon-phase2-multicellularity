"""Per-branch gene-content change from hierarchical orthologous groups (HOGs).

Reads an OrthoXML file (Ancestra ``hogs.orthoxml`` or ``.orthoxml.gz``; any
OrthoXML with nested orthologGroup/paralogGroup works) together with the rooted
guide species tree, reconstructs the *ancestral genes* present at every node of
the species tree, links each ancestral gene to its ancestor at the parent node,
and counts gains, losses and duplications on every branch. The file is read
with ``lxml.etree.iterparse`` one top-level group (gene family) at a time, so
memory scales with the number of genes, not with the file size, and all work is
linear in (genes x tree depth).

How Ancestra encodes HOGs (read from ``ancestra/orthoxml.py`` and ``merge.py``)
-------------------------------------------------------------------------------
* An ``orthologGroup`` is a speciation node. Its ``TaxRange`` property names the
  guide-tree node of the speciation; Ancestra always sets it to the merge node,
  which equals the MRCA of the member species.
* A ``paralogGroup`` is a duplication. It carries no TaxRange; its ``taxonId``
  attribute (OrthoXML 0.5, resolved through the ``<taxonomy>`` block) names the
  node *directly below* the branch on which the duplication happened. Without a
  usable taxonId the duplication is placed at the MRCA of its members
  (latest-possible placement, as pyham and Ancestra do).
* Groups with one child are collapsed and genes in no group (species-specific
  genes / unmatched singletons) appear only in the ``<species>`` block.
* A gene lineage that is not merged at a level "passes up" unchanged: e.g. an
  orthologGroup at the root with children (gene of species A, gene of species D)
  means that the ancestral gene existed at every node between the root and A
  and between the root and D, and was lost in all side lineages. Levels where a
  HOG is absent are therefore implicit, never written.
* Multifurcating nodes are merged k-way; an orthologGroup at a k-way node may
  contain members from only some of the children (the others lost the gene).
* Ancestra prunes guide-tree leaves without a proteome and collapses the unary
  nodes this creates; unnamed internal nodes are called ``node<k>``. Use the
  same tree (``stages/input/guide_tree.nwk``) or the embedded ``<taxonomy>``.

Definitions
-----------
The nested groups of one top-level group form a reconciled gene tree. Each node
has a *position* on the species tree: an orthologGroup or a gene sits *at* a
species-tree node (its level, or its species); a paralogGroup sits *on the
branch above* its level. An **ancestral gene at species-tree node X** is a
gene-tree edge (u -> v) that crosses X, i.e. v is at or below X and u is
strictly above X (for a top-level orthologGroup at X, the virtual edge above it
crosses X). Equivalently: the set of extant genes descending from one gene of
the ancestral genome X. At a leaf the ancestral genes are the extant genes.
The **ancestor** of an ancestral gene at X is the ancestral gene at parent(X)
on the same gene-tree path; if none exists the gene *originated* on the branch
parent(X) -> X. HOG identifiers are ``"<level>|<OrthoXML id of v>"`` (v is a
group id, or a gene id when a single gene lineage passes through a level).

For each branch p -> c:

* **gain** -- a gene family (top-level group or ungrouped gene) whose oldest
  position lies on p -> c, i.e. ancestral genes at c without an ancestor at p.
  A top-level paralogGroup on p -> c is *one* gain plus k-1 duplications
  (``dups_at_origin``). Gains on a leaf branch are species-specific genes
  (``leaf_specific_genes`` counts the genes). Families present at the root are
  reported as gains on the root "stem" (``branch_parent == ""``).
* **loss** -- an ancestral gene at p with no descendant ancestral gene at c.
* **duplication** -- an ancestral gene at p with k > 1 descendants at c
  contributes k-1 (column ``dups``; ``dups_at_origin`` is reported separately).

Conservation law, checked for every branch:
``n_child = n_parent - losses + dups + dups_at_origin + gains``.

Every gain row also carries the number of species containing the family and
the number of child clades of c it is present in (``n_child_clades_present``):
gains supported by >= 2 descendant lineages are less likely to be annotation
or orthology-inference artefacts. Deeper Dollo or "strict" variants can be
derived from the HOG and event tables.

Outputs (:class:`HogResult`)
----------------------------
``branches``  one row per branch (+ root stem): counts and the conservation check.
``hogs``      one row per ancestral gene at an internal node: id, level, family,
              ancestor id (``None`` = originated on the branch above), gene and
              species counts, representative genes.
``genes``     one row per extant gene: species, family, ancestral gene at the
              parent of its species (``None`` = species-specific gene).
``events``    one row per (branch, HOG, event) for enrichment analyses.
:meth:`HogResult.members` and :meth:`HogResult.copy_number` give, for a level X,
the extant genes of each ancestral gene at X and the per-species copy numbers.
"""
from __future__ import annotations

import gzip
import heapq
import itertools
import os
from dataclasses import dataclass, field
from typing import IO, Callable, Iterable

import numpy as np
import pandas as pd

__all__ = ["SpeciesTree", "HogResult", "parse_hogs", "load_species_tree"]

GENE, ORTHO, PARA = 0, 1, 2
_GROUP_TAGS = {"orthologGroup": ORTHO, "paralogGroup": PARA}


# ----------------------------------------------------------------------------
# species tree
# ----------------------------------------------------------------------------

class SpeciesTree:
    """Rooted species tree indexed in pre-order (node id = pre-order rank).

    ``is_anc(a, b)`` (a ancestor of b or a == b) and MRCA queries are O(1) /
    O(depth) using pre-order intervals.
    """

    def __init__(self, root):
        order = list(root.preorder())
        self.names: list[str] = [nd.name for nd in order]
        if any(not nm for nm in self.names):
            raise ValueError("all species-tree nodes must be named")
        self.index: dict[str, int] = {}
        for i, nm in enumerate(self.names):
            if nm in self.index:
                raise ValueError(f"duplicate species-tree node name {nm!r}")
            self.index[nm] = i
        n = len(order)
        self.parent: list[int] = [-1] * n
        self.children: list[list[int]] = [[] for _ in range(n)]
        for i, nd in enumerate(order):
            if nd.parent is not None:
                p = self.index[nd.parent.name]
                self.parent[i] = p
                self.children[p].append(i)
        self.end: list[int] = list(range(n))
        for i in range(n - 1, -1, -1):
            if self.children[i]:
                self.end[i] = max(self.end[c] for c in self.children[i])
        self.depth: list[int] = [0] * n
        for i in range(1, n):
            self.depth[i] = self.depth[self.parent[i]] + 1
        self.is_leaf: list[bool] = [not c for c in self.children]
        self.leaves: list[int] = [i for i in range(n) if self.is_leaf[i]]
        self.internal: list[int] = [i for i in range(n) if not self.is_leaf[i]]
        # species bit (for species-set masks) and clade masks
        self.leaf_bit: dict[int, int] = {s: k for k, s in enumerate(self.leaves)}
        self.clade_mask: list[int] = [0] * n
        for i in range(n - 1, -1, -1):
            if self.is_leaf[i]:
                self.clade_mask[i] = 1 << self.leaf_bit[i]
            else:
                m = 0
                for c in self.children[i]:
                    m |= self.clade_mask[c]
                self.clade_mask[i] = m

    def __len__(self) -> int:
        return len(self.names)

    def is_anc(self, a: int, b: int) -> bool:
        """True if node ``a`` is an ancestor of ``b`` or equal to it."""
        return a <= b <= self.end[a]

    def mrca(self, lo: int, hi: int) -> int:
        """MRCA of nodes whose pre-order ids span [lo, hi] (hi = max ``end``)."""
        x = lo
        while not (x <= lo and self.end[x] >= hi):
            x = self.parent[x]
        return x

    def child_towards(self, a: int, b: int) -> int:
        """Child of ``a`` on the path to its strict descendant ``b``."""
        for c in self.children[a]:
            if c <= b <= self.end[c]:
                return c
        raise ValueError(f"{self.names[b]} is not below {self.names[a]}")

    def clades(self) -> dict[str, frozenset]:
        return {self.names[i]: frozenset(self.names[s] for s in self.leaves
                                         if self.is_anc(i, s)) for i in range(len(self))}


def load_species_tree(newick: str | None = None, species: Iterable[str] | None = None,
                      taxonomy_root=None) -> SpeciesTree:
    """Build the species tree used for the analysis.

    With ``newick`` (path or string) the tree is pruned to ``species`` and unary
    nodes are collapsed exactly as Ancestra does (``ancestra.io.load_guide_tree``),
    so internal names match the ``TaxRange`` values. If ``taxonomy_root`` (the
    OrthoXML ``<taxonomy>``) is also given, both trees must have the same named
    clades. Without ``newick`` the taxonomy is used.
    """
    if newick is not None:
        from ancestra.io import load_guide_tree
        root = load_guide_tree(newick, sorted(species) if species is not None else [])
        tree = SpeciesTree(root)
        if taxonomy_root is not None:
            other = SpeciesTree(taxonomy_root).clades()
            mine = tree.clades()
            if mine != other:
                diff = sorted(set(mine.items()) ^ set(other.items()))[:5]
                raise ValueError("guide tree and the OrthoXML <taxonomy> disagree "
                                 f"(first differing (name, clade) entries: {diff})")
        return tree
    if taxonomy_root is None:
        raise ValueError("no species tree: give a Newick tree or use an OrthoXML with <taxonomy>")
    return SpeciesTree(taxonomy_root)


# ----------------------------------------------------------------------------
# result container
# ----------------------------------------------------------------------------

@dataclass
class HogResult:
    tree: SpeciesTree
    branches: pd.DataFrame
    hogs: pd.DataFrame
    genes: pd.DataFrame
    events: pd.DataFrame
    stats: dict = field(default_factory=dict)

    # -- derived tables -----------------------------------------------------
    def _gene_hog_at(self, level: str) -> np.ndarray:
        """For every gene, the row index in ``hogs`` of its ancestral gene at
        ``level`` (-1 if the gene's species is not below ``level`` or the gene
        has no ancestor at ``level``)."""
        X = self.tree.index[level]
        hl = self.hogs["level_idx"].to_numpy()
        hp = self.hogs["parent_idx"].to_numpy()
        cur = self.genes["parent_idx"].to_numpy().copy()
        out = np.full(len(cur), -1, dtype=np.int64)
        for _ in range(max(self.tree.depth) + 1):
            ok = cur >= 0
            if not ok.any():
                break
            safe = np.where(ok, cur, 0)
            hit = ok & (hl[safe] == X)
            out[hit] = cur[hit]
            cur = np.where(ok & ~hit, hp[safe], -1)
        return out

    def members(self, level: str) -> pd.DataFrame:
        """Long table (hog_id, gene, species) of the extant genes of every
        ancestral gene at ``level`` (internal node)."""
        if self.tree.is_leaf[self.tree.index[level]]:
            raise ValueError("members() is for internal levels; leaf ancestral genes are the genes")
        idx = self._gene_hog_at(level)
        k = idx >= 0
        return pd.DataFrame({"level": level,
                             "hog_id": self.hogs["hog_id"].to_numpy()[idx[k]],
                             "gene": self.genes["gene"].to_numpy()[k],
                             "species": self.genes["species"].to_numpy()[k]})

    def copy_number(self, level: str) -> pd.DataFrame:
        """Long/sparse copy-number matrix at ``level``: (hog_id, species, n_copies)
        for every ancestral gene at ``level`` and every species with >= 1 copy."""
        m = self.members(level)
        cn = m.groupby(["hog_id", "species"], sort=True, observed=True).size()
        cn = cn.rename("n_copies").reset_index()
        cn.insert(0, "level", level)
        return cn

    def check_conservation(self) -> pd.DataFrame:
        """Rows of ``branches`` that violate the conservation law (should be empty)."""
        b = self.branches
        lhs = b["n_child"]
        rhs = b["n_parent"] - b["losses"] + b["dups"] + b["dups_at_origin"] + b["gains"]
        return b[(lhs != rhs) | (b[["gains", "losses", "dups", "dups_at_origin"]] < 0).any(axis=1)]

    def write(self, out_dir: str, fmt: str = "parquet") -> dict[str, str]:
        """Write branches (TSV) and hogs/genes/events (parquet or TSV)."""
        os.makedirs(out_dir, exist_ok=True)
        paths = {"branches": os.path.join(out_dir, "branch_summary.tsv")}
        self.branches.to_csv(paths["branches"], sep="\t", index=False)
        for name in ("hogs", "genes", "events"):
            df = getattr(self, name)
            if fmt == "parquet":
                p = os.path.join(out_dir, f"{name}.parquet")
                df.to_parquet(p, index=False)
            else:
                p = os.path.join(out_dir, f"{name}.tsv.gz")
                df.to_csv(p, sep="\t", index=False)
            paths[name] = p
        return paths


# ----------------------------------------------------------------------------
# parsing
# ----------------------------------------------------------------------------

def _local(tag) -> str:
    if not isinstance(tag, str):
        return ""
    k = tag.rfind("}")
    return tag[k + 1:] if k >= 0 else tag


class _Family:
    """Nested groups of one top-level group, in creation (pre-)order."""
    __slots__ = ("kind", "parent", "children", "gene", "xid", "hint", "taxid")

    def __init__(self):
        self.kind: list[int] = []
        self.parent: list[int] = []
        self.children: list[list[int]] = []
        self.gene: list[int] = []
        self.xid: list[str | None] = []
        self.hint: list[str | None] = []      # TaxRange property
        self.taxid: list[str | None] = []     # taxonId attribute

    def add(self, kind: int, parent: int, gene: int = -1, xid=None, taxid=None) -> int:
        i = len(self.kind)
        self.kind.append(kind)
        self.parent.append(parent)
        self.children.append([])
        self.gene.append(gene)
        self.xid.append(xid)
        self.hint.append(None)
        self.taxid.append(taxid)
        if parent >= 0:
            self.children[parent].append(i)
        return i


class _Builder:
    """Accumulates the output tables family by family."""

    def __init__(self, tree: SpeciesTree, gene_names: list[str], gene_sp: list[int],
                 taxon_name: dict[str, str], n_rep: int):
        self.T = tree
        self.gene_names = gene_names
        self.gene_sp = gene_sp                  # gene -> species-tree node id
        self.taxon_name = taxon_name
        self.n_rep = n_rep
        n = len(tree)
        self.n_anc = [0] * n
        self.gains = [0] * n
        self.losses = [0] * n
        self.dups = [0] * n
        self.dups_origin = [0] * n
        self.leaf_specific = [0] * n
        self.gains_multi = [0] * n
        self.gene_parent = [-1] * len(gene_names)
        self.gene_family: list[str | None] = [None] * len(gene_names)
        self.seen = np.zeros(len(gene_names), dtype=bool)
        self.H = {k: [] for k in ("hog_id", "level_idx", "family_id", "parent_idx", "n_genes",
                                  "n_species", "rep_genes", "n_child_clades_present")}
        self.E = {k: [] for k in ("branch_parent_idx", "branch_child_idx", "event", "n_events",
                                  "at_origin", "hog_id", "family_id", "child_hog_ids",
                                  "n_extant_genes", "n_genes_in_child", "n_species", "rep_genes",
                                  "n_child_clades_present", "n_child_clades_total")}
        self.st = dict(n_families=0, n_singleton_genes=0, n_groups=0, n_level_hint_unknown=0,
                       n_level_hint_inconsistent=0, n_collapsed_inconsistent=0,
                       n_single_species_orthologgroups=0, n_empty_groups=0)

    # -- one family -----------------------------------------------------------
    def family(self, F: _Family, fam_index: int) -> None:
        T, st = self.T, self.st
        par_sp, end = T.parent, T.end
        kind, parent, children, gene = F.kind, F.parent, F.children, F.gene
        n = len(kind)
        st["n_families"] += 1
        st["n_groups"] += sum(1 for k in kind if k != GENE)
        gnames = self.gene_names
        # ---- bottom-up: species span, gene count, species mask, representatives
        lo = [0] * n
        hi = [0] * n
        ng = [0] * n
        mask = [0] * n
        rep: list = [None] * n
        R = self.n_rep
        for i in range(n - 1, -1, -1):
            if kind[i] == GENE:
                s = self.gene_sp[gene[i]]
                lo[i] = hi[i] = s
                ng[i] = 1
                mask[i] = T.clade_mask[s]
                rep[i] = [gnames[gene[i]]]
            else:
                ch = [c for c in children[i] if ng[c] > 0]
                if not ch:
                    st["n_empty_groups"] += 1
                    continue
                lo[i] = min(lo[c] for c in ch)
                hi[i] = max(hi[c] for c in ch)
                ng[i] = sum(ng[c] for c in ch)
                m = 0
                for c in ch:
                    m |= mask[c]
                mask[i] = m
                if len(ch) == 1:
                    rep[i] = rep[ch[0]]
                else:
                    rep[i] = heapq.nsmallest(R, itertools.chain.from_iterable(rep[c] for c in ch))
        root = 0
        if ng[root] == 0:
            return
        fam_id = F.xid[root] or f"FAM{fam_index:08d}"
        # ---- top-down: positions (node, above_branch) with consistency repair
        pn = [-1] * n
        pa = [False] * n
        is_anc = T.is_anc
        for i in range(n):
            if ng[i] == 0:
                continue
            m = T.mrca(lo[i], hi[i])
            k = kind[i]
            if k == GENE:
                cands = [(m, False)]
            else:
                above = (k == PARA)
                if k == ORTHO and T.is_leaf[m]:
                    # orthologGroup with a single species: only in-paralogs -> duplication
                    st["n_single_species_orthologgroups"] += 1
                    above = True
                h = None
                name = F.hint[i] if (k == ORTHO and F.hint[i] is not None) else None
                if name is None and F.taxid[i] is not None:
                    name = self.taxon_name.get(F.taxid[i])
                if name is not None:
                    h = T.index.get(name)
                    if h is None:
                        st["n_level_hint_unknown"] += 1
                    elif not is_anc(h, m):
                        st["n_level_hint_inconsistent"] += 1
                        h = None
                cands = [(h, above)] if (h is not None and h != m) else []
                cands.append((m, above))
            p = parent[i]
            chosen = None
            if p < 0:
                chosen = cands[0]
            else:
                qn, qa = pn[p], pa[p]
                for cn, ca in cands:
                    ok = (cn != qn and is_anc(qn, cn)) if not qa else is_anc(qn, cn)
                    if ok:
                        chosen = (cn, ca)
                        break
                if chosen is None:          # cannot be placed below its parent: merge into it
                    st["n_collapsed_inconsistent"] += 1
                    chosen = (qn, qa)
            pn[i], pa[i] = chosen
        # ---- crossings: ancestral genes = gene-tree edges crossing species nodes
        cross_base = [-1] * n
        cross_len = [0] * n
        cz: list[int] = []          # species node of each crossing
        cv: list[int] = []          # gene-tree node (lower end of the edge)
        for i in range(n):
            if ng[i] == 0:
                continue
            zn, za = pn[i], pa[i]
            Z = par_sp[zn] if za else zn
            p = parent[i]
            start = len(cz)
            if p < 0:
                if not za:
                    cz.append(zn)
                    cv.append(i)
            else:
                qn, qa = pn[p], pa[p]
                if not qa:
                    while Z >= 0 and Z != qn:
                        cz.append(Z)
                        cv.append(i)
                        Z = par_sp[Z]
                else:
                    while Z >= 0 and is_anc(qn, Z):
                        cz.append(Z)
                        cv.append(i)
                        Z = par_sp[Z]
            if len(cz) > start:
                cross_base[i] = start
                cross_len[i] = len(cz) - start
        # ancestor mapping
        anc_cross = [-1] * n
        for i in range(n):
            if cross_len[i]:
                anc_cross[i] = i
            elif parent[i] >= 0:
                anc_cross[i] = anc_cross[parent[i]]
        nc = len(cz)
        cpar = [-1] * nc
        origin: list[int] = []
        for i in range(n):
            L = cross_len[i]
            if not L:
                continue
            b = cross_base[i]
            for k in range(L - 1):
                cpar[b + k] = b + k + 1
            top = b + L - 1
            P = par_sp[cz[top]]
            w = anc_cross[parent[i]] if parent[i] >= 0 else -1
            if w >= 0:
                j = cross_base[w]
                if P < 0 or cz[j] != P:
                    raise AssertionError(f"family {fam_id}: inconsistent crossing at "
                                         f"{T.names[cz[top]]} (bug)")
                cpar[top] = j
            else:
                origin.append(top)
        # ---- ids, HOG rows, counts
        names = T.names
        cid: list[str] = [""] * nc
        gidx = [-1] * nc
        H = self.H
        for j in range(nc):
            Z, i = cz[j], cv[j]
            self.n_anc[Z] += 1
            if kind[i] == GENE:
                xid = gnames[gene[i]]
            else:
                xid = F.xid[i] or f"{fam_id}.n{i}"
            if T.is_leaf[Z]:
                cid[j] = xid
            else:
                cid[j] = f"{names[Z]}|{xid}"
                gidx[j] = len(H["hog_id"])
                H["hog_id"].append(cid[j])
                H["level_idx"].append(Z)
                H["family_id"].append(fam_id)
                H["n_genes"].append(ng[i])
                H["n_species"].append(mask[i].bit_count())
                H["rep_genes"].append(";".join(rep[i]))
                H["parent_idx"].append(-1)      # filled below
                H["n_child_clades_present"].append(
                    sum(1 for c in T.children[Z] if mask[i] & T.clade_mask[c]))
        for j in range(nc):
            if gidx[j] >= 0 and cpar[j] >= 0:
                H["parent_idx"][gidx[j]] = gidx[cpar[j]]
        # genes
        for i in range(n):
            if kind[i] == GENE:
                g = gene[i]
                b = cross_base[i]
                if b < 0 or not T.is_leaf[cz[b]]:
                    raise AssertionError(f"gene {gnames[g]} not placed at its species (bug)")
                self.gene_family[g] = fam_id
                self.gene_parent[g] = gidx[cpar[b]] if cpar[b] >= 0 else -1
        # gains (one per family) and duplications at the origin
        E = self.E
        if origin:
            Zs = {cz[j] for j in origin}
            if len(Zs) != 1:
                raise AssertionError(f"family {fam_id}: several origins (bug)")
            Z = Zs.pop()
            P = par_sp[Z]
            self.gains[Z] += 1
            self.dups_origin[Z] += len(origin) - 1
            fmask = 0
            for j in origin:
                fmask |= mask[cv[j]]
            n_fam = sum(ng[cv[j]] for j in origin)
            present = sum(1 for c in T.children[Z] if fmask & T.clade_mask[c])
            if T.is_leaf[Z]:
                self.leaf_specific[Z] += n_fam
            if present >= 2:
                self.gains_multi[Z] += 1
            reps = ";".join(heapq.nsmallest(R, itertools.chain.from_iterable(rep[cv[j]] for j in origin)))
            for ev, nev in (("gain", 1), ("dup", len(origin) - 1)):
                if nev <= 0:
                    continue
                E["branch_parent_idx"].append(P)
                E["branch_child_idx"].append(Z)
                E["event"].append(ev)
                E["n_events"].append(nev)
                E["at_origin"].append(True)
                E["hog_id"].append(cid[origin[0]])
                E["family_id"].append(fam_id)
                E["child_hog_ids"].append(";".join(cid[j] for j in origin))
                E["n_extant_genes"].append(n_fam)
                E["n_genes_in_child"].append(n_fam)
                E["n_species"].append(fmask.bit_count())
                E["rep_genes"].append(reps)
                E["n_child_clades_present"].append(present)
                E["n_child_clades_total"].append(len(T.children[Z]))
        # losses and duplications of ancestral genes at internal nodes
        kids_of: dict[int, list[int]] = {}
        for j in range(nc):
            if cpar[j] >= 0:
                kids_of.setdefault(cpar[j], []).append(j)
        for j in range(nc):
            Pz = cz[j]
            if T.is_leaf[Pz]:
                continue
            i = cv[j]
            by_child: dict[int, list[int]] = {}
            for k in kids_of.get(j, ()):
                by_child.setdefault(cz[k], []).append(k)
            for c in T.children[Pz]:
                ks = by_child.get(c)
                if ks is None:
                    self.losses[c] += 1
                    ev, nev, ch_ids, n_in = "loss", 1, "", 0
                elif len(ks) > 1:
                    self.dups[c] += len(ks) - 1
                    ev, nev = "dup", len(ks) - 1
                    ch_ids = ";".join(cid[k] for k in ks)
                    n_in = sum(ng[cv[k]] for k in ks)
                else:
                    continue
                E["branch_parent_idx"].append(Pz)
                E["branch_child_idx"].append(c)
                E["event"].append(ev)
                E["n_events"].append(nev)
                E["at_origin"].append(False)
                E["hog_id"].append(cid[j])
                E["family_id"].append(fam_id)
                E["child_hog_ids"].append(ch_ids)
                E["n_extant_genes"].append(ng[i])
                E["n_genes_in_child"].append(n_in)
                E["n_species"].append(mask[i].bit_count())
                E["rep_genes"].append(";".join(rep[i]))
                E["n_child_clades_present"].append(-1)
                E["n_child_clades_total"].append(len(T.children[Pz]))

    # -- finish ------------------------------------------------------------------
    def finish(self) -> HogResult:
        T = self.T
        # ungrouped genes: species-specific singletons (one gain each on the leaf branch)
        k = 0
        for g in range(len(self.gene_names)):
            if not self.seen[g]:
                F = _Family()
                F.add(GENE, -1, gene=g, xid=self.gene_names[g])
                self.family(F, -1)
                k += 1
        self.st["n_singleton_genes"] = k
        names = T.names
        rows = []
        for c in range(len(T)):
            p = T.parent[c]
            n_par = self.n_anc[p] if p >= 0 else 0
            rows.append(dict(branch_parent=names[p] if p >= 0 else "", branch_child=names[c],
                             child_is_leaf=T.is_leaf[c], depth_child=T.depth[c],
                             n_parent=n_par, n_child=self.n_anc[c], gains=self.gains[c],
                             losses=self.losses[c], dups=self.dups[c],
                             dups_at_origin=self.dups_origin[c],
                             leaf_specific_genes=self.leaf_specific[c],
                             gains_ge2_child_clades=self.gains_multi[c]))
        br = pd.DataFrame(rows)
        br["conservation_ok"] = br["n_child"] == (br["n_parent"] - br["losses"] + br["dups"]
                                                  + br["dups_at_origin"] + br["gains"])
        lvl_cat = pd.Categorical.from_codes(np.asarray(self.H["level_idx"], dtype=np.int32),
                                            categories=names)
        hogs = pd.DataFrame({
            "hog_id": self.H["hog_id"], "level": lvl_cat,
            "level_idx": np.asarray(self.H["level_idx"], dtype=np.int32),
            "family_id": self.H["family_id"],
            "parent_idx": np.asarray(self.H["parent_idx"], dtype=np.int64),
            "n_genes": np.asarray(self.H["n_genes"], dtype=np.int64),
            "n_species": np.asarray(self.H["n_species"], dtype=np.int32),
            "n_child_clades_present": np.asarray(self.H["n_child_clades_present"], dtype=np.int32),
            "rep_genes": self.H["rep_genes"]})
        pidx = hogs["parent_idx"].to_numpy()
        hid = hogs["hog_id"].to_numpy(dtype=object)
        hogs.insert(4, "parent_hog_id", [hid[x] if x >= 0 else None for x in pidx])
        hogs.insert(5, "origin", pidx < 0)
        gp = np.asarray(self.gene_parent, dtype=np.int64)
        sp_names = [names[s] for s in self.gene_sp]
        genes = pd.DataFrame({
            "gene": self.gene_names,
            "species": pd.Categorical(sp_names, categories=[names[s] for s in T.leaves]),
            "family_id": self.gene_family, "parent_idx": gp,
            "parent_hog_id": [hid[x] if x >= 0 else None for x in gp],
            "species_specific": gp < 0})
        E = self.E
        bp = np.asarray(E["branch_parent_idx"], dtype=np.int64)
        bc = np.asarray(E["branch_child_idx"], dtype=np.int64)
        nm = np.asarray(names + [""], dtype=object)
        events = pd.DataFrame({
            "branch_parent": nm[bp], "branch_child": nm[bc],
            "child_is_leaf": np.asarray([T.is_leaf[c] for c in E["branch_child_idx"]], dtype=bool),
            "event": E["event"], "n_events": np.asarray(E["n_events"], dtype=np.int64),
            "at_origin": np.asarray(E["at_origin"], dtype=bool),
            "hog_id": E["hog_id"], "family_id": E["family_id"], "child_hog_ids": E["child_hog_ids"],
            "n_extant_genes": np.asarray(E["n_extant_genes"], dtype=np.int64),
            "n_genes_in_child": np.asarray(E["n_genes_in_child"], dtype=np.int64),
            "n_species": np.asarray(E["n_species"], dtype=np.int32),
            "n_child_clades_present": np.asarray(E["n_child_clades_present"], dtype=np.int32),
            "n_child_clades_total": np.asarray(E["n_child_clades_total"], dtype=np.int32),
            "rep_genes": E["rep_genes"]})
        st = dict(self.st)
        st.update(n_genes=len(self.gene_names), n_species=len(T.leaves),
                  n_internal_nodes=len(T.internal), n_hogs_internal=len(hogs),
                  n_event_rows=len(events), conservation_violations=int((~br["conservation_ok"]).sum()))
        return HogResult(tree=T, branches=br, hogs=hogs, genes=genes, events=events, stats=st)


def _open(src):
    if isinstance(src, (str, os.PathLike)):
        src = os.fspath(src)
        return (gzip.open(src, "rb") if src.endswith(".gz") else open(src, "rb")), True
    return src, False


def parse_hogs(src: str | IO[bytes], newick: str | None = None, n_rep: int = 5,
               progress: Callable[[int], None] | None = None) -> HogResult:
    """Stream an OrthoXML file and compute ancestral genes and branch events.

    Parameters
    ----------
    src : path (``.gz`` allowed) or binary file object.
    newick : path or Newick string of the guide tree used for the run. If
        omitted, the file's ``<taxonomy>`` is used. If both are present they
        must agree (same named clades after Ancestra-style pruning).
    n_rep : number of representative gene ids kept per HOG / event row.
    progress : optional callback called with the number of families processed
        (every 10,000 families).
    """
    from lxml import etree

    fh, own = _open(src)
    gene_names: list[str] = []
    gene_species: list[str] = []
    internal: dict[str, int] = {}
    species_order: list[str] = []
    species_taxon: dict[str, str] = {}
    taxa_nodes: dict[str, object] = {}
    tax_stack: list = []
    tax_root = None
    cur_species = None
    builder: _Builder | None = None
    fam: _Family | None = None
    gstack: list[int] = []
    n_fam = 0

    def start_groups():
        nonlocal builder
        from ancestra.newick import Node  # noqa: F401  (taxonomy nodes already built)
        if tax_root is not None:
            rev = {v: k for k, v in species_taxon.items()}
            for tid, nd in taxa_nodes.items():
                if not nd.children and tid in rev:
                    nd.name = rev[tid]
        tree = load_species_tree(newick, species_order, tax_root)
        missing = [s for s in species_order if s not in tree.index or not tree.is_leaf[tree.index[s]]]
        if missing:
            raise ValueError(f"{len(missing)} OrthoXML species are not leaves of the species tree, "
                             f"e.g. {missing[:5]}")
        taxon_name = {tid: nd.name for tid, nd in taxa_nodes.items()}
        gsp = [tree.index[s] for s in gene_species]
        builder = _Builder(tree, gene_names, gsp, taxon_name, n_rep)

    try:
        for event, el in etree.iterparse(fh, events=("start", "end"), huge_tree=True):
            tag = _local(el.tag)
            if event == "start":
                if tag in _GROUP_TAGS:
                    if builder is None:
                        start_groups()
                    k = _GROUP_TAGS[tag]
                    if not gstack:
                        fam = _Family()
                        i = fam.add(k, -1, xid=el.get("id"), taxid=el.get("taxonId"))
                    else:
                        i = fam.add(k, gstack[-1], xid=el.get("id"), taxid=el.get("taxonId"))
                    gstack.append(i)
                elif tag == "species":
                    cur_species = el.get("name")
                    if cur_species in species_taxon or cur_species in species_order:
                        raise ValueError(f"species {cur_species!r} declared twice")
                    species_order.append(cur_species)
                    if el.get("taxonId") is not None:
                        species_taxon[cur_species] = el.get("taxonId")
                elif tag == "taxon":
                    from ancestra.newick import Node
                    nd = Node(el.get("name") or el.get("id"))
                    taxa_nodes[el.get("id")] = nd
                    if tax_stack:
                        tax_stack[-1].add_child(nd)
                    else:
                        tax_root = nd
                    tax_stack.append(nd)
                continue
            # end events
            if tag == "geneRef":
                if gstack:
                    g = internal.get(el.get("id"))
                    if g is None:
                        raise KeyError(f"geneRef id={el.get('id')!r} not declared")
                    if builder.seen[g]:
                        raise ValueError(f"gene {gene_names[g]!r} referenced twice")
                    builder.seen[g] = True
                    fam.add(GENE, gstack[-1], gene=g)
            elif tag == "property":
                if gstack and el.get("name") == "TaxRange":
                    fam.hint[gstack[-1]] = el.get("value")
            elif tag in _GROUP_TAGS:
                gstack.pop()
                if not gstack:
                    builder.family(fam, n_fam)
                    n_fam += 1
                    fam = None
                    if progress is not None and n_fam % 10000 == 0:
                        progress(n_fam)
                    el.clear()
                    while el.getprevious() is not None:
                        del el.getparent()[0]
            elif tag == "gene":
                pid = el.get("protId") or el.get("geneId") or el.get("id")
                if el.get("id") in internal:
                    raise ValueError(f"duplicate gene id {el.get('id')!r}")
                internal[el.get("id")] = len(gene_names)
                gene_names.append(pid)
                gene_species.append(cur_species)
                el.clear()
                while el.getprevious() is not None:
                    del el.getparent()[0]
            elif tag == "species":
                cur_species = None
                el.clear()
            elif tag == "taxon":
                tax_stack.pop()
    finally:
        if own:
            fh.close()
    if builder is None:          # file without any group
        start_groups()
    if len(set(gene_names)) != len(gene_names):
        raise ValueError("gene identifiers (protId) are not unique")
    res = builder.finish()
    res.stats["source"] = os.fspath(src) if isinstance(src, (str, os.PathLike)) else "<stream>"
    return res
