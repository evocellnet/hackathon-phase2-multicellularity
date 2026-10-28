"""True per-branch gene-content events from ``ancestra.simulate`` gene trees.

This is the reference against which :mod:`mcell.hogs` is validated. It uses a
*different algorithm and a different input*: it walks the simulator's true
reconciled gene trees (``SimResult.gene_trees``; nodes carry
``data["kind"]`` in {S, D, T, leaf} and ``data["species"]``) and enumerates
events directly, instead of parsing OrthoXML and counting crossings.

The simulator records every event (``events.tsv``), including events on
lineages that later went extinct. Those are invisible in extant genomes, so the
reference is the *observable* history: the gene trees pruned to surviving genes
(as the simulator stores them, unary nodes collapsed). Under this convention a
family surviving only in a subclade "originates" at the MRCA of its survivors,
and a family transferred horizontally (subtree marked ``xfer``) is a separate
family originating where it was received, as in the simulator's true HOGs.

Rules (positions as in :mod:`mcell.hogs`):

* family root S at Y / D at c / leaf at s -> one gain on the branch into Y / c / s;
* D node at c with m children -> m-1 duplications on the branch into c
  (``dups_at_origin`` if the D node is the family root or only preceded by
  other D nodes on the same branch);
* S node at Y -> each child clade of Y without a gene-tree child: one loss;
* a gene lineage that passes a species node Z without speciating (because the
  other descendants of Z lost it) -> one loss on every other child branch of Z;
* ancestral genes at X = S nodes at X + lineages passing X (+ genes at a leaf).
"""
from __future__ import annotations

import pandas as pd

from .hogs import SpeciesTree

__all__ = ["true_branch_events", "true_ancestral_genes"]


def _effective_children(nd):
    """Children of ``nd`` within its own HOG family.

    Transferred subtrees (``xfer``) are cut; T (transfer-donor) nodes are
    spliced out; S/D nodes left with a single child after the cut are unary and
    are replaced by that child (the simulator's true HOGs collapse them too).
    """
    out = []
    stack = list(reversed(nd.children))
    while stack:
        c = stack.pop()
        if "xfer" in c.data:
            continue
        if c.data.get("kind") == "T":
            stack.extend(reversed(c.children))
            continue
        if c.children:
            kids = _effective_children(c)
            if len(kids) == 1:
                stack.append(kids[0])
                continue
        out.append(c)
    return out


def _family_roots(root):
    """Roots of the HOG families inside one simulated gene tree (unary chains
    created by cutting transfers are descended)."""
    roots = [root] + [nd for nd in root.preorder() if nd is not root and "xfer" in nd.data]
    out = []
    for r in roots:
        while r.children:
            kids = _effective_children(r)
            if len(kids) != 1:
                break
            r = kids[0]
        out.append(r)
    return out


def true_ancestral_genes(gene_trees: dict, species_tree) -> dict[str, set]:
    """Internal species node -> set of frozensets (extant genes of each true ancestral gene)."""
    rec: dict = {}
    true_branch_events(gene_trees, species_tree, record=rec)
    return rec


def true_branch_events(gene_trees: dict, species_tree, record: dict | None = None) -> pd.DataFrame:
    """Per-branch counts (same columns as ``HogResult.branches``) from true gene trees.

    Parameters
    ----------
    gene_trees : family -> root :class:`ancestra.newick.Node` (``SimResult.gene_trees``).
    species_tree : species-tree root Node (``SimResult.species_tree``).
    record : if a dict is given, it is filled with internal node name -> set of
        frozensets of extant gene names (one per true ancestral gene).
    """
    T = SpeciesTree(species_tree)
    n = len(T)
    n_anc = [0] * n
    gains = [0] * n
    losses = [0] * n
    dups = [0] * n
    dups_origin = [0] * n
    leaf_specific = [0] * n
    gains_multi = [0] * n

    def loc(nd):
        """(species node, kind) of a gene-tree node; kind in {S, D, L}."""
        k = nd.data.get("kind", "leaf")
        s = T.index[nd.data["species"]]
        if not nd.children or k == "leaf":
            return s, "L"
        return s, k

    def family_leaves(r):
        """Extant genes of the family rooted at ``r`` (transferred subtrees excluded)."""
        out, stack = [], [r]
        while stack:
            x = stack.pop()
            if not x.children:
                out.append(x)
            else:
                stack.extend(_effective_children(x))
        return out

    def note(Z: int, v) -> None:
        if record is not None and not T.is_leaf[Z]:
            record.setdefault(T.names[Z], set()).add(frozenset(x.name for x in family_leaves(v)))

    def walk(Z: int, v) -> None:
        """Lineage enters species node Z (or the branch above it) heading to v."""
        vs, vk = loc(v)
        while Z != vs:
            # passes the speciation at Z without splitting: loss in the other children
            n_anc[Z] += 1
            note(Z, v)
            nxt = T.child_towards(Z, vs)
            for c in T.children[Z]:
                if c != nxt:
                    losses[c] += 1
            Z = nxt

    for fam, groot in gene_trees.items():
        for r in _family_roots(groot):
            rs, rk = loc(r)
            gains[rs] += 1
            leaves = family_leaves(r)
            if T.is_leaf[rs]:
                leaf_specific[rs] += len(leaves)
            else:   # gains supported by >= 2 child clades of the origin node
                fm = 0
                for x in leaves:
                    fm |= T.clade_mask[T.index[x.data["species"]]]
                if sum(1 for c in T.children[rs] if fm & T.clade_mask[c]) >= 2:
                    gains_multi[rs] += 1
            # traverse the family
            stack = [(r, True)]         # (node, still on the origin branch)
            while stack:
                u, at_origin = stack.pop()
                us, uk = loc(u)
                if uk == "L":
                    n_anc[us] += 1
                    continue
                kids = _effective_children(u)
                if uk == "S":
                    n_anc[us] += 1
                    note(us, u)
                    covered = set()
                    for v in kids:
                        vs, _ = loc(v)
                        c = T.child_towards(us, vs)
                        covered.add(c)
                        walk(c, v)
                        stack.append((v, False))
                    for c in T.children[us]:
                        if c not in covered:
                            losses[c] += 1
                else:  # D at us: duplication on the branch into us
                    if at_origin:
                        dups_origin[us] += len(kids) - 1
                    else:
                        dups[us] += len(kids) - 1
                    for v in kids:
                        walk(us, v)
                        vs, vk = loc(v)
                        stack.append((v, at_origin and vk == "D" and vs == us))

    rows = []
    for c in range(n):
        p = T.parent[c]
        rows.append(dict(branch_parent=T.names[p] if p >= 0 else "", branch_child=T.names[c],
                         child_is_leaf=T.is_leaf[c], n_parent=n_anc[p] if p >= 0 else 0,
                         n_child=n_anc[c], gains=gains[c], losses=losses[c], dups=dups[c],
                         dups_at_origin=dups_origin[c], leaf_specific_genes=leaf_specific[c],
                         gains_ge2_child_clades=gains_multi[c]))
    return pd.DataFrame(rows)
