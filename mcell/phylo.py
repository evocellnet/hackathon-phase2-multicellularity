"""Tree utilities for the fixed-topology branch-length fits (DESIGN.md §3.1, §4.5).

The guide trees are rooted, have named internal nodes and may contain
polytomies. IQ-TREE ``-te`` keeps the user topology (multifurcations included;
checked on IQ-TREE 3.1.4) but returns an *unrooted* tree, so lengths are
transferred back to the named rooted guide tree by matching bipartitions:

* every non-root guide edge is identified by the leaf set below it;
* if the root is bifurcating its two edges form one unrooted edge, whose length
  is split equally between them (the root position is not identifiable from
  the data; root-child branches are excluded from every control pool);
* the unrooted tree must not contain any other bipartition (that would mean the
  topology was changed); this is asserted.

Pruning keeps a node of the guide tree iff >= 2 of its children have present
tips; a kept node therefore carries the guide name of the MRCA of its present
tips, and the pruned edge into it spans the guide path from that node up to
(but excluding) the next kept ancestor.
"""
from __future__ import annotations

import ete3

__all__ = ["read_tree", "prune_named", "transfer_branch_lengths", "write_named", "path_nodes",
           "assert_same_topology"]


def read_tree(path_or_str) -> ete3.Tree:
    txt = path_or_str
    try:
        with open(path_or_str) as fh:
            txt = fh.read()
    except (OSError, TypeError):
        pass
    return ete3.Tree(txt.strip(), format=1)


def prune_named(tree: ete3.Tree, keep) -> ete3.Tree:
    """Copy of ``tree`` restricted to the leaves in ``keep``; unary nodes removed.

    Branch lengths are summed along collapsed paths; the kept internal node of a
    collapsed path is the lower one (the MRCA of the present tips), so names are
    always guide names of MRCAs. The new root is the MRCA of ``keep``.
    """
    keep = set(keep)

    def rec(n):
        if n.is_leaf():
            if n.name in keep:
                c = ete3.TreeNode(name=n.name, dist=n.dist)
                return c
            return None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].dist += n.dist
            return kids[0]
        new = ete3.TreeNode(name=n.name, dist=n.dist)
        for k in kids:
            new.add_child(k)
        return new

    r = rec(tree)
    if r is None or r.is_leaf():
        raise ValueError("fewer than 2 kept leaves")
    r.dist = 0.0
    return r


def transfer_branch_lengths(guide: ete3.Tree, unrooted: ete3.Tree, tol: float = 1e-5) -> dict:
    """Map lengths of ``unrooted`` (same leaf set) onto the named nodes of ``guide``.

    Returns {node name: length} for every non-root node of ``guide``. Raises if
    a guide bipartition is missing from ``unrooted`` or if ``unrooted`` has an
    extra non-trivial bipartition (topology changed).
    """
    leaves = set(guide.get_leaf_names())
    if set(unrooted.get_leaf_names()) != leaves:
        raise ValueError("leaf sets differ")
    ref = sorted(leaves)[0]

    def canon(side):
        side = frozenset(side)
        return side if ref not in side else frozenset(leaves - side)

    ulen = {}
    for n in unrooted.traverse():
        if n.up is None:
            continue
        k = canon(n.get_leaf_names())
        ulen[k] = ulen.get(k, 0.0) + n.dist  # a bifurcating root in `unrooted` gives 2 edges with one key
    out, used = {}, set()
    root_kids = guide.children
    for n in guide.traverse():
        if n.up is None:
            continue
        k = canon(n.get_leaf_names())
        if k not in ulen:
            raise ValueError(f"guide edge into {n.name} not found in the fitted tree (topology changed?)")
        L = ulen[k]
        if n.up is guide and len(root_kids) == 2:
            L = L / 2.0
        out[n.name] = L
        used.add(k)
    extra = [k for k in ulen if k not in used and 1 < len(k) < len(leaves) - 1]
    if extra:
        raise ValueError(f"fitted tree has {len(extra)} bipartitions absent from the guide tree")
    return out


def assert_same_topology(a: ete3.Tree, b: ete3.Tree) -> None:
    """Assert identical rooted topology and identical node names."""
    def sig(t):
        return {n.name: frozenset(n.get_leaf_names()) for n in t.traverse()}
    sa, sb = sig(a), sig(b)
    if sa != sb:
        diff = [k for k in set(sa) | set(sb) if sa.get(k) != sb.get(k)]
        raise AssertionError(f"topologies/names differ at {diff[:5]}")


def write_named(tree: ete3.Tree, lengths: dict | None = None) -> str:
    """Newick with internal names and lengths (format 1)."""
    t = tree.copy()
    if lengths is not None:
        for n in t.traverse():
            if n.up is not None:
                n.dist = float(lengths[n.name])
    return t.write(format=1, format_root_node=True)


def path_nodes(guide: ete3.Tree, child: str, upper: str) -> list:
    """Guide nodes whose edges make up the pruned edge upper -> child (child first)."""
    nodes = {n.name: n for n in guide.traverse()}
    n, out = nodes[child], []
    while n is not None and n.name != upper:
        out.append(n.name)
        n = n.up
    if n is None:
        raise ValueError(f"{upper} is not an ancestor of {child}")
    return out
