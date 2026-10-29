"""Shared helpers for the EukaryOMA cross-validation (scripts/eukaryoma/).

1. :func:`prune_tree` -- restrict our guide tree to a species subset. A node of the original
   tree is *retained* iff it is a kept leaf or >= 2 of its children contain kept species; every
   retained node keeps its original name, so a pruned node is literally the original node that is
   the MRCA of its kept descendants. Unary nodes disappear (their branches merge).
2. :func:`filter_orthoxml` -- stream an OrthoXML (Ancestra or FastOMA/EukaryOMA), keep only the
   given species (renamed to our codes), optionally rename gene protIds (EukaryOMA ids -> our ids),
   drop geneRefs of removed genes and groups that become empty, and (``strip_levels``) remove
   every level annotation (TaxRange properties, taxonId attributes on groups) so that
   mcell.hogs places every group purely topologically: an orthologGroup at the MRCA of its
   members, a paralogGroup on the branch above that MRCA (latest-possible placement). This is the
   same rule for both methods, on the same tree.
3. :func:`branch_map` -- how an original branch (or path stem) is represented on the pruned
   tree: status 'exact' (same edge set), 'approx' (pruned branch covers additional original
   edges, e.g. the clade's first split is unsampled) or 'missing'.
"""
from __future__ import annotations

import gzip
import os
from xml.sax.saxutils import quoteattr

from ete3 import Tree


# ----------------------------------------------------------------------------- tree

def load_tree(path):
    return Tree(open(path).read().strip(), format=1)


def prune_tree(tree: Tree, keep: set) -> Tree:
    """Return a new ete3 tree on ``keep`` leaves with retained-node names (see module doc)."""
    keep = set(keep)
    missing = keep - set(tree.get_leaf_names())
    if missing:
        raise ValueError(f"species not in tree: {sorted(missing)[:5]}")

    def build(nd):
        if nd.is_leaf():
            if nd.name in keep:
                t = Tree(name=nd.name)
                return t
            return None
        subs = [s for s in (build(c) for c in nd.children) if s is not None]
        if not subs:
            return None
        if len(subs) == 1:
            return subs[0]
        t = Tree(name=nd.name)
        for s in subs:
            t.add_child(s)
        return t

    out = build(tree)
    if out is None:
        raise ValueError("no species kept")
    return out


def write_newick(t: Tree, path):
    s = t.write(format=8)          # ete3 omits the root name
    if t.name and not s.rstrip(";").endswith(")" + t.name):
        s = s.rstrip().rstrip(";") + t.name + ";"
    with open(path, "w") as fh:
        fh.write(s + "\n")


def edges_on_path(tree: Tree, lower: str, upper: str) -> list:
    """Original edges (parent, child) from node ``lower`` up to ancestor ``upper``."""
    nd = tree & lower
    out = []
    while nd.name != upper:
        if nd.up is None:
            raise ValueError(f"{upper} is not an ancestor of {lower}")
        out.append((nd.up.name, nd.name))
        nd = nd.up
    return out


def branch_map(orig: Tree, pruned: Tree, top: str, bottom: str) -> dict:
    """Represent the original path top -> ... -> bottom (one edge or a path stem) on the pruned tree.

    bottom' = retained MRCA of the kept species below ``bottom``; top' = retained MRCA of the kept
    species below ``top`` -- if that equals bottom' (no kept species outside ``bottom`` within
    ``top``) the pruned parent of bottom' is used. The pruned path top' -> bottom' covers a set of
    original edges; status 'exact' iff it equals the original path's edge set.
    """
    kept = set(pruned.get_leaf_names())
    pn = {n.name: n for n in pruned.traverse()}
    sb = [s for s in (orig & bottom).get_leaf_names() if s in kept]
    if not sb:
        return dict(status="missing", reason=f"no kept species below {bottom}")
    st = [s for s in (orig & top).get_leaf_names() if s in kept]
    b2 = pruned.get_common_ancestor(sb).name if len(sb) > 1 else sb[0]
    t2 = pruned.get_common_ancestor(st).name if len(st) > 1 else st[0]
    note = ""
    if t2 == b2:
        par = pn[b2].up
        if par is None:
            return dict(status="missing", reason="no kept outgroup")
        t2 = par.name
        note = f"no kept species in {top} outside {bottom}; parent taken as {t2}"
    # pruned path edges (by node names)
    pp, nd = [], pn[b2]
    while nd.name != t2:
        pp.append((nd.up.name, nd.name))
        nd = nd.up
    covered = set(edges_on_path(orig, b2, t2))
    want = set(edges_on_path(orig, bottom, top))
    status = "exact" if covered == want else "approx"
    extra = sorted(covered - want)
    lost = sorted(want - covered)
    return dict(status=status, top=t2, bottom=b2, pruned_edges=";".join(f"{a}>{b}" for a, b in reversed(pp)),
                n_pruned_edges=len(pp), extra_orig_edges=";".join(f"{a}>{b}" for a, b in extra),
                missing_orig_edges=";".join(f"{a}>{b}" for a, b in lost), n_kept_below=len(sb),
                reason=note)


# ----------------------------------------------------------------------------- OrthoXML filter

_GROUPS = ("orthologGroup", "paralogGroup")


def _local(tag):
    if not isinstance(tag, str):
        return ""
    k = tag.rfind("}")
    return tag[k + 1:] if k >= 0 else tag


def _prune_group(el, keep_ids, strip_levels) -> int:
    """Remove dropped geneRefs / empty subgroups / scores (/ levels); return #kept children."""
    n = 0
    for ch in list(el):
        t = _local(ch.tag)
        if t == "geneRef":
            if ch.get("id") in keep_ids:
                n += 1
            else:
                el.remove(ch)
        elif t in _GROUPS:
            if _prune_group(ch, keep_ids, strip_levels):
                n += 1
            else:
                el.remove(ch)
        elif t == "property" and not strip_levels and ch.get("name") == "TaxRange":
            continue
        else:   # score, notes, other properties, comments
            el.remove(ch)
    if strip_levels:
        for a in ("taxonId",):
            if a in el.attrib:
                del el.attrib[a]
    return n


def filter_orthoxml(src, dst, species_rename: dict, gene_rename: dict | None = None,
                    gene_keep=None, strip_levels: bool = True, log=print) -> dict:
    """See module doc. ``species_rename``: OrthoXML species name -> new name (others dropped).
    ``gene_rename``: protId -> new protId. ``gene_keep``: optional predicate(new_protId) -> bool to
    drop further genes (e.g. only genes mapped in both datasets)."""
    from lxml import etree
    gene_rename = gene_rename or {}
    fh = gzip.open(src, "rb") if str(src).endswith(".gz") else open(src, "rb")
    out = gzip.open(dst, "wt", compresslevel=3) if str(dst).endswith(".gz") else open(dst, "w")
    out.write('<?xml version="1.0" encoding="UTF-8"?>\n'
              '<orthoXML xmlns="http://orthoXML.org/2011/" origin="filtered" originVersion="1" version="0.5">\n')
    keep_ids: set = set()
    cur_sp, cur_keep, cur_genes = None, False, []
    depth = 0
    st = dict(species_kept=0, genes_kept=0, genes_renamed=0, groups_in=0, groups_out=0)
    in_groups = False
    seen_species = set()
    for ev, el in etree.iterparse(fh, events=("start", "end"), huge_tree=True):
        tag = _local(el.tag)
        if ev == "start":
            if tag == "species":
                cur_sp = el.get("name")
                cur_keep = cur_sp in species_rename
                cur_genes = []
            elif tag in _GROUPS:
                depth += 1
            elif tag == "groups" and not in_groups:
                in_groups = True
                out.write("  <groups>\n")
            continue
        if tag == "gene" and cur_sp is not None:
            if cur_keep:
                pid = el.get("protId") or el.get("geneId") or el.get("id")
                new = gene_rename.get(pid, pid)
                if gene_keep is None or gene_keep(new):
                    cur_genes.append((el.get("id"), new))
                    st["genes_renamed"] += new != pid
            el.clear()
            while el.getprevious() is not None:
                del el.getparent()[0]
        elif tag == "species":
            if cur_keep:
                seen_species.add(cur_sp)
                st["species_kept"] += 1
                out.write(f'  <species name={quoteattr(species_rename[cur_sp])} NCBITaxId="0">\n'
                          '    <database name="db" version="1">\n      <genes>\n')
                for gid, pid in cur_genes:
                    keep_ids.add(gid)
                    out.write(f"        <gene id={quoteattr(gid)} protId={quoteattr(pid)}/>\n")
                out.write("      </genes>\n    </database>\n  </species>\n")
                st["genes_kept"] += len(cur_genes)
            cur_sp, cur_keep, cur_genes = None, False, []
            el.clear()
            while el.getprevious() is not None:
                del el.getparent()[0]
        elif tag in _GROUPS:
            depth -= 1
            if depth == 0:
                st["groups_in"] += 1
                if _prune_group(el, keep_ids, strip_levels):
                    st["groups_out"] += 1
                    s = etree.tostring(el, encoding="unicode", with_tail=False)
                    out.write("    " + s.replace(' xmlns="http://orthoXML.org/2011/"', "", 1) + "\n")
                el.clear()
                while el.getprevious() is not None:
                    del el.getparent()[0]
                if st["groups_in"] % 200000 == 0:
                    log(f"  {st['groups_in']:,} top-level groups, {st['groups_out']:,} kept")
    if in_groups:
        out.write("  </groups>\n")
    out.write("</orthoXML>\n")
    out.close()
    fh.close()
    st["species_missing"] = sorted(set(species_rename) - seen_species)
    return st
