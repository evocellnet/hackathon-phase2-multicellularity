"""Known-answer tests for the EukaryOMA cross-validation helpers (scripts/eukaryoma/xval_lib.py).

Toy tree  (((A,B)AB,C)ABC,(D,E)DE)R ; we keep A,B,C,D (drop E), so DE disappears (unary) and
D attaches directly to R.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "eukaryoma"))

xl = pytest.importorskip("xval_lib")

NWK = "(((A,B)AB,C)ABC,(D,E)DE)R;"

XML = """<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.5">
  <species name="spA" taxonId="1"><database name="d" version="1"><genes>
    <gene id="1" protId="a1"/><gene id="2" protId="a2"/><gene id="3" protId="a3"/></genes></database></species>
  <species name="spB" taxonId="2"><database name="d" version="1"><genes>
    <gene id="4" protId="b1"/><gene id="5" protId="b2"/></genes></database></species>
  <species name="spC" taxonId="3"><database name="d" version="1"><genes>
    <gene id="6" protId="c1"/></genes></database></species>
  <species name="spD" taxonId="4"><database name="d" version="1"><genes>
    <gene id="7" protId="d1"/></genes></database></species>
  <species name="spE" taxonId="5"><database name="d" version="1"><genes>
    <gene id="8" protId="e1"/><gene id="9" protId="e2"/></genes></database></species>
  <groups>
    <orthologGroup id="F1" taxonId="99"><score id="s" value="1"/><property name="TaxRange" value="Root"/>
      <orthologGroup id="F1a"><property name="TaxRange" value="ABC"/>
        <paralogGroup id="F1p" taxonId="7"><geneRef id="1"/><geneRef id="2"/></paralogGroup>
        <geneRef id="4"/><geneRef id="6"/>
      </orthologGroup>
      <geneRef id="8"/>
    </orthologGroup>
    <orthologGroup id="F2"><property name="TaxRange" value="DE"/><geneRef id="8"/><geneRef id="9"/></orthologGroup>
    <orthologGroup id="F3"><property name="TaxRange" value="AB"/><geneRef id="3"/><geneRef id="5"/></orthologGroup>
  </groups>
</orthoXML>
"""


def test_prune_tree_names_and_topology():
    t = xl.prune_tree(xl.Tree(NWK, format=1), {"A", "B", "C", "D"})
    assert t.name == "R"
    assert sorted(c.name for c in t.children) == ["ABC", "D"]
    assert sorted(n.name for n in t.traverse() if not n.is_leaf()) == ["AB", "ABC", "R"]


def test_branch_map_status():
    o = xl.Tree(NWK, format=1)
    p = xl.prune_tree(o, {"A", "B", "C", "D"})
    assert xl.branch_map(o, p, "ABC", "AB")["status"] == "exact"
    m = xl.branch_map(o, p, "DE", "D")      # DE collapsed: pruned branch R->D covers R>DE and DE>D
    assert m["status"] == "approx" and m["top"] == "R" and m["extra_orig_edges"] == "R>DE"
    assert xl.branch_map(o, p, "R", "ABC")["status"] == "exact"
    assert xl.branch_map(o, p, "DE", "E")["status"] == "missing"
    p2 = xl.prune_tree(o, {"A", "C", "D"})   # AB unsampled: ABC->AB maps to ABC->A
    m2 = xl.branch_map(o, p2, "ABC", "AB")
    assert m2["status"] == "approx" and m2["bottom"] == "A"


def test_filter_and_parse(tmp_path):
    from mcell.hogs import parse_hogs
    src = tmp_path / "in.xml"
    src.write_text(XML)
    dst = tmp_path / "out.xml"
    st = xl.filter_orthoxml(str(src), str(dst), {"spA": "A", "spB": "B", "spC": "C", "spD": "D"},
                            gene_rename={"a1": "A_1"}, strip_levels=True, log=lambda m: None)
    assert st["species_kept"] == 4 and st["genes_kept"] == 7 and st["groups_out"] == 2  # F2 only had E genes
    txt = dst.read_text()
    assert "TaxRange" not in txt and "taxonId" not in txt and "score" not in txt and 'protId="A_1"' in txt
    t = xl.prune_tree(xl.Tree(NWK, format=1), {"A", "B", "C", "D"})
    nwk = tmp_path / "t.nwk"
    xl.write_newick(t, str(nwk))
    res = parse_hogs(str(dst), newick=str(nwk))
    assert res.check_conservation().empty
    br = res.branches.set_index("branch_child")
    # F1 (after removing E): genes A_1,a2,b1,c1 -> MRCA ABC: one gain on R->ABC; the paralogGroup (a1,a2)
    # sits on the branch above A (MRCA = leaf A) -> 1 dup on AB->A. F3 (a3,b2): gain on ABC->AB.
    assert br.loc["ABC", "gains"] == 1
    assert br.loc["AB", "gains"] == 1
    assert br.loc["A", "dups"] == 1
    # F1 lineage passes AB (A and B members) and is absent in D-side: not present at R
    assert br.loc["D", "losses"] == 0
    # gene D1 is ungrouped -> species-specific gain on R->D
    assert br.loc["D", "gains"] == 1
    hogs = res.hogs
    assert set(hogs.loc[hogs["origin"], "level"].astype(str)) == {"ABC", "AB"}


DST_XML = XML.replace('<orthologGroup id="F3"><property name="TaxRange" value="AB"/><geneRef id="3"/><geneRef id="5"/>',
                      '<orthologGroup id="F3"><geneRef id="3"/><geneRef id="5"/><geneRef id="7"/>')


def _parsed(tmp_path, tag, xml):
    from mcell.hogs import parse_hogs
    src = tmp_path / f"{tag}.xml"
    src.write_text(xml)
    dst = tmp_path / f"{tag}.f.xml"
    xl.filter_orthoxml(str(src), str(dst), {"spA": "A", "spB": "B", "spC": "C", "spD": "D"}, log=lambda m: None)
    nwk = tmp_path / "t.nwk"
    xl.write_newick(xl.prune_tree(xl.Tree(NWK, format=1), {"A", "B", "C", "D"}), str(nwk))
    res = parse_hogs(str(dst), newick=str(nwk))
    out = tmp_path / tag
    res.write(str(out))
    return res.tree, str(out)


def test_classify_gain_and_dup(tmp_path):
    xa = pytest.importorskip("xval_analysis")
    tree, da = _parsed(tmp_path, "a", XML)
    _, de = _parsed(tmp_path, "e", DST_XML)
    A = xa.M("a", tree, d=da)
    E = xa.M("e", tree, d=de)
    ix = tree.index
    famF1 = A.genes.set_index("gene").loc["a2", "family_id"]
    famF3 = A.genes.set_index("gene").loc["a3", "family_id"]
    # F1 originates on R->ABC in both
    assert xa.classify_gain(A, E, famF1, ix["ABC"], {ix["ABC"]}, ix["R"])["cls"] == "same"
    # F3: Ancestra-like src gains it on ABC->AB; dst has (a3,b2,d1) rooted at R -> deeper by 2 nodes
    r = xa.classify_gain(A, E, famF3, ix["AB"], {ix["AB"]}, ix["ABC"])
    assert r["cls"] == "deeper" and r["nodes_deeper"] == 2
    # reverse: dst family originates at R; src has it only at AB -> shallower
    famE3 = E.genes.set_index("gene").loc["d1", "family_id"]
    assert xa.classify_gain(E, A, famE3, ix["R"], {ix["ABC"], ix["D"], ix["R"]}, ix["R"])["cls"] in ("shallower",)
    # dup a1/a2 on AB->A: same in both
    ev = A.events[(A.events["event"] == "dup") & (A.events["branch_child"] == "A")].iloc[0]
    ids = ev["child_hog_ids"].split(";")
    assert xa.classify_dup(A, E, ids, ix["A"], ix["AB"])["cls"] == "same"
