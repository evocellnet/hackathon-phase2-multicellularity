"""Hand-made OrthoXML cases with known gains, losses and duplications.

Species tree (multifurcation at DEF)::

    R
    |-- ABC
    |   |-- AB
    |   |   |-- A
    |   |   `-- B
    |   `-- C
    `-- DEF
        |-- D
        |-- E
        `-- F

Every case is written with Ancestra's own OrthoXML writer (same code path as
``ancestra run``: normalisation, TaxRange for orthologGroups, taxonId for
paralogGroups), parsed by :func:`mcell.hogs.parse_hogs`, and compared with
counts derived by hand.
"""
import gzip
import io

import pandas as pd
import pytest

from ancestra.newick import parse_newick
from ancestra.orthoxml import GeneEntry, Group, OrthoXML, write_orthoxml
from mcell.hogs import parse_hogs

TREE = "(((A,B)AB,C)ABC,(D,E,F)DEF)R;"
SPECIES = list("ABCDEF")
EVENTS = ("gains", "losses", "dups", "dups_at_origin")


def O(level, *kids):
    return Group("ortholog", list(kids), level=level)


def P(level, *kids):
    return Group("paralog", list(kids), level=level)


def genes_of(x):
    return x.genes() if isinstance(x, Group) else [x]


def build(groups, extra_genes=(), taxonomy=True):
    """OrthoXML bytes; gene ``x...`` belongs to species ``X`` (first letter, upper case)."""
    names = []
    for g in groups:
        names += genes_of(g)
    names += list(extra_genes)
    genes = {n: GeneEntry(n, n[0].upper()) for n in names}
    # species without genes are absent from the <species> block (as in Ancestra output);
    # with the embedded <taxonomy> they remain leaves of the species tree
    ox = OrthoXML(genes=genes, groups=list(groups), species_tree=parse_newick(TREE))
    buf = io.StringIO()
    write_orthoxml(ox, buf, include_taxonomy=taxonomy)
    return buf.getvalue().encode()


def run(groups, extra_genes=(), taxonomy=True):
    data = build(groups, extra_genes, taxonomy)
    return parse_hogs(io.BytesIO(data), newick=None if taxonomy else TREE)


def full(tag):
    """A family present in one copy in every species (no events below the root)."""
    return O("R", O("ABC", O("AB", f"a{tag}", f"b{tag}"), f"c{tag}"),
             O("DEF", f"d{tag}", f"e{tag}", f"f{tag}"))


def counts(res):
    return res.branches.set_index("branch_child")[list(EVENTS)]


def expect(res, exp: dict):
    """exp: branch_child -> {event: n}; all other counts must be 0."""
    got = counts(res)
    want = pd.DataFrame(0, index=got.index, columns=list(EVENTS))
    for child, d in exp.items():
        for k, v in d.items():
            want.loc[child, k] = v
    pd.testing.assert_frame_equal(got, want, check_dtype=False)
    assert res.check_conservation().empty
    assert res.branches["conservation_ok"].all()


# ----------------------------------------------------------------------------
# single-family cases
# ----------------------------------------------------------------------------

def test_full_family_only_root_gain():
    res = run([full(1)])
    expect(res, {"R": {"gains": 1}})
    assert (res.branches["n_child"] == 1).all()


def test_gene_lost_in_one_child():
    res = run([O("R", O("ABC", "a1", "c1"), O("DEF", "d1", "e1", "f1"))])
    expect(res, {"R": {"gains": 1}, "B": {"losses": 1}})
    # the ancestral gene at AB contains only the A gene (implicit, level skipped in the file)
    m = res.members("AB")
    assert m["gene"].tolist() == ["a1"]
    ev = res.events.query("event == 'loss'").iloc[0]
    assert (ev.branch_parent, ev.branch_child, ev.hog_id) == ("AB", "B", "AB|a1")


def test_duplication_on_internal_branch():
    fam = O("R", O("ABC", P("AB", O("AB", "a1", "b1"), O("AB", "a2", "b2")), "c1"),
            O("DEF", "d1", "e1", "f1"))
    res = run([fam])
    expect(res, {"R": {"gains": 1}, "AB": {"dups": 1}})
    nb = res.branches.set_index("branch_child")["n_child"]
    assert nb["AB"] == 2 and nb["A"] == 2 and nb["ABC"] == 1
    dup = res.events.query("event == 'dup'").iloc[0]
    assert (dup.branch_parent, dup.branch_child, dup.n_events) == ("ABC", "AB", 1)
    assert len(dup.child_hog_ids.split(";")) == 2
    # each AB copy maps to the single ABC ancestral gene
    h = res.hogs.set_index("hog_id")
    ab = h[h["level"] == "AB"]
    assert len(ab) == 2 and set(ab["parent_hog_id"]) == {h.index[h["level"] == "ABC"][0]}


def test_inparalogs_in_one_species():
    fam = O("R", O("ABC", O("AB", "a1", "b1"), "c1"), O("DEF", P("D", "d1", "d2"), "e1", "f1"))
    res = run([fam])
    expect(res, {"R": {"gains": 1}, "D": {"dups": 1}})
    cn = res.copy_number("R").set_index("species")["n_copies"]
    assert cn.to_dict() == {"A": 1, "B": 1, "C": 1, "D": 2, "E": 1, "F": 1}


def test_species_specific_gene_and_species_specific_paralogs():
    res = run([full(1), P("E", "e8", "e9")], extra_genes=["c7"])
    expect(res, {"R": {"gains": 1}, "C": {"gains": 1}, "E": {"gains": 1, "dups_at_origin": 1}})
    b = res.branches.set_index("branch_child")
    assert b.loc["C", "leaf_specific_genes"] == 1 and b.loc["E", "leaf_specific_genes"] == 2
    g = res.genes.set_index("gene")
    assert g.loc["c7", "species_specific"] and g.loc["e8", "species_specific"]
    assert not g.loc["c1", "species_specific"]


def test_multifurcation_partial_presence():
    res = run([O("DEF", "d1", "f1")])
    expect(res, {"DEF": {"gains": 1}, "E": {"losses": 1}})
    gain = res.events.query("event == 'gain'").iloc[0]
    assert (gain.branch_parent, gain.branch_child) == ("R", "DEF")
    assert gain.n_child_clades_present == 2 and gain.n_child_clades_total == 3
    assert gain.n_species == 2


def test_two_distant_species_origin_at_mrca_and_intermediate_losses():
    res = run([O("R", "a1", "d1")])
    expect(res, {"R": {"gains": 1}, "C": {"losses": 1}, "B": {"losses": 1},
                 "E": {"losses": 1}, "F": {"losses": 1}})
    nb = res.branches.set_index("branch_child")["n_child"]
    assert nb[["R", "ABC", "AB", "A", "DEF", "D"]].tolist() == [1] * 6
    assert nb[["B", "C", "E", "F"]].tolist() == [0] * 4
    assert set(res.hogs["hog_id"]) == {"R|" + res.hogs["family_id"][0], "ABC|a1", "AB|a1", "DEF|d1"}


def test_paralog_without_taxonid_is_placed_at_member_mrca():
    # level=None -> no taxonId attribute -> latest-possible placement (branch ABC -> AB)
    fam = O("R", P(None, O("AB", "a1", "b1"), O("AB", "a2", "b2")), "d1")
    res = run([fam])
    expect(res, {"R": {"gains": 1}, "AB": {"dups": 1}, "C": {"losses": 1},
                 "E": {"losses": 1}, "F": {"losses": 1}})


def test_paralog_with_earlier_taxonid_and_copy_loss():
    # duplication on branch R -> ABC, before the ABC speciation; copy 2 later lost in C
    fam = O("R", P("ABC", O("ABC", O("AB", "a1", "b1"), "c1"), O("AB", "a2", "b2")), "d1")
    res = run([fam])
    expect(res, {"R": {"gains": 1}, "ABC": {"dups": 1}, "C": {"losses": 1},
                 "E": {"losses": 1}, "F": {"losses": 1}})
    nb = res.branches.set_index("branch_child")["n_child"]
    assert nb["ABC"] == 2 and nb["AB"] == 2 and nb["C"] == 1


def test_paralog_taxonid_older_than_member_mrca_is_respected():
    # taxonId ABC = duplication on R -> ABC although both copies survive only in AB:
    # two ancestral genes at ABC, both lost in C, no duplication on ABC -> AB.
    # (Latest-possible placement would instead give 1 gene at ABC and a dup on ABC -> AB.)
    fam = O("R", P("ABC", O("AB", "a1", "b1"), O("AB", "a2", "b2")), "d1")
    res = run([fam])
    expect(res, {"R": {"gains": 1}, "ABC": {"dups": 1}, "C": {"losses": 2},
                 "E": {"losses": 1}, "F": {"losses": 1}})


def test_toplevel_paralog_on_internal_branch():
    fam = P("DEF", O("DEF", "d1", "e1", "f1"), O("DEF", "d2", "e2", "f2"))
    res = run([fam])
    expect(res, {"DEF": {"gains": 1, "dups_at_origin": 1}})
    assert res.branches.set_index("branch_child").loc["DEF", "n_child"] == 2


def test_all_cases_together_add_up():
    cases = [
        ([full(1)], {"R": {"gains": 1}}),
        ([O("R", O("ABC", "a2", "c2"), O("DEF", "d2", "e2", "f2"))], {"R": {"gains": 1}, "B": {"losses": 1}}),
        ([O("DEF", "d3", "f3")], {"DEF": {"gains": 1}, "E": {"losses": 1}}),
        ([O("R", "a4", "d4")], {"R": {"gains": 1}, "C": {"losses": 1}, "B": {"losses": 1},
                                "E": {"losses": 1}, "F": {"losses": 1}}),
        ([O("R", O("ABC", P("AB", O("AB", "a5", "b5"), O("AB", "a6", "b6")), "c5"),
            O("DEF", P("D", "d5", "d6"), "e5", "f5"))], {"R": {"gains": 1}, "AB": {"dups": 1},
                                                          "D": {"dups": 1}}),
    ]
    groups, tot = [], {}
    for gs, exp in cases:
        groups += gs
        for b, d in exp.items():
            for k, v in d.items():
                tot.setdefault(b, {}).setdefault(k, 0)
                tot[b][k] += v
    res = run(groups, extra_genes=["f9"])
    tot.setdefault("F", {})["gains"] = 1
    expect(res, tot)


# ----------------------------------------------------------------------------
# input handling
# ----------------------------------------------------------------------------

def test_newick_instead_of_taxonomy_and_gzip(tmp_path):
    fam = O("R", P("ABC", O("ABC", O("AB", "a1", "b1"), "c1"), O("AB", "a2", "b2")), "d1")
    # without <taxonomy> the paralog taxonId is not written -> latest-possible placement at the
    # MRCA of its members, which is ABC here, so the result equals the taxonomy-based one
    data = build([fam], extra_genes=["e9", "f9"], taxonomy=False)   # all species declared
    p = tmp_path / "x.orthoxml.gz"
    with gzip.open(p, "wb") as fh:
        fh.write(data)
    res = parse_hogs(str(p), newick=TREE)
    expect(res, {"R": {"gains": 1}, "ABC": {"dups": 1}, "C": {"losses": 1},
                 "E": {"losses": 1, "gains": 1}, "F": {"losses": 1, "gains": 1}})


def test_guide_tree_is_pruned_like_ancestra():
    # leaf G has no genes in the file -> pruned; DEFG becomes unary and is collapsed
    data = build([full(1)], taxonomy=False)
    res = parse_hogs(io.BytesIO(data), newick="(((A,B)AB,C)ABC,((D,E,F)DEF,G)DEFG)R;")
    assert "G" not in res.tree.index and "DEFG" not in res.tree.index
    expect(res, {"R": {"gains": 1}})


def test_tree_mismatch_with_taxonomy_raises():
    data = build([full(1)])
    with pytest.raises(ValueError, match="disagree"):
        parse_hogs(io.BytesIO(data), newick="(((A,C)AB,B)ABC,(D,E,F)DEF)R;")


def test_handwritten_xml_quirks():
    """Raw XML: TaxRange on a paralogGroup is ignored, an unknown TaxRange and a
    paralog taxonId that does not contain its members fall back to the MRCA."""
    xml = b"""<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" origin="hand" originVersion="1" version="0.4">
  <species name="A" NCBITaxId="0"><database name="x" version="1"><genes>
    <gene id="1" protId="a1"/><gene id="2" protId="a2"/></genes></database></species>
  <species name="B" NCBITaxId="0"><database name="x" version="1"><genes>
    <gene id="3" protId="b1"/><gene id="4" protId="b2"/></genes></database></species>
  <species name="C" NCBITaxId="0"><database name="x" version="1"><genes>
    <gene id="5" protId="c1"/></genes></database></species>
  <species name="D" NCBITaxId="0"><database name="x" version="1"><genes>
    <gene id="6" protId="d1"/></genes></database></species>
  <species name="E" NCBITaxId="0"><database name="x" version="1"><genes>
    <gene id="7" protId="e1"/></genes></database></species>
  <species name="F" NCBITaxId="0"><database name="x" version="1"><genes>
    <gene id="8" protId="f1"/></genes></database></species>
  <groups>
    <orthologGroup id="G1">
      <property name="TaxRange" value="Nowhere"/>
      <orthologGroup id="G1.0"><property name="TaxRange" value="ABC"/>
        <paralogGroup id="G1.1">
          <property name="TaxRange" value="R"/>
          <orthologGroup id="G1.2"><property name="TaxRange" value="AB"/>
            <geneRef id="1"/><geneRef id="3"/></orthologGroup>
          <orthologGroup id="G1.3"><property name="TaxRange" value="AB"/>
            <geneRef id="2"/><geneRef id="4"/></orthologGroup>
        </paralogGroup>
        <geneRef id="5"/>
      </orthologGroup>
      <orthologGroup id="G1.4"><property name="TaxRange" value="DEF"/>
        <geneRef id="6"/><geneRef id="7"/><geneRef id="8"/></orthologGroup>
    </orthologGroup>
  </groups>
</orthoXML>"""
    res = parse_hogs(io.BytesIO(xml), newick=TREE)
    assert res.stats["n_level_hint_unknown"] == 1
    # G1 -> R (MRCA); paralog -> latest placement on ABC -> AB
    expect(res, {"R": {"gains": 1}, "AB": {"dups": 1}})
    assert set(res.members("AB")["hog_id"]) == {"AB|G1.2", "AB|G1.3"}


def test_two_children_in_same_child_clade_imply_duplication():
    """Not produced by Ancestra, but valid OrthoXML: an orthologGroup at R with two
    children inside the ABC clade means two ancestral genes at ABC (one duplication on
    R -> ABC); each copy is then lost in one of the ABC lineages."""
    fam = O("R", O("AB", "a1", "b1"), "c1", O("DEF", "d1", "e1", "f1"))
    res = run([fam])
    expect(res, {"R": {"gains": 1}, "ABC": {"dups": 1}, "AB": {"losses": 1}, "C": {"losses": 1}})


def test_duplicate_gene_reference_rejected():
    xml = build([O("R", "a1", "d1")]).replace(b'<geneRef id="2"/>', b'<geneRef id="1"/>')
    with pytest.raises(ValueError, match="referenced twice"):
        parse_hogs(io.BytesIO(xml))
