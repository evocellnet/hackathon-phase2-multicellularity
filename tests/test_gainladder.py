"""Tests of :mod:`mcell.gainladder` (gain-evidence ladder, DESIGN §3.2).

Toy cases use hand-made OrthoXML (Ancestra's writer) and hand-made DIAMOND hit
tables with known answers; the slow test runs DIAMOND on a small simulation
with the TRUE HOGs, where every non-HGT internal gain must be sequence-novel.

Species tree::

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
"""
import io
import os
import shutil
import sys

import numpy as np
import pandas as pd
import pytest

from ancestra.newick import parse_newick
from ancestra.orthoxml import GeneEntry, Group, OrthoXML, write_orthoxml
from mcell import gainladder as gl
from mcell.hogs import parse_hogs

TREE = "(((A,B)AB,C)ABC,(D,E,F)DEF)R;"


def O(level, *kids):
    return Group("ortholog", list(kids), level=level)


def _genes(x):
    return x.genes() if isinstance(x, Group) else [x]


def parse(groups):
    names = [n for g in groups for n in _genes(g)]
    genes = {n: GeneEntry(n, n[0].upper()) for n in names}
    buf = io.StringIO()
    write_orthoxml(OrthoXML(genes=genes, groups=list(groups), species_tree=parse_newick(TREE)), buf)
    return parse_hogs(io.BytesIO(buf.getvalue().encode()))


def full(tag):
    return O("R", O("ABC", O("AB", f"a{tag}", f"b{tag}"), f"c{tag}"), O("DEF", f"d{tag}", f"e{tag}", f"f{tag}"))


# X: old family present at R, passes ABC through c1, LOST on ABC->AB
# G : split piece of X gained on ABC->AB (reciprocal best of X)  -> split_repaired
# G2: second piece hitting X more weakly                          -> homologous (not reciprocal)
# N : sequence-novel gain on ABC->AB
# H : gain on R->DEF whose only outside hit is a prokaryote       -> homologous, other-domain only
GROUPS = [full("9"),
          O("R", "c1", O("DEF", "d1", "e1", "f1")),
          O("AB", "a1", "b1"),
          O("AB", "a2", "b2"),
          O("AB", "a3", "b3"),
          O("DEF", "d4", "e4")]
HITS = [  # qseqid sseqid pident length evalue bitscore qcovhsp
    ("a1", "a1", 100, 300, 1e-100, 600, 100), ("a1", "b1", 80, 290, 1e-90, 500, 97),
    ("a1", "c1", 50, 280, 1e-60, 200, 95), ("a1", "d1", 40, 270, 1e-40, 150, 90),
    ("b1", "c1", 48, 200, 1e-50, 180, 95),
    ("a2", "c1", 35, 150, 1e-20, 100, 60), ("a2", "a1", 70, 200, 1e-70, 300, 90),
    ("a3", "a3", 100, 250, 1e-120, 500, 100), ("a3", "b3", 90, 250, 1e-110, 450, 100),
    ("b3", "a3", 90, 250, 1e-110, 450, 100),
    ("d4", "e4", 90, 200, 1e-100, 400, 100), ("d4", "Zpro_1", 40, 150, 1e-30, 120, 70),
]


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    res = parse(GROUPS)
    genes = res.genes["gene"].tolist()
    lengths = {g: 200 for g in genes}
    lengths.update(a1=300, b1=250, a2=210)
    targets = pd.DataFrame({"seqid": genes + ["Zpro_1"],
                            "species": [g[0].upper() for g in genes] + ["Zpro"],
                            "domain": ["euk"] * len(genes) + ["prok"],
                            "length": [lengths[g] for g in genes] + [180]})
    gains = gl.internal_gains(res)
    reps = gl.choose_representatives(res, gains, targets.set_index("seqid")["length"], n_rep=3, seed=11)
    p = tmp_path_factory.mktemp("toy") / "hits.tsv"
    pd.DataFrame(HITS).to_csv(p, sep="\t", header=False, index=False)
    ts = gl.TargetSpace.build(targets, res)
    pq, top = gl.reduce_hits(str(p), reps, gains, res.tree, ts, chunksize=3)  # tiny chunks: test merging
    cls = gl.classify_gains(res, gains, reps, pq, top, ts).set_index("family_id")
    return res, gains, reps, pq, cls


def fam_of(res, gene):
    return res.genes.set_index("gene").loc[gene, "family_id"]


def test_internal_gains_scope(toy):
    res, gains, *_ = toy
    assert set(gains["branch_child"]) == {"AB", "DEF"}
    assert len(gains) == 4                       # G, G2, N on AB; H on DEF (root stem and leaves excluded)


def test_representatives(toy):
    res, gains, reps, *_ = toy
    r = reps[reps["family_id"] == fam_of(res, "a1")].sort_values("rep_rank")
    assert r["gene"].tolist() == ["a1", "b1"]   # longest first; 2 species only
    assert reps.groupby("family_id").size().max() <= 3


def test_reduce_hits_clade_filter(toy):
    res, gains, reps, pq, cls = toy
    q = pq.set_index("gene")
    assert q.loc["a1", "n_hits"] == 4 and q.loc["a1", "n_in_clade"] == 2 and q.loc["a1", "n_outside"] == 2
    assert q.loc["a1", "best_out_target"] == "c1" and q.loc["a1", "n_outside_species"] == 2
    assert q.loc["a3", "n_outside"] == 0 and q.loc["b3", "n_outside"] == 0
    assert q.loc["d4", "n_outside_other_domain"] == 1 and q.loc["d4", "best_out_domain"] == "prok"


def test_labels(toy):
    res, gains, reps, pq, cls = toy
    lab = {g: cls.loc[fam_of(res, g), "ladder"] for g in ("a1", "a2", "a3", "d4")}
    assert lab == {"a1": "split_repaired", "a2": "homologous", "a3": "seq_novel", "d4": "homologous"}
    x = res.genes.set_index("gene").loc["c1", "parent_hog_id"]       # c1's ancestor at ABC
    assert cls.loc[fam_of(res, "a1"), "repair_hog"] == x
    assert cls.loc[fam_of(res, "a2"), "best_parent_hog"] == x        # same X, but not reciprocal
    assert bool(cls.loc[fam_of(res, "d4"), "outside_other_domain_only"])


def test_branch_ladder_counts(toy):
    res, gains, reps, pq, cls = toy
    br = gl.branch_ladder(res, cls.reset_index()).set_index("branch_child")
    ab = br.loc["AB"]
    assert (ab.hog_gains, ab.split_repaired, ab.seq_novel, ab.homologous_unrepaired) == (3, 1, 1, 1)
    assert ab.gains_after_repair == 2 and ab.losses_after_repair == ab.losses - 1
    assert ab.losses == 1                         # X lost on ABC->AB
    assert br.loc["DEF", "seq_novel"] == 0 and "ABC" in br.index and br.loc["ABC", "hog_gains"] == 0


def test_representatives_span_root_split():
    # family gained at R->DEF with several genes per species; d-genes longest
    res = parse([full("9"), O("DEF", O("DEF", "d1", "e1", "f1"), O("DEF", "d2", "e2"))])
    gains = gl.internal_gains(res)
    L = {g: 100 for g in res.genes["gene"]}
    L.update(d1=500, d2=400, e1=300, e2=350, f1=50)
    reps = gl.choose_representatives(res, gains, pd.Series(L), n_rep=3, seed=11)
    r = reps.sort_values("rep_rank")
    assert r["gene"].iloc[0] == "d1"
    assert len(set(r["species"])) == 3           # D, E, F: one per child clade
    assert r[r.species == "E"]["gene"].item() == "e2"   # longest gene of the chosen species
    again = gl.choose_representatives(res, gains, pd.Series(L), n_rep=3, seed=11)
    assert again.equals(reps)                    # deterministic


def test_architectures():
    d = pd.DataFrame({"protein": ["p1", "p1", "p1", "p1", "p2"],
                      "pfam_acc": ["PF1.2", "PF1.2", "PF9.1", "PF2.3", "PF2.3"],
                      "ali_from": [1, 50, 55, 200, 10], "ali_to": [40, 90, 95, 260, 70],
                      "dom_score": [30, 25, 20, 40, 50]})
    r = gl._resolve_overlaps(d)
    assert "PF9.1" not in r["pfam_acc"].tolist()  # overlaps the better PF1 hit
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, "X.pfam.tsv.gz")
        d.to_csv(f, sep="\t", index=False)
        a = gl.pfam_architectures({"X": f}).set_index("protein")["arch"]
    assert a["p1"] == "PF1-PF2" and a["p2"] == "PF2"   # tandem repeat collapsed


@pytest.mark.slow
def test_true_hogs_internal_gains_are_novel(tmp_path):
    """End to end with DIAMOND on a small simulation using the TRUE HOGs: families are
    independent random sequences, so every true internal gain that is not an HGT
    receipt must be sequence-novel, and nothing is split-repaired."""
    exe = shutil.which("diamond") or os.path.join(os.path.dirname(sys.executable), "diamond")
    if not os.path.exists(exe):
        pytest.skip("diamond not available")
    from ancestra.simulate import SimParams, simulate, write_simulation
    sim = simulate(SimParams(n_species=8, n_genes=300, loss_rate=0.4), 5)
    write_simulation(sim, str(tmp_path / "sim"))
    buf = io.StringIO()
    write_orthoxml(sim.orthoxml, buf)
    res = parse_hogs(io.BytesIO(buf.getvalue().encode()))
    db = gl.make_diamond_db({"sim": str(tmp_path / "sim" / "proteomes")}, str(tmp_path / "db" / "sim"),
                            threads=2, diamond=exe)
    targets = pd.read_parquet(db["targets"])
    gains = gl.internal_gains(res)
    assert len(gains) > 5
    reps = gl.choose_representatives(res, gains, targets.set_index("seqid")["length"])
    gl.write_queries(reps, [db["faa"]], str(tmp_path / "q.faa"))
    gl.run_diamond(str(tmp_path / "q.faa"), db["dmnd"], str(tmp_path / "hits.tsv.gz"), threads=2,
                   sensitivity="very-sensitive", diamond=exe, extra=["--dbsize", "1489260785"])
    ts = gl.TargetSpace.build(targets, res)
    pq, top = gl.reduce_hits(str(tmp_path / "hits.tsv.gz"), reps, gains, res.tree, ts)
    cls = gl.classify_gains(res, gains, reps, pq, top, ts)
    assert (pq["n_in_clade"] > 0).all()           # every representative finds its own family
    assert cls["seq_novel"].all(), cls[~cls["seq_novel"]].to_string()
    assert not cls["split_repaired"].any()
