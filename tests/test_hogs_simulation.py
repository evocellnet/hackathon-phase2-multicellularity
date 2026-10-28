"""Validation of :mod:`mcell.hogs` on ``ancestra.simulate`` data with known history.

(a) The simulator's TRUE HOGs, written with Ancestra's OrthoXML writer and read
    back by :func:`mcell.hogs.parse_hogs`, must give exactly the per-branch
    counts and ancestral gene sets that :mod:`mcell.simtruth` derives
    independently from the true reconciled gene trees.
(b) (slow) ``ancestra run --preset fast`` on a small simulation must parse
    cleanly and satisfy the conservation law; accuracy against the truth is
    measured by ``scripts/validate_hogparser_sim.py`` (not asserted here).
"""
import io
import os
import shutil
import subprocess
import sys

import pytest

from ancestra.orthoxml import write_orthoxml
from ancestra.simulate import SimParams, simulate, write_simulation
from mcell.hogs import parse_hogs
from mcell.simtruth import true_ancestral_genes, true_branch_events

COLS = ["n_parent", "n_child", "gains", "losses", "dups", "dups_at_origin",
        "leaf_specific_genes", "gains_ge2_child_clades"]

MULTI_TREE = ("((sp0:0.6,sp1:0.6,sp2:0.6)X:0.5,(sp3:0.7,(sp4:0.5,sp5:0.5,sp6:0.5)Y:0.2)Z:0.4,"
              "sp7:1.1)R;")

CONFIGS = {
    "default_hgt": (11, dict(n_species=12, n_genes=600, hgt_rate=0.01)),
    "high_turnover_wgd_hgt": (5, dict(n_species=10, n_genes=400, dup_rate=0.4, loss_rate=0.4,
                                      wgd_branch="random", hgt_rate=0.03)),
    "multifurcating_tree": (3, dict(newick=MULTI_TREE, n_genes=400, dup_rate=0.3, loss_rate=0.5)),
    "many_species": (7, dict(n_species=40, n_genes=150, dup_rate=0.3, loss_rate=0.3)),
}


def _parse_truth(res):
    buf = io.StringIO()
    write_orthoxml(res.orthoxml, buf)
    return parse_hogs(io.BytesIO(buf.getvalue().encode()))


@pytest.fixture(scope="module", params=sorted(CONFIGS))
def sim(request):
    seed, kw = CONFIGS[request.param]
    res = simulate(SimParams(**kw), seed)
    return res, _parse_truth(res)


def test_true_hogs_reproduce_true_branch_events(sim):
    res, parsed = sim
    truth = true_branch_events(res.gene_trees, res.species_tree).set_index("branch_child")[COLS]
    got = parsed.branches.set_index("branch_child")[COLS].loc[truth.index]
    diff = (truth - got).abs()
    assert diff.values.sum() == 0, "\n" + truth.join(got, rsuffix="_parser").to_string()
    # events actually happened (the test is not vacuous)
    assert truth["losses"].sum() > 0 and truth["dups"].sum() > 0
    assert truth.drop(index=parsed.tree.names[0])["gains"].sum() > 0


def test_true_hogs_reproduce_ancestral_gene_sets(sim):
    res, parsed = sim
    truth = true_ancestral_genes(res.gene_trees, res.species_tree)
    for level in [parsed.tree.names[i] for i in parsed.tree.internal]:
        m = parsed.members(level)
        got = {frozenset(x[1]) for x in m.groupby("hog_id")["gene"]}
        assert got == truth.get(level, set()), level


def test_conservation_and_invariants(sim):
    res, parsed = sim
    b = parsed.branches
    assert parsed.check_conservation().empty
    assert (b[["gains", "losses", "dups", "dups_at_origin", "n_child"]] >= 0).all().all()
    # each gene appears once; ancestral genes at leaves are the extant genes
    g = parsed.genes
    assert g["gene"].is_unique and len(g) == sum(len(v) for v in res.proteomes.values())
    leaf = b[b["child_is_leaf"]].set_index("branch_child")["n_child"]
    assert leaf.to_dict() == g.groupby("species", observed=True).size().to_dict()
    # event table sums equal the branch table
    ev = parsed.events
    for name, col in (("gain", "gains"), ("loss", "losses")):
        s = ev[ev.event == name].groupby("branch_child")["n_events"].sum()
        assert s.to_dict() == {k: v for k, v in b.set_index("branch_child")[col].items() if v}
    dup = ev[(ev.event == "dup") & ~ev.at_origin].groupby("branch_child")["n_events"].sum()
    assert dup.to_dict() == {k: v for k, v in b.set_index("branch_child")["dups"].items() if v}
    # HOG parent pointers go exactly one level up
    h = parsed.hogs
    T = parsed.tree
    has = h["parent_idx"] >= 0
    lv = h["level_idx"].to_numpy()
    assert all(T.parent[a] == lv[p] for a, p in zip(lv[has], h["parent_idx"][has]))
    # copy numbers at the root sum to the genes below root-level HOGs
    root = T.names[0]
    cn = parsed.copy_number(root)
    assert cn["n_copies"].sum() == h.loc[h["level"] == root, "n_genes"].sum()


@pytest.mark.slow
def test_ancestra_inferred_hogs_parse(tmp_path):
    if shutil.which("ancestra") is None and not os.path.exists(
            os.path.join(os.path.dirname(sys.executable), "ancestra")):
        pytest.skip("ancestra CLI not available")
    exe = shutil.which("ancestra") or os.path.join(os.path.dirname(sys.executable), "ancestra")
    res = simulate(SimParams(n_species=6, n_genes=200), 13)
    write_simulation(res, str(tmp_path / "sim"))
    env = dict(os.environ, PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""),
               OMP_NUM_THREADS="2", NUMBA_NUM_THREADS="2")
    subprocess.run([exe, "run", "--proteomes", str(tmp_path / "sim" / "proteomes"),
                    "--tree", str(tmp_path / "sim" / "species_tree.nwk"), "--out", str(tmp_path / "out"),
                    "--threads", "2", "--preset", "fast"], check=True, env=env, capture_output=True)
    parsed = parse_hogs(str(tmp_path / "out" / "hogs.orthoxml"),
                        newick=str(tmp_path / "sim" / "species_tree.nwk"))
    assert parsed.check_conservation().empty
    assert len(parsed.genes) == sum(len(v) for v in res.proteomes.values())
    assert parsed.stats["n_collapsed_inconsistent"] == 0
    assert parsed.stats["n_level_hint_inconsistent"] == 0
