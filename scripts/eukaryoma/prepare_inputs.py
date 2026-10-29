#!/usr/bin/env python
"""Project both HOG sets onto the shared species set and our guide-tree topology, then parse them
with the same parser (mcell.hogs) and the same level rule.

Steps (one ``--step`` per process so they can run in parallel):
  tree      prune runs/euk_fast/stages/input/guide_tree.nwk to the species mapped to EukaryOMA
            (results/eukaryoma_xval/species_map.tsv) -> runs/eukaryoma_xval/tree_pruned.nwk
  ancestra  filter runs/euk_fast/hogs.orthoxml to those species, strip level annotations, parse
            -> runs/eukaryoma_xval/parsed/ancestra/
  ancestra_hints  same but keep Ancestra's TaxRange/taxonId (sensitivity; names absent from the
            pruned tree fall back to the MRCA rule) -> parsed/ancestra_hints/
  eukaryoma filter data/raw/eukaryoma/oma-hogs.orthoXML.gz (species renamed to our codes, gene ids
            kept as EukaryOMA ids), strip levels, parse -> parsed/eukaryoma/
Level rule after stripping (both methods): orthologGroup at the MRCA of its members on the pruned
tree, paralogGroup on the branch above that MRCA; an ancestral gene is present at node X if its
gene-tree edge crosses X (mcell.hogs definitions).
"""
import argparse, json, os, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd  # noqa: E402
import xval_lib as xl  # noqa: E402

W = os.path.join(ROOT, "runs/eukaryoma_xval")
TREE = os.path.join(W, "tree_pruned.nwk")


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def species_map():
    return pd.read_csv(os.path.join(ROOT, "results/eukaryoma_xval/species_map.tsv"), sep="\t")


def step_tree():
    sm = species_map()
    orig = xl.load_tree(os.path.join(ROOT, "runs/euk_fast/stages/input/guide_tree.nwk"))
    t = xl.prune_tree(orig, set(sm["code"]))
    xl.write_newick(t, TREE)
    log(f"pruned tree: {len(t.get_leaves())} leaves, {sum(1 for n in t.traverse() if not n.is_leaf())} internal nodes "
        f"(original {len(orig.get_leaves())} / {sum(1 for n in orig.traverse() if not n.is_leaf())})")


def parse_and_write(xml, name):
    from mcell.hogs import parse_hogs
    t0 = time.time()
    res = parse_hogs(xml, newick=TREE, progress=lambda n: log(f"{name}: {n:,} families") if n % 200000 == 0 else None)
    bad = res.check_conservation()
    out = os.path.join(W, "parsed", name)
    res.write(out, fmt="parquet")
    json.dump(dict(stats=res.stats, conservation_violations=len(bad), wall_s=time.time() - t0, xml=xml),
              open(os.path.join(out, "run_info.json"), "w"), indent=1, default=str)
    log(f"{name}: parsed {res.stats['n_genes']:,} genes, {res.stats['n_families']:,} families, "
        f"{len(res.hogs):,} ancestral genes; conservation violations {len(bad)}; stats {res.stats}")


def step_ancestra(strip=True):
    sm = species_map()
    name = "ancestra" if strip else "ancestra_hints"
    xml = os.path.join(W, f"{name}.orthoxml.gz")
    if not os.path.exists(xml):
        st = xl.filter_orthoxml(os.path.join(ROOT, "runs/euk_fast/hogs.orthoxml"), xml + ".tmp.gz",
                                {c: c for c in sm["code"]}, strip_levels=strip, log=log)
        os.replace(xml + ".tmp.gz", xml)
        log(f"{name} filter: {st}")
    parse_and_write(xml, name)


def step_eukaryoma():
    sm = species_map()
    xml = os.path.join(W, "eukaryoma.orthoxml.gz")
    if not os.path.exists(xml):
        st = xl.filter_orthoxml(os.path.join(ROOT, "data/raw/eukaryoma/oma-hogs.orthoXML.gz"), xml + ".tmp.gz",
                                dict(zip(sm["eukaryoma_name"], sm["code"])), strip_levels=True, log=log)
        os.replace(xml + ".tmp.gz", xml)
        log(f"eukaryoma filter: {st}")
    parse_and_write(xml, "eukaryoma")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True, choices=["tree", "ancestra", "ancestra_hints", "eukaryoma"])
    a = ap.parse_args()
    {"tree": step_tree, "ancestra": step_ancestra, "ancestra_hints": lambda: step_ancestra(False),
     "eukaryoma": step_eukaryoma}[a.step]()
