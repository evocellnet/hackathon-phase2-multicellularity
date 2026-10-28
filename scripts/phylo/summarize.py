#!/usr/bin/env python
"""Write results/phylo/summary.md (+ TSVs) from the outputs of steps 1-3 (euk and prok).

Reports per set: markers kept, alignment lengths, supermatrix size, total tree
length, the 10 longest terminal branches of the concatenated tree, tips whose
terminal branch is extreme (log-scale robust z > 3 within the set, a sign of
contamination / misplacement), and from the gene trees: per taxon, how often
its terminal branch is the longest of the gene tree, plus genes with any branch
> 2 subst./site.
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from mcell.phylo import read_tree  # noqa: E402

OUT = os.path.join(ROOT, "results", "phylo")


def species_info(set_name):
    f = os.path.join(ROOT, "data", "species_final.tsv" if set_name == "euk" else "species_final_prok.tsv")
    d = pd.read_csv(f, sep="\t", dtype=str).set_index("code")
    return d[["species", "data_type_final", "low_quality", "busco_C", "busco_F"]]


def jload(p):
    return json.load(open(p)) if os.path.exists(p) else None


def section(set_name):
    base = os.path.join(ROOT, "runs", "phylo", set_name)
    lines = [f"## {'Eukaryotes' if set_name == 'euk' else 'Prokaryotes'} ({set_name})", ""]
    bm = jload(os.path.join(base, "build_markers.DONE.json"))
    st = jload(os.path.join(base, "species_tree_bl.DONE.json"))
    gt = jload(os.path.join(base, "gene_trees.DONE.json"))
    if bm is None:
        return lines + ["Step 1 not finished.", ""]
    mk = pd.read_csv(os.path.join(base, "markers", "markers.tsv"), sep="\t")
    info = species_info(set_name)
    lines += [f"**Step 1 (markers).** {bm['n_markers']} markers present (single-copy) in ≥ {bm['min_frac']:.0%} of "
              f"{bm['n_species']} species" + (f"; kinds: {bm['kinds']}" if set_name == "prok" else "") +
              f". MAFFT alignment length median {mk['aln_len'].median():.0f} (range {mk['aln_len'].min()}–"
              f"{mk['aln_len'].max()}); after trimAl -automated1 median {mk['trim_len'].median():.0f} "
              f"(range {mk['trim_len'].min()}–{mk['trim_len'].max()}), total {mk['trim_len'].sum():,} columns; "
              f"species per marker median {mk['n_species'].median():.0f}. Wall time {bm['seconds'] / 60:.0f} min.", ""]
    if set_name == "prok" and os.path.exists(os.path.join(base, "universal_pairs.tsv")):
        up = pd.read_csv(os.path.join(base, "universal_pairs.tsv"), sep="\t")
        lines += [f"Universal (bacteria+archaea) markers: {len(up)} reciprocal majority best-hit marker pairs found by "
                  f"DIAMOND, {int(up['kept'].sum())} pass the ≥ 50 % per-domain presence filter and are merged into one "
                  f"alignment each (`runs/phylo/prok/universal_pairs.tsv`); they alone inform the Archaea–Bacteria path.", ""]
    pres = pd.read_csv(os.path.join(base, "markers", "species_presence.tsv"), sep="\t").set_index("species")
    low = pres[pres["frac_markers"] < 0.5]
    lines += [f"Species with single-copy sequences in < 50 % of the kept markers: {len(low)}" +
              (": " + ", ".join(f"{s} ({r.frac_markers:.0%})" for s, r in low.sort_values('frac_markers').iterrows())
               if len(low) else "") + ".", ""]
    if st is None:
        return lines + ["Step 2 not finished.", ""]
    t = read_tree(os.path.join(ROOT, st["output"]))
    lines += [f"**Step 2 (concatenated tree).** {st['n_markers']} markers, {st['n_sites']:,} sites, {st['n_taxa']} taxa; "
              f"IQ-TREE 3 LG+G4 on the fixed guide topology (`-te`, polytomies kept), logL {st['logL']}. "
              f"Total tree length **{st['total_tree_length']:.2f}** subst./site. Wall time {st['seconds'] / 60:.0f} min "
              f"({st['threads']} threads). Output `{st['output']}` (guide internal node names preserved; topology "
              f"asserted identical to the guide tree).", ""]
    term = pd.Series({lf.name: lf.dist for lf in t.get_leaves()})
    rtt = pd.Series({lf.name: lf.get_distance(t) for lf in t.get_leaves()})
    lt = np.log(term.clip(lower=1e-6))
    med, mad = lt.median(), 1.4826 * np.median(np.abs(lt - lt.median()))
    z = (lt - med) / mad
    occ = pd.read_csv(os.path.join(base, "concat", "taxon_occupancy.tsv"), sep="\t").set_index("taxon")["frac_sites_present"]
    tab = pd.DataFrame({"terminal_bl": term, "root_to_tip": rtt, "robust_z_log_terminal": z, "occupancy": occ}).join(info)
    tab = tab.sort_values("terminal_bl", ascending=False)
    os.makedirs(OUT, exist_ok=True)
    tab.to_csv(os.path.join(OUT, f"terminal_branches_{set_name}.tsv"), sep="\t", float_format="%.5g")
    lines += [f"Terminal branch lengths: median {term.median():.3f}, max {term.max():.3f}. Root-to-tip: median "
              f"{rtt.median():.3f}, range {rtt.min():.3f}–{rtt.max():.3f}.", "",
              "10 longest terminal branches:", "",
              "| tip | species | data type | terminal BL | robust z (log) | root-to-tip | supermatrix occupancy | BUSCO C/F |",
              "|---|---|---|---|---|---|---|---|"]
    for s, r in tab.head(10).iterrows():
        lines.append(f"| {s} | *{r['species']}* | {r['data_type_final']}{' (low quality)' if r['low_quality'] == '1' else ''} | "
                     f"{r['terminal_bl']:.3f} | {r['robust_z_log_terminal']:.1f} | {r['root_to_tip']:.3f} | "
                     f"{r['occupancy']:.0%} | {r['busco_C']}/{r['busco_F']} |")
    ext = tab[tab["robust_z_log_terminal"] > 3]
    lines += ["", f"Tips with extreme terminal branches (robust z of log length > 3): "
              + (", ".join(f"{s} ({r.terminal_bl:.2f})" for s, r in ext.iterrows()) if len(ext) else "none") + ".", ""]
    il = pd.Series({n.name: n.dist for n in t.traverse() if not n.is_leaf() and n.up is not None})
    lines += [f"Internal branches: median {il.median():.4f}; {int((il < 1e-5).sum())} of {len(il)} are < 1e-5 "
              f"(effectively zero; usually short/unsupported guide-tree splits).", ""]
    if gt is None:
        return lines + ["Step 3 not finished.", ""]
    genes = pd.read_csv(os.path.join(base, "genetrees", "genes.tsv"), sep="\t")
    bl = pd.read_csv(os.path.join(base, "genetrees", "branch_lengths.tsv"), sep="\t")
    ok = genes[genes["status"] == "ok"]
    ndrop = ok["dropped"].fillna("").map(lambda s: len([x for x in str(s).split(",") if x]))
    lines += [f"**Step 3 (gene trees).** {gt['n_ok']}/{gt['n_markers']} markers fitted (IQ-TREE `-te` LG+G4, pruned guide "
              f"tree, 1 thread each, {gt['jobs']} in parallel); {gt['n_failed']} failed/skipped"
              + (f" ({', '.join(genes.loc[genes.status != 'ok', 'status'].value_counts().index)})" if gt['n_failed'] else "")
              + f". Taxa per gene median {ok['n_taxa'].median():.0f}; taxa dropped for < 50 % site coverage: "
              f"median {ndrop.median():.0f} per gene, total {int(ndrop.sum())}. Gene tree length median "
              f"{ok['tree_length'].median():.2f} (range {ok['tree_length'].min():.2f}–{ok['tree_length'].max():.2f}). "
              f"Genes with any branch > 2 subst./site: {int((ok['max_branch'] > 2).sum())}. Wall time "
              f"{gt['seconds'] / 60:.0f} min, {gt['cpu_seconds'] / 3600:.1f} CPU-h.", ""]
    # per taxon: fraction of genes in which its terminal branch is the longest terminal; median gene-wise terminal z
    leaf = bl[bl["n_tips"] == 1].copy()
    leaf["rank_in_gene"] = leaf.groupby("gene")["length"].rank(ascending=False)
    leaf["rel"] = leaf["length"] / leaf.groupby("gene")["length"].transform("median")
    per = leaf.groupby("node").agg(n_genes=("gene", "size"), frac_longest=("rank_in_gene", lambda r: np.mean(r == 1)),
                                   median_rel_terminal=("rel", "median"),
                                   frac_rel_gt5=("rel", lambda r: np.mean(r > 5)))
    per = per.join(info).sort_values("median_rel_terminal", ascending=False)
    per.to_csv(os.path.join(OUT, f"genetree_terminal_by_taxon_{set_name}.tsv"), sep="\t", float_format="%.4g")
    lines += ["Gene-tree terminal branches by taxon (relative = terminal length / median terminal length of the gene); "
              "top 10 by median relative terminal length:", "",
              "| tip | species | data type | genes | median relative terminal | frac genes > 5× | frac genes longest |",
              "|---|---|---|---|---|---|---|"]
    for s, r in per.head(10).iterrows():
        lines.append(f"| {s} | *{r['species']}* | {r['data_type_final']} | {r['n_genes']} | {r['median_rel_terminal']:.2f} | "
                     f"{r['frac_rel_gt5']:.0%} | {r['frac_longest']:.0%} |")
    lines.append("")
    return lines


def main():
    lines = ["# Phylogenomic branch lengths — summary", "",
             "Generated by `scripts/phylo/summarize.py` from `runs/phylo/{euk,prok}` (steps: "
             "`scripts/phylo/build_markers.py` → `species_tree_bl.py` → `gene_trees.py`). Topology is always the "
             "guide tree (DESIGN.md §3.1); only branch lengths are estimated. Tables: "
             "`results/phylo/terminal_branches_<set>.tsv`, `results/phylo/genetree_terminal_by_taxon_<set>.tsv`.", ""]
    for s in ("euk", "prok"):
        lines += section(s)
    notes = os.path.join(OUT, "notes.md")
    if os.path.exists(notes):
        lines += open(notes).read().splitlines()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
