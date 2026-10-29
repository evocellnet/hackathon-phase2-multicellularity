#!/usr/bin/env python
"""W3 per panel: Pfam domains, architectures and functional categories of the sequence-novel gained families and
the duplicated HOGs on every laddered branch (events and controls; the controls feed the recurrence null).

    w3_panel.py R [--cpus 12]

Representatives per family: up to 3 genes; for sequence-novel gains the ladder representatives; for duplications
(non-origin 'dup' events: the parent-level HOG duplicated on the branch) one gene from each of up to 3 daughter
copies. Domains: round-1 full Pfam tables (Pfam-A 38.2, GA) where they exist; new panel species are annotated here
with pyhmmer hmmsearch (Pfam-A 38.2, gathering cutoffs), only for the representative proteins. Family domain set =
union over representatives, transposon mask (data/r2_sets/transposon_mask.tsv) removed; architecture = ordered
domains of the longest annotated representative (overlaps resolved as mcell.gainladder._resolve_overlaps).
Outputs: results/r2/panels/<R>/w3_families.tsv.gz (all laddered branches), w3_event_families.tsv (simple->complex
events), gene_domains.tsv.gz.
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, paths, guard, log, kingdom, panel_codes  # noqa: E402
from panel import load_res, transitions  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import gainladder as gl  # noqa: E402

PFAM_HMM = os.path.join(ROOT, "data", "raw", "pfam", "Pfam-A.hmm")
ANNOT_PY = os.path.join(ROOT, "conda_envs", "annot", "bin", "python")


def category_tables():
    tk = pd.read_csv(os.path.join(ROOT, "data", "toolkit_pfam.tsv"), sep="\t")
    cat = {}
    for r in tk.itertuples():
        cat.setdefault(r.pfam_acc, set()).add(r.category)
    for f, name in (("elm_motif_binding.tsv", "elm_motif_binding"), ("cell_death.tsv", "cell_death"),
                    ("tf_celltype.tsv", "tf_celltype")):
        t = pd.read_csv(os.path.join(ROOT, "data", "r2_sets", f), sep="\t")
        for acc in t["pfam"]:
            cat.setdefault(acc, set()).add(name)
    go = {}
    for line in open(os.path.join(ROOT, "data", "raw", "go", "pfam2go")):
        if line.startswith("!"):
            continue
        acc = line.split()[0].split(":")[1]
        term = line.split(">", 1)[1].split(";")[0].strip().replace("GO:", "", 1)
        go.setdefault(acc, []).append(term)
    mask = set(pd.read_csv(os.path.join(ROOT, "data", "r2_sets", "transposon_mask.tsv"), sep="\t")["pfam"])
    return cat, go, mask


def annotate_new(fasta, out_tsv, cpus):
    code = f"""
import gzip, pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile
abc = Alphabet.amino()
with SequenceFile({fasta!r}, digital=True, alphabet=abc) as sf:
    seqs = sf.read_block()
with HMMFile({PFAM_HMM!r}) as hf:
    hmms = list(hf)
with gzip.open({out_tsv!r} + '.tmp', 'wt') as out:
    out.write('protein\\tpfam_acc\\tpfam_name\\tali_from\\tali_to\\tenv_from\\tenv_to\\tdom_score\\n')
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus={cpus}, bit_cutoffs='gathering'):
        acc = hits.query.accession; acc = (acc.decode() if isinstance(acc, bytes) else acc).split('.')[0]
        name = hits.query.name; name = name.decode() if isinstance(name, bytes) else name
        for hit in hits.included:
            pid = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            for d in hit.domains.included:
                a = d.alignment
                out.write(f'{{pid}}\\t{{acc}}\\t{{name}}\\t{{a.target_from}}\\t{{a.target_to}}\\t{{d.env_from}}\\t{{d.env_to}}\\t{{d.score:.1f}}\\n')
import os; os.replace({out_tsv!r} + '.tmp', {out_tsv!r})
"""
    subprocess.run(["nice", "-n", "10", ANNOT_PY, "-c", code], check=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("R")
    ap.add_argument("--cpus", type=int, default=12)
    a = ap.parse_args(argv)
    R = a.R
    p = paths(R)
    dst = os.path.join(p["out"], "w3_families.tsv.gz")
    guard(dst)
    os.makedirs(p["pfam"], exist_ok=True)
    res = load_res(p)
    T = res.tree
    bt = pd.read_csv(os.path.join(p["out"], "branch_table.tsv"), sep="\t").set_index("branch_child")
    # internal laddered branches + event tips (non-event tips are left out: their orphan-rich gains would need
    # Pfam for >100k new-species proteins; tip events are therefore described but not in the calibrated null)
    lad = set(bt.index[bt["laddered"] & (~bt["is_leaf"] | (bt["events"].fillna("") != ""))])
    cls = pd.read_parquet(os.path.join(p["ladder"], "classified.parquet"))
    reps = pd.read_parquet(os.path.join(p["ladder"], "reps.parquet"))
    g2sp = dict(zip(res.genes["gene"], res.genes["species"].astype(str)))
    # ---- novel families
    nov = cls[cls["seq_novel"] & cls["branch_child"].isin(lad)][["family_id", "branch_child", "n_genes", "n_species"]]
    rep_of = reps.groupby("family_id")["gene"].apply(list)
    fam_rows = [dict(kind="seq_novel", branch_child=r.branch_child, unit_id=r.family_id, family_id=r.family_id,
                     n_genes=r.n_genes, n_dup_events=0, reps=rep_of.get(r.family_id, [])) for r in nov.itertuples()]
    # ---- duplications (non-origin)
    ev = res.events
    d = ev[(ev["event"] == "dup") & ~ev["at_origin"].astype(bool) & ev["branch_child"].isin(lad)]
    hrep = dict(zip(res.hogs["hog_id"], res.hogs["rep_genes"]))
    for r in d.itertuples():
        kids = str(r.child_hog_ids).split(";")
        genes = []
        for k in kids:
            if k in g2sp:
                genes.append(k)
            elif k in hrep and isinstance(hrep[k], str):
                genes.append(hrep[k].split(";")[0])
            if len(genes) >= 3:
                break
        fam_rows.append(dict(kind="dup", branch_child=r.branch_child, unit_id=r.hog_id, family_id=r.family_id,
                             n_genes=r.n_extant_genes, n_dup_events=r.n_events, reps=genes))
    fams = pd.DataFrame(fam_rows)
    genes = sorted({g for gs in fams["reps"] for g in gs})
    log(f"{len(fams):,} units ({(fams['kind'] == 'seq_novel').sum():,} novel, {(fams['kind'] == 'dup').sum():,} dup); "
        f"{len(genes):,} representative proteins")
    # ---- domains: round-1 tables + new annotation
    r1dir = os.path.join(ROOT, "data", "annotations", "pfam", kingdom(R))
    by_sp = pd.Series(genes).groupby(pd.Series([g2sp[g] for g in genes])).apply(set)
    hits, new_genes = [], []
    for sp, gs in by_sp.items():
        f = os.path.join(r1dir, f"{sp}.pfam.tsv.gz")
        if os.path.exists(f):
            t = pd.read_csv(f, sep="\t", usecols=["protein", "pfam_acc", "pfam_name", "ali_from", "ali_to", "env_from",
                                                  "env_to", "dom_score"])
            hits.append(t[t["protein"].isin(gs)])
        else:
            new_genes += sorted(gs)
    newtsv = os.path.join(p["pfam"], "new_species_reps.pfam.tsv.gz")
    if new_genes and not os.path.exists(newtsv):
        fa = os.path.join(p["pfam"], "new_species_reps.faa")
        want = set(new_genes)
        with open(fa, "w") as fh:
            for f in sorted(os.listdir(p["proteomes"])):
                for h, s in gl.read_fasta(os.path.join(p["proteomes"], f), keep=want):
                    fh.write(f">{h}\n{s.replace('*', '')}\n")
        t0 = time.time()
        annotate_new(fa, newtsv, a.cpus)
        log(f"pyhmmer Pfam-A on {len(new_genes):,} new-species proteins in {time.time() - t0:.0f}s")
    if new_genes:
        hits.append(pd.read_csv(newtsv, sep="\t"))
    H = pd.concat(hits, ignore_index=True) if hits else pd.DataFrame(columns=["protein", "pfam_acc"])
    H["pfam_acc"] = H["pfam_acc"].astype(str).str.split(".").str[0]
    H = gl._resolve_overlaps(H.rename(columns={})) if len(H) else H
    H.to_csv(os.path.join(p["out"], "gene_domains.tsv.gz"), sep="\t", index=False)
    cat, go, mask = category_tables()
    name_of = dict(zip(H["pfam_acc"], H["pfam_name"]))
    dom_of = H[~H["pfam_acc"].isin(mask)].sort_values(["protein", "ali_from"]).groupby("protein")["pfam_acc"].apply(list)
    masked_of = set(H.loc[H["pfam_acc"].isin(mask), "protein"])
    out = []
    for r in fams.itertuples():
        reps_ = list(r.reps)
        ann = [g for g in reps_ if g in dom_of.index]
        pset = sorted({x for g in ann for x in dom_of[g]})
        arch = ""
        if ann:
            g = max(ann, key=lambda x: len(dom_of[x]))
            arch = "+".join(dom_of[g])
        cats = sorted({c for x in pset for c in cat.get(x, ())})
        gos = sorted({t for x in pset for t in go.get(x, [])})
        out.append(dict(panel=R, kind=r.kind, branch_child=r.branch_child, unit_id=r.unit_id, family_id=r.family_id,
                        n_genes=r.n_genes, n_dup_events=r.n_dup_events, reps=";".join(reps_),
                        rep_species=";".join(sorted({g2sp[g] for g in reps_})),
                        transposon_domain=any(g in masked_of for g in reps_),
                        pfam=";".join(pset), pfam_names=";".join(name_of.get(x, x) for x in pset),
                        architecture=arch, categories=";".join(cats), go="; ".join(gos[:6])))
    W = pd.DataFrame(out)
    W.to_csv(dst, sep="\t", index=False)
    # ---- simple->complex events
    tr = transitions(R)
    sc = tr[(tr["from_state"] == "simple") & (tr["to_state"] == "complex")]
    E = W[W["branch_child"].isin(set(sc["M"]))].merge(sc[["M", "event_id", "kind"]].rename(
        columns={"M": "branch_child", "kind": "event_kind"}), on="branch_child")
    E.to_csv(os.path.join(p["out"], "w3_event_families.tsv"), sep="\t", index=False)
    log(f"wrote {dst} ({len(W):,} units; {int((W['pfam'] != '').sum()):,} with Pfam) and {len(E):,} event units")


if __name__ == "__main__":
    main()
