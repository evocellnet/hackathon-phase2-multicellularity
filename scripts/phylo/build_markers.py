#!/usr/bin/env python
"""Step 1 of the phylogenomic pipeline: BUSCO single-copy markers -> trimmed alignments.

DESIGN.md §3.1 / §4.5; reviews/round2/methods_compute.md §5-6.

For every BUSCO marker we collect the *single-copy* protein of each species
(``busco_sequences/single_copy_busco_sequences/<marker>.faa``; duplicated and
fragmented BUSCOs are not used), rename the sequence to the species code, and
keep markers present in >= ``--min-frac`` (default 0.5) of the species of the
set. Each kept marker is aligned with ``mafft --auto`` and trimmed with
``trimal -automated1``.

Sets
----
euk   eukaryota_odb10, all tips of data/guide_tree_euk.nwk.
prok  bacteria_odb10 (bacteria) and archaea_odb10 (archaea) have disjoint
      marker IDs. We (i) keep domain-specific markers present in >= min-frac of
      the species *of that domain* (they inform within-domain branch lengths
      only), and (ii) build *universal* markers by pairing a bacterial with an
      archaeal marker when they are reciprocal majority best-hit markers in an
      all-vs-all DIAMOND search of the raw single-copy sequences (> 50 % of the
      sequences of each marker have their best cross-domain hit in the other).
      The two sequence sets of a pair are merged into one marker (ID
      ``U_<bac>__<arc>``) before alignment; the merged marker must also pass
      min-frac in each domain. Universal markers are the only ones that inform
      the Archaea-Bacteria root path. A domain marker used in a pair is not
      also kept on its own.

Outputs (runs/phylo/<set>/markers/): ``<m>.raw.fa`` (unaligned), ``<m>.aln.fa``
(MAFFT), ``<m>.trim.fa`` (trimAl), ``markers.tsv`` (per-marker stats),
``species_presence.tsv`` and ``../build_markers.DONE.json``.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BUSCO = os.path.join(ROOT, "data", "busco_runs")
ENV_BIN = os.path.join(ROOT, "conda_envs", "ancestra", "bin")
os.environ["PATH"] = ENV_BIN + os.pathsep + os.environ.get("PATH", "")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), 80):
                fh.write(v[i:i + 80] + "\n")


def tree_tips(path):
    import ete3
    return ete3.Tree(open(path).read().strip(), format=1).get_leaf_names()


def species_lineages(set_name):
    """species code -> BUSCO lineage dataset name."""
    tsv = os.path.join(ROOT, "data", "species_final.tsv" if set_name == "euk" else "species_final_prok.tsv")
    out = {}
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["code"]] = r["busco_lineage"]
    return out


def collect(codes, lineage_of):
    """lineage -> marker -> {species: seq} (single-copy only, '*' and stop removed)."""
    out = defaultdict(lambda: defaultdict(dict))
    for c in codes:
        lin = lineage_of[c]
        d = os.path.join(BUSCO, c, f"run_{lin}", "busco_sequences", "single_copy_busco_sequences")
        if not os.path.isdir(d):
            print(f"WARNING: no single-copy dir for {c}", file=sys.stderr)
            continue
        for fn in os.listdir(d):
            if not fn.endswith(".faa"):
                continue
            seqs = read_fasta(os.path.join(d, fn))
            if len(seqs) != 1:
                print(f"WARNING: {c}/{fn} has {len(seqs)} sequences; skipped", file=sys.stderr)
                continue
            s = next(iter(seqs.values())).replace("*", "").upper()
            out[lin][fn[:-4]][c] = s
    return out


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed ({r.returncode}): {r.stderr[-2000:]}")
    return r


def align_one(marker, outdir):
    raw = os.path.join(outdir, f"{marker}.raw.fa")
    aln = os.path.join(outdir, f"{marker}.aln.fa")
    trim = os.path.join(outdir, f"{marker}.trim.fa")
    t0 = time.time()
    if not (os.path.exists(trim) and os.path.getsize(trim) > 0):
        r = run(["mafft", "--auto", "--thread", "1", "--quiet", raw])
        tmp = aln + ".tmp"
        with open(tmp, "w") as fh:
            fh.write(r.stdout)
        os.replace(tmp, aln)
        run(["trimal", "-in", aln, "-out", trim + ".tmp", "-automated1"])
        os.replace(trim + ".tmp", trim)
    a = read_fasta(aln)
    tr = read_fasta(trim) if os.path.getsize(trim) > 0 else {}
    tlen = len(next(iter(tr.values()))) if tr else 0
    cov = {k: (sum(ch not in "-X?" for ch in v) / tlen if tlen else 0.0) for k, v in tr.items()}
    return dict(marker=marker, n_species=len(a), aln_len=len(next(iter(a.values()))),
                trim_len=tlen, n_species_trim=len(tr),
                n_cov_lt50=sum(c < 0.5 for c in cov.values()),
                median_cov=round(sorted(cov.values())[len(cov) // 2], 3) if cov else 0.0,
                seconds=round(time.time() - t0, 1))


def diamond_pairs(bac, arc, workdir, threads):
    """Reciprocal majority best-hit pairing of bacterial and archaeal markers."""
    os.makedirs(workdir, exist_ok=True)
    qa, qb = os.path.join(workdir, "arc.fa"), os.path.join(workdir, "bac.fa")
    write_fasta(qa, {f"{m}|{sp}": s for m, d in arc.items() for sp, s in d.items()})
    write_fasta(qb, {f"{m}|{sp}": s for m, d in bac.items() for sp, s in d.items()})
    best = {}
    for q, db in ((qb, qa), (qa, qb)):
        dbp = db[:-3]
        run(["diamond", "makedb", "--in", db, "-d", dbp, "--quiet"])
        out = q[:-3] + ".vs." + os.path.basename(dbp) + ".tsv"
        run(["diamond", "blastp", "-q", q, "-d", dbp, "-o", out, "--more-sensitive", "-e", "1e-5",
             "-k", "5", "--threads", str(threads), "--quiet", "-f", "6", "qseqid", "sseqid", "bitscore"])
        bs = {}
        with open(out) as fh:
            for line in fh:
                qid, sid, b = line.split("\t")
                if qid not in bs or float(b) > bs[qid][1]:
                    bs[qid] = (sid.split("|")[0], float(b))
        best.update({k: v[0] for k, v in bs.items()})
    votes_b = defaultdict(Counter)
    votes_a = defaultdict(Counter)
    for qid, tm in best.items():
        m = qid.split("|")[0]
        (votes_b if m in bac else votes_a)[m][tm] += 1
    pairs = []
    for b, cnt in votes_b.items():
        a, n = cnt.most_common(1)[0]
        if n > 0.5 * len(bac[b]) and a in votes_a:
            b2, n2 = votes_a[a].most_common(1)[0]
            if b2 == b and n2 > 0.5 * len(arc[a]):
                pairs.append((b, a, n / len(bac[b]), n2 / len(arc[a])))
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=["euk", "prok"], required=True)
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--jobs", type=int, default=4, help="parallel single-thread MAFFT/trimAl jobs")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t0 = time.time()
    outroot = a.out or os.path.join(ROOT, "runs", "phylo", a.set)
    mdir = os.path.join(outroot, "markers")
    os.makedirs(mdir, exist_ok=True)

    codes = tree_tips(os.path.join(ROOT, "data", f"guide_tree_{a.set}.nwk"))
    lin = species_lineages(a.set)
    data = collect(codes, lin)
    final = {}  # marker -> {sp: seq}
    meta = {}
    if a.set == "euk":
        n = len(codes)
        for m, d in data["eukaryota_odb10"].items():
            if len(d) >= a.min_frac * n:
                final[m] = d
                meta[m] = dict(kind="euk", n_dom1=len(d), n_dom2="")
    else:
        bac_sp = [c for c in codes if lin[c] == "bacteria_odb10"]
        arc_sp = [c for c in codes if lin[c] == "archaea_odb10"]
        bac, arc = data["bacteria_odb10"], data["archaea_odb10"]
        pairs = diamond_pairs(bac, arc, os.path.join(outroot, "diamond_pairing"), a.jobs)
        with open(os.path.join(outroot, "universal_pairs.tsv"), "w") as fh:
            fh.write("bacteria_marker\tarchaea_marker\tfrac_bac_besthit\tfrac_arc_besthit\tkept\n")
            used = set()
            for b, ar, fb, fa in sorted(pairs):
                ok = len(bac[b]) >= a.min_frac * len(bac_sp) and len(arc[ar]) >= a.min_frac * len(arc_sp)
                fh.write(f"{b}\t{ar}\t{fb:.3f}\t{fa:.3f}\t{int(ok)}\n")
                if ok:
                    m = f"U_{b}__{ar}"
                    final[m] = {**bac[b], **arc[ar]}
                    meta[m] = dict(kind="universal", n_dom1=len(bac[b]), n_dom2=len(arc[ar]))
                    used.update([b, ar])
        for dom, dd, nsp in (("bacteria", bac, len(bac_sp)), ("archaea", arc, len(arc_sp))):
            for m, d in dd.items():
                if m not in used and len(d) >= a.min_frac * nsp:
                    final[m] = d
                    meta[m] = dict(kind=dom, n_dom1=len(d) if dom == "bacteria" else "",
                                   n_dom2=len(d) if dom == "archaea" else "")
    for m, d in final.items():
        raw = os.path.join(mdir, f"{m}.raw.fa")
        write_fasta(raw + ".tmp", dict(sorted(d.items())))
        os.replace(raw + ".tmp", raw)
    print(f"[{a.set}] {len(codes)} species, {len(final)} markers kept; aligning with {a.jobs} jobs", flush=True)

    rows = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(align_one, m, mdir): m for m in sorted(final)}
        for i, f in enumerate(cf.as_completed(futs), 1):
            r = f.result()
            r.update(meta[r["marker"]])
            rows.append(r)
            if i % 25 == 0:
                print(f"  {i}/{len(futs)} aligned ({time.time() - t0:.0f} s)", flush=True)
    rows.sort(key=lambda r: r["marker"])
    cols = ["marker", "kind", "n_species", "n_dom1", "n_dom2", "aln_len", "trim_len", "n_species_trim",
            "n_cov_lt50", "median_cov", "seconds"]
    with open(os.path.join(mdir, "markers.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    pres = Counter(sp for d in final.values() for sp in d)
    with open(os.path.join(mdir, "species_presence.tsv"), "w") as fh:
        fh.write("species\tbusco_lineage\tn_markers\tfrac_markers\n")
        for c in codes:
            fh.write(f"{c}\t{lin[c]}\t{pres[c]}\t{pres[c] / max(len(final), 1):.3f}\n")
    empty = [r["marker"] for r in rows if r["trim_len"] == 0]
    done = dict(step="build_markers", set=a.set, n_species=len(codes), n_markers=len(final),
                n_empty_after_trim=len(empty), min_frac=a.min_frac,
                kinds=dict(Counter(meta[m]["kind"] for m in final)),
                total_trim_len=sum(r["trim_len"] for r in rows),
                seconds=round(time.time() - t0, 1), finished=dt.datetime.now().isoformat(timespec="seconds"),
                versions=dict(mafft=subprocess.run(["mafft", "--version"], capture_output=True, text=True).stderr.strip(),
                              trimal=subprocess.run(["trimal", "--version"], capture_output=True, text=True).stdout.strip()))
    with open(os.path.join(outroot, "build_markers.DONE.json"), "w") as fh:
        json.dump(done, fh, indent=1)
    print(json.dumps(done, indent=1))


if __name__ == "__main__":
    main()
