#!/usr/bin/env python
"""Second-pass protein mapping for species whose exact-sequence (MD5) mapping is < 90 %:
reciprocal best DIAMOND hits between our not-yet-mapped proteins and the not-yet-mapped
EukaryOMA proteins of the same species (sequences from fetch_unmapped_seqs.py).
A pair is accepted if it is reciprocal best (by bitscore), identity >= 95 % and the alignment
covers >= 50 % of both proteins (same gene from a different assembly / gene-model version;
not orthology).
Output: runs/eukaryoma_xval/protein_map_rbh.tsv.gz (code, our_id, oma_id, pident, qcov, scov),
        results/eukaryoma_xval/protein_mapping_stats_rbh.tsv
"""
import os, subprocess, sys, tempfile
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from map_proteins import read_fasta  # noqa: E402

DIAMOND = os.path.join(ROOT, "conda_envs/ancestra/bin/diamond")
THREADS = 12
MIN_ID, MIN_COV = 95.0, 0.5


def best(path):
    cols = ["q", "s", "pident", "length", "qlen", "slen", "bitscore"]
    d = pd.read_csv(path, sep="\t", names=cols)
    if d.empty:
        return d
    d = d.sort_values(["q", "bitscore"], ascending=[True, False]).drop_duplicates("q")
    d["qcov"] = d["length"] / d["qlen"]
    d["scov"] = d["length"] / d["slen"]
    return d


def main():
    st = pd.read_csv(os.path.join(ROOT, "results/eukaryoma_xval/protein_mapping_stats.tsv"), sep="\t")
    pm = pd.read_csv(os.path.join(ROOT, "runs/eukaryoma_xval/protein_map.tsv.gz"), sep="\t")
    mapped_ours = set(pm["our_id"])
    out, stats = [], []
    for _, r in st.iterrows():
        code = r["code"]
        sq = os.path.join(ROOT, f"runs/eukaryoma_xval/seqs/{code}.fa")
        if not os.path.exists(sq):
            continue
        with tempfile.TemporaryDirectory() as td:
            a = os.path.join(td, "ours.fa")
            n_un = 0
            with open(a, "w") as fh:
                for pid, s in read_fasta(os.path.join(ROOT, f"data/proteomes/euk/{code}.fa")):
                    if pid not in mapped_ours:
                        fh.write(f">{pid}\n{s.rstrip('*')}\n")
                        n_un += 1
            res = {}
            for tag, q, db in (("f", a, sq), ("r", sq, a)):
                subprocess.run([DIAMOND, "makedb", "--in", db, "-d", os.path.join(td, f"db{tag}"), "--quiet"], check=True)
                o = os.path.join(td, f"{tag}.tsv")
                subprocess.run([DIAMOND, "blastp", "-q", q, "-d", os.path.join(td, f"db{tag}"), "-o", o, "--quiet",
                                "--threads", str(THREADS), "--id", str(MIN_ID - 5), "--max-target-seqs", "5",
                                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore"],
                               check=True)
                res[tag] = best(o)
        f, rv = res["f"], res["r"]
        if f.empty or rv.empty:
            stats.append(dict(code=code, n_ours_unmapped=n_un, n_rbh=0))
            continue
        back = dict(zip(rv["q"], rv["s"]))
        f = f[[back.get(s) == q for q, s in zip(f["q"], f["s"])]]
        f = f[(f["pident"] >= MIN_ID) & (f["qcov"] >= MIN_COV) & (f["scov"] >= MIN_COV)]
        for x in f.itertuples():
            out.append((code, x.q, x.s, x.pident, round(x.qcov, 3), round(x.scov, 3)))
        stats.append(dict(code=code, n_ours_unmapped=n_un, n_rbh=len(f)))
        print(code, n_un, len(f), flush=True)
    o = pd.DataFrame(out, columns=["code", "our_id", "oma_id", "pident", "qcov", "scov"])
    o.to_csv(os.path.join(ROOT, "runs/eukaryoma_xval/protein_map_rbh.tsv.gz"), sep="\t", index=False)
    s = pd.DataFrame(stats).merge(st[["code", "n_ours", "n_oma", "n_mapped"]], on="code")
    s["frac_ours_mapped_total"] = ((s["n_mapped"] + s["n_rbh"]) / s["n_ours"]).round(4)
    s.to_csv(os.path.join(ROOT, "results/eukaryoma_xval/protein_mapping_stats_rbh.tsv"), sep="\t", index=False)
    print(s.to_string())


if __name__ == "__main__":
    main()
