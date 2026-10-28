#!/usr/bin/env python
"""Pfam-A domain annotation of every proteome with pyhmmer (hmmsearch, gathering cutoffs).

One output per species: <out>/<code>.pfam.tsv.gz with one row per domain hit that passes the
per-domain gathering threshold (columns: protein, pfam_acc, pfam_name, ali_from, ali_to, env_from,
env_to, dom_score, dom_ievalue, seq_score). A species is skipped when its DONE marker exists, so the
script can be re-launched after an interruption.

Usage: run_pfam.py --proteomes DIR --hmm Pfam-A.hmm --out DIR --cpus 12 [--species code1,code2]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import time

import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def annotate(fasta: str, hmms: list, out_tsv: str, cpus: int) -> int:
    abc = Alphabet.amino()
    with SequenceFile(fasta, digital=True, alphabet=abc) as sf:
        seqs = sf.read_block()
    n = 0
    tmp = out_tsv + ".tmp"
    with gzip.open(tmp, "wt") as out:
        out.write("protein\tpfam_acc\tpfam_name\tali_from\tali_to\tenv_from\tenv_to\tdom_score\tdom_ievalue\tseq_score\n")
        for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=cpus, bit_cutoffs="gathering"):
            acc = hits.query.accession
            acc = (acc.decode() if isinstance(acc, bytes) else acc).split(".")[0]
            name = hits.query.name
            name = name.decode() if isinstance(name, bytes) else name
            for hit in hits.included:
                pid = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                for d in hit.domains.included:
                    a = d.alignment
                    out.write(f"{pid}\t{acc}\t{name}\t{a.target_from}\t{a.target_to}\t{d.env_from}\t{d.env_to}"
                              f"\t{d.score:.1f}\t{d.i_evalue:.2e}\t{hit.score:.1f}\n")
                    n += 1
    os.replace(tmp, out_tsv)
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cpus", type=int, default=12)
    ap.add_argument("--species", default=None, help="comma-separated species codes (default: all)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    files = sorted(f for f in os.listdir(args.proteomes) if f.endswith(".fa"))
    if args.species:
        keep = set(args.species.split(","))
        files = [f for f in files if f[:-3] in keep]
    t0 = time.time()
    with HMMFile(args.hmm) as hf:
        hmms = list(hf)
    print(f"loaded {len(hmms)} HMMs in {time.time() - t0:.0f}s", flush=True)
    for f in files:
        code = f[:-3]
        done = os.path.join(args.out, f"{code}.DONE.json")
        if os.path.exists(done):
            continue
        fasta = os.path.join(args.proteomes, f)
        t = time.time()
        n = annotate(fasta, hmms, os.path.join(args.out, f"{code}.pfam.tsv.gz"), args.cpus)
        rec = {"species": code, "fasta_md5": md5(fasta), "n_domain_hits": n, "seconds": round(time.time() - t, 1),
               "pyhmmer": pyhmmer.__version__, "hmm": os.path.abspath(args.hmm), "cutoff": "gathering"}
        with open(done, "w") as fh:
            json.dump(rec, fh)
        print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
