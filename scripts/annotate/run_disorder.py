#!/usr/bin/env python
"""Per-protein length and predicted intrinsic disorder (metapredict v3) for every proteome.

One output per species: <out>/<code>.disorder.tsv.gz with columns protein, length, frac_disordered
(fraction of residues with metapredict score >= 0.5), mean_disorder. Resumable via DONE markers.

Usage: run_disorder.py --proteomes DIR --out DIR [--threads 8] [--species code1,code2]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time

import numpy as np


def read_fasta(path: str):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--species", default=None)
    ap.add_argument("--batch", type=int, default=2000)
    args = ap.parse_args()

    import torch
    torch.set_num_threads(args.threads)
    import metapredict as meta

    os.makedirs(args.out, exist_ok=True)
    files = sorted(f for f in os.listdir(args.proteomes) if f.endswith(".fa"))
    if args.species:
        keep = set(args.species.split(","))
        files = [f for f in files if f[:-3] in keep]
    for f in files:
        code = f[:-3]
        done = os.path.join(args.out, f"{code}.DONE.json")
        if os.path.exists(done):
            continue
        t = time.time()
        recs = list(read_fasta(os.path.join(args.proteomes, f)))
        tmp = os.path.join(args.out, f"{code}.disorder.tsv.gz.tmp")
        with gzip.open(tmp, "wt") as out:
            out.write("protein\tlength\tfrac_disordered\tmean_disorder\n")
            for i in range(0, len(recs), args.batch):
                chunk = recs[i:i + args.batch]
                # metapredict rejects non-canonical letters; X/U/B/Z/O are mapped to A only for prediction
                clean = {n: s.translate(str.maketrans("XUBZOJ", "AAAAAA")) for n, s in chunk}
                res = meta.predict_disorder(clean)
                for n, s in chunk:
                    sc = np.asarray(res[n][1] if isinstance(res[n], (tuple, list)) else res[n], dtype=float)
                    out.write(f"{n}\t{len(s)}\t{(sc >= 0.5).mean():.4f}\t{sc.mean():.4f}\n")
        os.replace(tmp, os.path.join(args.out, f"{code}.disorder.tsv.gz"))
        rec = {"species": code, "n": len(recs), "seconds": round(time.time() - t, 1), "metapredict": meta.__version__,
               "threshold": 0.5}
        with open(done, "w") as fh:
            json.dump(rec, fh)
        print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
