#!/usr/bin/env python
"""Concatenate per-species ProstT5 3Di FASTAs for the proteins of an Ancestra run, in the run's gene order,
with integrity checks (one record per protein, 3Di length == protein length, alphabet).

Usage: make_3di_fasta.py <ancestra_run_dir> <out.fa>
Gene table: <run>/stages/input/genes.tsv (Ancestra's own protein ids). 3Di: data/structure3di/all/<sp>.3di.fa
(only species with a .done marker). Proteins without 3Di (e.g. listed in .skipped) are omitted and counted;
Ancestra's rescue skips genes without 3Di. Refuses to overwrite <out.fa>. Writes <out.fa>.report.json.
"""
import json, os, sys, time
import pandas as pd

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
D3 = os.path.join(ROOT, "data/structure3di/all")
ALPHA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fa(p):
    name, buf = None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main(run, out):
    if os.path.exists(out):
        sys.exit(f"refusing to overwrite {out}")
    g = pd.read_csv(os.path.join(run, "stages/input/genes.tsv"), sep="\t")
    lens = dict(zip(g["protein"], g["length"]))
    rep = dict(n_genes=len(g), per_species={}, started=time.strftime("%Y-%m-%dT%H:%M:%S"))
    tmp = out + ".part"
    n_bad = 0
    with open(tmp, "w") as fo:
        for sp, sub in g.groupby("species", sort=True):
            f = os.path.join(D3, f"{sp}.3di.fa")
            if not os.path.exists(f + ".done"):
                sys.exit(f"3Di not complete for {sp}")
            want = set(sub["protein"])
            seen, mism, dup, alpha_bad, extra = {}, 0, 0, 0, 0
            for pid, s in read_fa(f):
                s = s.upper()
                if pid not in want:
                    extra += 1
                    continue
                if pid in seen:
                    dup += 1
                    continue
                if len(s) != lens[pid]:
                    mism += 1
                    continue
                if set(s) - ALPHA:
                    alpha_bad += 1
                    continue
                seen[pid] = s
            for pid in sub["protein"]:      # gene order of the run
                if pid in seen:
                    fo.write(f">{pid}\n{seen[pid]}\n")
            r = dict(n=len(want), with_3di=len(seen), missing=len(want) - len(seen), length_mismatch=mism,
                     duplicates=dup, bad_alphabet=alpha_bad, extra_ids=extra)
            n_bad += mism + dup + alpha_bad
            rep["per_species"][sp] = r
    os.replace(tmp, out)
    tot = pd.DataFrame(rep["per_species"]).T.sum()
    rep["totals"] = {k: int(v) for k, v in tot.items()}
    rep["coverage"] = round(tot["with_3di"] / tot["n"], 6)
    json.dump(rep, open(out + ".report.json", "w"), indent=1)
    print(json.dumps(rep["totals"]), "coverage", rep["coverage"])
    if n_bad:
        print("WARNING: integrity problems (see report)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
