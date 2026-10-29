#!/usr/bin/env python3
"""Integrity check of ProstT5 3Di FASTAs. Usage: check_3di.py [codes...] (default: all with .done). Writes TSV report."""
import sys, os, glob
P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
D = P + "/data/structure3di/all"
OUT = P + "/results/r2/w5_3di_validation/integrity_report.tsv"
ALPHA = set("ACDEFGHIKLMNPQRSTVWY")

def fasta(path):
    n = None; buf = []
    with open(path) as fh:
        for l in fh:
            l = l.rstrip("\n")
            if l.startswith(">"):
                if n is not None: yield n, "".join(buf)
                n = l[1:].split()[0]; buf = []
            else: buf.append(l.strip())
    if n is not None: yield n, "".join(buf)

def check(code):
    src = next((f"{P}/data/proteomes/{k}/{code}.fa" for k in ("euk", "prok", "dense") if os.path.exists(f"{P}/data/proteomes/{k}/{code}.fa")), None)
    seqs = {}; dup_in = 0
    for n, s in fasta(src):
        if n in seqs: dup_in += 1
        seqs[n] = s
    skipped = set()
    sk = f"{D}/{code}.3di.fa.skipped"
    if os.path.exists(sk):
        skipped = {l.split("\t")[0].strip() for l in open(sk) if l.strip()}
    seen = {}; dup = 0; badalpha = 0; badlen = 0; extra = 0; badchars = set()
    for n, s in fasta(f"{D}/{code}.3di.fa"):
        if n in seen: dup += 1; continue
        seen[n] = 1
        if n not in seqs: extra += 1; continue
        if len(s) != len(seqs[n].rstrip("*")) and len(s) != len(seqs[n]): badlen += 1
        bc = set(s) - ALPHA
        if bc: badalpha += 1; badchars |= bc
    missing = [n for n in seqs if n not in seen and n not in skipped]
    skipped_with_rec = sum(1 for n in skipped if n in seen)
    ok = not (dup or badalpha or badlen or extra or missing or dup_in)
    return [code, len(seqs), len(seen), len(skipped & set(seqs)), len(missing), dup, extra, badlen, badalpha, "".join(sorted(badchars)), skipped_with_rec, dup_in, "PASS" if ok else "FAIL"]

if __name__ == "__main__":
    codes = sys.argv[1:] or sorted(os.path.basename(f)[:-len(".3di.fa.done")] for f in glob.glob(D + "/*.3di.fa.done"))
    hdr = ["code", "n_proteome", "n_3di_records", "n_skipped", "n_missing", "n_dup_ids", "n_extra_ids", "n_length_mismatch", "n_bad_alphabet", "bad_chars", "skipped_but_present", "dup_ids_in_proteome", "status"]
    with open(OUT, "w") as o:
        o.write("\t".join(hdr) + "\n")
        for c in codes:
            r = check(c); o.write("\t".join(map(str, r)) + "\n"); print(*r, sep="\t")
