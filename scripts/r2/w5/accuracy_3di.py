#!/usr/bin/env python3
"""ProstT5 3Di (ours) vs 3Di derived by Foldseek from AFDB models, for proteins in completed species.
Sources: (pdb) AFDB PDB models in runs/structure/{query_models,afdb_ref/*}; (sp) Foldseek AFDB-Swiss-Prot DB (sp_ss; no pLDDT).
Usage: accuracy_3di.py [--max 5000] [--seed 11]. CPU only (foldseek, 8 threads, run under nice)."""
import os, re, gzip, glob, random, subprocess, sys, argparse, collections, shutil
import numpy as np
P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
FS = P + "/conda_envs/ancestra/bin/foldseek"
OD = P + "/results/r2/w5_3di_validation"; W = OD + "/work"
S3 = P + "/data/structure3di/all"; R = P + "/runs/structure"
ap = argparse.ArgumentParser(); ap.add_argument("--max", type=int, default=5000); ap.add_argument("--seed", type=int, default=11)
ap.add_argument("--species", help="file with species codes (frozen list)"); a = ap.parse_args()
os.makedirs(W, exist_ok=True)

def fasta(path, op=open):
    n = None; b = []
    with op(path, "rt") as fh:
        for l in fh:
            l = l.rstrip("\n")
            if l.startswith(">"):
                if n is not None: yield n, "".join(b)
                n = l[1:].split()[0]; b = []
            else: b.append(l.strip())
    if n is not None: yield n, "".join(b)

def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode: sys.exit("FAILED: " + cmd + "\n" + r.stderr[-2000:])

codes = sorted(os.path.basename(f)[:-len(".3di.fa.done")] for f in glob.glob(S3 + "/*.3di.fa.done"))
if a.species: codes = open(a.species).read().split()
models = {}
for d in [R + "/query_models"] + [x for x in glob.glob(R + "/afdb_ref/*") if os.path.isdir(x)]:
    for f in os.listdir(d):
        m = re.match(r"AF-(.+)-F1-model", f)
        if m: models[m.group(1)] = d + "/" + f
spkey = {}
for l in open(R + "/afdb_sp/sp.lookup"):
    p = l.split("\t"); m = re.match(r"AF-(.+)-F1", p[1])
    if m: spkey[m.group(1)] = p[0]
cand = []  # (code,newid,acc,src)
for c in codes:
    if not os.path.exists(f"{P}/data/proteomes/idmaps/{c}.tsv.gz"): print("no idmap, skipped", c); continue
    for l in gzip.open(f"{P}/data/proteomes/idmaps/{c}.tsv.gz", "rt"):
        p = l.rstrip("\n").split("\t")
        if p[0] == "original_id" or p[-1] != "kept": continue
        if p[0] in models: cand.append((c, p[1], p[0], "pdb"))
        elif p[0] in spkey: cand.append((c, p[1], p[0], "sp"))
print("completed species", len(codes), "candidates with AFDB model", len(cand), collections.Counter(x[3] for x in cand))
random.seed(a.seed); cand.sort()
if len(cand) > a.max: cand = random.sample(cand, a.max)
bysrc = collections.defaultdict(list)
for x in cand: bysrc[x[3]].append(x)

afdb = {}  # acc -> (aa, ss, plddt or None)
if bysrc["pdb"]:
    d = W + "/pdbs"; shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    for _, _, acc, _ in bysrc["pdb"]: os.symlink(models[acc], f"{d}/{acc}.pdb.gz")
    for f in glob.glob(W + "/pdbdb*"): os.remove(f)
    run(f"nice -n 10 {FS} createdb {d} {W}/pdbdb --threads 8")
    run(f"{FS} lndb {W}/pdbdb_h {W}/pdbdb_ss_h")
    run(f"{FS} convert2fasta {W}/pdbdb {W}/pdb_aa.fa"); run(f"{FS} convert2fasta {W}/pdbdb_ss {W}/pdb_ss.fa")
    aa = {n.split(".pdb")[0]: s for n, s in fasta(W + "/pdb_aa.fa")}
    ss = {n.split(".pdb")[0]: s for n, s in fasta(W + "/pdb_ss.fa")}
    for _, _, acc, _ in bysrc["pdb"]:
        pl = []
        with gzip.open(models[acc], "rt") as fh:
            for l in fh:
                if l.startswith("ATOM") and l[12:16].strip() == "CA": pl.append(float(l[60:66]))
        if acc in aa and acc in ss: afdb[acc] = (aa[acc], ss[acc], float(np.mean(pl)) if pl else None)
if bysrc["sp"]:
    keys = W + "/spkeys.txt"; open(keys, "w").write("\n".join(spkey[x[2]] for x in bysrc["sp"]) + "\n")
    for f in glob.glob(W + "/spsub*"): os.remove(f)
    for s in ("", "_ss", "_h"): run(f"{FS} createsubdb {keys} {R}/afdb_sp/sp{s} {W}/spsub{s}")
    run(f"{FS} lndb {W}/spsub_h {W}/spsub_ss_h")
    run(f"{FS} convert2fasta {W}/spsub {W}/sp_aa.fa"); run(f"{FS} convert2fasta {W}/spsub_ss {W}/sp_ss.fa")
    g = lambda n: re.match(r"AF-(.+)-F1", n).group(1)
    aa = {g(n): s for n, s in fasta(W + "/sp_aa.fa")}; ss = {g(n): s for n, s in fasta(W + "/sp_ss.fa")}
    for _, _, acc, _ in bysrc["sp"]:
        if acc in aa and acc in ss: afdb[acc] = (aa[acc], ss[acc], None)

# optional pLDDT cache for sp-sourced proteins (scripts/r2/w5/fetch_plddt.py; AFDB API mean pLDDT)
cf = OD + "/plddt_cache.tsv"
if os.path.exists(cf):
    pc = {}
    for l in open(cf):
        k, v = l.rstrip("\n").split("\t")
        if v: pc[k] = float(v)
    for k, (x, y, z) in list(afdb.items()):
        if z is None and k in pc: afdb[k] = (x, y, pc[k])
# our data
ours = {}; seqs = {}
for c in sorted({x[0] for x in cand}):
    want = {x[1] for x in cand if x[0] == c}
    for n, s in fasta(f"{S3}/{c}.3di.fa"):
        if n in want: ours[n] = s
    src = next(f"{P}/data/proteomes/{k}/{c}.fa" for k in ("euk", "prok", "dense") if os.path.exists(f"{P}/data/proteomes/{k}/{c}.fa"))
    for n, s in fasta(src):
        if n in want: seqs[n] = s.rstrip("*")
dis = {}
for c in sorted({x[0] for x in cand}):
    f = f"{P}/data/annotations/disorder/euk/{c}.disorder.tsv.gz"
    if os.path.exists(f):
        for l in gzip.open(f, "rt"):
            p = l.split("\t")
            if p[0] != "protein": dis[p[0]] = float(p[2])

def lbin(n): return "<=1000 (unwindowed)" if n <= 1000 else ("1000-3000 (windowed)" if n <= 3000 else ">3000 (not predicted)")
def pbin(x): return "NA" if x is None else ("<50" if x < 50 else "50-70" if x < 70 else "70-90" if x < 90 else ">=90")
def dbin(x): return "NA" if x is None else ("<0.1" if x < 0.1 else "0.1-0.4" if x < 0.4 else ">=0.4")
rows = []; n_mis = n_no3 = 0; mism = []
for c, nid, acc, src in cand:
    if acc not in afdb: continue
    aas, ss, pl = afdb[acc]
    if nid not in ours: n_no3 += 1; continue
    if aas != seqs[nid]: n_mis += 1; mism.append((nid, acc, len(aas), len(seqs[nid]))); continue
    o = ours[nid]
    if len(o) != len(ss): n_mis += 1; continue
    ok = np.frombuffer(ss.encode(), np.uint8) != ord("X")
    eq = (np.frombuffer(o.encode(), np.uint8) == np.frombuffer(ss.encode(), np.uint8)) & ok
    rows.append(dict(protein=nid, acc=acc, src=src, L=len(o), nvalid=int(ok.sum()), nagree=int(eq.sum()), acc_p=eq.sum() / max(ok.sum(), 1),
                     plddt=pl, disfrac=dis.get(nid), lbin=lbin(len(o)), pbin=pbin(pl), dbin=dbin(dis.get(nid))))
with open(OD + "/per_protein_accuracy.tsv", "w") as f:
    ks = list(rows[0].keys()); f.write("\t".join(ks) + "\n")
    for r in rows: f.write("\t".join("" if r[k] is None else (f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])) for k in ks) + "\n")
open(OD + "/length_or_sequence_mismatches.tsv", "w").write("protein\tacc\tafdb_len\tour_len\n" + "".join("\t".join(map(str, m)) + "\n" for m in mism))
# baseline: agreement expected if ProstT5 letters were independent of AFDB letters (same marginals)
fo = collections.Counter(); fa = collections.Counter()
for r in rows:
    fo.update(ours[r["protein"]]); fa.update(afdb[r["acc"]][1])
tot = sum(fo.values()); base = sum(fo[k] / tot * fa[k] / sum(fa.values()) for k in fo)
top = max(fa.values()) / sum(fa.values())
def summ(name, key):
    out = []
    g = collections.defaultdict(list)
    for r in rows: g[r[key]].append(r)
    for k in sorted(g):
        v = g[k]; N = sum(r["nvalid"] for r in v); A = sum(r["nagree"] for r in v)
        ap_ = np.array([r["acc_p"] for r in v])
        out.append(f"{name}\t{k}\t{len(v)}\t{N}\t{A/N:.4f}\t{np.mean(ap_):.4f}\t{np.median(ap_):.4f}")
    return out
lines = ["stratum\tbin\tn_proteins\tn_residues\tpooled_identity\tmean_per_protein\tmedian_per_protein"]
N = sum(r["nvalid"] for r in rows); A = sum(r["nagree"] for r in rows)
lines.append(f"all\tall\t{len(rows)}\t{N}\t{A/N:.4f}\t{np.mean([r['acc_p'] for r in rows]):.4f}\t{np.median([r['acc_p'] for r in rows]):.4f}")
for nm, k in (("length", "lbin"), ("plddt", "pbin"), ("disorder", "dbin"), ("source", "src")): lines += summ(nm, k)
open(OD + "/accuracy_by_stratum.tsv", "w").write("\n".join(lines) + "\n")
open(OD + "/accuracy_meta.txt", "w").write(f"completed_species\t{len(codes)}\ncandidates_sampled\t{len(cand)}\nwith_afdb_ss\t{len(afdb)}\nseq_or_len_mismatch_skipped\t{n_mis}\nno_3di_record\t{n_no3}\nevaluated\t{len(rows)}\nindependent_baseline_identity\t{base:.4f}\nmost_common_state_frequency\t{top:.4f}\nspecies\t{','.join(codes)}\n")
print("\n".join(lines)); print(open(OD + "/accuracy_meta.txt").read())
