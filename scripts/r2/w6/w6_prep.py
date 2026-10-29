#!/usr/bin/env python
"""W6 (exploratory): structural (3Di) vs sequence (AA) divergence per branch.

Steps (per panel R; outputs under results/r2/w6/<R>/work):
    w6_prep.py select R   # strictly single-copy root HOGs in >= 80 % of panel species; <= 300 (seed 11)
    w6_prep.py pfam   R   # Pfam-A 38.2 hmmsearch (gathering cutoffs) on the selected proteins (pyhmmer)
    w6_prep.py align  R   # MAFFT (L-INS-i via --auto) on AA; 3Di projected onto the same columns;
                          # residues outside Pfam envelopes masked; columns kept if >= 50 % of rows are in-domain;
                          # taxa with < 50 % in-domain residues over kept columns dropped; >= 50 columns, >= 8 taxa
    w6_prep.py model  R   # panel 3Di exchangeability matrix by pair counting (see cmd_model; IQ-TREE GTR20 failed)
    w6_prep.py fit    R   # per gene, fixed pruned topology: AA LG+G4; 3Di <panel 3Di matrix>+G4; shuffled-3Di (column-wise
                          # permutation among taxa) <panel 3Di matrix>+G4; 3Di Poisson+FO+G4 (sensitivity); AA split-half (random column halves) LG+G4 x2
Seed 11 throughout. nice -n 10; <= 14 threads.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
sys.path.insert(0, ROOT)
from mcell.phylo import read_tree, prune_named, transfer_branch_lengths  # noqa: E402

OUT = os.path.join(ROOT, "results", "r2", "w6")
ANNOT_PY = os.path.join(ROOT, "conda_envs", "annot", "bin", "python")
BIN = os.path.join(ROOT, "conda_envs", "ancestra", "bin")
PFAM_HMM = os.path.join(ROOT, "data", "raw", "pfam", "Pfam-A.hmm")
DI_DIR = os.path.join(ROOT, "data", "structure3di", "all")
SEED = 11
MAX_GENES = 300
OCC = 0.8
MAX_MULTI = 0.1
AA = "ACDEFGHIKLMNPQRSTVWY"


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", file=sys.stderr, flush=True)


def P(R):
    w = os.path.join(OUT, R, "work")
    return dict(work=w, genes=os.path.join(w, "genes"), sel=os.path.join(w, "selected.tsv"),
                pfam=os.path.join(w, "pfam_hits.tsv.gz"), cache=os.path.join(ROOT, "results", "r2", "panels", R, "work", "cache"),
                guide=os.path.join(ROOT, "data", "dense", R, "guide_tree.nwk"),
                proteomes=os.path.join(ROOT, "data", "dense", R, "proteomes"),
                model=os.path.join(w, "model3di"))


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[:3]} failed ({r.returncode}): {r.stderr[-1500:]}")
    return r


def read_fasta(path, keep=None):
    out, h, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None and (keep is None or h in keep):
                    out[h] = "".join(buf)
                h, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if h is not None and (keep is None or h in keep):
        out[h] = "".join(buf)
    return out


def write_fasta(path, d):
    with open(path, "w") as fh:
        for k, v in d.items():
            fh.write(f">{k}\n{v}\n")


# ---------------------------------------------------------------------------
def cmd_select(a):
    p = P(a.R)
    os.makedirs(p["genes"], exist_ok=True)
    codes = sorted(f[:-3] for f in os.listdir(p["proteomes"]) if f.endswith(".fa"))
    missing = [c for c in codes if not os.path.exists(os.path.join(DI_DIR, f"{c}.3di.fa.done"))]
    if missing:
        raise SystemExit(f"3Di incomplete for {missing}")
    N = len(codes)
    mem = pd.read_parquet(os.path.join(p["cache"], "members.parquet"))
    genes = pd.read_parquet(os.path.join(p["cache"], "genes.parquet"), columns=["gene", "species"])
    hogs = pd.read_parquet(os.path.join(p["cache"], "hogs.parquet"), columns=["hog_id"])
    root = read_tree(p["guide"]).name
    m = mem[mem["level"] == root]
    m = pd.DataFrame({"hog_id": hogs["hog_id"].to_numpy()[m["hog_row"].to_numpy()],
                      "gene": genes["gene"].to_numpy()[m["gene_row"].to_numpy()],
                      "species": genes["species"].astype(str).to_numpy()[m["gene_row"].to_numpy()]})
    cn = m.groupby(["hog_id", "species"]).size().rename("n").reset_index()
    agg = cn.groupby("hog_id").agg(n_sp=("n", "size"), n_multi=("n", lambda x: int((x > 1).sum())))
    agg["n_sc"] = cn[cn["n"] == 1].groupby("hog_id").size().reindex(agg.index, fill_value=0)
    # one gene per species: species with >1 copy are dropped from the gene (allowed in <= 10 % of species);
    # >= 80 % of panel species must be single-copy (same rule as the panel phylogenomic markers)
    cand = agg[(agg["n_multi"] <= MAX_MULTI * N) & (agg["n_sc"] >= OCC * N)]
    log(f"{a.R}: {len(agg):,} root HOGs, {len(cand)} one-gene-per-species (<=10% multi-copy dropped) with >= {OCC:.0%} of {N} species")
    n_cand = len(cand)
    rng = np.random.default_rng(SEED)
    cand = cand.assign(r=rng.random(len(cand))).sort_values("r")
    # take up to 1.5x MAX_GENES so that the Pfam/alignment filters can still leave MAX_GENES
    cand = cand.head(int(MAX_GENES * 1.5))
    sc = cn[(cn["n"] == 1) & cn["hog_id"].isin(cand.index)][["hog_id", "species"]]
    mm = m[m["hog_id"].isin(cand.index)].merge(sc, on=["hog_id", "species"])
    want = set(mm["gene"])
    aa, di = {}, {}
    for c in codes:
        aa.update(read_fasta(os.path.join(p["proteomes"], f"{c}.fa"), keep=want))
        di.update(read_fasta(os.path.join(DI_DIR, f"{c}.3di.fa"), keep=want))
    rows, bad = [], 0
    for k, (hog, g) in enumerate(mm.groupby("hog_id", sort=False)):
        seqs, dis = {}, {}
        for r in g.itertuples():
            s = aa.get(r.gene, "")
            d = di.get(r.gene)
            if d is None:
                bad += 1
                continue
            s2 = s.rstrip("*")
            if len(s2) != len(d):
                if len(s) == len(d):
                    s2 = s
                else:
                    bad += 1
                    continue
            seqs[r.species] = (r.gene, re.sub(r"[^A-Z]", "X", s2.upper()).replace("U", "X").replace("*", "X"), d.upper())
        nm = f"g{k:03d}"
        with open(os.path.join(p["genes"], f"{nm}.raw.fa"), "w") as fa, \
             open(os.path.join(p["genes"], f"{nm}.raw3di.fa"), "w") as fd:
            for sp, (gid, s, d) in sorted(seqs.items()):
                fa.write(f">{sp}\n{s}\n")
                fd.write(f">{sp}\n{d}\n")
        rows.append(dict(gene=nm, hog_id=hog, n_species=len(seqs),
                         members=";".join(f"{sp}={gid}" for sp, (gid, _, _) in sorted(seqs.items()))))
    pd.DataFrame(rows).to_csv(p["sel"], sep="\t", index=False)
    json.dump(dict(panel=a.R, n_species=N, n_root_hogs=int(len(agg)), n_candidates=int(n_cand), n_strict_single_copy=int(((agg["n_multi"] == 0) & (agg["n_sp"] >= OCC * N)).sum()),
                   n_drawn=len(rows), n_member_seq_problems=bad, seed=SEED),
              open(os.path.join(p["work"], "select_info.json"), "w"), indent=1)
    log(f"{a.R}: wrote {len(rows)} genes ({bad} member sequences without matching 3Di)")


def cmd_pfam(a):
    p = P(a.R)
    sel = pd.read_csv(p["sel"], sep="\t")
    fa = os.path.join(p["work"], "pfam_input.fa")
    with open(fa, "w") as out:
        for g in sel["gene"]:
            for sp, s in read_fasta(os.path.join(p["genes"], f"{g}.raw.fa")).items():
                out.write(f">{g}|{sp}\n{s}\n")
    code = f"""
import gzip, os, pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile
abc = Alphabet.amino()
with SequenceFile({fa!r}, digital=True, alphabet=abc) as sf:
    seqs = sf.read_block()
with HMMFile({PFAM_HMM!r}) as hf:
    hmms = list(hf)
with gzip.open({p['pfam']!r} + '.tmp', 'wt') as out:
    out.write('protein\\tpfam_acc\\tpfam_name\\tenv_from\\tenv_to\\tdom_score\\n')
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus={a.threads}, bit_cutoffs='gathering'):
        acc = hits.query.accession; acc = (acc.decode() if isinstance(acc, bytes) else acc).split('.')[0]
        name = hits.query.name; name = name.decode() if isinstance(name, bytes) else name
        for hit in hits.included:
            pid = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            for d in hit.domains.included:
                out.write(f'{{pid}}\\t{{acc}}\\t{{name}}\\t{{d.env_from}}\\t{{d.env_to}}\\t{{d.score:.1f}}\\n')
os.replace({p['pfam']!r} + '.tmp', {p['pfam']!r})
"""
    t0 = time.time()
    subprocess.run(["nice", "-n", "10", ANNOT_PY, "-c", code], check=True)
    log(f"{a.R}: pfam done in {time.time() - t0:.0f}s")


def _align_one(args):
    g, p = args
    raw = os.path.join(p["genes"], f"{g}.raw.fa")
    aln = os.path.join(p["genes"], f"{g}.mafft.fa")
    if not (os.path.exists(aln) and os.path.getsize(aln) > 0):
        r = run(["nice", "-n", "10", os.path.join(BIN, "mafft"), "--auto", "--thread", "1", "--quiet", raw])
        with open(aln + ".tmp", "w") as fh:
            fh.write(r.stdout)
        os.replace(aln + ".tmp", aln)
    return aln


def cmd_align(a):
    p = P(a.R)
    sel = pd.read_csv(p["sel"], sep="\t")
    dom = pd.read_csv(p["pfam"], sep="\t")
    dom[["gene", "species"]] = dom["protein"].str.split("|", expand=True)
    with ThreadPoolExecutor(a.threads) as ex:
        list(ex.map(_align_one, [(g, p) for g in sel["gene"]]))
    rng = np.random.default_rng(SEED)
    info = []
    for g in sel["gene"]:
        al = read_fasta(os.path.join(p["genes"], f"{g}.mafft.fa"))
        di = read_fasta(os.path.join(p["genes"], f"{g}.raw3di.fa"))
        dg = dom[dom["gene"] == g]
        sps = sorted(al)
        L = len(al[sps[0]])
        A = np.array([list(al[s]) for s in sps])
        D = np.full(A.shape, "-", dtype="<U1")
        IN = np.zeros(A.shape, dtype=bool)
        pfams = set()
        for i, s in enumerate(sps):
            res_pos = np.where(A[i] != "-")[0]
            assert len(res_pos) == len(di[s]), (g, s)
            D[i, res_pos] = list(di[s])
            ind = np.zeros(len(res_pos), dtype=bool)
            for r in dg[dg["species"] == s].itertuples():
                ind[r.env_from - 1:r.env_to] = True
                pfams.add(r.pfam_name)
            IN[i, res_pos[ind]] = True
        keepc = IN.mean(axis=0) >= 0.5
        nk = int(keepc.sum())
        row = dict(gene=g, hog_id=sel.loc[sel["gene"] == g, "hog_id"].iloc[0], n_raw=len(sps), aln_len=L, dom_cols=nk,
                   pfams=";".join(sorted(pfams)))
        if nk < 50:
            info.append({**row, "status": "few_domain_columns"})
            continue
        A2, D2, IN2 = A[:, keepc].copy(), D[:, keepc].copy(), IN[:, keepc]
        A2[~IN2] = "-"   # residues outside the taxon's own Pfam envelopes are masked in both alphabets
        D2[~IN2] = "-"
        cov = (A2 != "-").mean(axis=1)
        keepr = cov >= 0.5
        if keepr.sum() < 8:
            info.append({**row, "status": "few_taxa"})
            continue
        sps2 = [s for s, k in zip(sps, keepr) if k]
        A2, D2 = A2[keepr], D2[keepr]
        # shuffled-3Di negative control: permute 3Di characters among taxa within each column (gap pattern of the
        # column is permuted with them), destroying the phylogenetic structure of the 3Di but keeping composition
        S2 = D2.copy()
        for j in range(S2.shape[1]):
            S2[:, j] = S2[rng.permutation(S2.shape[0]), j]
        # AA split-half null: random halves of the columns
        half = rng.permutation(A2.shape[1]) < A2.shape[1] // 2
        out = os.path.join(p["genes"], g)
        write_fasta(out + ".aa.fa", {s: "".join(r) for s, r in zip(sps2, A2)})
        write_fasta(out + ".3di.fa", {s: "".join(r) for s, r in zip(sps2, D2)})
        write_fasta(out + ".shuf3di.fa", {s: "".join(r) for s, r in zip(sps2, S2)})
        write_fasta(out + ".aaH1.fa", {s: "".join(r) for s, r in zip(sps2, A2[:, half])})
        write_fasta(out + ".aaH2.fa", {s: "".join(r) for s, r in zip(sps2, A2[:, ~half])})
        info.append({**row, "status": "ok", "n_taxa": len(sps2), "n_cols": int(A2.shape[1]),
                     "aa_identity_proxy": float(np.mean([(A2[:, j] != "-").mean() for j in range(A2.shape[1])]))})
    df = pd.DataFrame(info)
    ok = df[df["status"] == "ok"]
    if len(ok) > MAX_GENES:
        df.loc[ok.index[MAX_GENES:], "status"] = "over_cap"
    df.to_csv(os.path.join(p["work"], "alignments.tsv"), sep="\t", index=False)
    log(f"{a.R}: {df['status'].value_counts().to_dict()}")


# ---------------------------------------------------------------------------
def write_topology(guide, sps, path):
    t = prune_named(guide, sps)
    with open(path, "w") as fh:
        fh.write(t.write(format=9) + "\n")
    return t


def parse_gtr20(iqtree_file):
    """Exchangeabilities (dict pair->rate) and state frequencies from an IQ-TREE .iqtree report."""
    txt = open(iqtree_file).read()
    rates, freqs = {}, {}
    m = re.search(r"Rate parameter R:\s*\n\n(.*?)\n\n", txt, re.S)
    for line in m.group(1).splitlines():
        k, v = line.split(":")
        a_, b_ = k.strip().split("-")
        rates[(a_, b_)] = float(v)
    m = re.search(r"State frequencies:.*?\n\n(.*?)\n\n", txt, re.S)
    for line in m.group(1).splitlines():
        k, v = line.split("=")
        freqs[k.strip().replace("pi(", "").replace(")", "")] = float(v)
    return rates, freqs


THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
         "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}
PAML_ORDER = "ARNDCQEGHILKMFPSTWYV"


def write_paml(rates, freqs, path):
    inv3 = {v: k for k, v in THREE.items()}

    def key(x):
        return inv3.get(x, x)
    R = {}
    for (x, y), v in rates.items():
        x, y = key(x), key(y)
        R[(x, y)] = R[(y, x)] = v
    F = {key(k): v for k, v in freqs.items()}
    with open(path, "w") as fh:
        for i in range(1, 20):
            fh.write(" ".join(f"{R[(PAML_ORDER[i], PAML_ORDER[j])]:.6f}" for j in range(i)) + "\n")
        fh.write("\n" + " ".join(f"{F[c]:.6f}" for c in PAML_ORDER) + "\n")


def cmd_model(a):
    """Panel 3Di exchangeability matrix by pair counting (Dayhoff-style, small-divergence approximation).

    IQ-TREE 3.1.4 GTR20+FO+G4 on the concatenated 3Di did not estimate the exchangeabilities (all stayed at the
    initial value 10 on R2 = a Poisson model) and aborted on a subset (assertion in optimizeAllBranches); see
    results/r2/w6/<R>/work/model3di/gtr20.log. Instead, for every gene and every pair of taxa whose domain-masked AA
    p-distance is 0.05-0.30 and 3Di p-distance <= 0.40, aligned 3Di letter pairs (both non-gap) are counted
    symmetrically: N_ij. With pi the panel 3Di letter frequencies, exchangeability r_ij = N_ij / (pi_i pi_j)
    (since N_ij ~ pi_i Q_ij t = pi_i pi_j r_ij t at small t), scaled to mean 1, written as a PAML file with pi.
    """
    p = P(a.R)
    os.makedirs(p["model"], exist_ok=True)
    al = pd.read_csv(os.path.join(p["work"], "alignments.tsv"), sep="\t")
    ok = al[al["status"] == "ok"]
    idx = {c: i for i, c in enumerate(PAML_ORDER)}
    N = np.zeros((20, 20))
    F = np.zeros(20)
    npairs = 0
    for g in ok["gene"]:
        A = read_fasta(os.path.join(p["genes"], f"{g}.aa.fa"))
        D = read_fasta(os.path.join(p["genes"], f"{g}.3di.fa"))
        sps = sorted(A)
        Aa = np.array([list(A[s]) for s in sps])
        Dd = np.array([list(D[s]) for s in sps])
        for ch in Dd[Dd != "-"]:
            if ch in idx:
                F[idx[ch]] += 1
        for i in range(len(sps)):
            for j in range(i + 1, len(sps)):
                m = (Aa[i] != "-") & (Aa[j] != "-") & (Dd[i] != "-") & (Dd[j] != "-")
                if m.sum() < 30:
                    continue
                pa = float((Aa[i][m] != Aa[j][m]).mean())
                pd_ = float((Dd[i][m] != Dd[j][m]).mean())
                if not (0.05 <= pa <= 0.30 and pd_ <= 0.40):
                    continue
                npairs += 1
                for x, y in zip(Dd[i][m], Dd[j][m]):
                    if x != y and x in idx and y in idx:
                        N[idx[x], idx[y]] += 1
                        N[idx[y], idx[x]] += 1
    pi = F / F.sum()
    R = N / np.outer(pi, pi)
    np.fill_diagonal(R, 0)
    R = R / R[np.triu_indices(20, 1)].mean()
    R[R == 0] = 1e-3
    with open(os.path.join(p["model"], "panel3di.paml"), "w") as fh:
        for i in range(1, 20):
            fh.write(" ".join(f"{R[i, j]:.6f}" for j in range(i)) + "\n")
        fh.write("\n" + " ".join(f"{x:.6f}" for x in pi) + "\n")
    pd.DataFrame(R, index=list(PAML_ORDER), columns=list(PAML_ORDER)).to_csv(os.path.join(p["model"], "exchangeabilities.tsv"), sep="\t")
    json.dump(dict(method="pair counting, AA p 0.05-0.30, 3Di p <= 0.40, r_ij = N_ij/(pi_i pi_j), mean 1",
                   n_pairs=npairs, n_substitution_pairs=int(N.sum() / 2), n_genes=int(len(ok)),
                   pi={c: round(float(x), 4) for c, x in zip(PAML_ORDER, pi)},
                   r_min=float(R[np.triu_indices(20, 1)].min()), r_max=float(R.max())),
              open(os.path.join(p["model"], "model_info.json"), "w"), indent=1)
    log(f"{a.R}: 3Di exchangeabilities from {npairs} pairs, {int(N.sum() / 2)} substitutions")


def _fit_one(args):
    g, p, guide = args
    gd = os.path.join(p["genes"])
    aa = read_fasta(os.path.join(gd, f"{g}.aa.fa"))
    sps = sorted(aa)
    te = os.path.join(gd, f"{g}.topo.nwk")
    t = write_topology(guide, sps, te)
    paml = os.path.join(p["model"], "panel3di.paml")
    jobs = [("aa", f"{g}.aa.fa", "LG+G4"), ("di", f"{g}.3di.fa", f"{paml}+G4"),
            ("shuf", f"{g}.shuf3di.fa", f"{paml}+G4"), ("dipois", f"{g}.3di.fa", "Poisson+FO+G4"), ("aaH1", f"{g}.aaH1.fa", "LG+G4"), ("aaH2", f"{g}.aaH2.fa", "LG+G4")]
    out = {}
    for tag, fa, model in jobs:
        prefix = os.path.join(gd, f"{g}.{tag}")
        if not os.path.exists(prefix + ".treefile"):
            run(["nice", "-n", "10", os.path.join(BIN, "iqtree3"), "-s", os.path.join(gd, fa), "-te", te, "-m", model,
                 "-T", "1", "-seed", str(SEED), "-keep-ident", "--prefix", prefix, "-redo", "-quiet"])
        out[tag] = transfer_branch_lengths(t, read_tree(prefix + ".treefile"))
    gpar = {n.name: (n.up.name if n.up else None) for n in guide.traverse()}
    rows = []
    for n in t.traverse():
        if n.up is None:
            continue
        rows.append(dict(gene=g, node=n.name, parent=n.up.name, single_edge=gpar[n.name] == n.up.name,
                         root_child=n.up.up is None, **{f"bl_{k}": v[n.name] for k, v in out.items()}))
    return rows


def cmd_fit(a):
    p = P(a.R)
    al = pd.read_csv(os.path.join(p["work"], "alignments.tsv"), sep="\t")
    genes = list(al.loc[al["status"] == "ok", "gene"])
    guide = read_tree(p["guide"])
    t0 = time.time()
    rows = []
    with ThreadPoolExecutor(a.threads) as ex:
        for k, r in enumerate(ex.map(_fit_one, [(g, p, guide) for g in genes])):
            rows.extend(r)
            if k % 50 == 0:
                log(f"{a.R}: {k + 1}/{len(genes)} genes fitted ({time.time() - t0:.0f}s)")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, a.R, "branch_lengths.tsv.gz"), sep="\t", index=False)
    log(f"{a.R}: fitted {len(genes)} genes in {time.time() - t0:.0f}s")


def _sim_one(args):
    """Constant-ratio parametric null for one gene: 3Di simulated on the gene's fitted AA tree scaled by the gene's
    total 3Di/AA tree-length ratio (so the true 3Di/AA ratio is identical on every branch), panel 3Di matrix + G4 with
    the gene's fitted 3Di gamma shape, the real 3Di gap pattern imposed, then refitted exactly like the real 3Di."""
    g, p, guide, seed = args
    gd = p["genes"]
    paml = os.path.join(p["model"], "panel3di.paml")
    pre = os.path.join(gd, f"{g}.sim")
    if not os.path.exists(pre + ".treefile"):
        taa = read_tree(os.path.join(gd, f"{g}.aa.treefile"))
        tdi = read_tree(os.path.join(gd, f"{g}.di.treefile"))
        c = sum(n.dist for n in tdi.traverse()) / max(sum(n.dist for n in taa.traverse()), 1e-9)
        for n in taa.traverse():
            n.dist *= c
        alpha = float(re.search(r"Gamma shape alpha: ([0-9.eE+-]+)", open(os.path.join(gd, f"{g}.di.iqtree")).read()).group(1))
        real = read_fasta(os.path.join(gd, f"{g}.3di.fa"))
        L = len(next(iter(real.values())))
        st = pre + ".simtree.nwk"
        with open(st, "w") as fh:
            fh.write(taa.write(format=5) + "\n")
        run(["nice", "-n", "10", os.path.join(BIN, "iqtree3"), "--alisim", pre + ".raw", "-m", f"{paml}+G4{{{alpha}}}",
             "--seqtype", "AA", "-t", st, "--length", str(L), "-af", "fasta", "-seed", str(seed), "-redo", "-quiet"])
        sim = read_fasta(pre + ".raw.fa")
        out = {}
        for sp, r in real.items():
            out[sp] = "".join("-" if a_ == "-" else b_ for a_, b_ in zip(r, sim[sp]))
        write_fasta(pre + ".fa", out)
        run(["nice", "-n", "10", os.path.join(BIN, "iqtree3"), "-s", pre + ".fa", "-te", os.path.join(gd, f"{g}.topo.nwk"),
             "-m", f"{paml}+G4", "-T", "1", "-seed", str(SEED), "-keep-ident", "--prefix", pre, "-redo", "-quiet"])
    t = prune_named(guide, sorted(read_fasta(os.path.join(gd, f"{g}.aa.fa"))))   # named internal nodes
    bl = transfer_branch_lengths(t, read_tree(pre + ".treefile"))
    return [dict(gene=g, node=k, bl_sim=v) for k, v in bl.items()]


def cmd_simnull(a):
    p = P(a.R)
    al = pd.read_csv(os.path.join(p["work"], "alignments.tsv"), sep="\t")
    genes = list(al.loc[al["status"] == "ok", "gene"])
    guide = read_tree(p["guide"])
    rows = []
    with ThreadPoolExecutor(a.threads) as ex:
        for r in ex.map(_sim_one, [(g, p, guide, SEED + k) for k, g in enumerate(genes)]):
            rows.extend(r)
    pd.DataFrame(rows).to_csv(os.path.join(OUT, a.R, "sim_lengths.tsv.gz"), sep="\t", index=False)
    log(f"{a.R}: constant-ratio simulation null for {len(genes)} genes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["select", "pfam", "align", "model", "fit", "simnull"])
    ap.add_argument("R")
    ap.add_argument("--threads", type=int, default=14)
    a = ap.parse_args()
    {"select": cmd_select, "pfam": cmd_pfam, "align": cmd_align, "model": cmd_model, "fit": cmd_fit, "simnull": cmd_simnull}[a.cmd](a)


if __name__ == "__main__":
    main()
