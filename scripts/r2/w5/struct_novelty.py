#!/usr/bin/env python
"""W5 task B: are sequence-novel gained families also structurally new?

For a dense panel R with a finished gain ladder (results/r2/panels/<R>/work/ladder, written by the panel
agent; read-only) and complete ProstT5 3Di for all its species:

  prep R    sample families per laddered branch (seed 11):
              novel  = up to N_NOVEL sequence-novel gains (ladder == 'seq_novel')
              old    = up to N_OLD sequence-old gains (ladder == 'homologous') whose DIAMOND outside-clade hit
                       (best_out_species of any representative) is a PANEL species  -> positive control
              shuf   = the representatives of the 'old' sample with AA and 3Di shuffled by one shared random
                       permutation of positions (composition and AA/3Di pairing kept)   -> negative control
            queries = the ladder's own representatives (<= 3 per family, per_query.parquet).
            Builds one Foldseek DB of all panel proteins (AA + ProstT5 3Di as <db>_ss; base:createdb of the
            AA and 3Di FASTAs with --shuffle 0 so keys agree, as ancestra.structure.build_db) and one query DB.
  search R  per branch: target = panel proteins of species OUTSIDE the branch's clade (createsubdb of DB,
            DB_ss, DB_h); foldseek search --alignment-type 2 (3Di+AA), -e 1e-3, default sensitivity,
            --max-seqs 1000; 12 threads.
  classify R  a family is 'structurally old' if any representative has a hit with E < 1e-3 whose alignment
            covers >= 50 % of the query ((qend - qstart + 1) / qlen); otherwise 'structurally new'.
Outputs: results/r2/w5/novelty/<R>/{sample.tsv, hits.tsv.gz, families.tsv, branches.tsv, DONE.json}.
"""
import json, os, subprocess, sys, time

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

FS = os.path.join(ROOT, "conda_envs/ancestra/bin/foldseek")
D3 = os.path.join(ROOT, "data/structure3di/all")
OUT = os.path.join(ROOT, "results/r2/w5/novelty")
WK = os.path.join(ROOT, "results/r2/w5/work/novelty")
N_NOVEL, N_OLD, SEED, THREADS = 200, 60, 11, 12
EVAL, QCOV = 1e-3, 0.5
FMT = "query,target,evalue,bits,qstart,qend,qlen,tstart,tend,tlen,alnlen,fident"


def log(m):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}", flush=True)


def run(cmd, logf):
    cmd = [str(c) for c in cmd]
    with open(logf, "a") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        r = subprocess.run(["nice", "-n", "10"] + cmd, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"failed ({r.returncode}): {' '.join(cmd)}")


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


def panel_species(R):
    return sorted(f.rsplit(".", 1)[0] for f in os.listdir(os.path.join(ROOT, "data/dense", R, "proteomes")))


#: species whose 3Di was never queued for ProstT5 (not in data/structure3di/order_r2b.txt): searched without them
#: (their proteins are neither queries nor targets); recorded in sample_info.json
ALLOW_MISSING = {"R4": {"Eudosp"}}


def three_di_ready(R):
    return [s for s in panel_species(R) if not os.path.exists(os.path.join(D3, f"{s}.3di.fa.done"))
            and s not in ALLOW_MISSING.get(R, set())]


def clades(R):
    from ete3 import Tree
    t = Tree(open(os.path.join(ROOT, "runs/r2", R, "stages/input/guide_tree.nwk")).read().strip(), format=1)
    return {n.name: set(n.get_leaf_names()) for n in t.traverse()}


def build_db(pairs, prefix, logf):
    """pairs: iterable of (id, aa, 3di). Returns ids in key order."""
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    ids = []
    with open(prefix + ".aa.fa", "w") as fa, open(prefix + ".3di.fa", "w") as fs:
        for i, a, s in pairs:
            if len(a) != len(s):
                raise ValueError(f"length mismatch {i}")
            fa.write(f">{i}\n{a}\n")
            fs.write(f">{i}\n{s}\n")
            ids.append(i)
    run([FS, "base:createdb", prefix + ".aa.fa", prefix, "--shuffle", "0", "-v", "1"], logf)
    run([FS, "base:createdb", prefix + ".3di.fa", prefix + "_sstmp", "--shuffle", "0", "-v", "1"], logf)
    for suf in ("", ".index", ".dbtype"):
        os.replace(prefix + "_sstmp" + suf, prefix + "_ss" + suf)
    for f in os.listdir(os.path.dirname(prefix)):
        if f.startswith(os.path.basename(prefix) + "_sstmp"):
            os.remove(os.path.join(os.path.dirname(prefix), f))
    os.remove(prefix + ".aa.fa")
    os.remove(prefix + ".3di.fa")
    # key order check against the lookup
    lk = pd.read_csv(prefix + ".lookup", sep="\t", header=None, usecols=[0, 1])
    assert list(lk[1]) == ids and list(lk[0]) == list(range(len(ids))), "key order mismatch"
    return ids


def cmd_prep(R):
    miss = three_di_ready(R)
    if miss:
        sys.exit(f"{R}: 3Di not complete for {len(miss)} species: {miss[:8]}")
    od, wd = os.path.join(OUT, R), os.path.join(WK, R)
    if os.path.exists(os.path.join(od, "sample.tsv")):
        sys.exit(f"refusing to overwrite {od}/sample.tsv")
    os.makedirs(od, exist_ok=True)
    os.makedirs(wd, exist_ok=True)
    logf = os.path.join(wd, "commands.log")
    L = os.path.join(ROOT, "results/r2/panels", R, "work/ladder")
    cls = pd.read_parquet(os.path.join(L, "classified.parquet"))
    pq = pd.read_parquet(os.path.join(L, "per_query.parquet"))
    sp = set(panel_species(R))
    cl = clades(R)
    rng = np.random.default_rng(SEED)
    rows = []
    for br, sub in cls.groupby("branch_child", sort=True):
        inside = cl[br]
        outside_panel = sp - inside
        nov = sub[sub["ladder"] == "seq_novel"]["family_id"].to_numpy()
        hom = sub[sub["ladder"] == "homologous"]["family_id"]
        q = pq[pq["family_id"].isin(set(hom))]
        okf = set(q.loc[q["best_out_species"].isin(outside_panel), "family_id"])
        hom = np.array(sorted(okf))
        for kind, arr, n in (("novel", nov, N_NOVEL), ("old", hom, N_OLD)):
            pick = np.sort(rng.choice(arr, size=min(n, len(arr)), replace=False)) if len(arr) else []
            for f in pick:
                rows.append(dict(branch_child=br, kind=kind, family_id=f, n_available=len(arr)))
    S = pd.DataFrame(rows)
    reps = pq[["family_id", "rep_rank", "gene", "species"]]
    S = S.merge(reps, on="family_id", how="left")
    # sequences
    need = set(S["gene"])
    aa, ss = {}, {}
    for s in sorted(sp):
        for pid, a in read_fa(os.path.join(ROOT, "data/dense", R, "proteomes", f"{s}.fa")):
            pid = pid.split("|")[1] if pid.startswith(("sp|", "tr|")) else pid
            aa[pid] = a.upper()
        if s in ALLOW_MISSING.get(R, set()):
            continue
        for pid, x in read_fa(os.path.join(D3, f"{s}.3di.fa")):
            ss[pid] = x.upper()
    log(f"{R}: {len(aa):,} proteins, {len(ss):,} with 3Di")
    S["has_3di"] = S["gene"].map(lambda g: g in ss and len(ss[g]) == len(aa.get(g, "")))
    # panel DB (all proteins with 3Di of matching length)
    pairs = [(p, aa[p], ss[p]) for p in aa if p in ss and len(ss[p]) == len(aa[p])]
    ids = build_db(pairs, os.path.join(wd, "db", "panel"), logf)
    json.dump(dict(n_proteins=len(aa), n_in_db=len(ids), species_without_3di=sorted(ALLOW_MISSING.get(R, set()))),
              open(os.path.join(wd, "db", "info.json"), "w"))
    json.dump(dict(species_without_3di=sorted(ALLOW_MISSING.get(R, set()))), open(os.path.join(od, "sample_info.json"), "w"))
    # query DB: real reps (Q|gene) + shuffled reps of the 'old' sample (S|gene)
    qp, seen = [], set()
    for r in S[S["has_3di"]].itertuples():
        if ("Q", r.gene) not in seen:
            seen.add(("Q", r.gene))
            qp.append((f"Q|{r.gene}", aa[r.gene], ss[r.gene]))
    srng = np.random.default_rng(SEED)
    for r in S[S["has_3di"] & (S["kind"] == "old")].itertuples():
        if ("S", r.gene) in seen:
            continue
        seen.add(("S", r.gene))
        perm = srng.permutation(len(aa[r.gene]))
        a, s3 = np.array(list(aa[r.gene])), np.array(list(ss[r.gene]))
        qp.append((f"S|{r.gene}", "".join(a[perm]), "".join(s3[perm])))
    build_db(qp, os.path.join(wd, "qdb", "q"), logf)
    shuf = S[S["kind"] == "old"].copy()
    shuf["kind"] = "shuffled"
    S = pd.concat([S, shuf], ignore_index=True)
    S["query"] = np.where(S["kind"] == "shuffled", "S|", "Q|") + S["gene"]
    S.to_csv(os.path.join(od, "sample.tsv"), sep="\t", index=False)
    log(f"{R}: sample {S.groupby('kind')['family_id'].nunique().to_dict()} families; {len(qp):,} query sequences")


def subdb(src, keys, dst, logf):
    kf = dst + ".keys"
    with open(kf, "w") as fh:
        fh.write("\n".join(map(str, keys)) + "\n")
    for suf in ("", "_ss", "_h"):
        run([FS, "createsubdb", kf, src + suf, dst + suf, "--subdb-mode", "1", "-v", "1"], logf)


def cmd_search(R):
    od, wd = os.path.join(OUT, R), os.path.join(WK, R)
    dst = os.path.join(od, "hits.tsv.gz")
    if os.path.exists(dst):
        sys.exit(f"refusing to overwrite {dst}")
    logf = os.path.join(wd, "commands.log")
    S = pd.read_csv(os.path.join(od, "sample.tsv"), sep="\t")
    S = S[S["has_3di"]]
    tl = pd.read_csv(os.path.join(wd, "db", "panel.lookup"), sep="\t", header=None, usecols=[0, 1], names=["k", "id"])
    tl["species"] = tl["id"].str.split("_").str[0]
    ql = pd.read_csv(os.path.join(wd, "qdb", "q.lookup"), sep="\t", header=None, usecols=[0, 1], names=["k", "id"])
    qkey = dict(zip(ql["id"], ql["k"]))
    cl = clades(R)
    sp = set(panel_species(R))
    # species code = proteome file stem; protein ids start with '<code>_' (checked)
    bad = set(tl["species"]) - sp
    if bad:
        raise SystemExit(f"protein ids not prefixed by species code: {sorted(bad)[:5]}")
    parts = []
    t_all = time.time()
    for br, sub in S.groupby("branch_child", sort=True):
        bd = os.path.join(wd, "search", br)
        os.makedirs(bd, exist_ok=True)
        out_s = sp - cl[br]
        tk = tl.loc[tl["species"].isin(out_s), "k"].tolist()
        qk = sorted({qkey[q] for q in sub["query"]})
        subdb(os.path.join(wd, "db", "panel"), tk, os.path.join(bd, "t"), logf)
        subdb(os.path.join(wd, "qdb", "q"), qk, os.path.join(bd, "q"), logf)
        t0 = time.time()
        run([FS, "search", os.path.join(bd, "q"), os.path.join(bd, "t"), os.path.join(bd, "aln"), os.path.join(bd, "tmp"),
             "--alignment-type", "2", "-e", EVAL, "--max-seqs", 1000, "--threads", THREADS, "-v", "1"], logf)
        run([FS, "convertalis", os.path.join(bd, "q"), os.path.join(bd, "t"), os.path.join(bd, "aln"),
             os.path.join(bd, "aln.m8"), "--format-output", FMT, "--threads", THREADS, "-v", "1"], logf)
        h = pd.read_csv(os.path.join(bd, "aln.m8"), sep="\t", header=None, names=FMT.split(","))
        h.insert(0, "branch_child", br)
        h["search_s"] = round(time.time() - t0, 1)
        parts.append(h)
        log(f"{R} {br}: {len(qk)} queries vs {len(tk):,} targets ({len(out_s)} species) -> {len(h):,} hits "
            f"in {time.time() - t0:.0f}s")
        for f in os.listdir(bd):
            if f.startswith(("aln", "tmp")):
                p = os.path.join(bd, f)
                (subprocess.run(["rm", "-rf", p]) if os.path.isdir(p) else os.remove(p))
    H = pd.concat(parts, ignore_index=True)
    H.to_csv(dst, sep="\t", index=False, compression="gzip")
    json.dump(dict(wall_s=round(time.time() - t_all, 1), n_hits=len(H)), open(os.path.join(od, "search_info.json"), "w"))


def cmd_classify(R):
    od = os.path.join(OUT, R)
    dst = os.path.join(od, "families.tsv")
    if os.path.exists(dst):
        sys.exit(f"refusing to overwrite {dst}")
    S = pd.read_csv(os.path.join(od, "sample.tsv"), sep="\t")
    H = pd.read_csv(os.path.join(od, "hits.tsv.gz"), sep="\t")
    H["qcov"] = (H["qend"] - H["qstart"] + 1) / H["qlen"]
    good = H[(H["evalue"] < EVAL) & (H["qcov"] >= QCOV)]
    best = H.sort_values("evalue").groupby(["branch_child", "query"]).head(1).set_index(["branch_child", "query"])
    hit_q = set(zip(good["branch_child"], good["query"]))
    S["struct_hit"] = [(b, q) in hit_q for b, q in zip(S["branch_child"], S["query"])]
    S["any_hit"] = [(b, q) in best.index for b, q in zip(S["branch_child"], S["query"])]
    F = (S[S["has_3di"]].groupby(["branch_child", "kind", "family_id"])
         .agg(n_reps=("query", "size"), struct_old=("struct_hit", "any"), any_hit=("any_hit", "any")).reset_index())
    F["struct_class"] = np.where(F["struct_old"], "structurally_old", "structurally_new")
    F.insert(0, "panel", R)
    F.to_csv(dst, sep="\t", index=False)
    B = F.groupby(["branch_child", "kind"]).agg(n=("family_id", "size"), n_struct_old=("struct_old", "sum")).reset_index()
    B["frac_struct_new"] = (1 - B["n_struct_old"] / B["n"]).round(4)
    B.insert(0, "panel", R)
    B.to_csv(os.path.join(od, "branches.tsv"), sep="\t", index=False)
    import hashlib
    md5 = {f: hashlib.md5(open(os.path.join(od, f), "rb").read()).hexdigest()
           for f in sorted(os.listdir(od)) if f != "DONE.json"}
    json.dump(dict(status="done", finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"), md5=md5),
              open(os.path.join(od, "DONE.json"), "w"), indent=1)
    print(B.to_string())


if __name__ == "__main__":
    cmd, R = sys.argv[1], sys.argv[2]
    {"prep": cmd_prep, "search": cmd_search, "classify": cmd_classify}[cmd](R)
