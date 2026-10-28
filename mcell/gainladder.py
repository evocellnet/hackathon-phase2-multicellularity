"""Gain-evidence ladder for HOG-level gains (DESIGN §3.2).

For every family that :mod:`mcell.hogs` places as a *gain* on an internal branch
p -> c (root stem and leaf branches excluded; stems and control branches are
treated identically) the ladder asks, in order:

1. **hog_gain**       the raw Ancestra HOG-level gain (upper bound).
2. **split_repaired** the gain is a split of an older HOG: its best-scoring
                      outside-clade HOG X at the parent level p is *lost* on the
                      same branch p -> c, and the gain is X's best-scoring gain
                      on that branch (reciprocal best hit at the HOG level). The
                      gain and the loss cancel (gain re-merged into X).
3. **seq_novel**      no representative of the gained family has a DIAMOND hit
                      (``-k 0 -e 1e-5 --query-cover 30``) in any proteome outside
                      the clade below c (all proteomes of both domains are
                      targets; in-clade hits are filtered afterwards).
4. **arch_novel**     (secondary, optional) the family's majority ordered Pfam
                      architecture occurs in no protein outside the clade.

``split_repaired`` and ``seq_novel`` are mutually exclusive (a repair needs an
outside-clade hit); gains that are neither are "homologous, not repaired".

Representatives (≤ 3 per gain, deterministic): the longest member; then the
longest member of one random species (``numpy`` PCG64 seeded with
``[seed, crc32(family_id)]``) in each *other* child clade of c, so that the
representatives span the clade's root split; then further random species; and
finally, when the family has fewer species than representatives, the next
longest genes.

Definitions used for split repair
---------------------------------
* score(G, X) = sum over representatives r of G of the best bitscore of r to any
  member of X (X = ancestral gene at level p; outside-clade hits only, top
  ``top_n`` per representative). Summing per-representative maxima (rather than
  all hits) keeps large families from winning by size.
* X*(G) = argmax_X score(G, X) over all HOGs at p hit by G (ties: smallest id).
* Accept if X*(G) is lost on p -> c **and** G = argmax_{G' gained on p -> c}
  score(G', X*(G)) (the reciprocal condition; each X absorbs ≤ 1 gain).
* Corrected counts on p -> c: gains - split_repaired and losses - split_repaired.

The DIAMOND output is reduced in chunks (:func:`reduce_hits`), so the full hit
table never has to fit in memory.
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import time
import zlib
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from .hogs import HogResult, SpeciesTree

__all__ = [
    "DIAMOND_COLUMNS", "read_fasta", "fasta_lengths", "build_target_table", "make_diamond_db",
    "internal_gains", "choose_representatives", "write_queries", "run_diamond",
    "TargetSpace", "reduce_hits", "classify_gains", "branch_ladder", "pfam_architectures",
    "architecture_novelty",
]

DIAMOND_COLUMNS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovhsp"]


# ----------------------------------------------------------------------------
# FASTA / target database
# ----------------------------------------------------------------------------

def _open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path: str, keep: set | None = None) -> Iterable[tuple[str, str]]:
    """Yield (id, sequence); ``id`` is the header up to the first whitespace.
    With ``keep``, only those ids are assembled (fast skip of the others)."""
    h, buf, take = None, [], False
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None and take:
                    yield h, "".join(buf)
                h = line[1:].split(None, 1)[0]
                take = keep is None or h in keep
                buf = []
            elif take:
                buf.append(line.strip())
    if h is not None and take:
        yield h, "".join(buf)


def fasta_lengths(paths: Iterable[str]) -> pd.Series:
    """Protein id -> length for all records of ``paths``."""
    ids, lens = [], []
    for p in paths:
        h, n = None, 0
        with _open_text(p) as fh:
            for line in fh:
                if line.startswith(">"):
                    if h is not None:
                        ids.append(h)
                        lens.append(n)
                    h, n = line[1:].split(None, 1)[0], 0
                else:
                    n += len(line.strip())
        if h is not None:
            ids.append(h)
            lens.append(n)
    s = pd.Series(np.asarray(lens, dtype=np.int32), index=pd.Index(ids, name="protein"), name="length")
    if not s.index.is_unique:
        raise ValueError("protein ids are not unique across the FASTA files")
    return s


def build_target_table(proteome_dirs: dict[str, str], suffix: str = ".fa") -> pd.DataFrame:
    """One row per target protein: (seqid, species, domain, length). Species =
    FASTA file stem (the project's species code), domain = key of ``proteome_dirs``."""
    rows = []
    for dom, d in proteome_dirs.items():
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(suffix):
                continue
            sp = fn[: -len(suffix)]
            L = fasta_lengths([os.path.join(d, fn)])
            rows.append(pd.DataFrame({"seqid": L.index.to_numpy(), "species": sp, "domain": dom,
                                      "length": L.to_numpy()}))
    t = pd.concat(rows, ignore_index=True)
    if not t["seqid"].is_unique:
        dup = t.loc[t["seqid"].duplicated(), "seqid"].head(5).tolist()
        raise ValueError(f"target ids are not unique, e.g. {dup}")
    t["species"] = t["species"].astype("category")
    t["domain"] = t["domain"].astype("category")
    return t


def make_diamond_db(proteome_dirs: dict[str, str], out_prefix: str, threads: int = 4,
                    diamond: str = "diamond", suffix: str = ".fa") -> dict:
    """Concatenate all proteomes (domains in the given order, files sorted) into
    ``<out_prefix>.faa``, build ``<out_prefix>.dmnd`` and ``<out_prefix>.targets.parquet``.
    Existing outputs are reused."""
    faa, dmnd, tab = out_prefix + ".faa", out_prefix + ".dmnd", out_prefix + ".targets.parquet"
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    if not os.path.exists(faa):
        with open(faa + ".tmp", "wb") as out:
            for d in proteome_dirs.values():
                for fn in sorted(os.listdir(d)):
                    if fn.endswith(suffix):
                        with open(os.path.join(d, fn), "rb") as fh:
                            shutil.copyfileobj(fh, out)
        os.replace(faa + ".tmp", faa)
    if not os.path.exists(dmnd):
        subprocess.run([diamond, "makedb", "--in", faa, "-d", out_prefix, "--threads", str(threads)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if not os.path.exists(tab):
        build_target_table(proteome_dirs, suffix).to_parquet(tab, index=False)
    return dict(faa=faa, dmnd=dmnd, targets=tab)


# ----------------------------------------------------------------------------
# gains and representatives
# ----------------------------------------------------------------------------

def internal_gains(res: HogResult) -> pd.DataFrame:
    """Gains (one per family) on internal branches, excluding the root stem."""
    ev = res.events
    g = ev[(ev["event"] == "gain") & (ev["branch_parent"] != "") & ~ev["child_is_leaf"]]
    g = g[["branch_parent", "branch_child", "family_id", "hog_id", "n_extant_genes", "n_species",
           "n_child_clades_present", "n_child_clades_total"]].rename(columns={"n_extant_genes": "n_genes"})
    if not g["family_id"].is_unique:
        raise AssertionError("more than one gain per family (parser invariant broken)")
    return g.reset_index(drop=True)


def _rng(seed: int, key: str) -> np.random.Generator:
    return np.random.default_rng([seed, zlib.crc32(key.encode())])


def choose_representatives(res: HogResult, gains: pd.DataFrame, lengths: pd.Series,
                           n_rep: int = 3, seed: int = 11) -> pd.DataFrame:
    """Up to ``n_rep`` representative proteins per gained family (see module doc).

    Returns (family_id, rep_rank, gene, species, length, child_clade)."""
    T = res.tree
    fams = set(gains["family_id"])
    G = res.genes
    sub = G[G["family_id"].isin(fams)][["gene", "species", "family_id"]].copy()
    sub["length"] = lengths.reindex(sub["gene"]).to_numpy()
    if sub["length"].isna().any():
        miss = sub.loc[sub["length"].isna(), "gene"].head(5).tolist()
        raise KeyError(f"no protein length for {int(sub['length'].isna().sum())} genes, e.g. {miss}")
    sub["sp_idx"] = [T.index[s] for s in sub["species"].astype(str)]
    # longest first, ties by gene id -> deterministic
    sub = sub.sort_values(["family_id", "length", "gene"], ascending=[True, False, True], kind="mergesort")
    child_of = dict(zip(gains["family_id"], gains["branch_child"]))
    out = []
    fam_arr = sub["family_id"].to_numpy()
    starts = np.flatnonzero(np.r_[True, fam_arr[1:] != fam_arr[:-1]])
    ends = np.r_[starts[1:], len(fam_arr)]
    gene_arr, sp_arr = sub["gene"].to_numpy(), sub["sp_idx"].to_numpy()
    len_arr = sub["length"].to_numpy()
    for a, b in zip(starts, ends):
        fam = fam_arr[a]
        c = T.index[child_of[fam]]
        clade_of = {}
        for s in set(sp_arr[a:b].tolist()):
            clade_of[s] = T.child_towards(c, s) if s != c else c
        # first (longest) gene per species, in length order
        first_gene: dict[int, int] = {}
        for k in range(a, b):
            first_gene.setdefault(sp_arr[k], k)
        rng = _rng(seed, fam)
        chosen = [a]
        used_sp = {sp_arr[a]}
        used_cl = {clade_of[sp_arr[a]]}
        # one random species from each other child clade (clades in random order)
        clades = sorted({clade_of[s] for s in first_gene} - used_cl)
        clades = [clades[i] for i in rng.permutation(len(clades))]
        for cl in clades:
            if len(chosen) >= n_rep:
                break
            sps = sorted(s for s in first_gene if clade_of[s] == cl)
            s = sps[int(rng.integers(len(sps)))]
            chosen.append(first_gene[s])
            used_sp.add(s)
        # further random species
        rest = sorted(s for s in first_gene if s not in used_sp)
        for i in rng.permutation(len(rest)):
            if len(chosen) >= n_rep:
                break
            chosen.append(first_gene[rest[i]])
            used_sp.add(rest[i])
        # fewer species than representatives: next longest genes
        k = a
        while len(chosen) < n_rep and k < b:
            if k not in chosen:
                chosen.append(k)
            k += 1
        for r, k in enumerate(chosen):
            out.append((fam, r, gene_arr[k], T.names[sp_arr[k]], int(len_arr[k]),
                        T.names[clade_of[sp_arr[k]]]))
    return pd.DataFrame(out, columns=["family_id", "rep_rank", "gene", "species", "length", "child_clade"])


def write_queries(reps: pd.DataFrame, fasta_paths: Iterable[str], out_faa: str) -> int:
    """Write the representative sequences (FASTA ids = gene ids)."""
    want = set(reps["gene"])
    n = 0
    with open(out_faa + ".tmp", "w") as out:
        for p in fasta_paths:
            for h, s in read_fasta(p, keep=want):
                out.write(f">{h}\n{s}\n")
                want.discard(h)
                n += 1
    if want:
        raise KeyError(f"{len(want)} representative sequences not found, e.g. {sorted(want)[:5]}")
    os.replace(out_faa + ".tmp", out_faa)
    return n


def run_diamond(query: str, db: str, out: str, threads: int = 4, sensitivity: str = "ultra-sensitive",
                evalue: float = 1e-5, query_cover: float = 30, max_hsps: int = 1,
                block_size: float | None = 4, index_chunks: int | None = 1, compress: bool = True,
                diamond: str = "diamond", extra: Iterable[str] = (), log: str | None = None) -> dict:
    """``diamond blastp -k 0`` with the ladder's settings; returns timing info."""
    cmd = [diamond, "blastp", "-q", query, "-d", db, "-o", out, f"--{sensitivity}", "-k", "0",
           "-e", str(evalue), "--query-cover", str(query_cover), "--threads", str(threads),
           "-f", "6", *DIAMOND_COLUMNS]
    if max_hsps:
        cmd += ["--max-hsps", str(max_hsps)]
    if block_size:
        cmd += ["-b", str(block_size)]
    if index_chunks:
        cmd += ["-c", str(index_chunks)]
    if compress:
        cmd += ["--compress", "1"]
    cmd += list(extra)
    t0 = time.time()
    with open(log or os.devnull, "w") as lh:
        subprocess.run(cmd, check=True, stdout=lh, stderr=subprocess.STDOUT)
    if compress and not out.endswith(".gz") and os.path.exists(out + ".gz"):
        os.replace(out + ".gz", out)      # DIAMOND appends .gz to compressed output names
    return dict(cmd=" ".join(cmd), wall_s=round(time.time() - t0, 1))


# ----------------------------------------------------------------------------
# hit reduction (clade filter)
# ----------------------------------------------------------------------------

@dataclass
class TargetSpace:
    """Target proteins with, for each, the pre-order index of its species in the
    run's species tree (-1 if the species is not in the run) and the row of the
    gene in ``res.genes`` (-1 if not a gene of the run)."""
    seqid: pd.Index
    species: np.ndarray        # species code (object)
    domain: np.ndarray         # domain label (object)
    sp_pre: np.ndarray         # int64, -1 = not in the run tree
    gene_row: np.ndarray       # int64, -1 = not a gene of the run
    transcriptome: np.ndarray  # bool

    @classmethod
    def build(cls, targets: pd.DataFrame, res: HogResult,
              transcriptome_species: Iterable[str] = ()) -> "TargetSpace":
        T = res.tree
        sp = targets["species"].astype(str).to_numpy(dtype=object)
        pre = np.array([T.index.get(s, -1) if s in T.index and T.is_leaf[T.index[s]] else -1
                        for s in pd.unique(sp)], dtype=np.int64)
        lut = dict(zip(pd.unique(sp), pre))
        sp_pre = np.fromiter((lut[s] for s in sp), dtype=np.int64, count=len(sp))
        idx = pd.Index(targets["seqid"].to_numpy(dtype=object))
        gene_row = idx.get_indexer(res.genes["gene"].to_numpy(dtype=object))
        if (gene_row < 0).any():
            miss = res.genes["gene"].to_numpy()[gene_row < 0][:5].tolist()
            raise KeyError(f"{int((gene_row < 0).sum())} run genes are not in the target DB, e.g. {miss}")
        grow = np.full(len(idx), -1, dtype=np.int64)
        grow[gene_row] = np.arange(len(gene_row))
        tx = set(transcriptome_species)
        return cls(seqid=idx, species=sp, domain=targets["domain"].astype(str).to_numpy(dtype=object),
                   sp_pre=sp_pre, gene_row=grow,
                   transcriptome=np.fromiter((s in tx for s in sp), dtype=bool, count=len(sp)))


def reduce_hits(hits_path: str, reps: pd.DataFrame, gains: pd.DataFrame, tree: SpeciesTree,
                ts: TargetSpace, top_n: int = 50, chunksize: int = 5_000_000,
                tx_pident: float = 90.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Clade-filter a DIAMOND tabular file.

    Returns
    -------
    per_query : one row per representative: n_hits, n_in_clade, n_outside,
        n_outside_robust (outside hits not from a transcriptome at ≥ ``tx_pident``%),
        n_outside_species, n_outside_other_domain, best outside hit (seqid, species,
        domain, pident, bitscore, evalue).
    top : the ``top_n`` best outside-clade hits per representative
        (query, target row, bitscore, pident).
    """
    q_index = pd.Index(reps["gene"].to_numpy(dtype=object))
    child = dict(zip(gains["family_id"], gains["branch_child"]))
    q_c = np.array([tree.index[child[f]] for f in reps["family_id"]], dtype=np.int64)
    q_end = np.array([tree.end[c] for c in q_c], dtype=np.int64)
    run_domain = None
    nq = len(q_index)
    n_hits = np.zeros(nq, np.int64)
    n_out = np.zeros(nq, np.int64)
    n_rob = np.zeros(nq, np.int64)
    n_od = np.zeros(nq, np.int64)
    best_bs = np.full(nq, -1.0)
    best_t = np.full(nq, -1, np.int64)
    best_pid = np.full(nq, np.nan)
    best_ev = np.full(nq, np.nan)
    sp_pairs = []
    tops = []
    # the run's domain = domain of the run's genes
    in_run = ts.gene_row >= 0
    doms = pd.unique(ts.domain[in_run])
    run_domain = doms[0] if len(doms) == 1 else None
    if os.path.getsize(hits_path) > 0:
        reader = pd.read_csv(hits_path, sep="\t", header=None, names=DIAMOND_COLUMNS, chunksize=chunksize,
                             dtype={"qseqid": object, "sseqid": object, "pident": np.float32,
                                    "length": np.int32, "evalue": np.float64, "bitscore": np.float32,
                                    "qcovhsp": np.float32})
        for ch in reader:
            q = q_index.get_indexer(ch["qseqid"].to_numpy(dtype=object))
            t = ts.seqid.get_indexer(ch["sseqid"].to_numpy(dtype=object))
            if (q < 0).any():
                raise KeyError("hit table contains queries that are not representatives")
            if (t < 0).any():
                raise KeyError("hit table contains targets missing from the target table")
            np.add.at(n_hits, q, 1)
            pre = ts.sp_pre[t]
            inside = (pre >= 0) & (pre >= q_c[q]) & (pre <= q_end[q])
            o = ~inside
            qo, to = q[o], t[o]
            bs = ch["bitscore"].to_numpy()[o]
            pid = ch["pident"].to_numpy()[o]
            ev = ch["evalue"].to_numpy()[o]
            np.add.at(n_out, qo, 1)
            rob = ~(ts.transcriptome[to] & (pid >= tx_pident))
            np.add.at(n_rob, qo[rob], 1)
            if run_domain is not None:
                np.add.at(n_od, qo[ts.domain[to] != run_domain], 1)
            if len(qo):
                df = pd.DataFrame({"q": qo, "t": to, "bitscore": bs, "pident": pid, "evalue": ev})
                df = df.sort_values(["q", "bitscore", "t"], ascending=[True, False, True], kind="mergesort")
                first = df.drop_duplicates("q")
                fq = first["q"].to_numpy()
                better = first["bitscore"].to_numpy() > best_bs[fq]
                fq = fq[better]
                best_bs[fq] = first["bitscore"].to_numpy()[better]
                best_t[fq] = first["t"].to_numpy()[better]
                best_pid[fq] = first["pident"].to_numpy()[better]
                best_ev[fq] = first["evalue"].to_numpy()[better]
                tops.append(df.groupby("q", sort=False).head(top_n)[["q", "t", "bitscore", "pident"]])
                sp_pairs.append(pd.DataFrame({"q": qo, "sp": ts.species[to]}).drop_duplicates())
    if tops:
        top = pd.concat(tops, ignore_index=True)
        top = top.sort_values(["q", "bitscore", "t"], ascending=[True, False, True], kind="mergesort")
        top = top.groupby("q", sort=False).head(top_n).reset_index(drop=True)
        spp = pd.concat(sp_pairs, ignore_index=True).drop_duplicates()
        n_sp = spp.groupby("q").size().reindex(range(nq), fill_value=0).to_numpy()
    else:
        top = pd.DataFrame({"q": pd.Series([], dtype=np.int64), "t": pd.Series([], dtype=np.int64),
                            "bitscore": pd.Series([], dtype=np.float32),
                            "pident": pd.Series([], dtype=np.float32)})
        n_sp = np.zeros(nq, np.int64)
    has = best_t >= 0
    pq = reps[["family_id", "rep_rank", "gene", "species"]].copy()
    pq["n_hits"] = n_hits
    pq["n_in_clade"] = n_hits - n_out
    pq["n_outside"] = n_out
    pq["n_outside_robust"] = n_rob
    pq["n_outside_species"] = n_sp
    pq["n_outside_other_domain"] = n_od if run_domain is not None else -1
    pq["best_out_target"] = np.where(has, ts.seqid.to_numpy()[np.where(has, best_t, 0)], None)
    pq["best_out_species"] = np.where(has, ts.species[np.where(has, best_t, 0)], None)
    pq["best_out_domain"] = np.where(has, ts.domain[np.where(has, best_t, 0)], None)
    pq["best_out_pident"] = best_pid
    pq["best_out_bitscore"] = np.where(has, best_bs, np.nan)
    pq["best_out_evalue"] = best_ev
    return pq.reset_index(drop=True), top


# ----------------------------------------------------------------------------
# classification
# ----------------------------------------------------------------------------

def _hog_at_level(res: HogResult, gene_rows: np.ndarray, level_idx: int) -> np.ndarray:
    """Row in ``res.hogs`` of the ancestral gene at ``level_idx`` for each gene
    row (-1 where there is none)."""
    hl = res.hogs["level_idx"].to_numpy()
    hp = res.hogs["parent_idx"].to_numpy()
    cur = res.genes["parent_idx"].to_numpy()[gene_rows].copy()
    out = np.full(len(cur), -1, dtype=np.int64)
    for _ in range(max(res.tree.depth) + 1):
        ok = cur >= 0
        if not ok.any():
            break
        safe = np.where(ok, cur, 0)
        hit = ok & (hl[safe] == level_idx)
        out[hit] = cur[hit]
        # levels only go up; stop following once above the target level
        cur = np.where(ok & ~hit, hp[safe], -1)
    return out


def classify_gains(res: HogResult, gains: pd.DataFrame, reps: pd.DataFrame, per_query: pd.DataFrame,
                   top: pd.DataFrame, ts: TargetSpace) -> pd.DataFrame:
    """Novelty and split-repair labels per gain (one row per gained family)."""
    T = res.tree
    agg = per_query.groupby("family_id").agg(
        n_reps=("gene", "size"), n_hits=("n_hits", "sum"), n_outside=("n_outside", "sum"),
        n_outside_robust=("n_outside_robust", "sum"),
        n_outside_species_max=("n_outside_species", "max"),
        n_outside_other_domain=("n_outside_other_domain", "sum"),
        best_out_bitscore=("best_out_bitscore", "max"))
    best = (per_query.dropna(subset=["best_out_bitscore"])
            .sort_values(["family_id", "best_out_bitscore", "gene"], ascending=[True, False, True])
            .drop_duplicates("family_id").set_index("family_id")
            [["best_out_target", "best_out_species", "best_out_domain", "best_out_pident"]])
    g = gains.set_index("family_id").join(agg).join(best)
    g["n_reps"] = g["n_reps"].fillna(0).astype(int)
    for c in ("n_hits", "n_outside", "n_outside_robust", "n_outside_other_domain"):
        g[c] = g[c].fillna(0).astype(np.int64)
    g["seq_novel"] = g["n_outside"] == 0
    g["seq_novel_tx_sens"] = g["n_outside_robust"] == 0
    g["outside_other_domain_only"] = (g["n_outside"] > 0) & (g["n_outside_other_domain"] == g["n_outside"])

    # ---- split repair -------------------------------------------------------
    g["repair_hog"] = None
    g["repair_score"] = np.nan
    g["best_parent_hog"] = None
    g["best_parent_hog_lost"] = False
    g["split_repaired"] = False
    if len(top):
        tp = top.copy()
        tp["family_id"] = reps["family_id"].to_numpy()[tp["q"].to_numpy()]
        tp["gene_row"] = ts.gene_row[tp["t"].to_numpy()]
        tp = tp[tp["gene_row"] >= 0]
        parent_of = dict(zip(gains["family_id"], gains["branch_parent"]))
        tp["p_idx"] = [T.index[parent_of[f]] for f in tp["family_id"]]
        tp["hog_row"] = -1
        for p, idx in tp.groupby("p_idx").groups.items():
            rows = tp.loc[idx, "gene_row"].to_numpy()
            tp.loc[idx, "hog_row"] = _hog_at_level(res, rows, int(p))
        tp = tp[tp["hog_row"] >= 0]
        # score(G, X) = sum over reps of best bitscore to X
        s = (tp.groupby(["family_id", "q", "hog_row"], sort=False)["bitscore"].max()
             .groupby(level=["family_id", "hog_row"]).sum().rename("score").reset_index())
        s["hog_id"] = res.hogs["hog_id"].to_numpy()[s["hog_row"].to_numpy()]
        s["branch_child"] = gains.set_index("family_id").loc[s["family_id"], "branch_child"].to_numpy()
        # X*(G)
        xs = s.sort_values(["family_id", "score", "hog_id"], ascending=[True, False, True]).drop_duplicates("family_id")
        # losses on the gain's branch
        ev = res.events
        lost = set(zip(ev.loc[ev["event"] == "loss", "branch_child"], ev.loc[ev["event"] == "loss", "hog_id"]))
        xs["lost"] = [(c, h) in lost for c, h in zip(xs["branch_child"], xs["hog_id"])]
        # G*(X) among gains on the same branch
        gs = s.sort_values(["branch_child", "hog_id", "score", "family_id"],
                           ascending=[True, True, False, True]).drop_duplicates(["branch_child", "hog_id"])
        gbest = dict(zip(zip(gs["branch_child"], gs["hog_id"]), gs["family_id"]))
        xs["reciprocal"] = [gbest.get((c, h)) == f for c, h, f in zip(xs["branch_child"], xs["hog_id"],
                                                                         xs["family_id"])]
        xs = xs.set_index("family_id")
        g.loc[xs.index, "best_parent_hog"] = xs["hog_id"]
        g.loc[xs.index, "best_parent_hog_lost"] = xs["lost"]
        ok = xs[xs["lost"] & xs["reciprocal"]]
        ok = ok[~g.loc[ok.index, "seq_novel"].to_numpy()]
        g.loc[ok.index, "split_repaired"] = True
        g.loc[ok.index, "repair_hog"] = ok["hog_id"]
        g.loc[ok.index, "repair_score"] = ok["score"]
    # why a homologous gain was not repaired
    st = np.full(len(g), "no_parent_level_hit", dtype=object)
    has = g["best_parent_hog"].notna().to_numpy()
    lostv = g["best_parent_hog_lost"].astype(bool).to_numpy()
    st[has & ~lostv] = "best_parent_hog_not_lost"
    st[has & lostv] = "not_reciprocal"
    st[g["split_repaired"].to_numpy()] = "repaired"
    st[g["seq_novel"].to_numpy()] = "novel"
    g["repair_status"] = st
    lab = np.where(g["seq_novel"], "seq_novel", np.where(g["split_repaired"], "split_repaired", "homologous"))
    g["ladder"] = lab
    return g.reset_index()


def branch_ladder(res: HogResult, cls: pd.DataFrame) -> pd.DataFrame:
    """Per internal branch: HOG-level gains, split-repaired, corrected gains and
    losses, homologous (not repaired), sequence-novel (+ architecture-novel if present)."""
    b = res.branches
    b = b[(b["branch_parent"] != "") & ~b["child_is_leaf"]][
        ["branch_parent", "branch_child", "depth_child", "n_parent", "n_child", "gains", "losses",
         "dups", "dups_at_origin"]].rename(columns={"gains": "hog_gains"}).set_index("branch_child")
    grp = cls.groupby("branch_child")
    out = b.join(pd.DataFrame({
        "split_repaired": grp["split_repaired"].sum(),
        "seq_novel": grp["seq_novel"].sum(),
        "seq_novel_tx_sens": grp["seq_novel_tx_sens"].sum(),
        "homologous_unrepaired": grp["ladder"].apply(lambda x: int((x == "homologous").sum())),
        "outside_other_domain_only": grp["outside_other_domain_only"].sum()}))
    for c in ("split_repaired", "seq_novel", "seq_novel_tx_sens", "homologous_unrepaired",
              "outside_other_domain_only"):
        out[c] = out[c].fillna(0).astype(np.int64)
    out["gains_after_repair"] = out["hog_gains"] - out["split_repaired"]
    out["losses_after_repair"] = out["losses"] - out["split_repaired"]
    if "arch_novel" in cls.columns:
        an = cls[cls["seq_novel"]].groupby("branch_child")["arch_novel"]
        out["seq_novel_with_arch"] = an.apply(lambda x: int(x.notna().sum())).reindex(out.index).fillna(0).astype(int)
        out["arch_novel"] = an.apply(lambda x: int((x == True).sum())).reindex(out.index).fillna(0).astype(int)  # noqa: E712
    out["frac_seq_novel"] = np.where(out["hog_gains"] > 0, out["seq_novel"] / out["hog_gains"].clip(lower=1), np.nan)
    return out.reset_index()[["branch_parent", "branch_child"] + [c for c in out.columns if c != "branch_parent"]]


# ----------------------------------------------------------------------------
# Pfam architectures (optional, re-runnable)
# ----------------------------------------------------------------------------

def _resolve_overlaps(df: pd.DataFrame, min_frac: float = 0.5) -> pd.DataFrame:
    """Greedy by domain score: drop a hit overlapping a kept hit of the same
    protein by > ``min_frac`` of the shorter one (clan competition proxy)."""
    df = df.sort_values(["protein", "dom_score"], ascending=[True, False], kind="mergesort")
    keep = np.ones(len(df), bool)
    prot = df["protein"].to_numpy()
    a, b = df["ali_from"].to_numpy(), df["ali_to"].to_numpy()
    starts = np.flatnonzero(np.r_[True, prot[1:] != prot[:-1]])
    ends = np.r_[starts[1:], len(prot)]
    for s, e in zip(starts, ends):
        if e - s == 1:
            continue
        kept = []
        for i in range(s, e):
            ok = True
            for j in kept:
                ov = min(b[i], b[j]) - max(a[i], a[j]) + 1
                if ov > min_frac * min(b[i] - a[i] + 1, b[j] - a[j] + 1):
                    ok = False
                    break
            if ok:
                kept.append(i)
            else:
                keep[i] = False
    return df[keep]


def pfam_architectures(pfam_files: dict[str, str]) -> pd.DataFrame:
    """(protein, species, arch) from per-species Pfam tables (``run_pfam.py``
    format). Architecture = Pfam accessions ordered by ``ali_from`` after overlap
    resolution, with tandem repeats of the same domain collapsed (``A-B-A``)."""
    rows = []
    for sp, path in pfam_files.items():
        d = pd.read_csv(path, sep="\t", usecols=["protein", "pfam_acc", "ali_from", "ali_to", "dom_score"])
        if d.empty:
            continue
        d["pfam_acc"] = d["pfam_acc"].str.split(".").str[0]
        d = _resolve_overlaps(d).sort_values(["protein", "ali_from"], kind="mergesort")
        prev = d.groupby("protein")["pfam_acc"].shift()
        d = d[d["pfam_acc"] != prev]
        a = d.groupby("protein", sort=False)["pfam_acc"].agg("-".join).rename("arch").reset_index()
        a["species"] = sp
        rows.append(a)
    if not rows:
        return pd.DataFrame(columns=["protein", "arch", "species"])
    return pd.concat(rows, ignore_index=True)


def architecture_novelty(res: HogResult, cls: pd.DataFrame, arch: pd.DataFrame,
                         pfam_species: Iterable[str], all_species: Iterable[str]) -> pd.DataFrame:
    """Add ``arch``, ``arch_novel`` (None if the family has no Pfam domain or no
    member species has Pfam output) and ``pfam_outside_coverage`` (fraction of
    outside-clade species with Pfam output) to ``cls``."""
    T = res.tree
    pf = set(pfam_species)
    allsp = list(all_species)
    a = arch.set_index("protein")["arch"]
    G = res.genes[res.genes["family_id"].isin(set(cls["family_id"]))][["gene", "family_id"]].copy()
    G["arch"] = a.reindex(G["gene"]).to_numpy()
    maj = (G.dropna(subset=["arch"]).groupby(["family_id", "arch"]).size().rename("n").reset_index()
           .sort_values(["family_id", "n", "arch"], ascending=[True, False, True]).drop_duplicates("family_id")
           .set_index("family_id")["arch"])
    arch_species = arch.groupby("arch")["species"].agg(lambda x: frozenset(x))
    out = cls.copy()
    out["arch"] = maj.reindex(out["family_id"]).to_numpy()
    novel, cov = [], []
    for fam, c, ar in zip(out["family_id"], out["branch_child"], out["arch"]):
        ci = T.index[c]
        clade = {T.names[s] for s in T.leaves if T.is_anc(ci, s)}
        outside = [s for s in allsp if s not in clade]
        cov.append(sum(s in pf for s in outside) / max(1, len(outside)))
        if not isinstance(ar, str):
            novel.append(None)
            continue
        hit = arch_species.get(ar, frozenset())
        novel.append(not any(s not in clade for s in hit))
    out["arch_novel"] = novel
    out["pfam_outside_coverage"] = np.round(cov, 4)
    return out
