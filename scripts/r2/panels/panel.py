#!/usr/bin/env python
"""Round-2 dense-panel analysis, per panel (W2/W3; DESIGN_ROUND2 §6; task brief steps 1-5).

    panel.py parse    R            # HOGs -> per-branch gains/losses/dups (mcell.hogs; analysis cache)
    panel.py phylo    R --threads 4  # near-single-copy root HOGs -> MAFFT/trimAl -> IQ-TREE LG+G4 on the FIXED topology
    panel.py ladder   R --threads 12 # sequence-novelty ladder (mcell.gainladder; round-1 DIAMOND settings)
    panel.py tables   R            # branch table (covariates, crown/control flags, depth classes) + event table
    panel.py prereg_pfam R --cpus 4  # pre-registered Pfam subset for new (non-round-1) panel species (controls)
    panel.py controls R            # re-score round-1 positive controls on the panel

Intermediates: results/r2/panels/<R>/work/; final per-panel tables: results/r2/panels/<R>/.
Every final output is refused if it exists (guard). Definitions are in METHODS_DEFINITIONS.txt.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (ROOT, PANELS, EXTRA, base, OUT, log, kingdom, paths, guard, md5, write_done, panel_codes, r1_codes,  # noqa: E402
                    place_round1)

import ete3  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mcell import gainladder as gl  # noqa: E402
from mcell.hogs import parse_hogs  # noqa: E402
from mcell.phylo import read_tree, transfer_branch_lengths, write_named, assert_same_topology  # noqa: E402

SEED = 11
N_REP = 3
SENS = "very-sensitive"          # round-1 setting (results/gainladder/euk/DONE.json)
ALL_TIPS_PANELS = {"R4", "R6"}   # P1 panels: all tips laddered for the tip-control sensitivity (DEVIATIONS)
EVENT_ROLES_EXCLUDED = None      # every coded event branch is excluded from controls


# ----------------------------------------------------------------------------
# shared loaders
# ----------------------------------------------------------------------------

def load_res(p):
    t0 = time.time()
    res = parse_hogs(p["orthoxml"], newick=p["tree"])
    if not res.check_conservation().empty:
        raise SystemExit("conservation law violated (parser bug)")
    log(f"parsed {res.stats['n_genes']:,} genes, {res.stats['n_families']:,} families in {time.time() - t0:.0f}s")
    return res


def species_table(R):
    st = pd.read_csv(paths(R)["species"], sep="\t")
    st["is_tx"] = ~st["data_type"].astype(str).str.lower().eq("genome")
    return st


def transitions(R):
    return pd.read_csv(paths(R)["transitions"], sep="\t")


def require_done(p):
    if not os.path.exists(p["done"]):
        raise SystemExit(f"{p['done']} missing: Ancestra run not finished")


# ----------------------------------------------------------------------------
# 1. parse
# ----------------------------------------------------------------------------

def cmd_parse(a):
    p = paths(a.R)
    require_done(p)
    os.makedirs(p["out"], exist_ok=True)
    dst = os.path.join(p["out"], "branch_events.tsv")
    guard(dst)
    res = load_res(p)
    b = res.branches.copy()
    b.to_csv(dst, sep="\t", index=False)
    with open(os.path.join(p["out"], "parse_info.json"), "w") as fh:
        json.dump(dict(orthoxml=p["orthoxml"], orthoxml_md5=md5(p["orthoxml"]), stats=res.stats,
                       ancestra_done=json.load(open(p["done"]))), fh, indent=1, default=str)
    log(f"wrote {dst}")


# ----------------------------------------------------------------------------
# 2. phylogenomic branch lengths on the fixed topology
# ----------------------------------------------------------------------------

def _run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed ({r.returncode}): {r.stderr[-2000:]}")
    return r


def _align(args):
    raw, aln, trim = args
    if not (os.path.exists(trim) and os.path.getsize(trim) > 0):
        r = _run(["nice", "-n", "10", "mafft", "--auto", "--thread", "1", "--quiet", raw])
        with open(aln, "w") as fh:
            fh.write(r.stdout)
        _run(["nice", "-n", "10", "trimal", "-in", aln, "-out", trim + ".tmp", "-automated1"])
        os.replace(trim + ".tmp", trim)
    return trim


def read_fa(path):
    out, h = {}, None
    for line in open(path):
        if line.startswith(">"):
            h = line[1:].split()[0]
            out[h] = []
        elif h:
            out[h].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def cmd_phylo(a):
    p = paths(a.R)
    require_done(p)
    dst = os.path.join(p["out"], "tree_bl.nwk")
    guard(dst)
    os.makedirs(p["phylo"], exist_ok=True)
    mdir = os.path.join(p["phylo"], "markers")
    os.makedirs(mdir, exist_ok=True)
    res = load_res(p)
    T = res.tree
    root = T.names[0]
    sps = [T.names[s] for s in T.leaves]
    N = len(sps)
    m = res.members(root)
    cn = m.groupby(["hog_id", "species"], observed=True).size().rename("n").reset_index()
    agg = cn.groupby("hog_id").agg(n_sc=("n", lambda x: int((x == 1).sum())), n_multi=("n", lambda x: int((x > 1).sum())))
    sel = None
    for occ in (0.8, 0.7, 0.6, 0.5):
        cand = agg[(agg["n_sc"] >= occ * N) & (agg["n_multi"] <= 0.1 * N)]
        if len(cand) >= 100 or occ == 0.5:
            sel = cand
            used_occ = occ
            break
    rng = np.random.default_rng(SEED)
    sel = sel.assign(r=rng.random(len(sel))).sort_values(["n_sc", "r"], ascending=[False, True]).head(a.max_markers)
    log(f"{len(agg):,} root HOGs; {len(sel)} markers (single-copy occupancy >= {used_occ}, <=10% multi-copy)")
    mm = m[m["hog_id"].isin(sel.index)]
    sc = cn[(cn["n"] == 1) & cn["hog_id"].isin(sel.index)][["hog_id", "species"]]
    mm = mm.merge(sc, on=["hog_id", "species"])
    want = set(mm["gene"])
    seqs = {}
    for f in sorted(glob.glob(os.path.join(p["proteomes"], "*.fa"))):
        for h, s in gl.read_fasta(f, keep=want):
            seqs[h] = s.replace("*", "").replace("U", "X")
    jobs, names = [], []
    for k, (hog, g) in enumerate(mm.groupby("hog_id", sort=True)):
        nm = f"m{k:03d}"
        raw = os.path.join(mdir, f"{nm}.raw.fa")
        with open(raw, "w") as fh:
            for _, r in g.iterrows():
                fh.write(f">{r['species']}\n{seqs[r['gene']]}\n")
        jobs.append((raw, os.path.join(mdir, f"{nm}.aln.fa"), os.path.join(mdir, f"{nm}.trim.fa")))
        names.append((nm, hog))
    t0 = time.time()
    with ThreadPoolExecutor(a.threads) as ex:
        trims = list(ex.map(_align, jobs))
    log(f"aligned {len(trims)} markers in {time.time() - t0:.0f}s")
    cat = {s: [] for s in sps}
    parts, pos, minfo = [], 1, []
    for (nm, hog), tr in zip(names, trims):
        al = read_fa(tr) if os.path.getsize(tr) > 0 else {}
        if not al:
            continue
        L = len(next(iter(al.values())))
        if L < 30:
            continue
        for s in sps:
            cat[s].append(al.get(s, "-" * L))
        parts.append((nm, pos, pos + L - 1))
        minfo.append(dict(marker=nm, hog_id=hog, n_species=len(al), trim_len=L))
        pos += L
    nsites = pos - 1
    sm = os.path.join(p["phylo"], "supermatrix.fa")
    with open(sm, "w") as fh:
        for s in sps:
            fh.write(f">{s}\n{''.join(cat[s])}\n")
    pd.DataFrame(minfo).to_csv(os.path.join(p["phylo"], "markers.tsv"), sep="\t", index=False)
    occ = {s: sum(ch not in "-X?" for ch in "".join(cat[s])) / nsites for s in sps}
    guide = read_tree(p["tree"])
    te = os.path.join(p["phylo"], "guide_topology.nwk")
    with open(te, "w") as fh:
        fh.write(guide.write(format=9) + "\n")
    prefix = os.path.join(p["phylo"], "concat")
    cmd = ["nice", "-n", "10", "iqtree3", "-s", sm, "-te", te, "-m", "LG+G4", "-T", str(a.threads), "-seed", str(SEED),
           "-keep-ident", "--prefix", prefix, "-redo", "-quiet"]
    t0 = time.time()
    _run(cmd)
    fitted = read_tree(prefix + ".treefile")
    lengths = transfer_branch_lengths(guide, fitted)
    nwk = write_named(guide, lengths)
    assert_same_topology(guide, read_tree(nwk))
    with open(dst, "w") as fh:
        fh.write(nwk + "\n")
    info = dict(n_taxa=N, n_markers=len(minfo), n_sites=nsites, occupancy_threshold=used_occ, model="LG+G4",
                fixed_topology=True, iqtree_s=round(time.time() - t0, 1), cmd=" ".join(cmd),
                min_taxon_occupancy=round(min(occ.values()), 3),
                low_occupancy_taxa={s: round(v, 3) for s, v in occ.items() if v < 0.3},
                total_tree_length=round(sum(lengths.values()), 4))
    with open(os.path.join(p["out"], "tree_bl_info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    log(f"wrote {dst}: {info}")


# ----------------------------------------------------------------------------
# 3. sequence-novelty ladder
# ----------------------------------------------------------------------------

def ladder_branches(R, res):
    """Internal non-root branches (every internal node has >= 2 descendants) + event tips
    (+ all tips in ALL_TIPS_PANELS)."""
    T = res.tree
    tr = transitions(R)
    sel = [T.names[i] for i in T.internal if i != 0]
    tips = set(tr.loc[tr["M_is_tip"] == 1, "M"].astype(str))
    if base(R) in ALL_TIPS_PANELS:
        tips |= {T.names[s] for s in T.leaves}
    sel += sorted(t for t in tips if t in T.index)
    return sel


def target_table(R):
    p = paths(R)
    pc = set(panel_codes(R))
    rows = []
    files = [(c, os.path.join(p["proteomes"], c + ".fa"), "panel") for c in sorted(pc)]
    files += [(c, os.path.join(p["r1_proteomes"], c + ".fa"), "round1") for c in r1_codes(R) if c not in pc]
    for c, f, src in files:
        L = gl.fasta_lengths([f])
        rows.append(pd.DataFrame({"seqid": L.index.to_numpy(), "species": c, "domain": kingdom(R), "source": src,
                                  "length": L.to_numpy()}))
    t = pd.concat(rows, ignore_index=True)
    if not t["seqid"].is_unique:
        raise SystemExit(f"target ids not unique, e.g. {t.loc[t['seqid'].duplicated(), 'seqid'].head().tolist()}")
    return t, files


def transcriptome_species(R):
    st = species_table(R)
    tx = set(st.loc[st["is_tx"], "code"])
    tr = pd.read_csv(os.path.join(ROOT, "data", "species_traits.tsv"), sep="\t")
    tx |= set(tr.loc[tr["data_type"].astype(str).str.contains("transcriptome"), "code"])
    return tx


def cmd_ladder(a):
    p = paths(a.R)
    require_done(p)
    tips = getattr(a, "tips_only", False)
    sfx = "_tips" if tips else ""
    w = p["ladder"] + sfx
    os.makedirs(w, exist_ok=True)
    dst = os.path.join(p["out"], "tip_gain_ladder.tsv" if tips else "branch_gain_ladder.tsv")
    guard(dst)
    res = load_res(p)
    T = res.tree
    # ---- DB (shared with the main ladder)
    db = os.path.join(p["ladder"], "db")
    tpq = db + ".targets.parquet"
    if not os.path.exists(db + ".dmnd"):
        t, files = target_table(a.R)
        with open(db + ".faa.tmp", "wb") as out:
            for _, f, _ in files:
                with open(f, "rb") as fh:
                    shutil.copyfileobj(fh, out)
        os.replace(db + ".faa.tmp", db + ".faa")
        subprocess.run(["nice", "-n", "10", "diamond", "makedb", "--in", db + ".faa", "-d", db, "--threads", str(a.threads)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        t.to_parquet(tpq, index=False)
        log(f"DB: {len(t):,} proteins, {t['species'].nunique()} species ({(t.groupby('species')['source'].first() == 'round1').sum()} round-1-only)")
    targets = pd.read_parquet(tpq)
    # ---- gains, representatives, queries
    sel = ladder_branches(a.R, res)
    if tips:   # add-on: every tip not laddered by the main run (tip-control sensitivity, DEVIATIONS D-P6)
        sel = [T.names[i] for i in T.leaves if T.names[i] not in set(sel)]
    pd.DataFrame({"branch": sel}).to_csv(os.path.join(w, "branches.tsv"), sep="\t", index=False)
    gpq, rpq, qfa = (os.path.join(w, x) for x in ("gains.parquet", "reps.parquet", "queries.faa"))
    if not os.path.exists(qfa):
        gains = gl.internal_gains(res, branches=sel)
        lengths = targets.set_index("seqid")["length"]
        reps = gl.choose_representatives(res, gains, lengths, n_rep=N_REP, seed=SEED)
        gains.to_parquet(gpq, index=False)
        reps.to_parquet(rpq, index=False)
        n = gl.write_queries(reps, [db + ".faa"], qfa)
        log(f"{len(gains):,} gains on {len(sel)} branches, {n:,} representatives")
    gains, reps = pd.read_parquet(gpq), pd.read_parquet(rpq)
    # ---- target space with round-1 placement
    ts = gl.TargetSpace.build(targets.drop(columns=["source"]), res, transcriptome_species(a.R))
    pl = place_round1(a.R, ete3.Tree(p["tree"], format=1))
    if not tips:
        pl.to_csv(os.path.join(p["out"], "round1_placement.tsv"), sep="\t", index=False)
    lut = {r.species: T.index[r.panel_node] for r in pl.itertuples() if r.panel_node}
    for s, idx in lut.items():
        ts.sp_pre[ts.species == s] = idx
    # ---- two-pass search (round-1 cmd_search logic)
    seqs = dict(gl.read_fasta(qfa))
    h1 = os.path.join(w, "hits_pass1.tsv.gz")
    q1 = reps[reps["rep_rank"] == 0]
    passes = []
    if not os.path.exists(h1):
        f1 = os.path.join(w, "queries_pass1.faa")
        with open(f1, "w") as fh:
            fh.writelines(f">{g}\n{seqs[g]}\n" for g in q1["gene"])
        i1 = gl.run_diamond(f1, db + ".dmnd", h1 + ".part", threads=a.threads, sensitivity=SENS, nice=10,
                            log=os.path.join(w, "diamond_pass1.log"))
        os.replace(h1 + ".part", h1)
        i1["n_queries"] = len(q1)
        passes.append(i1)
        log(f"pass 1 done ({i1['wall_s']:.0f}s)")
    pq1, _ = gl.reduce_hits(h1, reps, gains, T, ts, top_n=1)
    hit_fams = set(pq1.loc[pq1["n_outside"] > 0, "family_id"])
    q2 = reps[(reps["rep_rank"] > 0) & ~reps["family_id"].isin(hit_fams)]
    h2 = os.path.join(w, "hits_pass2.tsv.gz")
    log(f"pass 1: {len(hit_fams):,}/{len(gains):,} gains have an outside hit; pass 2: {len(q2):,} queries")
    if len(q2) and not os.path.exists(h2):
        f2 = os.path.join(w, "queries_pass2.faa")
        with open(f2, "w") as fh:
            fh.writelines(f">{g}\n{seqs[g]}\n" for g in q2["gene"])
        i2 = gl.run_diamond(f2, db + ".dmnd", h2 + ".part", threads=a.threads, sensitivity=SENS, nice=10,
                            log=os.path.join(w, "diamond_pass2.log"))
        os.replace(h2 + ".part", h2)
        i2["n_queries"] = len(q2)
        passes.append(i2)
        log(f"pass 2 done ({i2['wall_s']:.0f}s)")
    hits = os.path.join(w, "hits.tsv.gz")
    with open(hits + ".part", "wb") as out:
        for h in (h1, h2):
            if os.path.exists(h):
                with open(h, "rb") as fh:
                    shutil.copyfileobj(fh, out)
    os.replace(hits + ".part", hits)
    # ---- classify
    pq, top = gl.reduce_hits(hits, reps, gains, T, ts, top_n=50)
    pq.to_parquet(os.path.join(w, "per_query.parquet"), index=False)
    cls = gl.classify_gains(res, gains, reps, pq, top, ts)
    cls["supported"] = cls["n_child_clades_present"] >= 2
    cls.to_parquet(os.path.join(w, "classified.parquet"), index=False)
    br = gl.branch_ladder(res, cls, branches=sel)
    g = cls.groupby("branch_child")
    br = br.set_index("branch_child")
    br["supported_gains"] = g["supported"].sum().reindex(br.index).fillna(0).astype(int)
    br["seq_novel_supported"] = cls[cls["seq_novel"]].groupby("branch_child")["supported"].sum().reindex(br.index).fillna(0).astype(int)
    br.loc[br["child_is_leaf"], ["supported_gains", "seq_novel_supported"]] = np.nan
    br = br.reset_index()
    br.to_csv(dst, sep="\t", index=False)
    keep = ["family_id", "branch_parent", "branch_child", "hog_id", "n_genes", "n_species", "n_child_clades_present",
            "supported", "n_reps", "n_hits", "n_outside", "n_outside_robust", "n_outside_transcriptome",
            "best_out_target", "best_out_species", "best_out_pident", "seq_novel", "seq_novel_tx_sens",
            "seq_novel_genomes_only", "split_repaired", "ladder"]
    cls[keep].to_csv(os.path.join(p["out"], f"gains_classified{sfx}.tsv.gz"), sep="\t", index=False)
    info = dict(n_branches=len(sel), n_gains=int(len(cls)), n_seq_novel=int(cls["seq_novel"].sum()),
                n_split_repaired=int(cls["split_repaired"].sum()), passes=passes, sensitivity=SENS, n_rep=N_REP,
                db_proteins=int(len(targets)), db_species=int(targets["species"].nunique()),
                round1_placed=lut and {k: T.names[v] for k, v in lut.items()})
    with open(os.path.join(p["out"], f"ladder{sfx}_info.json"), "w") as fh:
        json.dump(info, fh, indent=1, default=str)
    log(f"wrote {dst}: {info['n_gains']:,} gains, {info['n_seq_novel']:,} sequence-novel")


# ----------------------------------------------------------------------------
# 4. branch table + event table
# ----------------------------------------------------------------------------

def node_states(R, T):
    """Primary-coding multicellular state of every node: 'mc' at/below the M of an origin,
    origin_complex or nested event, unless at/below the M of a reversal to unicellular
    (deepest coded event wins). Sensitivity / T_event codings are not used."""
    tr = transitions(R)
    st = {}
    marks = []
    for _, r in tr.iterrows():
        if r["M"] not in T.index:
            raise SystemExit(f"event {r['event_id']}: M {r['M']} not in the tree")
        k = r["kind"]
        if k in ("origin", "origin_complex", "nested"):
            marks.append((T.index[r["M"]], "mc"))
        elif k == "reversal" and r["to_state"] == "unicellular":
            marks.append((T.index[r["M"]], "uni"))
    for i in range(len(T.names)):
        best, d = "uni", -1
        for m, s in marks:
            if T.is_anc(m, i) and T.depth[m] > d:
                best, d = s, T.depth[m]
        st[T.names[i]] = best
    return st


def cmd_tables(a):
    p = paths(a.R)
    dst_b = os.path.join(p["out"], "branch_table.tsv")
    dst_e = os.path.join(p["out"], "event_table.tsv")
    guard(dst_b)
    guard(dst_e)
    res = load_res(p)
    T = res.tree
    bl_tree = read_tree(os.path.join(p["out"], "tree_bl.nwk"))
    nodes = {n.name: n for n in bl_tree.traverse()}
    st = species_table(a.R).set_index("code")
    tr = transitions(a.R)
    lad = pd.read_csv(os.path.join(p["out"], "branch_gain_ladder.tsv"), sep="\t").set_index("branch_child")
    be = res.branches.set_index("branch_child")
    state = node_states(a.R, T)
    ev_of = {}
    for _, r in tr.iterrows():
        ev_of.setdefault(r["M"], []).append(r["event_id"])
    leaves_all = [T.names[s] for s in T.leaves]
    rows = []
    for i in range(1, len(T.names)):
        c = T.names[i]
        par = T.names[T.parent[i]]
        tips = [T.names[s] for s in T.leaves if T.is_anc(i, s)]
        n = nodes[c]
        out_leaves = [x for x in leaves_all if x not in set(tips)]
        dist_out = min(bl_tree.get_distance(n, nodes[x]) for x in out_leaves)
        sub = st.loc[tips]
        rows.append(dict(
            panel=a.R, branch_child=c, branch_parent=par, is_leaf=T.is_leaf[i], depth=T.depth[i],
            root_dist=bl_tree.get_distance(bl_tree, n), n_desc=len(tips), bl=n.dist,
            frac_tx=float(sub["is_tx"].mean()), mean_busco=float(sub["busco_C"].mean()),
            dist_outgroup=dist_out, root_child=T.parent[i] == 0,
            state_child=state[c], state_parent=state[par],
            inside_crown=(state[c] == "mc" and state[par] == "mc"),
            events=";".join(ev_of.get(c, [])),
            hog_gains=be.loc[c, "gains"], losses=be.loc[c, "losses"], dups=be.loc[c, "dups"],
            dups_at_origin=be.loc[c, "dups_at_origin"],
        ))
    bt = pd.DataFrame(rows).set_index("branch_child")
    for col in ("supported_gains", "seq_novel", "seq_novel_supported", "seq_novel_tx_sens", "seq_novel_genomes_only",
                "split_repaired", "gains_after_repair"):
        bt[col] = lad[col].reindex(bt.index)
    bt["laddered"] = bt.index.isin(lad.index)
    bt["is_event"] = bt["events"] != ""
    bt["control"] = (~bt["is_leaf"]) & (~bt["root_child"]) & (~bt["is_event"]) & (~bt["inside_crown"]) & bt["laddered"]
    # tip controls (sensitivity, laddered panels only): non-event tips of non-crown (unicellular-state) species
    bt["tip_control"] = bt["is_leaf"] & (~bt["is_event"]) & (bt["state_child"] == "uni") & bt["laddered"]
    # complex-crown-only definition (sensitivity): exclude only branches inside complex crowns
    bt.reset_index().to_csv(dst_b, sep="\t", index=False)
    # ---- events
    erows = []
    ctrl = bt[bt["control"]]
    for _, r in tr.iterrows():
        c = r["M"]
        if c not in bt.index:
            continue
        b = bt.loc[c]
        P_ok = (b["branch_parent"] == r["P"])
        erows.append(dict(panel=a.R, event_id=r["event_id"], kind=r["kind"], from_state=r["from_state"],
                          to_state=r["to_state"], confidence=r["confidence"], M=c, M_is_tip=int(r["M_is_tip"]),
                          P=r["P"], P_is_parent=P_ok, depth=b["depth"], root_dist=b["root_dist"], n_desc=b["n_desc"],
                          bl=b["bl"], frac_tx=b["frac_tx"], mean_busco=b["mean_busco"], dist_outgroup=b["dist_outgroup"],
                          hog_gains=b["hog_gains"], supported_gains=b["supported_gains"], seq_novel=b["seq_novel"],
                          seq_novel_supported=b["seq_novel_supported"], dups=b["dups"], losses=b["losses"],
                          n_controls_panel=int(len(ctrl))))
    ev = pd.DataFrame(erows)
    # depth SMD (events vs in-panel controls), topological depth and root distance
    for col in ("depth", "root_dist"):
        sd = np.sqrt((ev[col].var(ddof=1) + ctrl[col].var(ddof=1)) / 2) if len(ev) > 1 and len(ctrl) > 1 else np.nan
        ev[f"smd_{col}_all_events_vs_controls"] = (ev[col].mean() - ctrl[col].mean()) / sd if sd and sd > 0 else np.nan
    ev.to_csv(dst_e, sep="\t", index=False)
    log(f"wrote {dst_b} ({int(bt['control'].sum())} controls) and {dst_e} ({len(ev)} events)")


# ----------------------------------------------------------------------------
# 5a. pre-registered Pfam subset for new species
# ----------------------------------------------------------------------------

def cmd_prereg_pfam(a):
    p = paths(a.R)
    d = os.path.join(p["work"], "pfam_prereg")
    os.makedirs(d, exist_ok=True)
    r1dir = os.path.join(ROOT, "data", "annotations", "pfam_prereg", kingdom(a.R))
    new = [c for c in panel_codes(a.R) if not os.path.exists(os.path.join(r1dir, f"{c}.pfam.tsv.gz"))]
    for c in panel_codes(a.R):
        src = os.path.join(r1dir, f"{c}.pfam.tsv.gz")
        dst = os.path.join(d, f"{c}.pfam.tsv.gz")
        if os.path.exists(src) and not os.path.exists(dst):
            os.symlink(src, dst)
    if new:
        cmd = ["nice", "-n", "10", os.path.join(ROOT, "conda_envs", "annot", "bin", "python"),
               os.path.join(ROOT, "scripts", "annotate", "run_pfam.py"), "--proteomes", p["proteomes"],
               "--hmm", os.path.join(ROOT, "data", "raw", "pfam", "Pfam-A.prereg_subset.hmm"), "--out", d,
               "--cpus", str(a.cpus), "--species", ",".join(new)]
        subprocess.run(cmd, check=True)
    log(f"prereg Pfam: {len(new)} new species annotated in {d}")


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["parse", "phylo", "ladder", "tables", "prereg_pfam", "controls"])
    ap.add_argument("R", choices=PANELS + EXTRA)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--max-markers", type=int, default=150)
    ap.add_argument("--tips-only", action="store_true", help="ladder: add-on run over all remaining tips")
    a = ap.parse_args(argv)
    if a.cmd == "controls":
        from controls_panel import cmd_controls
        cmd_controls(a)
    else:
        globals()[f"cmd_{a.cmd}"](a)


if __name__ == "__main__":
    main()
