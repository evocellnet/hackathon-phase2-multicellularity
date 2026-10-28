"""Output-identical, faster drop-in for ``ancestra.families.leiden_split`` (freeze-v1).

Original algorithm (families.py L400-431): stack (LIFO) of node lists; for each popped node list ``ns``
that is heavier than ``gmax`` it (1) rescans the WHOLE component edge list to collect the edges inside
``ns`` (in edge-list order, re-indexed by position in ``ns``), (2) builds an igraph graph, (3) runs Leiden
RBConfiguration (fixed seed) up the resolution ladder until >1 community, (4) pushes the children
``[[ns[i] for i in c] for c in partition]``.

This version computes exactly the same thing, faster:

* Edge distribution instead of rescans. A child's edge list in the original is
  ``filter(component_edges, both ends in child)``. Because child ⊂ parent, that equals
  ``filter(parent_edges, both ends in child)`` (a subsequence of a subsequence keeps the relative order),
  so each split hands every parent edge whose ends fall in the same child to that child, in parent order
  (stable grouping), and re-indexes it to the child's local positions (= position in
  ``[ns[i] for i in c]``, exactly the original ``idx``). Implemented on int32/float64 numpy arrays.
* Parallelism. The Leiden call for a node list depends only on (its node order, its edge sequence, its
  weights, seed) -- leidenalg seeds its own RNG per ``find_partition`` call -- so sibling subtrees are
  independent. The root split runs in the calling process; large children are farmed out to a
  ``forkserver`` process pool (fresh small workers, no copy-on-write inheritance of the ~30 GB parent
  heap, so per-process RSS stays honest for the memory watchdog); a worker splits a task once and
  recurses inline into children with few edges, returning big children as new tasks.
* Same output ORDER as the original: the original stack-DFS emits, for a split ``X -> [c1..ck]``,
  ``leaves(ck) + ... + leaves(c1)``; results here are assembled in that order, so the returned list is
  identical element-by-element (same node lists, same node order inside each list, same list order).
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np

LADDER = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0)   # == families.leiden_split

# set by the wrapper
CONFIG = {"workers": 12, "inline_edges": 50_000, "progress_log": None, "parallel": True}

_W = None   # per-worker node weights (global ids -> python numbers), set by the pool initializer


def _init_worker(weights):
    global _W
    _W = weights


def _leiden_once(n, src, dst, w, seed):
    """One node set: smallest ladder resolution giving >1 community. Returns list of member lists or None."""
    import igraph as ig
    import leidenalg
    g = ig.Graph(n=n, edges=list(zip(src.tolist(), dst.tolist())))
    g.es["weight"] = w.tolist()
    for gamma in LADDER:
        p = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight",
                                     resolution_parameter=gamma, seed=seed)
        if len(p) > 1:
            return [list(c) for c in p]
    return None


def _distribute(parts, src, dst, w, n):
    """Child edge arrays (parent order preserved, re-indexed to child positions)."""
    memb = np.empty(n, dtype=np.int64)
    pos = np.empty(n, dtype=np.int64)
    for k, c in enumerate(parts):
        ca = np.asarray(c, dtype=np.int64)
        memb[ca] = k
        pos[ca] = np.arange(len(ca))
    ms, md = memb[src], memb[dst]
    keep = np.nonzero(ms == md)[0]                      # ascending -> parent order
    ck = ms[keep]
    order = keep[np.argsort(ck, kind="stable")]         # group by child, stable within child
    counts = np.bincount(ck, minlength=len(parts))
    bounds = np.concatenate(([0], np.cumsum(counts)))
    out = []
    for k in range(len(parts)):
        sel = order[bounds[k]:bounds[k + 1]]
        out.append((pos[src[sel]].astype(np.int32), pos[dst[sel]].astype(np.int32), w[sel]))
    return out


def _needs_split(gids, W, gmax):
    return not (len(gids) == 1 or sum(W[i] for i in gids.tolist()) <= gmax)


def _process(gids, src, dst, w, gmax, seed, W, inline_edges, inline_all=False):
    """Process one node set that needs splitting. Returns DFS-ordered items:
    ('L', gids) leaf  |  ('T', (gids, src, dst, w)) subtree left to the caller."""
    items = []
    stack = [("X", (gids, src, dst, w))]      # explicit stack reproducing the original DFS order
    while stack:
        kind, payload = stack.pop()
        if kind == "L":
            items.append(payload)
            continue
        g_, s_, d_, w_ = payload
        parts = _leiden_once(len(g_), s_, d_, w_, seed)
        if parts is None:               # cannot split further: keep as is
            items.append(("L", g_))
            continue
        child_edges = _distribute(parts, s_, d_, w_, len(g_))
        # original pushes children c1..ck, pops ck first -> emit in reverse order; we push onto our own
        # stack in forward order so that ck is processed (and emitted) first.
        for c, (cs, cd, cw) in zip(parts, child_edges):
            cg = g_[np.asarray(c, dtype=np.int64)]
            if not _needs_split(cg, W, gmax):
                stack.append(("L", ("L", cg)))
            elif inline_all or len(cs) <= inline_edges:
                stack.append(("X", (cg, cs, cd, cw)))
            else:
                stack.append(("L", ("T", (cg, cs, cd, cw))))
    return items


def _task(gids, src, dst, w, gmax, seed, inline_edges):
    t = time.time()
    items = _process(gids, src, dst, w, gmax, seed, _W, inline_edges)
    return items, time.time() - t


class _Progress:
    def __init__(self, path, total_genes, total_nodes):
        self.path, self.total, self.total_nodes = path, total_genes, total_nodes
        self.leaves = self.leaf_genes = 0
        self.t0 = time.time()
        self.last = 0.0

    def log(self, msg, force=False):
        now = time.time()
        if not force and now - self.last < 30:
            return
        self.last = now
        line = (f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] +{now - self.t0:8.0f}s leaves={self.leaves:,} "
                f"genes_done={self.leaf_genes:,}/{self.total:,} remaining_genes={self.total - self.leaf_genes:,} "
                f"{msg}")
        if self.path:
            with open(self.path, "a") as fh:
                fh.write(line + "\n")
        print("[fastleiden] " + line, file=sys.stderr, flush=True)


def leiden_split_fast(nodes: list, edges: list, node_w: dict, gmax: int, seed: int = 11):
    """Drop-in replacement for ancestra.families.leiden_split (identical return value)."""
    ns = list(nodes)
    if len(ns) == 1 or sum(node_w[n] for n in ns) <= gmax:
        return [ns]
    # global ids: position in ns; weights in python objects (sum semantics unchanged)
    W = [node_w[n] for n in ns]
    idx = {n: i for i, n in enumerate(ns)}
    es = [(idx[a], idx[b], w) for a, b, w in edges if a in idx and b in idx]   # the one full scan
    src = np.fromiter((a for a, _, _ in es), dtype=np.int32, count=len(es))
    dst = np.fromiter((b for _, b, _ in es), dtype=np.int32, count=len(es))
    wts = np.fromiter((x for _, _, x in es), dtype=np.float64, count=len(es))
    del es, idx
    gids = np.arange(len(ns), dtype=np.int64)
    total_genes = sum(W)
    prog = _Progress(CONFIG.get("progress_log"), total_genes, len(ns))
    prog.log(f"start component nodes={len(ns):,} edges={len(src):,} gmax={gmax} workers={CONFIG['workers']}",
             force=True)
    inline_edges = CONFIG["inline_edges"]

    def leafcount(items):
        for kind, x in items:
            if kind == "L":
                prog.leaves += 1
                prog.leaf_genes += sum(W[i] for i in x.tolist())

    if not CONFIG.get("parallel", True) or CONFIG["workers"] <= 1:
        items = _process(gids, src, dst, wts, gmax, seed, W, inline_edges, inline_all=True)
        leafcount(items)
        prog.log("done (serial)", force=True)
        return [[ns[i] for i in x.tolist()] for _, x in items]

    # root split, serial in this process
    t = time.time()
    root_items = _process(gids, src, dst, wts, gmax, seed, W, inline_edges=-1)   # every big child -> task
    leafcount(root_items)
    prog.log(f"root split done in {time.time() - t:.0f}s: "
             f"{sum(k == 'T' for k, _ in root_items)} subtasks, {sum(k == 'L' for k, _ in root_items)} leaves",
             force=True)
    results = {}           # task id -> items
    tree = {0: root_items}
    import multiprocessing as mp
    ctx = mp.get_context("forkserver")
    ctx.set_forkserver_preload(["fastleiden"])
    next_id = [1]
    pending = {}
    n_tasks_done = 0
    busy = 0.0
    with ProcessPoolExecutor(max_workers=CONFIG["workers"], mp_context=ctx, initializer=_init_worker,
                             initargs=(W,)) as ex:
        def submit_all(items):
            out = []
            for kind, x in items:
                if kind == "T":
                    tid = next_id[0]; next_id[0] += 1
                    cg, cs, cd, cw = x
                    fut = ex.submit(_task, cg, cs, cd, cw, gmax, seed, inline_edges)
                    pending[fut] = (tid, len(cg), len(cs))
                    out.append(("T", tid))
                else:
                    out.append((kind, x))
            return out

        tree[0] = submit_all(root_items)
        while pending:
            done, _ = wait(list(pending), timeout=30, return_when=FIRST_COMPLETED)
            for fut in done:
                tid, nn, ne = pending.pop(fut)
                items, dt = fut.result()
                busy += dt
                n_tasks_done += 1
                leafcount(items)
                tree[tid] = submit_all(items)
            big = max((v[1] for v in pending.values()), default=0)
            prog.log(f"tasks_done={n_tasks_done} tasks_pending={len(pending)} largest_pending_nodes={big:,} "
                     f"worker_busy_s={busy:.0f}", force=not pending)
    # flatten in original DFS order
    out = []
    stack = [(0, 0)]
    while stack:
        tid, i = stack.pop()
        items = tree[tid]
        while i < len(items):
            kind, x = items[i]
            i += 1
            if kind == "T":
                stack.append((tid, i))
                stack.append((x, 0))
                break
            out.append([ns[j] for j in x.tolist()])
    assert sum(len(p) for p in out) == len(ns), "fastleiden: partition does not cover the node set"
    prog.log(f"done: {len(out):,} parts", force=True)
    return out
