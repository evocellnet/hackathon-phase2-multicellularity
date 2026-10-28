#!/usr/bin/env python
"""Identity test: fastleiden.leiden_split_fast vs the frozen ancestra.families.leiden_split.

Random weighted graphs (planted communities + noise, several sizes, small gmax -> deep recursion).
Compares the full return value (list of node lists, including list order and node order) for
serial mode and for the process-pool mode with tiny inline thresholds (forces many pool tasks).
Usage: test_fastleiden.py [n_graphs_per_size]
"""
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastleiden  # noqa: E402
from ancestra.families import leiden_split  # noqa: E402


def make_graph(n, k, p_in, n_noise, rng):
    nodes = [f"c{i:07d}" for i in rng.sample(range(10 * n), n)]
    comm = {v: rng.randrange(k) for v in nodes}
    by = {}
    for v in nodes:
        by.setdefault(comm[v], []).append(v)
    best = {}
    for c, vs in by.items():
        for _ in range(int(p_in * len(vs) * 3)):
            a, b = rng.sample(vs, 2) if len(vs) > 1 else (vs[0], vs[0])
            if a != b:
                best[(min(a, b), max(a, b))] = round(rng.uniform(30, 500), 1)
    for _ in range(n_noise):
        a, b = rng.sample(nodes, 2)
        best[(min(a, b), max(a, b))] = round(rng.uniform(30, 80), 1)
    # spanning chain so the node set is connected
    order = sorted(nodes)
    for a, b in zip(order, order[1:]):
        best.setdefault((a, b), 31.0)
    edges = sorted((a, b, w) for (a, b), w in best.items())
    node_w = {v: rng.randint(1, 20) for v in nodes}
    return sorted(nodes), edges, node_w


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    rng = random.Random(12345)
    cases = []
    for n in (30, 200, 1000, 5000, 20000):
        for r in range(reps):
            k = max(2, n // rng.choice((5, 20, 50)))
            cases.append((n, make_graph(n, k, rng.choice((1.0, 2.0, 4.0)), n // rng.choice((2, 5, 20)), rng)))
    # one clique-ish case that cannot be split at some level + edges outside the node set
    n_ok = 0
    for i, (n, (nodes, edges, node_w)) in enumerate(cases):
        tot = sum(node_w.values())
        for gmax in (max(20, tot // 50), max(10, tot // 400)):
            t = time.time()
            ref = leiden_split(nodes, edges, node_w, gmax, 11)
            t_ref = time.time() - t
            fastleiden.CONFIG.update(parallel=False, progress_log=None)
            t = time.time()
            ser = leiden_split_fast(nodes, edges, node_w, gmax)
            t_ser = time.time() - t
            fastleiden.CONFIG.update(parallel=True, workers=4, inline_edges=rng.choice((0, 20, 200)))
            t = time.time()
            par = leiden_split_fast(nodes, edges, node_w, gmax)
            t_par = time.time() - t
            same = (ref == ser) and (ref == par)
            same_set = {frozenset(p) for p in ref} == {frozenset(p) for p in par}
            n_ok += same
            print(f"case {i:2d} n={n:6d} edges={len(edges):7d} gmax={gmax:6d} parts={len(ref):5d} "
                  f"identical={same} (set-identical={same_set}) t_orig={t_ref:7.2f}s t_serial={t_ser:6.2f}s "
                  f"t_pool={t_par:6.2f}s", flush=True)
            if not same:
                print("MISMATCH", file=sys.stderr)
                sys.exit(1)
    # subset call with edges touching nodes outside the node set (generic contract of leiden_split)
    nodes, edges, node_w = make_graph(3000, 60, 2.0, 1000, rng)
    sub = sorted(rng.sample(nodes, 1500))
    gmax = sum(node_w[v] for v in sub) // 100
    ref = leiden_split(sub, edges, node_w, gmax, 7)
    fastleiden.CONFIG.update(parallel=True, workers=3, inline_edges=10)
    par = leiden_split_fast(sub, edges, node_w, gmax, 7)
    print(f"subset case parts={len(ref)} identical={ref == par}")
    assert ref == par
    print(f"ALL {n_ok + 1} CASES IDENTICAL")


leiden_split_fast = fastleiden.leiden_split_fast

if __name__ == "__main__":
    main()
