#!/usr/bin/env python
"""Run the frozen Ancestra pipeline with two runtime monkeypatches (nothing under external/ is modified).

1. ``ancestra.families.leiden_split`` -> ``fastleiden.leiden_split_fast``: identical return value
   (same node lists, same order), but edges are distributed parent->children instead of rescanning the
   whole component for every part, and independent sub-parts run in a forkserver process pool.
2. ``ancestra.families._run`` -> DIAMOND commands whose outputs are available in a cache directory
   (copied from an earlier run) are skipped and the cached output is hard-linked into place, but ONLY if
   (a) the command, with the s0 directory replaced by a placeholder, equals the cached run's command
   and (b) the md5 of every input file the command reads equals the md5 recorded for the cached run.
   Otherwise (strict mode, default) the run aborts.

Usage:
  run_ancestra_fastleiden.py [--s0-cache DIR] [--workers N] [--progress-log F] [--provenance F] \
      -- <ancestra run arguments, e.g. --proteomes .. --tree .. --out .. --threads 12 --preset fast>

The cache directory holds stage1_clusters.tsv / reps.dmnd / reps_hits.tsv (any subset), orig_md5.txt
(md5sum lines for the ORIGINAL run's canonical_input.fa, reps.fa, reps.dmnd, ...) and orig_commands.txt
(the original ``$ diamond ...`` lines from commands.log).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import fastleiden  # noqa: E402  (imported as top-level module name, as the pool workers do)


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


PROV = {"events": []}


def _event(**kw):
    kw["time"] = datetime.datetime.now().isoformat(timespec="seconds")
    PROV["events"].append(kw)
    _write_prov()


def _write_prov():
    p = PROV.get("path")
    if p:
        tmp = p + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(PROV, fh, indent=1, default=str)
        os.replace(tmp, p)


def make_run_patch(orig_run, cache_dir, strict=True):
    """Patched families._run: skip DIAMOND commands whose verified outputs are in ``cache_dir``."""
    orig_md5 = {}
    with open(os.path.join(cache_dir, "orig_md5.txt")) as fh:
        for line in fh:
            if line.strip():
                h, name = line.split(None, 1)
                orig_md5[name.strip()] = h
    orig_cmds = {}
    with open(os.path.join(cache_dir, "orig_commands.txt")) as fh:
        for line in fh:
            if line.startswith("$ diamond"):
                toks = line[2:].split()
                orig_cmds[toks[1]] = toks
    # s0 dir of the cached run = dirname of the -o / -d output of the cluster command
    def _norm(toks, s0):
        return [t.replace(s0 + "/", "<S0>/") if isinstance(t, str) else t for t in toks]

    def opt(toks, flag):
        return toks[toks.index(flag) + 1]

    spec = {   # subcommand -> (output flag, [input flags])
        "cluster": ("-o", ["-d"]),
        "makedb": ("-d", ["--in"]),
        "blastp": ("-o", ["-q", "-d"]),
    }

    def patched(cmd, log, cwd=None):
        toks = [str(c) for c in cmd]
        sub = toks[1] if len(toks) > 1 else None
        if toks[0] != "diamond" or sub not in spec or sub not in orig_cmds:
            return orig_run(cmd, log, cwd)
        oflag, iflags = spec[sub]
        out = opt(toks, oflag)
        out_file = out + ".dmnd" if sub == "makedb" and not out.endswith(".dmnd") else out
        name = os.path.basename(out_file)
        cached = os.path.join(cache_dir, name)
        if not os.path.exists(cached):
            _event(action="run", cmd=toks, reason=f"no cached {name}")
            return orig_run(cmd, log, cwd)
        s0_new = os.path.dirname(out)
        oc = orig_cmds[sub]
        s0_old = os.path.dirname(opt(oc, oflag))
        problems = []
        if _norm(toks, s0_new) != _norm(oc, s0_old):
            problems.append(f"command differs: {_norm(toks, s0_new)} vs {_norm(oc, s0_old)}")
        checks = {}
        for f in iflags:
            p = opt(toks, f)
            if sub == "blastp" and f == "-d" and not p.endswith(".dmnd"):
                p += ".dmnd"
            b = os.path.basename(p)
            if b not in orig_md5:
                problems.append(f"no recorded md5 for input {b}")
                continue
            h = md5_file(p)
            checks[b] = {"md5_new": h, "md5_orig": orig_md5[b], "match": h == orig_md5[b]}
            if h != orig_md5[b]:
                problems.append(f"input {b} md5 {h} != original {orig_md5[b]}")
        if problems:
            _event(action="mismatch", cmd=toks, problems=problems, checks=checks)
            if strict:
                raise RuntimeError("fastleiden cache: refusing to reuse / " + "; ".join(problems))
            return orig_run(cmd, log, cwd)
        if os.path.exists(out_file):
            os.remove(out_file)
        try:
            os.link(cached, out_file)
            how = "hardlink"
        except OSError:
            import shutil
            shutil.copyfile(cached, out_file)
            how = "copy"
        log.write("# [fastleiden wrapper] SKIPPED (verified cache hit, " + how + " of " + cached + "): $ "
                  + " ".join(toks) + "\n")
        log.write("# [fastleiden wrapper] input md5: " + json.dumps(checks) + "\n")
        log.flush()
        _event(action="skipped_cached", cmd=toks, cached=cached, placed=how, input_checks=checks,
               cached_output_md5=orig_md5.get(name))

    return patched


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--" not in argv:
        sys.exit("usage: run_ancestra_fastleiden.py [opts] -- <ancestra run args>")
    k = argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--s0-cache")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--inline-edges", type=int, default=50_000)
    ap.add_argument("--progress-log")
    ap.add_argument("--provenance", help="sidecar JSON (default <out>/fastleiden_patch.json)")
    ap.add_argument("--no-strict", action="store_true")
    ap.add_argument("--stop-after", help="stop after this pipeline stage (e.g. families; for tests)")
    a = ap.parse_args(argv[:k])
    run_args = argv[k + 1:]
    out = run_args[run_args.index("--out") + 1]
    os.makedirs(out, exist_ok=True)

    import ancestra
    import ancestra.families as fam
    import ancestra.pipeline as pipe

    fastleiden.CONFIG.update(workers=max(1, min(12, a.workers)), inline_edges=a.inline_edges,
                             progress_log=a.progress_log, parallel=a.workers > 1)
    orig_leiden, orig_run = fam.leiden_split, fam._run
    fam.leiden_split = fastleiden.leiden_split_fast
    if a.s0_cache:
        fam._run = make_run_patch(orig_run, os.path.abspath(a.s0_cache), strict=not a.no_strict)

    git = None
    try:
        root = os.path.dirname(os.path.dirname(os.path.abspath(ancestra.__file__)))
        git = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        pass
    PROV.update(
        path=a.provenance or os.path.join(out, "fastleiden_patch.json"),
        wrapper=os.path.abspath(__file__), wrapper_sha256=sha256_file(os.path.abspath(__file__)),
        fastleiden=os.path.abspath(fastleiden.__file__), fastleiden_sha256=sha256_file(fastleiden.__file__),
        ancestra_families=fam.__file__, ancestra_families_sha256=sha256_file(fam.__file__),
        ancestra_git=git, started=datetime.datetime.now().isoformat(timespec="seconds"), pid=os.getpid(),
        argv=sys.argv, ancestra_run_args=run_args,
        patches={"ancestra.families.leiden_split": "fastleiden.leiden_split_fast (identical output; edge "
                                                   "distribution + forkserver pool)",
                 "ancestra.families._run": (f"skip DIAMOND commands with verified cached outputs from "
                                            f"{os.path.abspath(a.s0_cache)}") if a.s0_cache else None},
        fastleiden_config=dict(fastleiden.CONFIG), s0_cache=os.path.abspath(a.s0_cache) if a.s0_cache else None)
    _write_prov()
    if a.stop_after:
        class _Stop(Exception):
            pass
        orig_stage = pipe.Runner.stage

        def stage(self, name, key, fn):
            r = orig_stage(self, name, key, fn)
            if name == a.stop_after:
                raise _Stop(name)
            return r
        pipe.Runner.stage = stage
        PROV["stop_after"] = a.stop_after
    # pipeline records sys.argv in its manifest: make it the equivalent `ancestra run ...` line + marker
    sys.argv = ["ancestra(run_ancestra_fastleiden.py)", "run"] + run_args
    rc = 1
    try:
        try:
            rc = pipe.main(run_args)
        except Exception as e:
            if type(e).__name__ != "_Stop":
                raise
            print(f"[fastleiden wrapper] stopped after stage {e}", file=sys.stderr)
            rc = 0
    finally:
        PROV["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
        PROV["rc"] = rc
        _write_prov()
    return rc


if __name__ == "__main__":
    sys.exit(main())
