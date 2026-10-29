"""Shared helpers for round-2 W4 (comparative per-family analyses)."""
import os
import sys

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = os.path.join(P, "results/r2/w4")
SEEDS = [11, 23, 37, 41, 53]
SEED = 11
RSCRIPT = os.path.join(P, "conda_envs/rstats/bin/Rscript")

# Primary eukaryotic origins (data/transitions.tsv: kingdom=euk, primary_eligible=1)
# and secondary nested simple->complex events (analysis_set=secondary, euk).
# These lists are read from the file; the constants below are only used for checks.


def out(rel):
    """Return an output path under results/r2/w4, refusing to overwrite."""
    path = os.path.join(OUT, rel)
    if os.path.exists(path):
        sys.exit(f"REFUSE: {path} exists; not overwriting")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def load_transitions():
    import pandas as pd
    t = pd.read_csv(os.path.join(P, "data/transitions.tsv"), sep="\t")
    t = t[t.kingdom == "euk"].copy()
    t["mc_list"] = t.mc_species.fillna("").map(lambda s: [x for x in s.split(",") if x])
    t["sis_list"] = t.sister_species.fillna("").map(lambda s: [x for x in s.split(",") if x])
    return t


def primary_origins(t):
    return t[t.primary_eligible == 1].copy()


def secondary_events(t):
    return t[(t.analysis_set == "secondary")].copy()
