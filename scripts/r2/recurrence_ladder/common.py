"""Shared helpers for the recurrence ladder (D17; rules in results/r2/recurrence_ladder/RULES.txt)."""
import os
import numpy as np

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = f"{P}/results/r2/recurrence_ladder"
R = f"{P}/results"

R1_ORIGINS = ["T_metazoa", "T_fungi_hyphae", "T_dictyostelia", "T_streptophyta", "T_rhodophyta", "T_px_filaments",
              "T_volvocine_tgv"]
R1_NESTED = {"T_fungi_pezizomycotina": "T_fungi_hyphae", "T_fungi_agaricomycetes": "T_fungi_hyphae",
             "T_embryophyta": "T_streptophyta"}
LABEL = {"T_metazoa": "animals", "T_fungi_hyphae": "fungal hyphae", "T_dictyostelia": "dictyostelids",
         "T_streptophyta": "streptophytes", "T_rhodophyta": "red algae", "T_px_filaments": "brown-algal (PX) filaments",
         "T_volvocine_tgv": "volvocines", "T_fungi_pezizomycotina": "Pezizomycotina", "T_fungi_agaricomycetes":
         "Agaricomycetes", "T_embryophyta": "land plants"}


def out(name):
    """Refuse to overwrite an existing output."""
    path = f"{OUT}/{name}"
    if os.path.exists(path):
        raise SystemExit(f"refusing to overwrite {path}")
    return path


def poibin_pmf(pis):
    """Exact Poisson-binomial pmf by dynamic programming."""
    pmf = np.zeros(len(pis) + 1)
    pmf[0] = 1.0
    for p in pis:
        pmf[1:] = pmf[1:] * (1 - p) + pmf[:-1] * p
        pmf[0] *= (1 - p)
    return pmf


def poibin_sf(pis, k):
    """P(X >= k)."""
    if k <= 0:
        return 1.0
    pmf = poibin_pmf(list(pis))
    return float(pmf[k:].sum()) if k < len(pmf) else 0.0


def holm(p):
    p = np.asarray(p, float)
    n = np.sum(~np.isnan(p))
    o = np.argsort(np.where(np.isnan(p), np.inf, p))
    adj = np.full_like(p, np.nan)
    run = 0.0
    for i, j in enumerate(o):
        if np.isnan(p[j]):
            continue
        run = max(run, min(1.0, (n - i) * p[j]))
        adj[j] = run
    return adj


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.minimum(r, 1.0)
    q[ok] = qq
    return q
