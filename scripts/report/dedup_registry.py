#!/usr/bin/env python
"""Write a de-duplicated view of results/analysis_registry.tsv (the append-only log is left untouched).

Rows that repeat the same (hypothesis, statistic, tier, scheme/note-free key, estimate, p) from re-runs are collapsed to the last
occurrence; Holm annotation rows (no p-value of their own) are dropped from the multiplicity count. Output:
results/analysis_registry_dedup.tsv and a tier count printed to stdout.
"""
import pandas as pd

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
r = pd.read_csv(f"{P}/results/analysis_registry.tsv", sep="\t", dtype=str)
key = [c for c in ["hypothesis", "statistic", "tier", "n", "k", "estimate", "p"] if c in r.columns]
d = r.drop_duplicates(subset=key, keep="last")
d = d[d["p"].notna() & (d["p"] != "")] if "p" in d.columns else d
d.to_csv(f"{P}/results/analysis_registry_dedup.tsv", sep="\t", index=False)
print(f"rows: {len(r)} raw -> {len(d)} distinct analyses with a p-value")
print(d["tier"].value_counts().to_string())
