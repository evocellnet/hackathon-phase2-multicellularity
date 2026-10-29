#!/usr/bin/env python3
"""Fetch mean pLDDT (AFDB API globalMetricValue) for accessions of the Swiss-Prot-sourced proteins (no PDB on disk). <=3 req/s, resumable.
Usage: fetch_plddt.py  -> results/r2/w5_3di_validation/plddt_cache.tsv (reads per_protein_accuracy.tsv)"""
import json, time, os, urllib.request
OD = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/w5_3di_validation"
cache = OD + "/plddt_cache.tsv"
have = {l.split("\t")[0] for l in open(cache)} if os.path.exists(cache) else set()
accs = [l.split("\t")[1] for l in list(open(OD + "/per_protein_accuracy.tsv"))[1:] if l.split("\t")[2] == "sp"]
with open(cache, "a") as o:
    for a in accs:
        if a in have: continue
        v = ""
        for t in range(3):
            try:
                d = json.load(urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{a}", timeout=30))
                v = d[0].get("globalMetricValue", ""); break
            except Exception as e: time.sleep(2)
        o.write(f"{a}\t{v}\n"); o.flush(); time.sleep(0.4)
