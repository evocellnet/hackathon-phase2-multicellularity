#!/usr/bin/env python3
"""Median BUSCO completeness (busco_C) per domain, for the report (Section 3). Output: results/qc/busco_medians.tsv"""
import csv, os, statistics
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rows = []
for dom, f in (("eukaryotes", "data/species_final.tsv"), ("prokaryotes", "data/species_final_prok.tsv")):
    v = [float(r["busco_C"]) for r in csv.DictReader(open(os.path.join(ROOT, f)), delimiter="\t") if r["busco_C"] not in ("", "NA")]
    rows.append(dict(domain=dom, n=len(v), median_busco_C=round(statistics.median(v), 1),
                     n_below_50=sum(x < 50 for x in v), n_below_80=sum(x < 80 for x in v)))
out = os.path.join(ROOT, "results/qc/busco_medians.tsv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
