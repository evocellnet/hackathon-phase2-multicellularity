#!/usr/bin/env python
"""Download per-protein metadata (OMA id, canonical id, sequence MD5, length, main-isoform flag)
for the mapped EukaryOMA genomes from the REST API (/api/genome/<code>/proteins/, 5,000 per page).

Polite: 3 concurrent requests, retries with back-off; resumable (one TSV.gz per genome).
Output: runs/eukaryoma_xval/md5/<code>.tsv.gz  (columns omaid, canonicalid, md5, length, main_isoform, oma_hog_id)
"""
import gzip, json, os, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
API = "https://eukaryoma.unil.ch/api/genome/{code}/proteins/?per_page=5000&page={page}"
OUT = os.path.join(ROOT, "runs/eukaryoma_xval/md5")


def get(url, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return json.load(r), int(r.headers.get("X-Total-Count", "0"))
        except Exception as e:  # noqa
            if k == tries - 1:
                raise
            time.sleep(5 * 2 ** k)


def fetch(code):
    dst = os.path.join(OUT, f"{code}.tsv.gz")
    if os.path.exists(dst):
        return code, "cached"
    rows, page, total = [], 1, None
    while True:
        d, tot = get(API.format(code=code, page=page))
        total = total or tot
        for p in d:
            rows.append((p["omaid"], p.get("canonicalid", ""), p["sequence_md5"], p["sequence_length"],
                         int(bool(p.get("is_main_isoform", True))), p.get("oma_hog_id", "")))
        if len(rows) >= total or not d:
            break
        page += 1
    if len(rows) != total:
        raise RuntimeError(f"{code}: got {len(rows)} of {total}")
    tmp = dst + ".tmp"
    with gzip.open(tmp, "wt") as fh:
        fh.write("omaid\tcanonicalid\tmd5\tlength\tmain_isoform\toma_hog_id\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    os.replace(tmp, dst)
    return code, len(rows)


def main():
    os.makedirs(OUT, exist_ok=True)
    sm = pd.read_csv(os.path.join(ROOT, "results/eukaryoma_xval/species_map.tsv"), sep="\t")
    with ThreadPoolExecutor(3) as ex:
        for code, n in ex.map(fetch, sm["eukaryoma_code"]):
            print(code, n, flush=True)


if __name__ == "__main__":
    main()
