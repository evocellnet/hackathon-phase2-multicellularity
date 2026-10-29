#!/usr/bin/env python
"""Fetch EukaryOMA sequences of proteins not mapped by exact MD5, for species whose MD5 mapping
covers < 90 % of our proteome (different proteome versions). POST /api/protein/bulk_retrieve/
(200 ids per request, 4 concurrent), resumable per species.
Output: runs/eukaryoma_xval/seqs/<code>.fa (EukaryOMA ids, sequences).
"""
import json, os, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
URL = "https://eukaryoma.unil.ch/api/protein/bulk_retrieve/"
OUT = os.path.join(ROOT, "runs/eukaryoma_xval/seqs")
THRESH = 0.9


def post(ids, tries=6):
    for k in range(tries):
        try:
            req = urllib.request.Request(URL, data=json.dumps({"ids": ids}).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.load(r)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(5 * 2 ** k)


def main():
    os.makedirs(OUT, exist_ok=True)
    st = pd.read_csv(os.path.join(ROOT, "results/eukaryoma_xval/protein_mapping_stats.tsv"), sep="\t")
    pm = pd.read_csv(os.path.join(ROOT, "runs/eukaryoma_xval/protein_map.tsv.gz"), sep="\t")
    mapped = set(pm["oma_id"])
    todo = st[st["frac_ours_mapped"] < THRESH]
    for _, r in todo.iterrows():
        dst = os.path.join(OUT, f"{r['code']}.fa")
        if os.path.exists(dst):
            continue
        oma = pd.read_csv(os.path.join(ROOT, f"runs/eukaryoma_xval/md5/{r['eukaryoma_code']}.tsv.gz"), sep="\t")
        ids = [x for x in oma["omaid"] if x not in mapped]
        chunks = [ids[i:i + 200] for i in range(0, len(ids), 200)]
        t0 = time.time()
        with ThreadPoolExecutor(4) as ex:
            res = list(ex.map(post, chunks))
        with open(dst + ".tmp", "w") as fh:
            n = 0
            for d in res:
                for e in d:
                    t = e.get("target") or {}
                    if t.get("sequence"):
                        fh.write(f">{t['omaid']}\n{t['sequence']}\n")
                        n += 1
        os.replace(dst + ".tmp", dst)
        print(r["code"], len(ids), n, f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
