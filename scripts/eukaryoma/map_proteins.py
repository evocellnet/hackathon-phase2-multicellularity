#!/usr/bin/env python
"""Protein ID mapping our proteomes <-> EukaryOMA proteins by exact sequence (MD5).

EukaryOMA's OrthoXML carries only internal OMA ids (protId="X259900001"), not UniProt/EukProt
accessions, so mapping is by identical amino-acid sequence: MD5 of our sequence (as stored, and
with a terminal '*' removed, upper-case) against the sequence_md5 that the EukaryOMA REST API
reports for every protein (fetch_eukaryoma_md5.py). When several proteins of one species share a
sequence on either side they are paired in sorted-id order (identical sequences are
interchangeable for orthology). Inputs: results/eukaryoma_xval/species_map.tsv,
data/proteomes/euk/<code>.fa, runs/eukaryoma_xval/md5/<X-code>.tsv.gz.
Outputs: runs/eukaryoma_xval/protein_map.tsv.gz (code, our_id, oma_id),
         results/eukaryoma_xval/protein_mapping_stats.tsv (per species).
"""
import gzip, hashlib, os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main():
    sm = pd.read_csv(os.path.join(ROOT, "results/eukaryoma_xval/species_map.tsv"), sep="\t")
    pairs, stats = [], []
    for _, r in sm.iterrows():
        code, x = r["code"], r["eukaryoma_code"]
        oma = pd.read_csv(os.path.join(ROOT, f"runs/eukaryoma_xval/md5/{x}.tsv.gz"), sep="\t",
                          keep_default_na=False)
        by_md5 = oma.groupby("md5")["omaid"].apply(lambda s: sorted(s)).to_dict()
        ours = {}
        n_ours = 0
        for pid, s in read_fasta(os.path.join(ROOT, f"data/proteomes/euk/{code}.fa")):
            n_ours += 1
            s = s.upper()
            h = hashlib.md5(s.encode()).hexdigest()
            if h not in by_md5 and s.endswith("*"):
                h = hashlib.md5(s.rstrip("*").encode()).hexdigest()
            ours.setdefault(h, []).append(pid)
        n_map = 0
        for h, ids in ours.items():
            if h in by_md5:
                for a, b in zip(sorted(ids), by_md5[h]):
                    pairs.append((code, a, b))
                    n_map += 1
        n_main = int(oma["main_isoform"].astype(int).sum())
        stats.append(dict(code=code, eukaryoma_code=x, species=r["species"], match=r["match"],
                          n_ours=n_ours, n_oma=len(oma), n_oma_main_isoform=n_main, n_mapped=n_map,
                          frac_ours_mapped=round(n_map / n_ours, 4), frac_oma_mapped=round(n_map / len(oma), 4)))
    pm = pd.DataFrame(pairs, columns=["code", "our_id", "oma_id"])
    assert pm["our_id"].is_unique and pm["oma_id"].is_unique
    pm.to_csv(os.path.join(ROOT, "runs/eukaryoma_xval/protein_map.tsv.gz"), sep="\t", index=False)
    st = pd.DataFrame(stats)
    st.to_csv(os.path.join(ROOT, "results/eukaryoma_xval/protein_mapping_stats.tsv"), sep="\t", index=False)
    print(st[["frac_ours_mapped", "frac_oma_mapped"]].describe().to_string())
    print(f"total mapped {len(pm):,} of {st['n_ours'].sum():,} ours / {st['n_oma'].sum():,} EukaryOMA")
    print(st.sort_values("frac_ours_mapped").head(25).to_string())


if __name__ == "__main__":
    main()
