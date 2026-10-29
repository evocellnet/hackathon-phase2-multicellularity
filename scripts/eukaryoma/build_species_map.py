#!/usr/bin/env python
"""Species mapping: our 187 eukaryote codes <-> EukaryOMA genome codes.

Starts from data/eukaryoma_species_overlap.tsv (scouting: NCBI taxid and/or lowercase name),
re-checks every row against data/raw/eukaryoma/oma-species.txt and adds same-species matches
that the scouting missed because of strain/forma-level taxids (listed in EXTRA, each checked by
hand: same binomial, different infraspecific taxid). Whether a matched genome is really the same
protein set is verified afterwards by exact-sequence (MD5) mapping (fetch_eukaryoma_md5.py,
map_proteins.py): rows with a low mapped fraction are reported there.

Output: results/eukaryoma_xval/species_map.tsv
"""
import os, sys
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# our code -> (EukaryOMA code, reason); binomial identical, taxid differs only at strain/forma level
EXTRA = {
    "Volvcar": ("X0135", "species-level: ours V. carteri f. nagariensis (3068), EukaryOMA V. carteri (3067)"),
    "Trypbru": ("X0035", "species-level: ours T. brucei brucei (185431), EukaryOMA T. brucei (5691)"),
    "Thectra": ("X2653", "species-level: ours T. trahens ATCC 50062 (461836), EukaryOMA T. trahens (529818)"),
}


def main():
    ov = pd.read_csv(os.path.join(ROOT, "data/eukaryoma_species_overlap.tsv"), sep="\t", dtype=str)
    oma = pd.read_csv(os.path.join(ROOT, "data/raw/eukaryoma/oma-species.txt"), sep="\t", comment="#",
                      header=None, dtype=str,
                      names=["oma_code", "oma_taxid", "ncbi_taxid", "gtdb", "name", "source", "release"])
    oma = oma[oma["oma_code"] != "OMA_Code"].set_index("oma_code")
    rows = []
    for _, r in ov.iterrows():
        code, m = r["code"], r["match"]
        ecode = r["eukaryoma_code"] if isinstance(r["eukaryoma_code"], str) else ""
        note = ""
        if code in EXTRA:
            ecode, note = EXTRA[code]
            m = "species_manual"
        if ecode:
            assert ecode in oma.index, ecode
            o = oma.loc[ecode]
            ours_bin = " ".join(r["species"].split()[:2]).lower()
            their_bin = " ".join(o["name"].split()[:2]).lower()
            if m == "taxid" and not note:
                note = f"taxid match; names differ ({r['species']} vs {o['name']}): synonym/strain, accepted"
            elif m == "name" and not note:
                note = f"binomial match; taxid differs ({r['ncbi_taxid']} vs {o['ncbi_taxid']}): strain vs species taxid, accepted"
            rows.append(dict(code=code, species=r["species"], ncbi_taxid=r["ncbi_taxid"],
                             eukaryoma_code=ecode, eukaryoma_name=o["name"], eukaryoma_taxid=o["ncbi_taxid"],
                             match=m, binomial_equal=ours_bin == their_bin, note=note,
                             multicellular=r["multicellular"], transition_roles=r["transition_roles"]))
    out = pd.DataFrame(rows)
    assert out["eukaryoma_code"].is_unique and out["code"].is_unique
    p = os.path.join(ROOT, "results/eukaryoma_xval/species_map.tsv")
    out.to_csv(p, sep="\t", index=False)
    print(f"{len(out)} species mapped -> {p}")
    print(out["match"].value_counts().to_string())


if __name__ == "__main__":
    main()
