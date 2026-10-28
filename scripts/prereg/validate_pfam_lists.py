#!/usr/bin/env python3
"""Validate Pfam accessions/names in the pre-registration tables against Pfam-A.hmm.

Checks (stdlib only):
  * every accession in data/positive_controls.tsv, data/toolkit_pfam.tsv and
    data/h12_sets.tsv exists in data/raw/pfam/Pfam-A.hmm (version-stripped ACC);
  * pfam_name in toolkit_pfam.tsv / h12_sets.tsv equals the HMM NAME (mismatches
    reported; exit code 1 if any);
  * (category, pfam_acc) / (module, pfam_acc) rows are unique;
  * positive_controls.tsv: controlled vocabularies for run/feature_type/
    expected_event, non-empty pass_rule, '@toolkit:<category>' references resolve
    to categories present in toolkit_pfam.tsv, hog_anchor rows name an anchor;
  * toolkit categories and h12 modules are within the allowed sets.
Usage: python3 scripts/prereg/validate_pfam_lists.py [--root PROJECT_ROOT]
"""
import argparse, csv, gzip, os, re, sys

CATS = {"adhesion_ecm", "signalling", "kinase_ptyr", "tf_dbd", "pcd", "cell_cycle", "c2h2_separate"}
MODULES = {"IFT-B", "BBSome", "axoneme", "transition_zone", "phagotrophy_exploratory"}
EVENTS = {"gain", "co-option_present_at_P", "loss", "duplication", "expansion", "presence", "null_no_expansion"}
FTYPES = {"pfam", "pfam_architecture", "hog_anchor"}
EXPECTED = {"loss", "no_loss", "not_applicable"}
ACC_RE = re.compile(r"^PF\d{5}$")


def load_pfam(path):
    """Return {acc_without_version: (full_acc, name)}."""
    out, name = {}, None
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC "):
                full = line.split(None, 1)[1].strip()
                out[full.split(".")[0]] = (full, name)
    return out


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return rd.fieldnames, list(rd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    a = ap.parse_args()
    D = os.path.join(a.root, "data")
    pfam = load_pfam(os.path.join(D, "raw", "pfam", "Pfam-A.hmm"))
    ver = ""
    vp = os.path.join(D, "raw", "pfam", "Pfam.version.gz")
    if os.path.exists(vp):
        ver = gzip.open(vp, "rt").read().splitlines()[0]
    print(f"Pfam-A.hmm: {len(pfam)} families; {ver}")
    errors, warnings = [], []

    def check_acc(acc, where):
        if not ACC_RE.match(acc):
            errors.append(f"{where}: malformed accession '{acc}'")
        elif acc not in pfam:
            errors.append(f"{where}: accession {acc} not found in Pfam-A.hmm")

    # --- toolkit
    cols, tk = read_tsv(os.path.join(D, "toolkit_pfam.tsv"))
    for c in ["category", "pfam_acc", "pfam_name", "rule_notes", "source"]:
        if c not in cols: errors.append(f"toolkit_pfam.tsv: missing column {c}")
    seen, cats = set(), set()
    for i, r in enumerate(tk, 2):
        where = f"toolkit_pfam.tsv:{i} ({r['pfam_acc']})"
        cats.add(r["category"])
        if r["category"] not in CATS: errors.append(f"{where}: unknown category {r['category']}")
        k = (r["category"], r["pfam_acc"])
        if k in seen: errors.append(f"{where}: duplicate row {k}")
        seen.add(k)
        check_acc(r["pfam_acc"], where)
        if r["pfam_acc"] in pfam and pfam[r["pfam_acc"]][1] != r["pfam_name"]:
            errors.append(f"{where}: name mismatch file='{r['pfam_name']}' HMM='{pfam[r['pfam_acc']][1]}'")
    missing_cats = CATS - cats
    if missing_cats: errors.append(f"toolkit_pfam.tsv: categories with no rows: {sorted(missing_cats)}")
    n_norm = sum("name normalised" in r["rule_notes"] for r in tk)
    print(f"toolkit_pfam.tsv: {len(tk)} rows, {len({r['pfam_acc'] for r in tk})} unique accessions, "
          f"{n_norm} rows with a reviewer-label -> HMM NAME normalisation logged in rule_notes")

    # --- h12 sets
    cols, hs = read_tsv(os.path.join(D, "h12_sets.tsv"))
    for c in ["module", "pfam_acc", "pfam_name", "notes"]:
        if c not in cols: errors.append(f"h12_sets.tsv: missing column {c}")
    seen = set()
    for i, r in enumerate(hs, 2):
        where = f"h12_sets.tsv:{i} ({r['pfam_acc']})"
        if r["module"] not in MODULES: errors.append(f"{where}: unknown module {r['module']}")
        k = (r["module"], r["pfam_acc"])
        if k in seen: errors.append(f"{where}: duplicate row {k}")
        seen.add(k)
        check_acc(r["pfam_acc"], where)
        if r["pfam_acc"] in pfam and pfam[r["pfam_acc"]][1] != r["pfam_name"]:
            errors.append(f"{where}: name mismatch file='{r['pfam_name']}' HMM='{pfam[r['pfam_acc']][1]}'")
    for m in MODULES:
        if not any(r["module"] == m and r.get("role", "member") == "member" for r in hs):
            errors.append(f"h12_sets.tsv: module {m} has no member rows")
    print(f"h12_sets.tsv: {len(hs)} rows ({sum(r.get('role') == 'member' for r in hs)} members, "
          f"{sum(r.get('role') == 'excluded' for r in hs)} excluded)")

    # --- h12 expected direction
    cols, ed = read_tsv(os.path.join(D, "h12_expected_direction.tsv"))
    for i, r in enumerate(ed, 2):
        if r["expected"] not in EXPECTED: errors.append(f"h12_expected_direction.tsv:{i}: bad expected '{r['expected']}'")
    print(f"h12_expected_direction.tsv: {len(ed)} rows")

    # --- positive controls
    cols, pc = read_tsv(os.path.join(D, "positive_controls.tsv"))
    need = ["control_id", "run", "feature_type", "pfam_accessions", "anchor_protein", "expected_node_description",
            "expected_event", "confidence", "source", "pass_rule"]
    for c in need:
        if c not in cols: errors.append(f"positive_controls.tsv: missing column {c}")
    ids = set()
    n_acc = 0
    for i, r in enumerate(pc, 2):
        where = f"positive_controls.tsv:{i} ({r['control_id']})"
        if r["control_id"] in ids: errors.append(f"{where}: duplicate control_id")
        ids.add(r["control_id"])
        if r["run"] not in ("euk", "prok"): errors.append(f"{where}: bad run")
        if r["feature_type"] not in FTYPES: errors.append(f"{where}: bad feature_type")
        for ev in r["expected_event"].split("+"):
            if ev not in EVENTS: errors.append(f"{where}: bad expected_event component '{ev}'")
        if not r["pass_rule"].strip(): errors.append(f"{where}: empty pass_rule")
        if r["feature_type"] == "hog_anchor" and not r["anchor_protein"].strip():
            errors.append(f"{where}: hog_anchor row without anchor_protein")
        for item in r["pfam_accessions"].split(","):
            item = item.strip()
            if not item: continue
            if item.startswith("@toolkit:"):
                if item.split(":", 1)[1] not in cats: errors.append(f"{where}: unknown toolkit ref {item}")
                continue
            for acc in item.split("+"):
                n_acc += 1
                check_acc(acc.strip(), where)
                # every accession a pass_rule names must be listed in the row
        for acc in set(re.findall(r"PF\d{5}", r["pass_rule"])):
            if acc not in r["pfam_accessions"]:
                warnings.append(f"{where}: pass_rule mentions {acc} not in pfam_accessions")
    n_euk = sum(r["run"] == "euk" for r in pc); n_prok = sum(r["run"] == "prok" for r in pc)
    print(f"positive_controls.tsv: {len(pc)} controls ({n_euk} euk, {n_prok} prok), {n_acc} accession mentions")

    for w in warnings: print("WARNING:", w)
    for e in errors: print("ERROR:", e)
    print(f"RESULT: {'FAIL' if errors else 'OK'} ({len(errors)} errors, {len(warnings)} warnings)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
