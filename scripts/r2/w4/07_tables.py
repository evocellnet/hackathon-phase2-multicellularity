"""Annotated per-origin top-family tables, recurrent-family table and frozen-set families by origin.
Inputs: results/r2/w4/scan/*, matrix/architectures.tsv.gz, data/raw/pfam/Pfam-A.clans.tsv.gz, data/raw/go/pfam2go,
        data/r2_sets/*.tsv, data/toolkit_pfam.tsv, data/h12_sets.tsv.
Outputs: results/r2/w4/scan/top_families_per_origin.tsv, top_families_per_origin.txt, recurrent_families.tsv,
         set_families_by_origin.tsv
"""
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import P, OUT, out  # noqa: E402

TOPN = 20
S = f"{OUT}/scan"


def load_tags():
    tags = defaultdict(list)
    e = pd.read_csv(f"{P}/data/r2_sets/elm_motif_binding.tsv", sep="\t")
    for _, r in e.iterrows():
        tags[r.pfam].append("ELM" if r.set == "primary" else "ELM_repeat")
    for fn, t in (("cell_death.tsv", "CELLDEATH"), ("tf_celltype.tsv", "TF")):
        for f in pd.read_csv(f"{P}/data/r2_sets/{fn}", sep="\t").pfam:
            tags[f].append(t)
    tk = pd.read_csv(f"{P}/data/toolkit_pfam.tsv", sep="\t")
    for _, r in tk[tk.role != "exclusion_marker"].iterrows():  # exclusion markers (NB-ARC etc.) are not toolkit members
        tags[r.pfam_acc].append(f"R1toolkit:{r.category}")
    h = pd.read_csv(f"{P}/data/h12_sets.tsv", sep="\t")
    for _, r in h.iterrows():
        tags[r.pfam_acc].append(f"R1h12:{r.module}")
    return {k: ";".join(dict.fromkeys(v)) for k, v in tags.items()}


def load_go():
    go = defaultdict(list)
    pat = re.compile(r"^Pfam:(PF\d+) \S+ > GO:(.*) ; (GO:\d+)")
    for line in open(f"{P}/data/raw/go/pfam2go"):
        m = pat.match(line)
        if m:
            go[m.group(1)].append(f"{m.group(2)} ({m.group(3)})")
    return go


def main():
    ev = pd.read_csv(f"{S}/events_used.tsv", sep="\t")
    rec = pd.read_csv(f"{S}/recurrence_observed.tsv", sep="\t").set_index("pfam")
    clans = pd.read_csv(f"{P}/data/raw/pfam/Pfam-A.clans.tsv.gz", sep="\t", header=None,
                        names=["pfam", "clan", "clan_name", "name", "description"]).set_index("pfam")
    go = load_go(); tags = load_tags()
    arch = pd.read_csv(f"{OUT}/matrix/architectures.tsv.gz", sep="\t")
    arch["accs"] = arch.accs.str.split(";")
    archx = arch.explode("accs")
    f_top, f_txt, f_rec, f_set = (out("scan/top_families_per_origin.tsv"), out("scan/top_families_per_origin.txt"),
                                  out("scan/recurrent_families.tsv"), out("scan/set_families_by_origin.tsv"))

    def top_arch(acc, species, k=3):
        sub = archx[(archx.accs == acc) & (archx.species.isin(species))]
        if not len(sub):
            return "", 0
        c = Counter(sub.arch); n = len(sub)
        return "; ".join(f"{a} [{v / n:.0%}]" for a, v in c.most_common(k)), n

    allrows, setrows, lfc_ash = [], [], {}
    lines = []
    for _, e in ev.iterrows():
        d = pd.read_csv(f"{S}/event_{e.event}.tsv", sep="\t").merge(pd.read_csv(f"{S}/ashr_{e.event}.tsv", sep="\t"), on="pfam")
        lfc_ash[e.event] = d.set_index("pfam").lfc_ash
        mc, sis = e.mc_species.split(","), e.sister_species.split(",")
        d["tags"] = d.pfam.map(lambda f: tags.get(f, ""))
        for _, r in d[d.tags != ""].iterrows():
            setrows.append(dict(event=e.event, tier=e.tier, **r.to_dict()))
        ups = d.sort_values("lfc_ash", ascending=False).head(TOPN).assign(rank_dir="expanded")
        dns = d.sort_values("lfc_ash", ascending=True).head(TOPN).assign(rank_dir="contracted")
        lines.append(f"\n=== {e.event} ({e.label}; {e.tier}) — {e.n_mc} multicellular vs {e.n_sister} sister species "
                     f"({e.n_transcriptome} transcriptomes); {e.n_families_tested} families tested; q<0.1: "
                     f"{e.n_sig_expanded} expanded, {e.n_sig_contracted} contracted")
        for part in (ups, dns):
            lines.append(f"--- top {TOPN} {part.rank_dir.iloc[0]} (ranked by ashr posterior mean log rate ratio)")
            for rk, (_, r) in enumerate(part.iterrows(), 1):
                a_mc, n_mc = top_arch(r.pfam, mc); a_s, n_s = top_arch(r.pfam, sis)
                nm = clans.name.get(r.pfam, ""); desc = clans.description.get(r.pfam, "")
                g = go.get(r.pfam, [])
                row = dict(event=e.event, tier=e.tier, direction=r.rank_dir, rank=rk, pfam=r.pfam, name=nm, description=desc,
                           clan=clans.clan.get(r.pfam, ""), mean_count_mc=r.mean_count_mc, mean_count_sister=r.mean_count_sister,
                           prev_mc=r.prev_mc, prev_sister=r.prev_sister, lfc=r.lfc, se=r.se, lfc_ash=r.lfc_ash, lfsr=r.lfsr,
                           p=r.p, q_bh=r.q_bh, rec_primary_expanded=rec.rec_primary_expanded.get(r.pfam, 0),
                           rec_primary_contracted=rec.rec_primary_contracted.get(r.pfam, 0), set_tags=r.tags,
                           go_terms=" | ".join(g[:6]), n_proteins_mc=n_mc, top_arch_mc=a_mc, n_proteins_sister=n_s, top_arch_sister=a_s)
                allrows.append(row)
                flag = f"  <<{r.tags}>>" if r.tags else ""
                lines.append(f"{rk:2d}. {r.pfam} {nm} — {desc}{flag}\n    mean copies MC {r.mean_count_mc:.1f} vs sister {r.mean_count_sister:.1f}; "
                             f"lfc {r.lfc:+.2f} (ash {r.lfc_ash:+.2f}, lfsr {r.lfsr:.2g}); q {r.q_bh:.2g}; recurrence(primary) "
                             f"exp {row['rec_primary_expanded']} / con {row['rec_primary_contracted']}\n    GO: {'; '.join(g[:3]) or '-'}\n"
                             f"    arch MC: {a_mc or '-'}\n    arch sister: {a_s or '-'}")
    pd.DataFrame(allrows).to_csv(f_top, sep="\t", index=False)
    open(f_txt, "w").write("Per-origin top expanded/contracted Pfam families (exploratory; task C). Transposon mask applied.\n"
                           "Architectures: ordered Pfam domain strings of proteins carrying the family, consecutive repeats collapsed; % of those proteins.\n"
                           + "\n".join(lines) + "\n")
    pd.DataFrame(setrows).to_csv(f_set, sep="\t", index=False)
    # recurrent families
    prim = ev.event[ev.tier == "primary"].tolist(); sec = ev.event[ev.tier == "secondary"].tolist()
    L = pd.DataFrame(lfc_ash)
    rr = rec[(rec.rec_primary_expanded >= 2) | (rec.rec_primary_contracted >= 2) | (rec.rec_secondary_expanded >= 2) | (rec.rec_secondary_contracted >= 2)].copy()
    rr = rr[[c for c in rr.columns if c.startswith(("rec_", "perm_p_", "n_primary", "sig_"))]]
    rr["name"] = [clans.name.get(f, "") for f in rr.index]
    rr["description"] = [clans.description.get(f, "") for f in rr.index]
    rr["mean_lfc_ash_primary"] = L.reindex(rr.index)[prim].mean(1)
    rr["n_primary_lfc_ash_pos"] = (L.reindex(rr.index)[prim] > 0).sum(1)
    rr["set_tags"] = [tags.get(f, "") for f in rr.index]
    rr["go_terms"] = [" | ".join(go.get(f, [])[:4]) for f in rr.index]
    rr = rr.sort_values(["rec_primary_expanded", "perm_p_primary_expanded"], ascending=[False, True])
    rr.to_csv(f_rec, sep="\t")
    print(rr.head(30)[["name", "rec_primary_expanded", "perm_p_primary_expanded", "rec_primary_contracted", "mean_lfc_ash_primary", "set_tags"]].to_string())


if __name__ == "__main__":
    main()
