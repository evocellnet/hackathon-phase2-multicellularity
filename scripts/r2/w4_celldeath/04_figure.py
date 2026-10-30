#!/usr/bin/env python
"""Figure 5 (report v2; was Figure 4b before renumbering): cell-death family copy numbers and NB-ARC/NACHT architectures per lineage
(multicellular clade vs sister), and cell-death families gained/duplicated on the dense-panel simple->complex events.

Reads results/r2/w4_celldeath/{celldeath_proteins.tsv.gz, stand_proteins.tsv.gz, panels_event_celldeath.tsv,
panels_celldeath_families.tsv} and results/r2/w4/scan/events_used.tsv. Style from scripts/r2/figures/common.py.
Writes report/figures_v2/fig5_celldeath_architectures.{png,pdf}; also the plotted matrices to
results/r2/w4_celldeath/fig4b_matrix_{families,stand}.tsv.
(The matrix files keep their historical fig4b_ prefix so that source references in the report stay valid.)
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/scripts/r2/figures")
from common import (P, MM, INK, INK2, GRID, COMPLEX, plt, panel_label, save)  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, PowerNorm  # noqa: E402

OUT = os.path.join(P, "results/r2/w4_celldeath")
EVENTS = [("T_metazoa", "Animals"), ("T_fungi_hyphae", "Fungi\n(hyphae)"),
          ("T_fungi_pezizomycotina", "Pezizo-\nmycotina"), ("T_fungi_agaricomycetes", "Agarico-\nmycetes"),
          ("T_dictyostelia", "Dictyo-\nstelids"), ("T_streptophyta", "Strepto-\nphytes"),
          ("T_embryophyta", "Land\nplants"), ("T_rhodophyta", "Red\nalgae"),
          ("T_px_filaments", "Brown algae\n& relatives"), ("T_volvocine_tgv", "Volvo-\ncines")]

# rows of panel a: (label, function group, predicate on the set of cell-death accessions of a protein)
FAM_ROWS = [
    ("NB-ARC", "receptors", lambda s: "PF00931" in s),
    ("NACHT", "receptors", lambda s: "PF05729" in s),
    ("TIR", "receptors", lambda s: bool(s & {"PF01582", "PF13676"})),
    ("HET", "fungal non-self", lambda s: "PF06985" in s),
    ("HeLo", "fungal non-self", lambda s: "PF14479" in s),
    ("Death", "adaptors", lambda s: "PF00531" in s),
    ("CARD", "adaptors", lambda s: "PF00619" in s),
    ("DED", "adaptors", lambda s: "PF01335" in s),
    ("Caspase", "proteases", lambda s: "PF00656" in s and "PF31181" not in s),
    ("Metacaspase", "proteases", lambda s: "PF31181" in s),
    ("Legumain (VPE)", "proteases", lambda s: "PF01650" in s),
    ("Bcl-2 / BH4", "regulators", lambda s: bool(s & {"PF00452", "PF02180"})),
    ("BIR (IAP)", "regulators", lambda s: "PF00653" in s),
    ("Bax inhibitor-1", "regulators", lambda s: "PF01027" in s),
    ("zf-LSD1", "regulators", lambda s: "PF06943" in s),
    ("AIF (C-term.)", "regulators", lambda s: "PF14721" in s),
    ("MACPF", "pore formers", lambda s: "PF01823" in s),
    ("Gasdermin", "pore formers", lambda s: bool(s & {"PF04598", "PF26164"})),
]
GROUP_LABEL = {"receptors": "immune-type\nreceptors", "fungal non-self": "fungal\nnon-self", "adaptors": "death-fold\nadaptors",
               "proteases": "proteases", "regulators": "regulators", "pore formers": "pore\nformers"}


def stand_row(r):
    ct = r.cterm
    if r.nlr_like:
        reps = ct.split("+")
        if len(reps) > 1:
            return "with mixed repeats"
        return {"LRR": "with LRR", "WD40": "with WD40", "TPR": "with TPR", "ANK": "with ankyrin"}.get(
            reps[0], "with mixed repeats")
    if r.nterm not in ("none", "other"):
        return "effector domain, no repeats"
    return "no repeats, no effector"


STAND_ROWS = ["with LRR", "with WD40", "with TPR", "with ankyrin", "with mixed repeats",
              "effector domain, no repeats", "no repeats, no effector"]


def main():
    ev = pd.read_csv(os.path.join(P, "results/r2/w4/scan/events_used.tsv"), sep="\t").set_index("event")
    prot = pd.read_csv(os.path.join(OUT, "celldeath_proteins.tsv.gz"), sep="\t")
    prot["s"] = prot.cd_accs.str.split(";").map(set)
    stand = pd.read_csv(os.path.join(OUT, "stand_proteins.tsv.gz"), sep="\t")
    stand["row"] = stand.apply(stand_row, axis=1)

    cols, Mf, Ms = [], [], []
    for e, lab in EVENTS:
        for grp, field in (("sister", "sister_species"), ("MC", "mc_species")):
            codes = ev.loc[e, field].split(",")
            sub = prot[prot.species.isin(codes)]
            Mf.append([sum(p(s) for s in sub.s) / len(codes) for _, _, p in FAM_ROWS])
            ss = stand[stand.species.isin(codes)]
            Ms.append([(ss.row == r).sum() / len(codes) for r in STAND_ROWS])
            cols.append((e, lab, grp, len(codes)))
    Mf, Ms = np.array(Mf).T, np.array(Ms).T
    cn = [f"{c[0]}|{c[2]}|n={c[3]}" for c in cols]
    pd.DataFrame(Mf, index=[r[0] for r in FAM_ROWS], columns=cn).to_csv(
        os.path.join(OUT, "fig4b_matrix_families.tsv"), sep="\t", float_format="%.3f")
    pd.DataFrame(Ms, index=STAND_ROWS, columns=cn).to_csv(
        os.path.join(OUT, "fig4b_matrix_stand.tsv"), sep="\t", float_format="%.3f")

    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f8fb", "#9cc3e0", COMPLEX, "#00365a"])
    fig = plt.figure(figsize=(180 * MM, 240 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[len(FAM_ROWS), len(STAND_ROWS) + 0.2, 12], hspace=0.36,
                          left=0.25, right=0.95, top=0.95, bottom=0.08)

    def heat(ax, M, rows, title, vmax):
        im = ax.imshow(M, aspect="auto", cmap=cmap, norm=PowerNorm(0.4, vmin=0, vmax=vmax), interpolation="nearest")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                if v == 0:
                    ax.text(j, i, "·", ha="center", va="center", fontsize=6, color=INK2)
                    continue
                txt = f"{v:.0f}" if v >= 9.5 else (f"{v:.1f}" if v >= 0.05 else "<0.1")
                frac = (v / vmax) ** 0.4
                ax.text(j, i, txt, ha="center", va="center", fontsize=4.6,
                        color="white" if frac > 0.55 else INK)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(["sis." if c[2] == "sister" else "MC" for c in cols], fontsize=5.5)
        ax.tick_params(length=0)
        for k in range(1, len(EVENTS)):
            ax.axvline(2 * k - 0.5, color="white", lw=2.2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        for k, (e, lab) in enumerate(EVENTS):
            ax.text(2 * k + 0.5, -0.9, lab, ha="center", va="bottom", fontsize=5.8, color=INK)
        ax.set_title(title, loc="left", pad=24, color=INK)
        return im

    ax = fig.add_subplot(gs[0])
    vmax = max(Mf.max(), 1)
    im = heat(ax, Mf, [r[0] for r in FAM_ROWS], "Proteins per species carrying each cell-death domain "
              "(MC = multicellular clade, sis. = its sister clade)", vmax)
    # function group brackets
    grp = [r[1] for r in FAM_ROWS]
    i = 0
    while i < len(grp):
        j = i
        while j + 1 < len(grp) and grp[j + 1] == grp[i]:
            j += 1
        x = -0.21
        ax.plot([x + 0.01, x + 0.01], [i - 0.35, j + 0.35], transform=ax.get_yaxis_transform(), color=INK2,
                lw=0.6, clip_on=False)
        ax.text(x, (i + j) / 2, GROUP_LABEL[grp[i]], transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=5.3, color=INK2)
        if j + 1 < len(grp):
            ax.axhline(j + 0.5, color="white", lw=1.2)
        i = j + 1
    panel_label(ax, "a", x=-0.30, y=1.07)

    ax2 = fig.add_subplot(gs[1])
    heat(ax2, Ms, STAND_ROWS, "NB-ARC / NACHT proteins per species, by domain architecture", vmax)
    panel_label(ax2, "b", x=-0.30, y=1.12)
    cb = fig.colorbar(im, ax=[ax, ax2], fraction=0.012, pad=0.01, aspect=40)
    cb.set_label("proteins per species", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)

    # panel c: dense-panel simple->complex events
    t = pd.read_csv(os.path.join(OUT, "panels_event_celldeath.tsv"), sep="\t")
    f = pd.read_csv(os.path.join(OUT, "panels_celldeath_families.tsv"), sep="\t")
    f["stand"] = f.celldeath_pfam.str.contains("PF00931|PF05729")
    t = t[t.event_kind.isin(["nested", "origin_complex"])]
    names = {"R1_metazoa": "Animals (origin)", "R2_pezizomycotina": "Pezizomycotina", "R2_neolecta": "Neolecta*",
             "R2_agaricomycetes": "Agaricomycetes", "R2_tremella": "Tremella*", "R2_endogone": "Endogone*",
             "R3_embryophyta": "Land plants", "R3_chara": "Chara*", "R4_astrephomene": "Astrephomene*",
             "R4_pleodorina": "Pleodorina*", "R4_volvox_carteri": "Volvox carteri*",
             "R4_volvox_africanus": "V. africanus/reticuliferus", "R5_florideophyceae": "Florideophyceae",
             "R6_fucophycidae": "Fucophycidae", "R6_dictyota": "Dictyota*", "R8_heterocyst": "Heterocystous cyanobact."}
    order = list(names)
    ax3 = fig.add_subplot(gs[2])
    y = np.arange(len(order))
    h = 0.42
    for off, kind, lab in ((-h / 2, "seq_novel", "newly arisen families"), (h / 2, "dup", "duplicated families")):
        st = [int(((f.event_id == e) & (f.kind == kind) & f.stand).sum()) for e in order]
        ot = [int(((f.event_id == e) & (f.kind == kind) & ~f.stand).sum()) for e in order]
        exp = [float((lambda r: r.n_pfam_annotated * r.bg_frac)(t[(t.event_id == e) & (t.kind == kind)].iloc[0]))
               for e in order]
        q = [t[(t.event_id == e) & (t.kind == kind)].q_bh_simple_to_complex.iloc[0] for e in order]
        hatch = "////" if kind == "seq_novel" else None
        ax3.barh(y + off, st, height=h * 0.92, color=INK2, edgecolor="white", lw=0.4, hatch=hatch)
        ax3.barh(y + off, ot, left=st, height=h * 0.92, color="#c9c9c9", edgecolor="white", lw=0.4, hatch=hatch)
        ax3.scatter(exp, y + off, marker="|", s=28, color=INK, lw=1.0, zorder=3)
        for yi, a, b, qq in zip(y + off, st, ot, q):
            if a + b > 0 or (qq == qq and qq < 0.1):
                star = " (q<0.1)" if qq == qq and qq < 0.1 else ""
                ax3.text(a + b + 0.3, yi, f"{a + b}{star}", va="center", fontsize=5.3, color=INK)
    ax3.set_yticks(y)
    ax3.set_yticklabels([names[e] for e in order], fontsize=5.8)
    ax3.invert_yaxis()
    ax3.set_xlabel("cell-death families on the event branch")
    ax3.xaxis.grid(True, color=GRID, lw=0.5)
    ax3.set_axisbelow(True)
    ax3.set_title("Cell-death families arising or duplicated on simple-to-complex branches (focused reconstructions)", loc="left", color=INK)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    ax3.legend(handles=[Patch(facecolor=INK2, label="NB-ARC / NACHT family"),
                        Patch(facecolor="#c9c9c9", label="other cell-death family"),
                        Patch(facecolor="white", edgecolor=INK2, hatch="////", label="newly arisen (upper bar)"),
                        Patch(facecolor="white", edgecolor=INK2, label="duplicated (lower bar)"),
                        Line2D([0], [0], marker="|", color=INK, lw=0, markersize=6,
                               label="expected from control branches")],
               loc="upper center", bbox_to_anchor=(0.45, -0.16), fontsize=5.5, ncol=3)
    panel_label(ax3, "c", x=-0.30, y=1.02)
    save(fig, "fig5_celldeath_architectures")


if __name__ == "__main__":
    main()
