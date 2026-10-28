#!/usr/bin/env python
"""Report figure 2: workflow schematic (data -> Ancestra HOGs -> per-branch events -> gain-evidence ladder ->
stem-vs-control tests; annotation; rates). Pure drawing; numbers in the boxes are copied from DATA.md,
PREREG.md and results/validation_hogparser/summary.md (see report/REPORT_draft_methods.md for sources).

Usage: python scripts/report/fig_workflow.py  -> figures/fig2_workflow.{png,svg}
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INK, INK2, MUTED = "#1b1b1b", "#52514e", "#8f8d88"
FILL = "#f4f3ef"
ACCENT = {"data": "#2a78d6", "recon": "#eb6834", "test": "#1baf7a"}
PREREG = "#4a3aa7"

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})

fig = plt.figure(figsize=(13, 7))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 13)
ax.axis("off")


def box(x, y, w, h, title, body, lane, dashed=False):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.08", facecolor=FILL,
                       edgecolor=INK2 if not dashed else MUTED, lw=1.0, ls="--" if dashed else "-", zorder=2)
    ax.add_patch(p)
    ax.add_patch(FancyBboxPatch((x, y), 0.07, h, boxstyle="square,pad=0", facecolor=ACCENT[lane],
                                edgecolor="none", zorder=3))
    ax.text(x + 0.18, y + h - 0.13, title, ha="left", va="top", fontsize=9.6, fontweight="bold", color=INK,
            zorder=4)
    ax.text(x + 0.18, y + h - 0.42, body, ha="left", va="top", fontsize=8.1, color=INK2, zorder=4,
            linespacing=1.38)
    return (x, y, w, h)


def arrow(a, b, side_a="r", side_b="l", rad=0.0, color=INK2, frac_a=0.5, frac_b=0.5, ls="-"):
    def pt(bx, side, frac):
        x, y, w, h = bx
        return {"r": (x + w, y + h * frac), "l": (x, y + h * frac), "t": (x + w * frac, y + h),
                "b": (x + w * frac, y)}[side]
    ar = FancyArrowPatch(pt(a, side_a, frac_a), pt(b, side_b, frac_b), arrowstyle="-|>", mutation_scale=11,
                         lw=1.1, color=color, connectionstyle=f"arc3,rad={rad}", zorder=1, linestyle=ls,
                         shrinkA=2, shrinkB=2)
    ax.add_patch(ar)


# lane headers
TOP = 6.55
for x, lab, lane in [(0.25, "1  Data", "data"), (4.55, "2  Reconstruction", "recon"), (8.85, "3  Tests", "test")]:
    ax.text(x, TOP + 0.33, lab, fontsize=12, fontweight="bold", color=INK, va="center")
    ax.plot([x, x + 3.95], [TOP + 0.13, TOP + 0.13], color=ACCENT[lane], lw=2.2)


def hgt(n_lines):
    return 0.5 + 0.155 * n_lines


W = 3.95
prot = box(0.25, TOP - hgt(5), W, hgt(5), "Proteomes",
           "187 eukaryotes (2.97 M proteins)\n82 prokaryotes (0.41 M proteins)\nUniProt, EukProt, LukProt, NCBI;\n"
           "one protein per gene; BUSCO completeness;\n36 eukaryotic sets from transcriptomes", "data")
y = prot[1] - 0.3
tree = box(0.25, y - hgt(3), W, hgt(3), "Guide species trees",
           "NCBI taxonomy, edited to current phylogenies\n(62 logged operations, each with a reference);\n"
           "topology fixed in every analysis", "data")
y = tree[1] - 0.3
trans = box(0.25, y - hgt(4), W, hgt(4), "Transitions table",
            "state changes P \u2192 M with an eligibility rule:\n7 primary eukaryotic origins; 3 nested\n"
            "simple \u2192 complex and 4 prokaryotic (secondary);\n6 reversals; the rest exploratory", "data")

anc = box(4.55, TOP - hgt(3), W, hgt(3), "Ancestra (preset 'fast')",
          "protein families \u2192 within-family similarity \u2192\nbottom-up merge along the guide tree\n"
          "\u2192 hierarchical orthologous groups (HOGs)", "recon")
y = anc[1] - 0.4
ev = box(4.55, y - hgt(2), W, hgt(2), "Per-branch events",
         "gains, losses and duplications on each branch\n(parser exact on simulated true HOGs)", "recon")
y = ev[1] - 0.4
lad = box(4.55, y - hgt(4), W, hgt(4), "Gain-evidence ladder",
          "1  HOG gain: upper bound (2.1\u20132.6\u00d7 too many\n    on simulations)\n"
          "2  sequence-novel: no DIAMOND hit outside clade\n3  architecture-novel: Pfam layout unseen outside", "recon")
y = lad[1] - 0.4
ann = box(4.55, y - hgt(3), 1.9, hgt(3), "Annotation", "Pfam-A 38.2 domains\nmetapredict disorder\nprotein length",
          "recon")
phy = box(6.6, y - hgt(3), 1.9, hgt(3), "Phylogenomics", "250 BUSCO markers\nIQ-TREE LG+G4 on\nthe fixed topology",
          "recon")

tests = box(8.85, TOP - hgt(8), 3.95, hgt(8), "Stems vs control branches",
            "7 primary stems vs a pool of 30 unicellular\ncontrol branches (20 nearest per stem)\n"
            "H1/H2  sequence-novel gains (NB residual)\nH4       transcription-factor duplications\n"
            "H5       adhesion/ECM + signalling duplications\nH12     flagellar-module losses (direction fixed)\n"
            "nulls: permutation; parametric bootstrap\n(H1/H2); domain classes from Pfam", "test")
y = tests[1] - 0.4
h9 = box(8.85, y - hgt(3), 2.75, hgt(3), "Rate test (H9)",
         "per-gene log rate ratio,\nP\u2192M vs P\u2192sister tips;\n10 control nodes", "test")
y = h9[1] - 0.4
holm = box(8.85, y - 0.5, 3.95, 0.5, "Holm correction across the 5 primary tests", "", "test")
y = holm[1] - 0.3
expl = box(8.85, y - hgt(2), 3.95, hgt(2), "Exploratory (labelled as such)",
           "structural novelty (AlphaFold DB + Foldseek), gene-\nlevel rates, nested and prokaryotic transitions", "test",
           dashed=True)

arrow(prot, anc, frac_a=0.72, frac_b=0.5)
arrow(tree, anc, frac_a=0.75, frac_b=0.15)
arrow(anc, ev, "b", "t")
arrow(ev, lad, "b", "t")
arrow(ann, lad, "t", "b", frac_a=0.5, frac_b=0.24)
arrow(trans, lad, frac_a=0.55, frac_b=0.5)
arrow(ev, tests, frac_a=0.5, frac_b=0.62)
arrow(lad, tests, frac_a=0.75, frac_b=0.15)
arrow(phy, h9, "r", "l", frac_a=0.75, frac_b=0.4)
arrow(h9, holm, "b", "t", frac_a=0.5, frac_b=0.35)
ax.add_patch(FancyArrowPatch((12.2, tests[1]), (12.2, holm[1] + holm[3]), arrowstyle="-|>", mutation_scale=11,
                             lw=1.1, color=INK2, shrinkA=2, shrinkB=2))

# pre-registration bracket around the frozen inputs
yb = trans[1] - 0.12
ax.add_patch(FancyBboxPatch((0.12, yb), 4.22, TOP + 0.05 - yb, boxstyle="round,pad=0,rounding_size=0.12",
                            facecolor="none", edgecolor=PREREG, lw=1.3, ls=(0, (4, 3)), zorder=1))
ax.text(0.2, yb - 0.12, "Species sets, trees and transitions were frozen and\nhashed (with the design, test specifications and\n"
        "power table) before any eukaryote output existed",
        fontsize=8.1, color=PREREG, va="top", ha="left", fontweight="bold", linespacing=1.3)
ax.set_ylim(min(expl[1], phy[1], yb - 0.75) - 0.2, TOP + 0.6)
fig.set_size_inches(13, TOP + 0.8 - min(expl[1], phy[1], yb - 0.75))

for ext in ("png", "svg"):
    fig.savefig(os.path.join(P, "figures", f"fig2_workflow.{ext}"), dpi=200, facecolor="white")
print("ok")
