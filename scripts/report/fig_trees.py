#!/usr/bin/env python
"""Report figures 1 and S1: circular cladograms of the eukaryote and prokaryote guide trees.

Tips are coloured (and shaped) by multicellularity state from data/species_traits.tsv; the stems of the
pre-registered primary eukaryotic transitions (Fig 1) or the descriptive prokaryotic transitions (Fig S1)
from data/transitions.tsv are drawn thick and numbered; proteomes derived from transcriptomes carry an
asterisk. Topology only (cladogram): branch lengths are not drawn.

Usage: python scripts/report/fig_trees.py   (writes figures/fig1_euk_tree.{png,svg}, figures/figS1_prok_tree.{png,svg})
"""
from __future__ import annotations

import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Wedge  # noqa: E402
from ete3 import Tree  # noqa: E402

P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(P, "figures")

# Palette: reference categorical slots 1-3 (validated all-pairs, light mode) + a neutral for the
# baseline state; every state also has its own marker shape (secondary encoding).
INK = "#1b1b1b"
INK2 = "#52514e"
TREE = "#a9a7a1"
WEDGE = "#efede8"
STATE_STYLE = {  # key: (label, colour, marker, size)
    "unicellular": ("Unicellular", "#8f8d88", "o", 13),
    "simple": ("Multicellular, simple", "#2a78d6", "o", 26),
    "complex": ("Multicellular, complex", "#2a78d6", "s", 24),
    "grey": ("Grey zone (colonial, coenocytic,\nplasmodial, facultative)", "#eb6834", "D", 20),
    "reversal": ("Reversal to unicellularity", "#1baf7a", "v", 30),
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "svg.fonttype": "none",
    "font.size": 8,
})


def short_name(name: str) -> str:
    w = name.replace("Candidatus ", "Ca. ").split()
    pre = ""
    if w[0] == "Ca.":
        pre, w = "Ca. ", w[1:]
    if len(w) == 1:
        return pre + w[0]
    if w[1] in ("sp", "sp."):
        return pre + f"{w[0]} sp."
    return pre + f"{w[0][0]}. {w[1]}"


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def layout(tree: Tree, gap_deg: float = 14.0):
    leaves = tree.get_leaves()
    n = len(leaves)
    span = 2 * math.pi * (1 - gap_deg / 360.0)
    start = math.pi / 2 - math.radians(gap_deg) / 2  # gap centred at the top
    ang = {}
    for i, lf in enumerate(leaves):
        ang[lf] = start - span * (i + 0.5) / n  # clockwise
    height = {}
    for nd in tree.traverse("postorder"):
        height[nd] = 0 if nd.is_leaf() else 1 + max(height[c] for c in nd.children)
        if not nd.is_leaf():
            a = [ang[c] for c in nd.children]
            ang[nd] = (max(a) + min(a)) / 2
    D = height[tree]
    rad = {nd: (D - height[nd]) / D for nd in tree.traverse()}
    return ang, rad, D


def pol(r, a):
    return r * math.cos(a), r * math.sin(a)


def draw_tree(tree_path, traits_path, transitions_path, kingdom, sets, out_stem, title, clades,
              figsize, font_tip, R_label_pad, clade_r, legend_loc, show_reversal=True):
    t = Tree(tree_path, format=1)
    traits = {r["code"]: r for r in read_tsv(traits_path) if r["kingdom"] == kingdom}
    trans = [r for r in read_tsv(transitions_path) if r["kingdom"] == kingdom and r["analysis_set"] in sets]
    ang, rad, D = layout(t)
    nodes = {nd.name: nd for nd in t.traverse()}

    strip = 1.9 if figsize[0] >= 12 else 1.6
    fig = plt.figure(figsize=(figsize[0], figsize[1] + strip))
    fh = figsize[1] / (figsize[1] + strip)
    ax = fig.add_axes([0, 1 - fh, 1, fh])
    ax.set_aspect("equal")
    ax.axis("off")
    lim = 1 + clade_r + 0.18
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    # light wedges behind each highlighted multicellular crown
    stems = []
    for i, tr in enumerate(trans, 1):
        M = nodes.get(tr["M"])
        if M is None:
            continue
        path = [x for x in tr["stem_path"].split(">") if x]
        edges = [nodes[x] for x in path[1:]] if len(path) > 1 else [M]
        stems.append((i, tr, M, edges))
        lv = M.get_leaves()
        a = [ang[x] for x in lv]
        half = (math.pi * 2 * (1 - 14 / 360) / len(t)) / 2
        r_in = rad[M] - 0.01
        w = Wedge((0, 0), 1.03, math.degrees(min(a) - half), math.degrees(max(a) + half),
                  width=1.03 - r_in, facecolor=WEDGE, edgecolor="none", zorder=0)
        ax.add_patch(w)

    # branches
    for nd in t.traverse():
        if nd.is_root():
            continue
        r0, r1, a = rad[nd.up], rad[nd], ang[nd]
        ax.plot(*zip(pol(r0, a), pol(r1, a)), color=TREE, lw=0.7, zorder=1, solid_capstyle="butt")
        if not nd.is_leaf():
            pass
    for nd in t.traverse():
        if nd.is_leaf():
            continue
        a = [ang[c] for c in nd.children]
        th = [min(a) + (max(a) - min(a)) * k / 60 for k in range(61)]
        ax.plot([rad[nd] * math.cos(x) for x in th], [rad[nd] * math.sin(x) for x in th],
                color=TREE, lw=0.7, zorder=1)

    # highlighted stems (edges into each node on the stem path)
    for i, tr, M, edges in stems:
        for e in edges:
            ax.plot(*zip(pol(rad[e.up], ang[e]), pol(rad[e], ang[e])), color=INK, lw=3.6, zorder=3,
                    solid_capstyle="butt")
        # badge beside the stem (inward of P along the parent's direction would hide other branches)
        r_mid = (rad[edges[0].up] + rad[edges[-1]]) / 2
        a_mid = ang[edges[len(edges) // 2]] + 0.075 / max(r_mid, 0.1)
        x, y = pol(r_mid, a_mid)
        ax.scatter([x], [y], s=135, facecolor="white", edgecolor=INK, linewidth=1.2, zorder=4)
        ax.text(x, y, str(i), ha="center", va="center", fontsize=7.5, fontweight="bold", color=INK, zorder=5)

    # tips
    n_tx = 0
    for lf in t.get_leaves():
        tr = traits.get(lf.name, {})
        state = "reversal" if tr.get("reversal") == "1" else tr.get("mc_state", "unicellular")
        lab, col, mk, sz = STATE_STYLE[state]
        x, y = pol(1.0, ang[lf])
        ax.scatter([x], [y], s=sz * (figsize[0] / 12) ** 0, marker=mk, color=col, edgecolor="white",
                   linewidth=0.4, zorder=6)
        name = short_name(tr.get("species", lf.name))
        if tr.get("data_type") == "transcriptome":
            name += "*"
            n_tx += 1
        a = ang[lf]
        deg = math.degrees(a)
        x, y = pol(1.0 + R_label_pad * 0.25, a)
        flip = math.cos(a) < 0
        ax.text(x, y, name, rotation=deg + (180 if flip else 0), rotation_mode="anchor",
                ha="right" if flip else "left", va="center", fontsize=font_tip, style="italic", color=INK, zorder=6)

    # clade arcs
    for cname, label in clades:
        nd = nodes.get(cname)
        if nd is None:
            continue
        a = [ang[x] for x in nd.get_leaves()]
        half = (math.pi * 2 * (1 - 14 / 360) / len(t)) / 2 * 0.7
        a0, a1 = min(a) - half, max(a) + half
        rr = 1 + clade_r
        th = [a0 + (a1 - a0) * k / 80 for k in range(81)]
        ax.plot([rr * math.cos(x) for x in th], [rr * math.sin(x) for x in th], color=INK2, lw=1.4,
                solid_capstyle="butt")
        am = (a0 + a1) / 2
        x, y = pol(rr + 0.045, am)
        deg = math.degrees(am) - 90
        if math.sin(am) < 0:
            deg += 180
        ax.text(x, y, label, rotation=deg, rotation_mode="anchor", ha="center", va="center", fontsize=8.5,
                fontweight="bold", color=INK2)

    # legends
    handles = []
    for key in ["unicellular", "simple", "complex", "grey", "reversal"]:
        if key == "reversal" and not show_reversal:
            continue
        lab, col, mk, sz = STATE_STYLE[key]
        handles.append(Line2D([0], [0], marker=mk, color="none", markerfacecolor=col, markeredgecolor="white",
                              markersize=math.sqrt(sz) * 1.35, label=lab))
    handles.append(Line2D([0], [0], color=INK, lw=3, label="Transition stem (numbered)"))
    if n_tx:
        handles.append(Line2D([0], [0], color="none", label="* proteome from a transcriptome"))
    leg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 1 - fh + 0.005), frameon=False,
                     fontsize=8, title="Tip state", title_fontsize=8.5, alignment="left", handlelength=1.6)
    leg.get_title().set_fontweight("bold")

    key_lines = [f"{i}  {short_label(tr)}" for i, tr, M, e in stems]
    fig.text(legend_loc[2][0], 1 - fh - 0.012, legend_loc[3], fontsize=8.5, fontweight="bold", va="top",
             ha="left", color=INK)
    fig.text(legend_loc[2][0], 1 - fh - 0.012 - 0.022 * 12 / (figsize[1] + strip), "\n".join(key_lines),
             fontsize=8, va="top", ha="left", color=INK, linespacing=1.5)
    fig.text(0.015, 0.992, title, fontsize=11, fontweight="bold", va="top", ha="left", color=INK)

    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{out_stem}.{ext}"), dpi=200, facecolor="white")
    plt.close(fig)
    return len(t), n_tx, len(stems)


SHORT = {
    "T_metazoa": "Animals (Metazoa)",
    "T_fungi_hyphae": "Fungal hyphae (path stem to Amastigomycota)",
    "T_dictyostelia": "Dictyostelids (social amoebae; aggregative)",
    "T_streptophyta": "Streptophyte algae + land plants",
    "T_rhodophyta": "Red algae (Bangiophyceae + Florideophyceae)",
    "T_px_filaments": "Brown-algal lineage (PX clade: brown, yellow-green algae and kin)",
    "T_volvocine_tgv": "Volvocine colonies (path stem to Gonium + Volvocaceae)",
    "T_myxobacteria": "Myxobacteria (aggregative fruiting bodies)",
    "T_streptomyces": "Streptomycetes (mycelium)",
    "T_magnetotactic_mmp": "Multicellular magnetotactic bacteria",
    "T_methanosarcina": "Methanosarcina (cell packets)",
}


def short_label(tr):
    return SHORT.get(tr["transition_id"], tr["label"])


def main():
    os.makedirs(OUT, exist_ok=True)
    euk_clades = [
        ("Metazoa", "Animals"), ("Choanoflagellata", "Choanoflagellates"), ("Teretosporea", "Ichthyosporeans + kin"),
        ("Fungi", "Fungi"), ("Amoebozoa", "Amoebozoa"), ("Rhodophyta", "Red algae"),
        ("Chlorophyta", "Chlorophyta"), ("Streptophyta", "Streptophyta"), ("Stramenopiles", "Stramenopiles"),
        ("Alveolata", "Alveolates"), ("Rhizaria", "Rhizaria"), ("Discoba", "Discoba"),
    ]
    n, ntx, ns = draw_tree(
        os.path.join(P, "data/species_tree_bl_euk.nwk"), os.path.join(P, "data/species_traits.tsv"),
        os.path.join(P, "data/transitions.tsv"), "euk", {"primary"}, "fig1_euk_tree",
        "Eukaryote guide tree: 187 proteomes and 7 primary transitions",
        euk_clades, (12, 12), 5.3, 0.62, 0.40,
        ("lower left", (0.01, 0.01), (0.42, 0), "Primary transition stems (P \u2192 M)"))
    print("euk", n, "tips", ntx, "transcriptomes", ns, "stems")
    prok_clades = [
        ("Archaea", "Archaea"), ("Cyanophyceae", "Cyanobacteria"), ("Actinomycetota", "Actinomycetota"),
        ("Myxococcota", "Myxococcota"), ("Thermodesulfobacteriota", "Desulfobacterota"),
        ("Pseudomonadota", "Pseudomonadota"),
    ]
    n, ntx, ns = draw_tree(
        os.path.join(P, "data/species_tree_bl_prok.nwk"), os.path.join(P, "data/species_traits.tsv"),
        os.path.join(P, "data/transitions.tsv"), "prok", {"secondary"}, "figS1_prok_tree",
        "Prokaryote guide tree: 82 proteomes",
        prok_clades, (9, 9), 6.0, 0.62, 0.56,
        ("lower left", (0.01, 0.01), (0.45, 0), "Transition stems (descriptive only)"), show_reversal=False)
    print("prok", n, "tips", ntx, "transcriptomes", ns, "stems")


if __name__ == "__main__":
    main()
