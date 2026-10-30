#!/usr/bin/env python
"""Draw the schematic diagrams and figure crops used by the slide deck.

Output: report/slides/img/*.png. Run with the slides env:
    conda_envs/slides/bin/python scripts/report/slides/make_diagrams.py
All numbers drawn here are copied from report/REPORT.md (cited per function).
Schematics marked "schematic" use illustrative, not real, data.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import (Circle, FancyBboxPatch, FancyArrowPatch, Rectangle,
                                Polygon, Wedge, Arc, Ellipse)
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "report" / "slides" / "img"
OUT.mkdir(parents=True, exist_ok=True)

# fonts: Carlito is metric-compatible with Calibri (deck body font)
for f in (Path.home() / ".local/share/fonts").glob("*.ttf"):
    font_manager.fontManager.addfont(str(f))
FONT = "Carlito" if any("Carlito" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
plt.rcParams.update({"font.family": FONT, "font.size": 15, "axes.edgecolor": "#5C6B6A",
                     "axes.labelcolor": "#1F2A2A", "xtick.color": "#5C6B6A", "ytick.color": "#5C6B6A",
                     "savefig.dpi": 200, "savefig.transparent": False})

# palette (same as deck THEME; state colours follow report/figures_v2)
DARK, TEAL, TEXT, MUTED = "#0E3B39", "#1B7F79", "#1F2A2A", "#5C6B6A"
LIGHT = "#E6F1EF"
UNI, SIMPLE, COMPLEX, REV, CONF = "#9AA5A3", "#009E73", "#0072B2", "#E69F00", "#CC79A7"
ACCENT = "#E4572E"        # warm accent
GAIN, LOSS, DUP = "#009E73", "#C8442F", "#0072B2"


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)
    print("wrote", name)


def clean(ax):
    ax.set_axis_off()
    ax.set_aspect("equal")


def star(ax, x, y, s=260, c=GAIN, z=6):
    ax.scatter([x], [y], marker="*", s=s, c=c, edgecolors="white", linewidths=0.8, zorder=z)


def xmark(ax, x, y, s=0.18, c=LOSS, lw=3.5, z=6):
    ax.plot([x - s, x + s], [y - s, y + s], c=c, lw=lw, zorder=z, solid_capstyle="round")
    ax.plot([x - s, x + s], [y + s, y - s], c=c, lw=lw, zorder=z, solid_capstyle="round")


def gene(ax, x, y, c, w=0.32, h=0.22, z=5, ec="white"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.05",
                                fc=c, ec=ec, lw=1, zorder=z))


# ----------------------------------------------------------------------------------------------
# motif: tree glyph in a circle
def tree_glyph(ax, cx, cy, r, fg="white", lw=None):
    lw = lw or r * 9
    s = r * 0.55
    pts = [(-1, 0), (-0.3, 0)]
    ax.plot([cx - s, cx - 0.2 * s], [cy, cy], c=fg, lw=lw, solid_capstyle="round")
    ax.plot([cx - 0.2 * s, cx - 0.2 * s], [cy - 0.75 * s, cy + 0.75 * s], c=fg, lw=lw, solid_capstyle="round")
    for yy in (-0.75, 0.75):
        ax.plot([cx - 0.2 * s, cx + 0.35 * s], [cy + yy * s, cy + yy * s], c=fg, lw=lw, solid_capstyle="round")
    ax.plot([cx + 0.35 * s, cx + 0.35 * s], [cy + 0.35 * s, cy + 1.1 * s], c=fg, lw=lw, solid_capstyle="round")
    for yy in (0.35, 1.1):
        ax.plot([cx + 0.35 * s, cx + 0.95 * s], [cy + yy * s, cy + yy * s], c=fg, lw=lw, solid_capstyle="round")
    ax.plot([cx + 0.35 * s, cx + 0.95 * s], [cy - 0.75 * s, cy - 0.75 * s], c=fg, lw=lw, solid_capstyle="round")
    for yy, col in ((1.1, fg), (0.35, fg), (-0.75, fg)):
        ax.add_patch(Circle((cx + 0.95 * s + 0.08 * s, cy + yy * s), 0.14 * s, fc=col, ec="none"))


def motif():
    for name, bg, fg in (("motif_teal.png", TEAL, "white"), ("motif_light.png", "#CFE7E3", DARK),
                         ("motif_white.png", "white", DARK)):
        fig, ax = plt.subplots(figsize=(2, 2))
        ax.add_patch(Circle((0, 0), 1, fc=bg, ec="none"))
        tree_glyph(ax, -0.05, -0.1, 1, fg=fg, lw=7)
        ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); clean(ax)
        fig.savefig(OUT / name, transparent=True, dpi=200, bbox_inches="tight", pad_inches=0)
        plt.close(fig)
        print("wrote", name)


# ----------------------------------------------------------------------------------------------
# mini tree used by every hypothesis slide
TIPS = {"T0": 0.0, "T1": 1.0, "T2": 2.0, "T3": 3.0, "T4": 4.0}


def elbow(ax, x0, y0, x1, y1, c, lw, ls="-", z=2):
    ax.plot([x0, x0, x1], [y0, y1, y1], c=c, lw=lw, ls=ls, solid_capstyle="round", zorder=z)


def mini_tree(ax, stem=COMPLEX, crown=COMPLEX, nested=False, pre=False, crown_long=False,
              stem_label=None, tipcolors=None):
    """Cladogram: root -> T0 (outgroup), N1 -> T1 (sister), stem N1->N2, crown N2 -> T2, N3 -> T3,T4."""
    g, lw = "#8C9796", 3
    xr, x1, x2, x3, xt = 0.0, 1.0, 2.4, 3.2, 4.4 if not crown_long else 5.0
    y_n3 = 3.5; y_n2 = (2 + y_n3) / 2; y_n1 = (1 + y_n2) / 2; y_r = (0 + y_n1) / 2
    elbow(ax, xr, y_r, xt, 0, g, lw)
    elbow(ax, xr, y_r, x1, y_n1, ACCENT if pre else g, lw + (2 if pre else 0))
    elbow(ax, x1, y_n1, xt, 1, g, lw)
    # stem
    if nested:
        ax.plot([x1, x1], [y_n1, y_n2], c=g, lw=lw)
        ax.plot([x1, (x1 + x2) / 2], [y_n2, y_n2], c=SIMPLE, lw=9, solid_capstyle="butt")
        ax.plot([(x1 + x2) / 2, x2], [y_n2, y_n2], c=COMPLEX, lw=9, solid_capstyle="butt")
    else:
        ax.plot([x1, x1], [y_n1, y_n2], c=g, lw=lw)
        ax.plot([x1, x2], [y_n2, y_n2], c=stem, lw=9, solid_capstyle="butt")
    cl = 4.5 if crown_long else lw
    elbow(ax, x2, y_n2, xt, 2, crown, cl)
    elbow(ax, x2, y_n2, x3, y_n3, crown, cl)
    elbow(ax, x3, y_n3, xt, 3, crown, cl)
    elbow(ax, x3, y_n3, xt, 4, crown, cl)
    cols = tipcolors or [UNI, UNI, crown, crown, crown]
    for i, c in enumerate(cols):
        ax.add_patch(Circle((xt + 0.22, i), 0.17, fc=c, ec="white", lw=1.2, zorder=5))
    if stem_label:
        ax.text((x1 + x2) / 2, y_n2 + 0.35, stem_label, ha="center", va="bottom", fontsize=13, color=TEXT)
    return dict(stem=((x1 + x2) / 2, y_n2), xt=xt, pre=((xr + x1) / 2, y_n1))


def glyph_circle(ax, cx=7.4, cy=2.0, r=1.95):
    ax.add_patch(Circle((cx, cy), r, fc=LIGHT, ec="none", zorder=0))
    return cx, cy


def hyp_canvas():
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    clean(ax)
    ax.set_xlim(-0.3, 9.5); ax.set_ylim(-0.7, 4.7)
    return fig, ax


# glyph library (centred at cx, cy, radius ~1.6)
def g_stars(ax, cx, cy, n=5):
    rng = np.random.default_rng(3)
    for i in range(n):
        a = 2 * np.pi * i / n + 0.3
        star(ax, cx + 0.9 * np.cos(a), cy + 0.9 * np.sin(a), s=420)
    gene(ax, cx, cy, UNI, w=0.6, h=0.4)


def g_genes_present(ax, cx, cy):
    for i, c in enumerate([TEAL, COMPLEX, SIMPLE, TEAL]):
        gene(ax, cx - 0.75 + 0.5 * i, cy + 0.5, c, w=0.4, h=0.3)
    ax.annotate("", (cx, cy - 0.45), (cx, cy + 0.2), arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=2))
    for i, c in enumerate([TEAL, COMPLEX, SIMPLE, TEAL]):
        gene(ax, cx - 0.75 + 0.5 * i, cy - 0.8, c, w=0.4, h=0.3)
    ax.text(cx, cy + 1.05, "already there", ha="center", fontsize=12, color=MUTED)
    ax.text(cx, cy - 1.45, "re-used", ha="center", fontsize=12, color=MUTED)


def g_tf(ax, cx, cy):
    t = np.linspace(-1.2, 1.2, 200)
    ax.plot(cx + t, cy - 0.35 + 0.25 * np.sin(5 * t), c=TEXT, lw=2.5)
    ax.plot(cx + t, cy - 0.35 - 0.25 * np.sin(5 * t), c=MUTED, lw=2.5)
    for dx in (-0.45, 0.45):
        ax.add_patch(Ellipse((cx + dx, cy + 0.35), 0.75, 0.55, fc=COMPLEX, ec="white", lw=1.5, zorder=5))
    ax.text(cx, cy + 0.95, "TF  ×2", ha="center", fontsize=15, weight="bold", color=COMPLEX)


def g_adhesion(ax, cx, cy):
    for dx in (-0.65, 0.65):
        ax.add_patch(Circle((cx + dx, cy), 0.6, fc="white", ec=TEAL, lw=3, zorder=3))
        ax.add_patch(Circle((cx + dx, cy), 0.2, fc=TEAL, ec="none", zorder=4))
    for dy in (-0.3, 0, 0.3):
        ax.plot([cx - 0.12, cx + 0.12], [cy + dy, cy + dy], c=ACCENT, lw=4, zorder=5)
    ax.text(cx, cy + 0.95, "glue + signal", ha="center", fontsize=13, color=MUTED)


def g_kinase(ax, cx, cy):
    ax.add_patch(Rectangle((cx - 1.3, cy - 0.15), 2.6, 0.3, fc="#F2D3B8", ec="none", zorder=2))
    ax.plot([cx, cx], [cy - 0.25, cy + 0.6], c=TEXT, lw=3, zorder=3)
    ax.plot([cx, cx - 0.35], [cy + 0.6, cy + 1.0], c=TEXT, lw=3, zorder=3)
    ax.plot([cx, cx + 0.35], [cy + 0.6, cy + 1.0], c=TEXT, lw=3, zorder=3)
    ax.add_patch(Ellipse((cx, cy - 0.65), 0.9, 0.6, fc=COMPLEX, ec="white", lw=1.5, zorder=4))
    ax.text(cx, cy - 0.68, "P", ha="center", va="center", color="white", fontsize=14, weight="bold", zorder=5)
    ax.text(cx + 0.9, cy + 0.85, "out", fontsize=11, color=MUTED)
    ax.text(cx + 0.9, cy - 0.85, "in", fontsize=11, color=MUTED)


def g_aggregate(ax, cx, cy):
    # left: aggregation (arrows inward); right: clonal chain
    pts = [(-1.15, 0.55), (-1.25, -0.4), (-0.45, 0.85), (-0.5, -0.85)]
    for x, y in pts:
        ax.add_patch(Circle((cx + x, cy + y), 0.15, fc=UNI, ec="none"))
        ax.annotate("", (cx - 0.7, cy), (cx + x * 0.85, cy + y * 0.85),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2, shrinkA=8, shrinkB=10))
    ax.add_patch(Circle((cx - 0.7, cy), 0.25, fc=SIMPLE, ec="none"))
    for i in range(4):
        ax.add_patch(Circle((cx + 0.35 + 0.3 * i, cy), 0.15, fc=COMPLEX, ec="white", lw=1))
    ax.text(cx - 0.75, cy - 1.35, "aggregative", ha="center", fontsize=11, color=MUTED)
    ax.text(cx + 0.8, cy - 1.35, "clonal", ha="center", fontsize=11, color=MUTED)
    ax.plot([cx, cx], [cy - 1.0, cy + 1.0], c="white", lw=3)


def g_rate(ax, cx, cy):
    ax.add_patch(Wedge((cx, cy - 0.3), 1.2, 0, 180, width=0.25, fc="#CFE0DD", ec="none"))
    ax.add_patch(Wedge((cx, cy - 0.3), 1.2, 0, 60, width=0.25, fc=ACCENT, ec="none"))
    a = np.deg2rad(40)
    ax.plot([cx, cx + 0.95 * np.cos(a)], [cy - 0.3, cy - 0.3 + 0.95 * np.sin(a)], c=TEXT, lw=4, solid_capstyle="round")
    ax.add_patch(Circle((cx, cy - 0.3), 0.12, fc=TEXT))
    ax.text(cx, cy - 0.85, "faster protein\nevolution?", ha="center", va="top", fontsize=12, color=MUTED)


def g_disorder(ax, cx, cy):
    ax.add_patch(FancyBboxPatch((cx - 1.1, cy - 0.3), 1.0, 0.6, boxstyle="round,pad=0.05", fc=COMPLEX, ec="none"))
    t = np.linspace(0, 1, 200)
    ax.plot(cx - 0.05 + 1.15 * t, cy + 0.35 * np.sin(18 * t) * np.exp(-0.5 * t), c=ACCENT, lw=3)
    ax.text(cx - 0.6, cy - 0.75, "folded", ha="center", fontsize=11, color=MUTED)
    ax.text(cx + 0.6, cy - 0.75, "disordered", ha="center", fontsize=11, color=MUTED)


def g_domains(ax, cx, cy):
    cols = [COMPLEX, SIMPLE, ACCENT]
    for i, c in enumerate(cols):
        gene(ax, cx - 0.9 + i * 0.9, cy + 0.6, c, w=0.6, h=0.35)
    ax.annotate("", (cx, cy - 0.2), (cx, cy + 0.3), arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=2))
    for i, c in enumerate([ACCENT, COMPLEX, SIMPLE]):
        gene(ax, cx - 0.6 + i * 0.6, cy - 0.55, c, w=0.6, h=0.35, ec="white")
    ax.text(cx, cy - 1.15, "new combination", ha="center", fontsize=12, color=MUTED)


def g_flagellum(ax, cx, cy):
    ax.add_patch(Ellipse((cx - 0.4, cy), 1.0, 0.75, fc=UNI, ec="none"))
    t = np.linspace(0, 1, 100)
    ax.plot(cx + 0.1 + 1.1 * t, cy + 0.2 * np.sin(14 * t), c=TEXT, lw=2.5)
    xmark(ax, cx + 0.65, cy, s=0.35, lw=5)
    ax.text(cx, cy - 0.9, "flagellum lost", ha="center", fontsize=12, color=MUTED)


def g_cellcycle(ax, cx, cy):
    ax.add_patch(Arc((cx, cy), 2.3, 2.3, theta1=20, theta2=330, color=TEAL, lw=3))
    ax.annotate("", (cx + 1.08, cy + 0.38), (cx + 1.13, cy + 0.15),
                arrowprops=dict(arrowstyle="-|>", color=TEAL, lw=3, mutation_scale=20))
    ax.add_patch(Circle((cx - 0.3, cy), 0.35, fc=COMPLEX, ec="white", lw=1.5))
    ax.add_patch(Circle((cx + 0.3, cy), 0.35, fc=COMPLEX, ec="white", lw=1.5))
    ax.plot([cx - 0.05, cx + 0.05], [cy, cy], c="white", lw=5)
    ax.text(cx, cy - 0.75, "RB / ESCRT", ha="center", fontsize=12, color=MUTED)


def g_hgt(ax, cx, cy):
    ax.add_patch(FancyBboxPatch((cx - 1.25, cy + 0.35), 0.9, 0.4, boxstyle="round,pad=0.02,rounding_size=0.2",
                                fc=UNI, ec="none"))
    ax.text(cx - 0.8, cy + 0.98, "bacterium", ha="center", fontsize=11, color=MUTED)
    ax.add_patch(Circle((cx + 0.6, cy + 0.5), 0.45, fc=COMPLEX, ec="white", lw=1.5))
    ax.annotate("", (cx + 0.12, cy + 0.5), (cx - 0.3, cy + 0.55),
                arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=2.5))
    # site convergence: two aligned rows with matching column
    for r, y in enumerate((cy - 0.45, cy - 0.85)):
        for i, ch in enumerate("AKLVE"):
            col = ACCENT if i == 2 else "#D8E6E3"
            ax.add_patch(Rectangle((cx - 1.0 + 0.4 * i, y - 0.17), 0.36, 0.34, fc=col, ec="none"))
            ax.text(cx - 0.82 + 0.4 * i, y, "W" if i == 2 else ch, ha="center", va="center", fontsize=11,
                    color="white" if i == 2 else TEXT)


def g_regulatory(ax, cx, cy):
    ax.plot([cx - 1.4, cx + 1.4], [cy, cy], c=TEXT, lw=2)
    ax.add_patch(Rectangle((cx - 1.2, cy - 0.15), 0.6, 0.3, fc=ACCENT, ec="none"))
    ax.add_patch(Polygon([[cx - 0.1, cy - 0.25], [cx + 0.9, cy - 0.25], [cx + 1.25, cy], [cx + 0.9, cy + 0.25],
                          [cx - 0.1, cy + 0.25]], fc=COMPLEX, ec="none"))
    ax.text(cx - 0.9, cy + 0.4, "switch", ha="center", fontsize=11, color=MUTED)
    ax.text(cx + 0.55, cy + 0.4, "gene", ha="center", fontsize=11, color=MUTED)
    ax.text(cx, cy - 0.9, "invisible in\nproteomes", ha="center", va="center", fontsize=12, color=MUTED)


def g_motif(ax, cx, cy):
    ax.add_patch(Wedge((cx - 0.2, cy), 0.85, 30, 330, fc=COMPLEX, ec="none"))
    xs = np.linspace(cx + 0.15, cx + 1.3, 5)
    ax.plot(xs, cy + 0.08 * np.sin(np.arange(5) * 2), c=ACCENT, lw=5, solid_capstyle="round")
    for x in xs:
        ax.add_patch(Circle((x, cy), 0.08, fc="white", ec=ACCENT, lw=1.5, zorder=5))
    ax.text(cx, cy - 1.15, "binds a short motif", ha="center", fontsize=12, color=MUTED)


def g_celldeath(ax, cx, cy):
    ax.add_patch(Circle((cx - 0.15, cy + 0.05), 0.75, fc="#F6D5CC", ec=LOSS, lw=2.5, ls="--"))
    for (dx, dy, r) in ((0.85, 0.65, 0.15), (0.95, -0.2, 0.12), (0.7, -0.75, 0.1)):
        ax.add_patch(Circle((cx + dx, cy + dy), r, fc="#F6D5CC", ec=LOSS, lw=1.5))
    ax.add_patch(Circle((cx - 0.15, cy + 0.05), 0.25, fc=LOSS, ec="none"))
    ax.text(cx, cy - 1.2, "programmed cell death", ha="center", fontsize=12, color=MUTED)


def g_structure(ax, cx, cy):
    t = np.linspace(0, 1, 300)
    ax.plot(cx - 1.1 + 2.2 * t, cy + 0.35 + 0.3 * np.sin(30 * t), c=COMPLEX, lw=4)
    ax.text(cx, cy - 0.45, "D V L P C Q", ha="center", fontsize=15, family="monospace", color=ACCENT, weight="bold")
    ax.text(cx, cy - 1.05, "3Di structure letters", ha="center", fontsize=11, color=MUTED)


def g_complex(ax, cx, cy):
    pos = [(-0.35, 0.3), (0.3, 0.35), (0, -0.25), (-0.6, -0.35), (0.6, -0.3)]
    cols = [COMPLEX, SIMPLE, TEAL, ACCENT, UNI]
    for (x, y), c in zip(pos, cols):
        ax.add_patch(Circle((cx + x, cy + y), 0.38, fc=c, ec="white", lw=2))
    ax.text(cx, cy - 1.15, "protein complex", ha="center", fontsize=12, color=MUTED)


def g_tip(ax, cx, cy):
    star(ax, cx - 0.5, cy + 0.4, s=380)
    star(ax, cx + 0.5, cy + 0.4, s=380, c=COMPLEX)
    ax.text(cx, cy - 0.35, "new genes at\nsimple → complex?", ha="center", va="top", fontsize=12, color=MUTED)


HYPS = {
    # name: (tree kwargs, glyph fn, stem marks)
    "h01_h02": (dict(stem=COMPLEX, crown=COMPLEX), g_stars, "stars"),
    "h03": (dict(stem=COMPLEX, crown=COMPLEX, pre=True), g_genes_present, "pre"),
    "h04": (dict(stem=COMPLEX, crown=COMPLEX), g_tf, "dup"),
    "h05": (dict(stem=SIMPLE, crown=SIMPLE), g_adhesion, "dup"),
    "h06": (dict(stem=COMPLEX, crown=COMPLEX), g_kinase, "dup"),
    "h07_h08": (dict(stem=SIMPLE, crown=SIMPLE), g_aggregate, None),
    "h09": (dict(stem=COMPLEX, crown=COMPLEX, crown_long=True), g_rate, None),
    "h10": (dict(stem=COMPLEX, crown=COMPLEX), g_disorder, None),
    "h11": (dict(stem=COMPLEX, crown=COMPLEX), g_domains, "stars"),
    "h12": (dict(stem=SIMPLE, crown=SIMPLE), g_flagellum, "loss"),
    "h13": (dict(stem=SIMPLE, crown=SIMPLE), g_cellcycle, None),
    "h14_h15": (dict(stem=COMPLEX, crown=COMPLEX), g_hgt, None),
    "h16": (dict(stem=COMPLEX, crown=COMPLEX), g_regulatory, None),
    "r2_p1": (dict(nested=True, crown=COMPLEX), g_tip, "stars2"),
    "r2_p2": (dict(stem=COMPLEX, crown=COMPLEX), g_motif, "dup"),
    "r2_p3": (dict(stem=COMPLEX, crown=COMPLEX), g_celldeath, "dup"),
    "r2_struct": (dict(stem=COMPLEX, crown=COMPLEX), g_structure, None),
    "r2_complex": (dict(stem=COMPLEX, crown=COMPLEX), g_complex, "dup"),
}


def hypotheses():
    for name, (tk, glyph, mark) in HYPS.items():
        fig, ax = hyp_canvas()
        pos = mini_tree(ax, **tk)
        sx, sy = pos["stem"]
        if mark == "stars":
            for dx in (-0.4, 0.0, 0.4):
                star(ax, sx + dx, sy + 0.35, s=180)
        elif mark == "stars2":
            star(ax, sx + 0.35, sy + 0.35, s=200)
        elif mark == "dup":
            ax.text(sx, sy + 0.3, "×2", ha="center", va="bottom", fontsize=16, weight="bold", color=DUP)
        elif mark == "loss":
            xmark(ax, sx, sy + 0.42, s=0.17)
        elif mark == "pre":
            px, py = pos["pre"]
            for dx in (-0.25, 0.25):
                gene(ax, px + 0.15 + dx, py + 0.32, ACCENT, w=0.3, h=0.2)
        cx, cy = glyph_circle(ax)
        glyph(ax, cx, cy)
        save(fig, f"hyp_{name}.png")


# ----------------------------------------------------------------------------------------------
# E1: what is ancestral gene content
def explain_ancestral():
    fig, ax = plt.subplots(figsize=(10, 5.2)); clean(ax)
    ax.set_xlim(-0.5, 12.5); ax.set_ylim(-1.3, 5.6)
    g = "#7E8B8A"
    # tree: root(0,2.25) -> A(1.8, 0.75) [t0 y0, t1 y1.5], B(1.8,3.75)->[t2 y3, C(3.6,4.5) -> t3 y4, t4 y5]
    xt = 6.0
    elbow(ax, 0, 2.6, 2.0, 0.75, g, 3); elbow(ax, 0, 2.6, 2.0, 4.0, g, 3)
    elbow(ax, 2.0, 0.75, xt, 0, g, 3); elbow(ax, 2.0, 0.75, xt, 1.5, g, 3)
    elbow(ax, 2.0, 4.0, xt, 3.0, g, 3); elbow(ax, 2.0, 4.0, 4.0, 5.0, g, 3)
    elbow(ax, 4.0, 5.0, xt, 4.5, g, 3); elbow(ax, 4.0, 5.0, xt, 5.5, g, 3)
    fams = [TEAL, COMPLEX, ACCENT, SIMPLE]
    tips = {0: [1, 1, 0, 0], 1.5: [1, 1, 0, 0], 3.0: [1, 2, 1, 0], 4.5: [1, 2, 1, 1], 5.5: [1, 2, 0, 1]}
    names = {0: "species A", 1.5: "species B", 3.0: "species C", 4.5: "species D", 5.5: "species E"}
    for y, counts in tips.items():
        x = xt + 0.4
        for f, n in enumerate(counts):
            for k in range(n):
                gene(ax, x, y, fams[f], w=0.32, h=0.26)
                x += 0.38
        ax.text(xt + 2.75, y, names[y], va="center", fontsize=14, color=MUTED)
    # ancestor box
    def anc(x, y, counts, label):
        ax.add_patch(FancyBboxPatch((x - 0.15, y - 0.28), 0.38 * sum(counts) + 0.12, 0.56,
                                    boxstyle="round,pad=0.03", fc="white", ec=DARK, lw=1.5, zorder=4))
        xx = x
        for f, n in enumerate(counts):
            for k in range(n):
                gene(ax, xx + 0.04, y, fams[f], w=0.32, h=0.26, z=6); xx += 0.38
        ax.text(x - 0.25, y, label, ha="right", va="center", fontsize=12, color=DARK)
    anc(0.25, 2.6, [1, 1, 0, 0], "root\nancestor")
    # events (family 2 duplicated and family 3 gained on the C-D-E stem; family 4 gained on D-E; family 3 lost in E)
    ax.text(0.9, 4.3, "×2", fontsize=18, color=DUP, weight="bold", ha="center", va="center")
    ax.text(0.8, 3.65, "duplication", fontsize=13, color=DUP, ha="center")
    ax.text(1.65, 4.3, "+", fontsize=26, color=ACCENT, weight="bold", ha="center", va="center")
    ax.text(1.65, 4.75, "gain", fontsize=13, color=ACCENT, ha="center")
    ax.text(3.0, 5.3, "+", fontsize=26, color=SIMPLE, weight="bold", ha="center", va="center")
    ax.text(3.0, 4.62, "gain", fontsize=13, color=SIMPLE, ha="center")
    xmark(ax, 5.0, 5.5, s=0.14)
    ax.text(5.0, 5.85, "loss", fontsize=13, color=LOSS, ha="center")
    # legend
    for i, (c, l) in enumerate(zip(fams, ["family 1", "family 2", "family 3", "family 4"])):
        gene(ax, 0.3 + 2.1 * i, -0.9, c); ax.text(0.6 + 2.1 * i, -0.9, l, va="center", fontsize=13, color=MUTED)
    save(fig, "explain_ancestral.png")


# E2: Ancestra simplified
def explain_ancestra():
    fig, ax = plt.subplots(figsize=(12, 3.6)); clean(ax)
    ax.set_xlim(0, 15); ax.set_ylim(0, 4.9)
    rng = np.random.default_rng(1)
    cols = [TEAL, COMPLEX, ACCENT]
    titles = ["1  Sort proteins\n    into families", "2  Compare only\n    inside a family", "3  Walk the tree\n    tips → root"]
    for p in range(3):
        x0 = 0.3 + p * 5.0
        ax.add_patch(FancyBboxPatch((x0, 0.3), 4.3, 3.2, boxstyle="round,pad=0.05,rounding_size=0.25",
                                    fc=LIGHT, ec="none"))
        ax.text(x0 + 0.2, 4.15, titles[p], fontsize=16, weight="bold", color=DARK, va="center")
        if p < 2:
            ax.annotate("", (x0 + 4.75, 1.9), (x0 + 4.4, 1.9),
                        arrowprops=dict(arrowstyle="-|>", color=DARK, lw=2.5, mutation_scale=22))
    # panel 1: scattered dots coloured by family in clusters
    for k, c in enumerate(cols):
        cx, cy = 1.3 + 1.2 * k, 1.9 + (0.6 if k == 1 else -0.3)
        pts = rng.normal(0, 0.28, (9, 2))
        ax.scatter(cx + pts[:, 0], cy + pts[:, 1], s=60, c=c, edgecolors="white", zorder=4)
        ax.add_patch(Circle((cx, cy), 0.75, fc="none", ec=c, lw=1.5, ls="--"))
    # panel 2: matrix inside one family
    x0 = 5.0
    for i in range(5):
        for j in range(5):
            v = rng.uniform(0.2, 1.0) if i != j else 1
            ax.add_patch(Rectangle((x0 + 0.5 + 0.5 * j, 0.75 + 0.5 * i), 0.46, 0.46, fc=COMPLEX, alpha=v * 0.85, ec="none"))
    ax.text(x0 + 3.2, 1.9, "no all-\nvs-all\nsearch", fontsize=14, color=MUTED, va="center")
    # panel 3: small tree with merges
    x0 = 10.6; g = "#7E8B8A"
    elbow(ax, x0 + 0.4, 1.9, x0 + 2.6, 0.8, g, 3); elbow(ax, x0 + 0.4, 1.9, x0 + 1.5, 2.6, g, 3)
    elbow(ax, x0 + 1.5, 2.6, x0 + 2.6, 2.1, g, 3); elbow(ax, x0 + 1.5, 2.6, x0 + 2.6, 3.1, g, 3)
    for y in (0.8, 2.1, 3.1):
        gene(ax, x0 + 2.85, y, TEAL)
    for (x, y) in ((x0 + 1.5, 2.6), (x0 + 0.4, 1.9)):
        ax.add_patch(Circle((x, y), 0.2, fc=ACCENT, ec="white", lw=1.5, zorder=6))
    ax.text(x0 + 3.2, 1.9, "merge\n→ HOG\nat each\nnode", fontsize=13, color=MUTED, va="center")
    save(fig, "explain_ancestra_steps.png")


# E4: novelty ladder
def explain_ladder():
    fig, ax = plt.subplots(figsize=(11, 4.6)); clean(ax)
    ax.set_xlim(0, 16); ax.set_ylim(0, 6.4)
    tiers = [("All HOG gains on a branch", "raw Ancestra output; inflated 2.1–2.6× in simulations", 15.0, "#8C9796"),
             ("Sequence-novel", "no DIAMOND hit outside the clade (E ≤ 1e-5, ≥ 30% coverage)\n41–80% of gains in the panels", 12.5, TEAL),
             ("No structural match", "Foldseek search: about half of sequence-novel families\nhave a structural relative outside the clade", 10.0, COMPLEX)]
    for i, (t, sub, w, c) in enumerate(tiers):
        y = 4.6 - 2.0 * i
        x0 = 0.5
        ax.add_patch(FancyBboxPatch((x0, y), w, 1.5, boxstyle="round,pad=0.02,rounding_size=0.2", fc=c, ec="none"))
        ax.text(x0 + 0.3, y + 1.05, t, fontsize=18, weight="bold", color="white", va="center")
        ax.text(x0 + 0.3, y + 0.42, sub, fontsize=13, color="white", va="center")
    save(fig, "explain_ladder.png")


# E4b: why gains inflate (family split)
def explain_split():
    fig, ax = plt.subplots(figsize=(6, 4.2)); clean(ax)
    ax.set_xlim(0, 8); ax.set_ylim(0, 5.6)
    g = "#7E8B8A"
    for k, (x0, title) in enumerate(((0.3, "truth: one old family"), (4.3, "method splits it"))):
        ax.text(x0, 5.2, title, fontsize=15, weight="bold", color=DARK)
        elbow(ax, x0, 2.5, x0 + 2.6, 0.8, g, 3); elbow(ax, x0, 2.5, x0 + 1.3, 3.5, g, 3)
        elbow(ax, x0 + 1.3, 3.5, x0 + 2.6, 2.8, g, 3); elbow(ax, x0 + 1.3, 3.5, x0 + 2.6, 4.2, g, 3)
        if k == 0:
            for y in (0.8, 2.8, 4.2):
                gene(ax, x0 + 2.95, y, TEAL)
            gene(ax, x0 + 0.1, 2.85, TEAL)
        else:
            gene(ax, x0 + 2.95, 0.8, TEAL)
            for y in (2.8, 4.2):
                gene(ax, x0 + 2.95, y, ACCENT)
            star(ax, x0 + 0.65, 3.85, s=500, c=ACCENT)
            ax.text(x0 + 0.65, 4.3, "false gain", ha="center", fontsize=13, color=ACCENT)
    ax.text(4.0, 0.0, "old genes look new  →  raw gains inflated", ha="center", fontsize=14, color=MUTED)
    save(fig, "explain_split.png")


# E5: transition vs control, tip artefact (numbers: REPORT.md 6.3, Figure 3b)
def explain_controls():
    fig, ax = plt.subplots(figsize=(6.2, 4.6)); clean(ax)
    ax.set_xlim(-0.3, 7.6); ax.set_ylim(-0.9, 6.4)
    g = "#7E8B8A"
    xt = 6.0
    elbow(ax, 0, 2.6, 1.5, 1.0, g, 3); elbow(ax, 0, 2.6, 1.5, 4.2, g, 3)
    elbow(ax, 1.5, 1.0, xt, 0.3, g, 3); elbow(ax, 1.5, 1.0, 3.2, 1.7, "#B9C3C2", 7)
    elbow(ax, 3.2, 1.7, xt, 1.3, g, 3); elbow(ax, 3.2, 1.7, xt, 2.1, g, 3)
    ax.plot([1.5, 1.5], [4.2, 3.4], c=g, lw=3); ax.plot([1.5, xt], [3.4, 3.4], c=g, lw=3)
    ax.plot([1.5, 1.5], [4.2, 5.0], c=g, lw=3); ax.plot([1.5, 3.4], [5.0, 5.0], c=COMPLEX, lw=8)
    elbow(ax, 3.4, 5.0, xt, 4.4, COMPLEX, 3); elbow(ax, 3.4, 5.0, xt, 5.8, COMPLEX, 3)
    for y, c in ((0.3, UNI), (1.3, UNI), (2.1, UNI), (3.4, UNI), (4.4, COMPLEX), (5.8, COMPLEX)):
        ax.add_patch(Circle((xt + 0.22, y), 0.17, fc=c, ec="white", zorder=5))
    ax.text(2.0, 5.5, "event branch", fontsize=13, color=COMPLEX, ha="center", weight="bold")
    ax.text(3.9, 0.62, "control branch", fontsize=13, color=MUTED, ha="center")
    ax.text(3.8, -0.75, "observed ÷ expected from controls = rate ratio", fontsize=13, color=TEXT, ha="center")
    save(fig, "explain_controls.png")


def explain_tips():
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    labels = ["internal control\n(reference)", "unicellular tips\n(no event)", "origin tips", "reversal tips",
              "simple→complex\ntips"]
    vals = [1.0, 8.1, 19.0, 20.7, 40.8]
    cols = ["#B9C3C2", UNI, SIMPLE, REV, COMPLEX]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=cols, height=0.62)
    ax.set_xscale("log"); ax.set_xlim(0.8, 80)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=14)
    for yy, v in zip(y, vals):
        ax.text(v * 1.12, yy, f"{v:g}×", va="center", fontsize=14, color=TEXT)
    ax.set_xlabel("rate ratio of 'new' families vs internal controls (log)", fontsize=13)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=12)
    ax.set_xticks([1, 3, 10, 30]); ax.set_xticklabels(["1", "3", "10", "30"])
    save(fig, "explain_tips.png")


# E6: PGLS schematic
def explain_pgls():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw=dict(width_ratios=[0.9, 1.2]))
    ax = axs[0]; clean(ax); ax.set_xlim(-0.2, 6); ax.set_ylim(-0.8, 7.6)
    g = "#7E8B8A"
    ys = np.arange(8)
    states = [UNI, UNI, UNI, COMPLEX, COMPLEX, UNI, SIMPLE, SIMPLE]
    # simple balanced tree
    def node(y0, y1, x): return (y0 + y1) / 2
    xt = 4.5
    pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
    mids = []
    for a, b in pairs:
        m = (a + b) / 2; mids.append(m)
        elbow(ax, 3.2, m, xt, a, g, 2.5); elbow(ax, 3.2, m, xt, b, g, 2.5)
    m1 = (mids[0] + mids[1]) / 2; m2 = (mids[2] + mids[3]) / 2
    for m in mids[:2]: elbow(ax, 1.8, m1, 3.2, m, g, 2.5)
    for m in mids[2:]: elbow(ax, 1.8, m2, 3.2, m, g, 2.5)
    r = (m1 + m2) / 2
    elbow(ax, 0.3, r, 1.8, m1, g, 2.5); elbow(ax, 0.3, r, 1.8, m2, g, 2.5)
    for y, c in zip(ys, states):
        ax.add_patch(Circle((xt + 0.3, y), 0.25, fc=c, ec="white", zorder=5))
    ax.text(2.6, -0.75, "relatives are not independent", ha="center", fontsize=13, color=MUTED)
    ax.set_title("1  shared ancestry", fontsize=16, color=DARK, loc="left", weight="bold")
    ax = axs[1]
    rng = np.random.default_rng(4)
    x = rng.uniform(8.5, 11, 30)
    yu = 0.9 * x + rng.normal(0, 0.25, 30)
    x2 = rng.uniform(9, 11.3, 20)
    ym = 0.9 * x2 + 0.6 + rng.normal(0, 0.25, 20)
    ax.scatter(x, yu, c=UNI, s=45, label="unicellular", edgecolors="white")
    ax.scatter(x2, ym, c=COMPLEX, s=45, label="multicellular", edgecolors="white")
    xx = np.linspace(8.4, 11.4, 10)
    ax.plot(xx, 0.9 * xx, c=MUTED, lw=2); ax.plot(xx, 0.9 * xx + 0.6, c=COMPLEX, lw=2)
    ax.annotate("", (11.15, 0.9 * 11.15 + 0.6), (11.15, 0.9 * 11.15),
                arrowprops=dict(arrowstyle="<->", color=ACCENT, lw=2.5))
    ax.text(11.25, 0.9 * 11.15 + 0.3, "effect", color=ACCENT, fontsize=14, va="center")
    ax.set_xlabel("proteome size (log)", fontsize=14); ax.set_ylabel("family copy number (log)", fontsize=14)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=13, loc="upper left")
    ax.set_xlim(8.3, 12.2)
    ax.set_title("2  compare at equal proteome size", fontsize=16, color=DARK, loc="left", weight="bold")
    fig.tight_layout()
    save(fig, "explain_pgls.png")


# E7: pre-registration + Holm (p-values: REPORT.md 6.2, 7.2, 7.4)
def explain_holm():
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    tests = ["cell death (P3)", "simple→complex\nnovelty (P1)", "motif-binding\ndomains (P2)"]
    p = [0.0047, 0.0504, 0.72]
    thr = [0.05 / 3, 0.05 / 2, 0.05]
    y = np.arange(3)[::-1]
    cols = [SIMPLE, "#8C9796", "#8C9796"]
    ax.scatter(p, y, s=220, c=cols, zorder=5, edgecolors="white")
    for yy, t, pp in zip(y, thr, p):
        ax.plot([t, t], [yy - 0.3, yy + 0.3], c=ACCENT, lw=3)
        ax.text(pp, yy + 0.2, f"p = {pp:g}", fontsize=14, color=TEXT, va="bottom", ha="center")
        ax.text(t, yy - 0.36, f"must beat {t:.3f}", fontsize=12, color=ACCENT, va="top", ha="center")
    ax.set_xscale("log"); ax.set_xlim(0.001, 2)
    ax.set_yticks(y); ax.set_yticklabels(tests, fontsize=14)
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlabel("one-sided p-value (log); orange bar = Holm threshold", fontsize=13)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    save(fig, "explain_holm.png")


def explain_prereg():
    fig, ax = plt.subplots(figsize=(7.2, 2.4)); clean(ax)
    ax.set_xlim(0, 12); ax.set_ylim(0, 3.2)
    steps = [("write\ndesign", TEAL), ("hash +\ngit tag", DARK), ("run\nanalyses", COMPLEX), ("log every\nchange", ACCENT)]
    for i, (t, c) in enumerate(steps):
        x = 0.2 + i * 3.0
        ax.add_patch(FancyBboxPatch((x, 0.8), 2.3, 1.6, boxstyle="round,pad=0.03,rounding_size=0.3", fc=c, ec="none"))
        ax.text(x + 1.15, 1.6, t, ha="center", va="center", fontsize=16, color="white", weight="bold")
        if i < 3:
            ax.annotate("", (x + 2.95, 1.6), (x + 2.4, 1.6), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=2))
    ax.text(3.45, 0.3, "frozen before any output exists", ha="center", fontsize=12, color=MUTED)
    ax.text(9.45, 0.3, "dated, with 'Seen:' line", ha="center", fontsize=12, color=MUTED)
    save(fig, "explain_prereg.png")


# E8: ProstT5 -> 3Di -> Foldseek (numbers: REPORT.md 8.1; RESEARCH_TREE R2-3)
def explain_prostt5():
    fig, ax = plt.subplots(figsize=(12, 3.0)); clean(ax)
    ax.set_xlim(0, 16); ax.set_ylim(0, 3.6)
    boxes = [("protein sequence", "M K T L L V A …", "#8C9796"), ("ProstT5 on GPU", "language model", TEAL),
             ("3Di structure letters", "d v l p c q v …", COMPLEX), ("Foldseek search", "find relatives by shape", DARK)]
    for i, (t, s, c) in enumerate(boxes):
        x = 0.2 + i * 4.0
        ax.add_patch(FancyBboxPatch((x, 0.7), 3.3, 2.2, boxstyle="round,pad=0.03,rounding_size=0.3", fc=c, ec="none"))
        ax.text(x + 1.65, 2.25, t, ha="center", va="center", fontsize=15, color="white", weight="bold")
        ax.text(x + 1.65, 1.35, s, ha="center", va="center", fontsize=15, color="white",
                family="monospace" if i in (0, 2) else FONT)
        if i < 3:
            ax.annotate("", (x + 3.95, 1.8), (x + 3.4, 1.8), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=2.5))
    save(fig, "explain_prostt5.png")


def explain_gpu_bars():
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.0))
    ax = axs[0]
    ax.barh([1, 0], [15, 48], color=["#B9C3C2", TEAL], height=0.6)
    ax.set_yticks([1, 0]); ax.set_yticklabels(["Foldseek GPU build\n(silently on CPU)", "own implementation\n(A100)"], fontsize=13)
    for y, v in ((1, 15), (0, 48)): ax.text(v + 1, y, f"{v}", va="center", fontsize=14)
    ax.set_xlim(0, 60); ax.set_xlabel("proteins per second", fontsize=13)
    ax = axs[1]
    vals = [16.7, 31.4, 65.5]
    ax.barh([2, 1, 0], vals, color=["#D8E0DF", "#B9C3C2", COMPLEX], height=0.6)
    ax.set_yticks([2, 1, 0]); ax.set_yticklabels(["chance", "commonest letter", "our prediction"], fontsize=13)
    for y, v in zip([2, 1, 0], vals): ax.text(v + 1, y, f"{v}%", va="center", fontsize=14)
    ax.set_xlim(0, 85); ax.set_xlabel("residues matching AlphaFold letters", fontsize=13)
    for a in axs:
        for s in ("top", "right"): a.spines[s].set_visible(False)
        a.tick_params(axis="x", labelsize=11)
    fig.tight_layout()
    save(fig, "explain_gpu_bars.png")


# E9: EukaryOMA cross-check (numbers: REPORT.md section 5)
def explain_xcheck():
    fig, ax = plt.subplots(figsize=(11, 4.4)); clean(ax)
    ax.set_xlim(0, 16); ax.set_ylim(-0.4, 6.0)
    ax.add_patch(FancyBboxPatch((0.0, 2.0), 3.2, 1.8, boxstyle="round,pad=0.03,rounding_size=0.3", fc="#8C9796", ec="none"))
    ax.text(1.6, 2.9, "same proteins\n125 shared\nspecies", ha="center", va="center", fontsize=14, color="white", weight="bold")
    for k, (name, c, y) in enumerate((("Ancestra", TEAL, 4.0), ("EukaryOMA", COMPLEX, 0.6))):
        ax.add_patch(FancyBboxPatch((4.4, y), 2.8, 1.3, boxstyle="round,pad=0.03,rounding_size=0.3", fc=c, ec="none"))
        ax.text(5.8, y + 0.65, name, ha="center", va="center", fontsize=17, color="white", weight="bold")
        ax.annotate("", (4.35, y + 0.65), (3.25, 2.9), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=2))
        # small tree
        x0, g = 8.3, "#7E8B8A"
        yy = y + 0.65
        elbow(ax, x0, yy, x0 + 2.4, yy - 0.6, g, 2.5); elbow(ax, x0, yy, x0 + 1.2, yy + 0.4, g, 2.5)
        elbow(ax, x0 + 1.2, yy + 0.4, x0 + 2.4, yy + 0.1, g, 2.5); elbow(ax, x0 + 1.2, yy + 0.4, x0 + 2.4, yy + 0.8, g, 2.5)
        if k == 0:
            star(ax, x0 + 1.8, yy + 0.4 + 0.42, s=350, c=ACCENT)
        else:
            star(ax, x0 + 0.6, yy + 0.4 + 0.07, s=350, c=ACCENT)
        ax.annotate("", (x0 - 0.1, yy), (7.25, yy), arrowprops=dict(arrowstyle="-|>", color=DARK, lw=2))
    ax.text(9.5, 2.85, "where did the family\narise?", ha="center", va="center", fontsize=14, color=MUTED)
    ax.add_patch(FancyBboxPatch((11.6, 1.6), 4.2, 2.6, boxstyle="round,pad=0.03,rounding_size=0.3", fc=LIGHT, ec="none"))
    ax.text(13.7, 3.45, "29%", ha="center", va="center", fontsize=40, color=ACCENT, weight="bold")
    ax.text(13.7, 2.3, "same branch\n(540 of 1,859 stem gains)", ha="center", va="center", fontsize=13, color=TEXT)
    save(fig, "explain_xcheck.png")


# ----------------------------------------------------------------------------------------------
# Gantt (times: RESEARCH_TREE.md milestones; START.txt, START_ROUND2.txt)
def gantt():
    import datetime as dt
    t0 = dt.datetime(2026, 10, 28, 12, 0)
    def h(d, H, M): return (dt.datetime(2026, 10, d, H, M) - t0).total_seconds() / 3600
    fig, ax = plt.subplots(figsize=(13, 4.4))
    r1 = (h(28, 12, 21), h(29, 1, 55)); r2 = (h(29, 15, 7), h(30, 8, 31))
    ax.barh(1, r1[1] - r1[0], left=r1[0], color=TEAL, height=0.5)
    ax.barh(0, r2[1] - r2[0], left=r2[0], color=COMPLEX, height=0.5)
    ax.barh(0.5, r2[0] - r1[1], left=r1[1], color="#E3E8E7", height=0.25)
    ax.text((r1[1] + r2[0]) / 2, 0.5, "pause · human reads report · round-2 prompt", ha="center", va="center", fontsize=12, color=MUTED)
    ax.text(r1[0] + 0.2, 1.0, "Round 1 · 13.6 h · 24 CPU, 150 GB", va="center", fontsize=14, color="white", weight="bold")
    ax.text(r2[0] + 0.2, 0.0, "Round 2 · 17.4 h · 48 CPU, 500 GB, A100", va="center", fontsize=14, color="white", weight="bold")
    m1 = [(28, 12, 45, "literature"), (28, 14, 25, "pre-registration"), (28, 22, 36, "Ancestra run done"),
          (29, 1, 30, "tests + report")]
    m2 = [(29, 15, 41, "review + plan"), (29, 17, 33, "8 panels built"), (29, 22, 3, "confirmatory tests"),
          (30, 5, 20, "structure analyses"), (30, 8, 15, "verified + report")]
    for i, (d, H, M, lab) in enumerate(m1):
        x = h(d, H, M)
        ax.scatter([x], [1.33], marker="v", s=90, c=ACCENT, zorder=5)
        ax.text(x, 1.45 + (0.2 if i % 2 else 0), lab, ha="center", fontsize=13, color=TEXT)
    for i, (d, H, M, lab) in enumerate(m2):
        x = h(d, H, M)
        ax.scatter([x], [-0.33], marker="^", s=90, c=ACCENT, zorder=5)
        ax.text(x, -0.48 - (0.2 if i % 2 else 0), lab, ha="center", va="top", fontsize=13, color=TEXT)
    ticks = [0, 12, 24, 36, 48]
    ax.set_xticks(ticks); ax.set_xticklabels(["28 Oct 12:00", "29 Oct 00:00", "29 Oct 12:00", "30 Oct 00:00", "30 Oct 12:00"], fontsize=12)
    ax.set_yticks([]); ax.set_ylim(-0.95, 1.9); ax.set_xlim(-0.5, 46)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    save(fig, "gantt.png")


# ----------------------------------------------------------------------------------------------
# crops of report figures (report/figures_v2, figures/) at full resolution
def crops():
    F2 = ROOT / "report" / "figures_v2"
    jobs = {
        "fig1_ab.png": (F2 / "fig1_overview.png", (0, 0, 2125, 600)),
        "fig1_volvox.png": (F2 / "fig1_overview.png", (1630, 600, 2125, 1490)),
        "fig2_ab.png": (F2 / "fig2_orthology_agreement.png", (0, 0, 2125, 700)),
        "fig2_c.png": (F2 / "fig2_orthology_agreement.png", (0, 760, 2125, 1535)),
        "fig3_a.png": (F2 / "fig3_novelty_events.png", (0, 0, 2125, 525)),
        "fig3_b.png": (F2 / "fig3_novelty_events.png", (0, 550, 2125, 1035)),
        "fig3_c.png": (F2 / "fig3_novelty_events.png", (0, 1065, 2125, 1580)),
        "fig3_d.png": (F2 / "fig3_novelty_events.png", (0, 1600, 2125, 2539)),
        "fig4_a.png": (F2 / "fig4_family_sets.png", (0, 0, 2125, 850)),
        "fig4_bc.png": (F2 / "fig4_family_sets.png", (0, 880, 2125, 1771)),
        "fig5_a.png": (F2 / "fig5_celldeath_architectures.png", (0, 0, 2125, 1140)),
        "fig5_b.png": (F2 / "fig5_celldeath_architectures.png", (0, 1180, 2125, 1790)),
        "fig7_a.png": (F2 / "fig7_structure.png", (0, 0, 1060, 840)),
        "fig7_b.png": (F2 / "fig7_structure.png", (1090, 0, 2125, 840)),
        "fig7_c.png": (F2 / "fig7_structure.png", (0, 880, 1060, 1771)),
        "fig7_d.png": (F2 / "fig7_structure.png", (1090, 880, 2125, 1771)),
    }
    for name, (src, box) in jobs.items():
        im = Image.open(src).convert("RGB").crop(box)
        im.save(OUT / name)
        print("wrote", name, im.size)


if __name__ == "__main__":
    motif()
    hypotheses()
    explain_ancestral(); explain_ancestra(); explain_ladder(); explain_split()
    explain_controls(); explain_tips(); explain_pgls(); explain_holm(); explain_prereg()
    explain_prostt5(); explain_gpu_bars(); explain_xcheck()
    gantt()
    crops()
