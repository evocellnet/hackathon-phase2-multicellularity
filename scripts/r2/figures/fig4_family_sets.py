"""Figure 4. Copy number of three pre-defined family sets in multicellular vs unicellular genomes (PGLS, 148 genomes).
Sources: results/r2/w4/pgls/pgls_results.tsv (primary, three-state, leave-one-origin-out);
results/r2/verification/w4/v2_pgls_output.txt (independent re-derivation; proteome-size sensitivities).
"""
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from common import P, MM, SIMPLE, COMPLEX, INK, INK2, GRID, save, fig_label

SETS = [("P3_celldeath", "cell-death families\n(12 families; confirmatory)"),
        ("P2_elm", "motif-binding domains\n(115 families; confirmatory)"),
        ("P4_tf", "cell-type transcription factors\n(11 families; exploratory)")]
SET_STYLE = {"P3_celldeath": dict(marker="o", mfc=INK, mec=INK), "P2_elm": dict(marker="s", mfc="#8a8a8a", mec="#8a8a8a"),
             "P4_tf": dict(marker="^", mfc="white", mec=INK)}
SET_SHORT = {"P3_celldeath": "cell death", "P2_elm": "motif-binding", "P4_tf": "cell-type transcription factors"}
ORIG = {"T_metazoa": "animals", "T_fungi_hyphae": "hyphal fungi", "T_dictyostelia": "dictyostelids",
        "T_streptophyta": "streptophytes", "T_rhodophyta": "red algae", "T_px_filaments": "brown-algal lineage",
        "T_volvocine_tgv": "volvocine algae"}


def ptxt(p):
    return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.4f}"


def style_x(ax, lo=-0.4, hi=0.45):
    ax.axvline(0, color=INK2, lw=0.7, ls=(0, (2, 2)), zorder=0)
    ax.set_xlim(lo, hi)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def main():
    r = pd.read_csv(f"{P}/results/r2/w4/pgls/pgls_results.tsv", sep="\t")
    r = r[r.famset == "primary"]
    v2 = open(f"{P}/results/r2/verification/w4/v2_pgls_output.txt").read()

    def v2get(label):
        m = re.search(re.escape(label) + r"\s+n=(\d+) b=([+-][\d.]+) SE=[\d.]+ 95%CI\(t\)=\[([-\d.]+),([-\d.]+)\] p1=([\d.]+)", v2)
        assert m, label
        return dict(n=int(m.group(1)), b=float(m.group(2)), lo=float(m.group(3)), hi=float(m.group(4)), p=float(m.group(5)))

    fig = plt.figure(figsize=(180 * MM, 150 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], width_ratios=[1, 1], hspace=0.55, wspace=1.05,
                          left=0.22, right=0.86, top=0.93, bottom=0.08)

    # ---------- (a) main effects + three-state ----------
    ax = fig.add_subplot(gs[0, :])
    pos = ax.get_position(); ax.set_position([0.25, pos.y0, 0.47, pos.height])
    y = 0; yt, yl = [], []
    for s, lab in SETS:
        g = r[r.set == s].set_index("variant")
        rows = [("multicellular vs unicellular", g.loc["main"], INK, 2.0, 5.5),
                ("simple vs unicellular", g.loc["3-state: simple vs uni"], SIMPLE, 1.2, 4.2),
                ("complex vs unicellular", g.loc["3-state: complex(D or T) vs uni"], COMPLEX, 1.2, 4.2)]
        for k, (rl, row, col, lw, ms) in enumerate(rows):
            yy = -(y + k * 0.85)
            ax.plot([row.ci_lo, row.ci_hi], [yy, yy], color=col, lw=lw, solid_capstyle="butt")
            ax.plot(row.estimate, yy, "o", color=col, ms=ms, mec="white", mew=0.6, zorder=3)
            ax.text(1.02, yy, f"{row.estimate:+.3f} [{row.ci_lo:+.3f}, {row.ci_hi:+.3f}]  {ptxt(row.p_one)}",
                    transform=ax.get_yaxis_transform(), fontsize=5.9, va="center",
                    color=INK if k == 0 else INK2, fontweight="bold" if k == 0 else "normal")
            yt.append(yy); yl.append(rl)
        yt.append(-(y - 1.3)); yl.append(lab)
        y += 4.6
    ax.set_yticks(yt); ax.set_yticklabels(yl)
    for t, l in zip(ax.get_yticklabels(), yl):
        if "families" in l:
            t.set_fontweight("bold"); t.set_color(INK); t.set_fontsize(6.6)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-(y - 1.2), 1.9)
    style_x(ax)
    ax.set_xlabel("difference in copy number (SD of standardised log copy number, PGLS, 148 genomes)")
    ax.text(1.02, 1.0, "estimate [95% CI], one-sided p", transform=ax.transAxes, fontsize=5.8, color=INK2, va="bottom")
    ax.set_title("More cell-death proteins in multicellular genomes, mostly in complex ones", loc="left", color=INK, pad=12)

    # ---------- (b) leave one origin out ----------
    ax = fig.add_subplot(gs[1, 0])
    origins = list(ORIG)
    yy = np.arange(len(origins) + 1)[::-1].astype(float)
    off = {"P3_celldeath": 0.22, "P2_elm": 0.0, "P4_tf": -0.22}
    for s, _ in SETS:
        g = r[r.set == s]
        main_ = g[g.variant == "main"].iloc[0]
        ests = [(main_.estimate, main_.ci_lo, main_.ci_hi)]
        for o in origins:
            row = g[g.variant.str.startswith(f"LOO -{o} ")].iloc[0]
            ests.append((row.estimate, row.ci_lo, row.ci_hi))
        for yi, (e, lo, hi) in zip(yy, ests):
            ax.plot([lo, hi], [yi + off[s]] * 2, color=SET_STYLE[s]["mec"] if s != "P2_elm" else "#8a8a8a", lw=0.8)
            ax.plot(e, yi + off[s], ls="", ms=3.6, mew=0.8, zorder=3, **SET_STYLE[s])
    ax.set_yticks(yy)
    ax.set_yticklabels(["all origins"] + [f"without {ORIG[o]}" for o in origins])
    ax.tick_params(axis="y", length=0)
    ax.axhline(yy[0] - 0.5, color=GRID, lw=0.8)
    style_x(ax)
    ax.set_xlabel("multicellular − unicellular (SD units)")
    ax.set_title("Leaving out one origin at a time", loc="left", color=INK, pad=14)
    ax.legend(handles=[Line2D([], [], ls="", ms=4, label=SET_SHORT[s], **SET_STYLE[s]) for s, _ in SETS],
              loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, fontsize=6, handletextpad=0.2, columnspacing=0.8)

    # ---------- (c) proteome-size sensitivity (cell death) ----------
    ax = fig.add_subplot(gs[1, 1])
    rows = [("primary (log proteome size\n+ disorder as covariates)", v2get("P3 primary (rule families)")),
            ("no covariates", v2get("P3 without covariates (y ~ multicellular)")),
            ("proteome size only", v2get("P3 log_size only")),
            ("10 largest proteomes removed", v2get("P3 excl. 10 largest, standardisation kept (full S)")),
            ("10 largest removed,\nre-standardised", v2get("P3 excl. 10 largest, re-standardised")),
            ("share of proteome\n(count / proteome size)", v2get("P3 using proportion log(count/size) instead (diagnostic)"))]
    yy = np.arange(len(rows))[::-1]
    for yi, (lab, d) in zip(yy, rows):
        col = INK if yi == yy[0] else INK2
        ax.plot([d["lo"], d["hi"]], [yi, yi], color=col, lw=1.6 if yi == yy[0] else 1.0, solid_capstyle="butt")
        ax.plot(d["b"], yi, "o", color=col, ms=4.5 if yi == yy[0] else 3.6, mec="white", mew=0.6, zorder=3)
        ax.text(1.03, yi, f"{d['b']:+.3f}\n{ptxt(d['p'])}, n = {d['n']}", transform=ax.get_yaxis_transform(),
                fontsize=5.6, va="center", color=INK)
    ax.set_yticks(yy); ax.set_yticklabels([x[0] for x in rows], fontsize=6)
    ax.tick_params(axis="y", length=0)
    style_x(ax, -0.15, 0.45)
    ax.set_xlabel("cell-death difference (SD units)")
    ax.set_title("Cell death: dependence on proteome size", loc="left", color=INK, pad=14)

    fig_label(fig, "a", 0.01, 0.985)
    fig_label(fig, "b", 0.01, 0.47)
    fig_label(fig, "c", 0.50, 0.47)
    save(fig, "fig4_family_sets")
    for lab, d in rows:
        print(lab.replace("\n", " "), d)
    print(r[r.variant.str.startswith("LOO")].groupby("set").apply(lambda g: (round(g.p_one.min(), 4), round(g.p_one.max(), 3), int((g.estimate > 0).sum()))))


if __name__ == "__main__":
    main()
