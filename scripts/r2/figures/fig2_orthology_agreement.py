"""Figure 2. Ancestra vs EukaryOMA: per-stem counts of new families and duplications, and branch-placement agreement.
Sources: results/eukaryoma_xval/{counts_per_branch.tsv, agreement_per_branch.tsv, agreement_stems_depth_adjusted.tsv}.
The six exact internal stems are those used for the pooled gain agreement (agreement_pooled.tsv, stems_exact, gain);
duplication agreement is pooled over the same six stems (as in summary.txt), i.e. tip stems are excluded.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from common import P, MM, ANC, EOMA, INK, INK2, GRID, save, fig_label

D = f"{P}/results/eukaryoma_xval"
STEMS = ["T_metazoa", "T_dictyostelia", "T_streptophyta", "T_rhodophyta", "T_florideophyceae", "T_px_filaments"]
LAB = {"T_metazoa": "animals", "T_dictyostelia": "dictyostelid slime moulds", "T_streptophyta": "streptophytes",
       "T_rhodophyta": "red algae (Eurhodophytina)", "T_florideophyceae": "florideophyte red algae",
       "T_px_filaments": "brown-algal lineage (PX clade)"}


def pooled(per, branches, direction, event):
    s = per[(per.branch_id.isin(branches)) & (per.direction == direction) & (per.event == event)]
    return s.same.sum() / s.n_mapped.sum(), s


def main():
    cnt = pd.read_csv(f"{D}/counts_per_branch.tsv", sep="\t")
    per = pd.read_csv(f"{D}/agreement_per_branch.tsv", sep="\t")
    adj = pd.read_csv(f"{D}/agreement_stems_depth_adjusted.tsv", sep="\t")
    ctrls = sorted(per[(per.role == "control") & (per.status == "exact")].branch_id.unique())
    assert len(ctrls) == 26

    fig = plt.figure(figsize=(180 * MM, 130 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], wspace=0.08, hspace=0.62, left=0.225, right=0.985, top=0.94, bottom=0.11)
    y = np.arange(len(STEMS))[::-1]
    h = 0.36
    out = {}
    for j, (col, title, xl) in enumerate([("gains_all", "New gene families", "new families on the stem (log scale)"),
                                           ("dups", "Duplications", "duplications on the stem (log scale)")]):
        ax = fig.add_subplot(gs[0, j])
        for yi, b in zip(y, STEMS):
            a = cnt[(cnt.branch_id == b) & (cnt.method == "ancestra")][col].iloc[0]
            e = cnt[(cnt.branch_id == b) & (cnt.method == "eukaryoma")][col].iloc[0]
            out[(b, col)] = (a, e)
            ax.barh(yi + h / 2, a, height=h, color=ANC, edgecolor="white", linewidth=0.6)
            ax.barh(yi - h / 2, e, height=h, color=EOMA, edgecolor="white", linewidth=0.6)
            ax.text(a * 1.12, yi + h / 2, f"{a:,}", va="center", fontsize=5.6, color=INK2)
            ax.text(e * 1.12, yi - h / 2, f"{e:,}", va="center", fontsize=5.6, color=INK2)
        ax.set_xscale("log")
        ax.set_xlim(30, 3e4)
        ax.set_yticks(y)
        ax.set_yticklabels([LAB[b] for b in STEMS] if j == 0 else [])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_xlabel(xl)
        ax.set_title(title, loc="left", color=INK)
        ax.grid(axis="x", color=GRID, lw=0.5, which="major")
        ax.set_axisbelow(True)
    a, e = out[("T_metazoa", "gains_all")]
    ad, ed = out[("T_metazoa", "dups")]
    fig.axes[0].legend(handles=[Patch(color=ANC, label="Ancestra"), Patch(color=EOMA, label="EukaryOMA")],
                       loc="lower right", bbox_to_anchor=(1.0, 0.0), handlelength=1.0)

    # ---- (c) agreement ----
    ax = fig.add_subplot(gs[1, :])
    pos = ax.get_position(); ax.set_position([0.09, pos.y0, 0.6, pos.height])
    groups = [("ancestra->eukaryoma", "gain", "new families,\nAncestra →\nEukaryOMA"),
              ("eukaryoma->ancestra", "gain", "new families,\nEukaryOMA →\nAncestra"),
              ("ancestra->eukaryoma", "dup", "duplications,\nAncestra →\nEukaryOMA"),
              ("eukaryoma->ancestra", "dup", "duplications,\nEukaryOMA →\nAncestra")]
    w = 0.26
    rows = []
    for gi, (dr, ev, lab) in enumerate(groups):
        fs, s_st = pooled(per, STEMS, dr, ev)
        fc, s_ct = pooled(per, ctrls, dr, ev)
        aa = adj[(adj.direction == dr) & (adj.event == ev) & (adj.branch_id.isin(STEMS))]
        aa = aa.merge(s_st[["branch_id", "n_mapped"]], on="branch_id", suffixes=("_adj", ""))
        fe = (aa.expected_frac_same * aa.n_mapped).sum() / aa.n_mapped.sum()
        rows.append((lab.replace("\n", " "), fs, fc, fe, int(s_st.same.sum()), int(s_st.n_mapped.sum())))
        x0 = gi
        ax.bar(x0 - w, 100 * fs, width=w, color=INK, edgecolor="white", linewidth=0.6)
        ax.bar(x0, 100 * fc, width=w, color="#b0b0b0", edgecolor="white", linewidth=0.6)
        ax.bar(x0 + w, 100 * fe, width=w, facecolor="white", edgecolor=INK, linewidth=0.8, hatch="//////")
        rng = np.random.default_rng(11)
        ax.scatter(x0 - w + rng.uniform(-0.06, 0.06, len(s_st)), 100 * s_st.frac_same, s=6, color="white",
                   edgecolor=INK, linewidth=0.5, zorder=3)
        ax.scatter(x0 + rng.uniform(-0.06, 0.06, len(s_ct)), 100 * s_ct.frac_same, s=4, color="#6f6f6f",
                   edgecolor="none", zorder=3, alpha=0.8)
        tops = [max(fs, s_st.frac_same.max()), max(fc, s_ct.frac_same.max()), fe]
        for xx, v, tp in [(x0 - w, fs, tops[0]), (x0, fc, tops[1]), (x0 + w, fe, tops[2])]:
            ax.text(xx, 100 * v + 0.8, f"{100 * v:.0f}", ha="center", va="bottom", fontsize=6, color=INK, zorder=5,
                    bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85))
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[2] for g in groups], fontsize=6.3)
    ax.set_ylabel("events placed on the same\nbranch by both methods (%)")
    ax.set_ylim(0, 60)
    ax.set_title("Agreement on where events happened", loc="left", color=INK)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color=INK, label="6 stems (pooled)"), Patch(color="#b0b0b0", label="26 control branches"),
                       Patch(facecolor="white", edgecolor=INK, hatch="//////", label="expected for stem depth and size"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=3, label="single stem"),
                       Line2D([], [], marker="o", ls="", color="#6f6f6f", ms=2.5, label="single control")],
              loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=6.3, handlelength=1.2)
    fig_label(fig, "a", 0.005, 0.985)
    fig_label(fig, "b", 0.58, 0.985)
    fig_label(fig, "c", 0.005, 0.47)
    save(fig, "fig2_orthology_agreement")
    print("animal stem: gains A/E", a, e, f"ratio {e / a:.2f}", "dups A/E", ad, ed, f"ratio {ad / ed:.1f}")
    for r in rows:
        print(f"{r[0]:45s} stems {100 * r[1]:.1f}% ({r[4]}/{r[5]})  controls {100 * r[2]:.1f}%  expected {100 * r[3]:.1f}%")


if __name__ == "__main__":
    main()
