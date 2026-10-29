"""Task D figures (PNG + PDF, 300 dpi) in results/r2/w4/figures/."""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from common import P, OUT, out  # noqa: E402

BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
STATE_COL = {"unicellular": "#d9d8d4", "grey": GREY, "simple": AQUA, "complex": BLUE}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(out(f"figures/{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def forest():
    r = pd.read_csv(f"{OUT}/pgls/pgls_results.tsv", sep="\t")
    r = r[r.famset == "primary"]
    titles = {"P2_elm": "R2-P2 motif-binding (ELM, 115 families)\nconfirmatory",
              "P3_celldeath": "R2-P3 cell death (12 widespread families)\nconfirmatory",
              "P4_tf": "P4 cell-type TFs (11 widespread families)\nexploratory"}
    keep = ["main", "BM (lambda=1)", "no disorder covariate", "grey-zone excluded", "all 187 proteomes",
            "all 187 + transcriptome & BUSCO covariates"]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.6), sharey=True)
    for ax, (s, t) in zip(axes, titles.items()):
        g = r[r.set == s]
        rows = [g[g.variant == "main"]] + [g[g.variant == v] for v in keep[1:]] + [g[g.variant.str.startswith("LOO")]]
        g2 = pd.concat(rows)
        labels = [v.replace(" (", "\n(") if not v.startswith("LOO") else v.split(" (")[0].replace("LOO -T_", "drop ") for v in g2.variant]
        y = np.arange(len(g2))[::-1]
        for yi, (_, row), lab in zip(y, g2.iterrows(), labels):
            main = row.variant == "main"
            loo = row.variant.startswith("LOO")
            col = BLUE if main else (ORANGE if loo else INK2)
            ax.plot([row.ci_lo, row.ci_hi], [yi, yi], color=col, lw=2 if main else 1.2, solid_capstyle="round")
            ax.plot(row.estimate, yi, "o", color=col, ms=6 if main else 4.5, mec="white", mew=0.8, zorder=3)
            ax.text(1.02, yi, f"p={row.p_one:.3f}" if row.p_one >= 0.001 else f"p={row.p_one:.1e}", transform=ax.get_yaxis_transform(),
                    va="center", fontsize=6.5, color=INK2)
        ax.axvline(0, color=INK2, lw=0.8, ls=(0, (2, 2)))
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.5)
        ax.set_title(t, fontsize=8, loc="left", color=INK)
        ax.grid(axis="x", color="#e6e5e0", lw=0.5)
    axes[0].plot([], [], "o-", color=BLUE, label="primary PGLS (λ), 148 genomes")
    axes[0].plot([], [], "o-", color=INK2, label="sensitivity")
    axes[0].plot([], [], "o-", color=ORANGE, label="leave one primary origin out")
    fig.supxlabel("multicellular − unicellular coefficient (SD units of the equal-weight mean of standardised log copy number)", fontsize=8, y=0.06)
    fig.legend(loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("State effect on per-set copy number (95% CI; one-sided p for multicellular > unicellular)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    save(fig, "forest_sets_state_effect")


def tree_order():
    nwk = open(f"{P}/data/species_tree_bl_euk.nwk").read()
    return re.findall(r"[(,]([A-Za-z0-9_]+):", nwk)


def heatmap(fams, title, name, counts, sp, labels):
    order = [c for c in tree_order() if c in counts.index]
    M = np.log2(counts.loc[order, fams].to_numpy(float) + 1)
    fig = plt.figure(figsize=(max(4.5, 0.17 * len(fams) + 2.4), 13))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.25, 0.25, 0.17 * len(fams)], wspace=0.02)
    a0, a1, a2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
    st = sp.loc[order, "mc_state"].map(list(STATE_COL).index).to_numpy()
    a0.imshow(st[:, None], aspect="auto", cmap=ListedColormap(list(STATE_COL.values())), vmin=0, vmax=3, interpolation="nearest")
    dt = (sp.loc[order, "data_type"] == "transcriptome").to_numpy().astype(int)
    a1.imshow(dt[:, None], aspect="auto", cmap=ListedColormap(["#ffffff", INK2]), vmin=0, vmax=1, interpolation="nearest")
    im = a2.imshow(M, aspect="auto", cmap="Blues", interpolation="nearest", vmin=0)
    a0.set_yticks(range(len(order)))
    a0.set_yticklabels(order, fontsize=3.6)
    for a in (a0, a1):
        a.set_xticks([0]); a.tick_params(axis="x", labelsize=6)
    a0.set_xticklabels(["state"], rotation=90); a1.set_xticklabels(["transcr."], rotation=90)
    a1.set_yticks([]); a2.set_yticks([])
    a2.set_xticks(range(len(fams)))
    a2.set_xticklabels(labels, rotation=90, fontsize=5.5)
    a2.xaxis.tick_top()
    for a in (a0, a1, a2):
        for s in a.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, ax=a2, fraction=0.03, pad=0.01)
    cb.set_label("log2(proteins with domain + 1)")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c, label=k) for k, c in STATE_COL.items()], loc="lower left", ncol=4, frameon=False,
               bbox_to_anchor=(0.02, 0.0), fontsize=7, title="multicellularity state (round-1 coding)")
    fig.suptitle(title + "\nspecies in round-1 tree order", x=0.02, ha="left", fontsize=9)
    save(fig, name)


def heatmaps():
    counts = pd.read_csv(f"{OUT}/matrix/counts.tsv.gz", sep="\t", index_col=0)
    sp = pd.read_csv(f"{OUT}/matrix/species.tsv", sep="\t").set_index("code")
    clans = pd.read_csv(f"{P}/data/raw/pfam/Pfam-A.clans.tsv.gz", sep="\t", header=None, names=["pfam", "clan", "cn", "name", "d"]).set_index("pfam")
    rec = pd.read_csv(f"{OUT}/scan/recurrent_families.tsv", sep="\t", index_col=0)
    up = rec[rec.rec_primary_expanded >= 2].sort_values(["rec_primary_expanded", "perm_p_primary_expanded"], ascending=[False, True]).head(30)
    dn = rec[rec.rec_primary_contracted >= 2].sort_values(["rec_primary_contracted", "perm_p_primary_contracted"], ascending=[False, True]).head(10)
    fams = list(up.index) + [f for f in dn.index if f not in up.index]
    labels = [f"{clans.name.get(f, f)} (+{int(rec.loc[f, 'rec_primary_expanded'])}/−{int(rec.loc[f, 'rec_primary_contracted'])})" for f in fams]
    heatmap(fams, "Most recurrent families of the scan (top 30 expanded, top 10 contracted; +k/−k = primary origins with q<0.1)",
            "heatmap_recurrent_families", counts, sp, labels)
    cd = pd.read_csv(f"{P}/data/r2_sets/cell_death.tsv", sep="\t")
    fs = pd.read_csv(f"{OUT}/pgls_input/set_families.tsv", sep="\t").set_index("pfam")
    labs = [f"{n.split(' (')[0]}{'' if fs.loc[f, 'widespread'].max() == 1 else ' †'}" for f, n in zip(cd.pfam, cd["name"])]
    heatmap(list(cd.pfam), "R2-P3 cell-death families († lineage-restricted, not in pooled response)", "heatmap_celldeath", counts, sp, labs)


def recurrence():
    n = pd.read_csv(f"{OUT}/scan/recurrence_null_summary.tsv", sep="\t")
    c = pd.read_csv(f"{OUT}/scan/recurrence_conditional_null_summary.tsv", sep="\t")
    obs = pd.read_csv(f"{OUT}/scan/recurrence_observed.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
    for ax, d in zip(axes, ("expanded", "contracted")):
        o = obs[f"rec_primary_{d}"].to_numpy()
        ks = np.arange(1, 6)
        ob = np.array([(o >= k).sum() for k in ks])
        nn = n[(n.tier == "primary") & (n.direction == d)].set_index("k").loc[ks]
        cc = c[c.direction == d].set_index("k").loc[ks]
        w = 0.26
        ax.bar(ks - w - 0.02, ob, w, color=BLUE, label="observed")
        for off, t, col, lab in ((0, nn, "#c3c2b7", "null 1: label permutation (pre-specified)"),
                                 (w + 0.02, cc, ORANGE, "null 2: family shuffle, per-origin counts kept (added)")):
            ax.bar(ks + off, t.null_mean, w, color=col, label=lab)
            ax.errorbar(ks + off, t.null_mean, yerr=[(t.null_mean - t.null_q025).clip(lower=0), (t.null_q975 - t.null_mean).clip(lower=0)], fmt="none", ecolor=INK2, lw=0.9, capsize=2)
        for k, b, pc, pn in zip(ks, ob, cc.p_perm, nn.p_perm):
            if b > 0:
                ax.text(k - w - 0.02, b * 1.25 + 0.3, f"{b}", ha="center", fontsize=6.5, color=INK)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(0, 6000)
        ax.set_xticks(ks); ax.set_xticklabels([f"≥{k}" for k in ks])
        ax.set_xlabel(f"primary origins (of 7) where the family is significantly {d} (BH q<0.1)")
        ax.set_title(d.capitalize() + " in the multicellular clade vs its sister", loc="left", fontsize=9)
        ax.grid(axis="y", color="#e6e5e0", lw=0.5)
        txt = "; ".join(f"≥{k}: p₂={p:.3f}" for k, p in zip(ks[1:4], cc.p_perm.to_numpy()[1:4]))
        ax.text(0.98, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color=INK2)
    axes[0].set_ylabel("number of Pfam families (symlog)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Recurrence of per-family copy-number changes across 7 primary origins vs two permutation nulls (B=1000; bars = mean, whiskers = 95%)",
                 x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, "recurrence_vs_null")


if __name__ == "__main__":
    which = sys.argv[1:] or ["forest", "heatmaps", "recurrence"]
    for w in which:
        globals()[w]()
