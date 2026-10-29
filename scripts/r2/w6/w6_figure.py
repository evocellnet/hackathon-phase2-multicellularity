#!/usr/bin/env python
"""W6 figure: (A) 3Di vs AA branch length (positive sanity); (B) per-event stem residual d with gene-bootstrap 95 % CI,
real vs null pipelines; (C) per-gene crown-minus-unicellular contrast by panel and pipeline."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/w6"
panels = sys.argv[1:] or ["R1", "R2"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
PIPES = ["real", "sim_constratio", "shuffled3di", "aa_splithalf"]
COL = {"real": "#2a78d6", "sim_constratio": "#1baf7a", "shuffled3di": "#eb6834", "aa_splithalf": "#4a3aa7"}
LAB = {"real": "3Di / AA (real)", "sim_constratio": "constant-ratio simulation (null)",
       "shuffled3di": "shuffled 3Di / AA (null)", "aa_splithalf": "AA half / AA half (null)"}
SHORT = {"real": "real", "sim_constratio": "sim", "shuffled3di": "shuf", "aa_splithalf": "split"}
CLS = {"uni": ("#2a78d6", "unicellular"), "crown": ("#eb6834", "inside multicellular crown"), "stem": ("#0b0b0b", "event stem")}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 3, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1, 1.5, 1]})

# A
a = ax[0]
for R in panels:
    df = pd.read_csv(os.path.join(OUT, R, "gene_branch_residuals.tsv.gz"), sep="\t")
    for c, (col, lab) in CLS.items():
        d = df[df["cls"] == c]
        a.scatter(np.log10(d["bl_aa"]), np.log10(np.maximum(d["bl_di"], 1e-6)), s=3, alpha=0.15, color=col, lw=0,
                  label=lab if R == panels[0] else None)
a.set_xlabel("AA branch length (log10 subst./site)")
a.set_ylabel("3Di branch length (log10 subst./site)")
a.set_title("A  3Di tracks AA length (3Di = 0 drawn at 1e-6)", loc="left", color=INK)
lg = a.legend(frameon=False, markerscale=4, fontsize=8)
for h in lg.legend_handles:
    h.set_alpha(1)
a.grid(color=GRID, lw=0.6)

# B
b = ax[1]
st = pd.read_csv(os.path.join(OUT, "stem_tests.tsv"), sep="\t")
st = st[st["analysis"] == "primary"]
ev = st[st["pipeline"] == "real"].sort_values(["panel", "d"])["event"].tolist()
y = {e: i for i, e in enumerate(ev)}
off = {"real": 0.27, "sim_constratio": 0.09, "shuffled3di": -0.09, "aa_splithalf": -0.27}
for pipe in PIPES:
    g = st[(st["pipeline"] == pipe) & st["event"].isin(y)]
    yy = np.array([y[e] for e in g["event"]]) + off[pipe]
    b.errorbar(g["d"], yy, xerr=[g["d"] - g["d_lo"], g["d_hi"] - g["d"]], fmt="o", ms=4, color=COL[pipe], lw=1.2,
               capsize=0, label=LAB[pipe])
b.axvline(0, color=MUTED, lw=0.8)
b.set_yticks(range(len(ev)))
b.set_yticklabels([e + (" (tip)" if bool(st[(st.event == e) & (st.pipeline == "real")]["tip"].iloc[0]) else "") for e in ev], fontsize=8)
b.set_xlabel("stem minus matched controls, log-ratio residual\n(>0: more 3Di change per AA change)")
b.set_xlim(-2.5, 1.6)
b.set_title("B  Event stems vs matched controls (95% gene-bootstrap CI)", loc="left", color=INK)
b.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
b.grid(axis="x", color=GRID, lw=0.6)

# C
c = ax[2]
ge = pd.read_csv(os.path.join(OUT, "gene_crown_effects.tsv.gz"), sep="\t")
pos, ticks, tl = 0, [], []
for R in panels:
    for pipe in PIPES:
        v = ge[(ge["panel"] == R) & (ge["pipeline"] == pipe)]["crown_minus_uni"].to_numpy()
        bp = c.boxplot(v, positions=[pos], widths=0.6, showfliers=False, patch_artist=True,
                       medianprops=dict(color=INK, lw=1.5), whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
        bp["boxes"][0].set(facecolor=COL[pipe], alpha=0.55, edgecolor=MUTED)
        ticks.append(pos)
        tl.append(f"{SHORT[pipe]}" + (f"\n{R}" if pipe == "sim_constratio" else ""))
        pos += 1
    pos += 0.6
c.axhline(0, color=MUTED, lw=0.8)
c.set_xticks(ticks)
c.set_xticklabels(tl, fontsize=8)
c.set_ylabel("per-gene crown minus unicellular residual")
c.set_title("C  Crown vs unicellular branches (genes)", loc="left", color=INK)
c.grid(axis="y", color=GRID, lw=0.6)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "w6_figure.png"), dpi=150)
fig.savefig(os.path.join(OUT, "w6_figure.pdf"))
