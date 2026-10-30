"""D17 figure: report/figures_v2/fig8_recurrence_ladder.{png,pdf}."""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
sys.path.insert(0, os.path.dirname(__file__))
from common import P, OUT, poibin_pmf  # noqa: E402

FIG = f"{P}/report/figures_v2"
for ext in ("png", "pdf"):
    if os.path.exists(f"{FIG}/fig8_recurrence_ladder.{ext}"):
        raise SystemExit("refusing to overwrite fig8")
INK, INK2, GRID = "#1a1a1a", "#555555", "#e5e5e5"
GROUPS = {"R1": ("#009E73", "Seven origins, round 1"), "NEST": ("#0072B2", "Simple-to-complex steps"),
          "PORIG": ("#009E73", "Origins, focused datasets"), "TIP": ("#E69F00", "Single-species branches"),
          "PROK": ("#666666", "Cyanobacteria")}
E = pd.read_csv(f"{OUT}/hypothesis_events.tsv", sep="\t")
# (ladder, hypothesis, criterion column suffix, label, group)
ROWS = [
    ("R1 origins", "H1H2", "0.05", "New gene families", "R1"),
    ("R1 origins", "H4", "0.05", "Transcription-factor duplications", "R1"),
    ("R1 origins", "H5", "0.05", "Adhesion and signalling gains", "R1"),
    ("R1 origins", "H12", "0.05", "Loss of flagellar genes", "R1"),
    ("R1 origins", "H11 new_arch_multi", "0.05", "New multi-domain architectures", "R1"),
    ("R1 origins", "H1 co-option pfam_expanded", "0.05", "Co-option of pre-existing families", "R1"),
    ("R1 origins", "H10 disorder_sc", "top", "More disordered new genes*", "R1"),
    ("R1 origins", "H9", "top", "Faster protein evolution*", "R1"),
    ("R1 origins", "H1H2 novelty (eukaryoma, 125 shared species)", "top", "New gene families, EukaryOMA*", "R1"),
    ("R1 origins", "W7 complex subunits among stem duplications", "0.05", "Duplicated complex subunits (nominal)", "R1"),
    ("R1 nested", "H1H2", "0.05", "New gene families (round 1)", "NEST"),
    ("Panels nested (internal)", "Novelty (panels)", "top", "New gene families*", "NEST"),
    ("Panels nested (internal)", "Structural novelty (unadjusted)", "top", "No structural match, unadjusted*", "NEST"),
    ("Panels nested (internal)", "Cell-death families among dup families", "0.05", "Cell-death duplications (nominal)", "NEST"),
    ("Panels origins (internal)", "Novelty (panels)", "top", "New gene families*", "PORIG"),
    ("Panels origins (internal)", "Structure/sequence change ratio", "top", "Structure change per sequence change*", "PORIG"),
    ("Panels tips", "Structure/sequence change ratio", "top", "Structure change per sequence change*", "TIP"),
    ("Panels tips", "NEG CONTROL: ratio with shuffled 3Di", "top", "  same, shuffled-letter negative control*", "TIP"),
    ("Panels tips", "Novelty (panels, tip vs tip controls)", "top", "New gene families vs single-species controls*",
     "TIP"),
    ("Panels prokaryotes (internal)", "Novelty (panels)", "top", "New gene families*", "PROK"),
]
FL = pd.read_csv(f"{OUT}/family_ladder.tsv", sep="\t")
FS = pd.read_csv(f"{OUT}/family_ladder_sensitivity.tsv", sep="\t")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2})
fig = plt.figure(figsize=(11.2, 6.6))
axL = fig.add_axes([0.235, 0.09, 0.30, 0.83])
axE = fig.add_axes([0.665, 0.56, 0.32, 0.36])
axC = fig.add_axes([0.665, 0.09, 0.32, 0.36])

# ---- left ----
y = 0
yt, ylab, prev = [], [], None
xmax = 0
for lad, hyp, cr, lab, grp in ROWS:
    g = E[(E.ladder == lad) & (E.hypothesis == hyp)]
    if prev is not None and grp != prev:
        y += 0.8
    if grp != prev:
        axL.text(-0.4, y - 0.55, GROUPS[grp][1], ha="left", va="center", fontsize=8, fontweight="bold",
                 color=GROUPS[grp][0] if grp != "PROK" else INK2, transform=axL.transData, clip_on=False)
        y += 0.6
    prev = grp
    pis = g[f"pi_{cr}"].to_numpy(float)
    obs = int(g[f"sig_{cr}"].sum())
    n = len(g)
    xmax = max(xmax, n)
    cdf = np.cumsum(poibin_pmf(pis))
    lo, hi = int(np.searchsorted(cdf, 0.025)), int(np.searchsorted(cdf, 0.975))
    axL.plot([0, n], [y, y], color=GRID, lw=0.8, zorder=0)
    axL.add_patch(plt.Rectangle((lo - 0.3, y - 0.28), hi - lo + 0.6, 0.56, color="#d9d9d9", lw=0, zorder=1))
    axL.plot([pis.sum()] * 2, [y - 0.28, y + 0.28], color=INK2, lw=1.2, zorder=2)
    col = GROUPS[grp][0]
    axL.plot(obs, y, "o", ms=7, color=col, mec="white", mew=1.2, zorder=3)
    axL.text(n + 0.35, y, f"{obs}/{n}", va="center", fontsize=7, color=INK2)
    yt.append(y)
    ylab.append(lab)
    y += 1
axL.set_yticks(yt)
axL.set_yticklabels(ylab, fontsize=7.5)
axL.set_ylim(y - 0.2, -1.2)
axL.set_xlim(-0.5, xmax + 1.6)
axL.set_xticks(range(0, xmax + 1, 2))
axL.set_xlabel("Number of independent events with a significant result")
for s in ("top", "right", "left"):
    axL.spines[s].set_visible(False)
axL.tick_params(axis="y", length=0)
axL.set_title("a   Hypotheses: events significant vs chance", loc="left", fontsize=9, fontweight="bold", color=INK,
              x=-0.72)
# legend
axL.plot([], [], "o", color=INK2, ms=6, label="observed")
axL.plot([], [], "|", color=INK2, ms=9, mew=1.2, label="expected by chance")
axL.add_patch(plt.Rectangle((0, 0), 0, 0, color="#d9d9d9", label="95% chance range"))
axL.legend(loc="lower right", fontsize=7, frameon=False, bbox_to_anchor=(1.02, -0.005))

# ---- right ----
for ax, d, col, title in ((axE, "expanded", "#0072B2", "b   Gene families expanded in the multicellular clade"),
                          (axC, "contracted", "#D55E00", "c   Gene families contracted in the multicellular clade")):
    t = FL[(FL.ladder == "7 origins") & (FL.direction == d)].sort_values("tier_k")
    st = FS[(FS.analysis == "all 7 origins") & (FS.direction == d) &
            FS.null.str.startswith("detectability")].set_index("tier_k")
    k = t.tier_k.to_numpy()
    ax.fill_between(k, t.null_q025, t.null_q975, color="#d9d9d9", lw=0, step=None, zorder=1)
    ax.plot(k, t.null_mean, color=INK2, lw=1.2, zorder=2)
    sm = st.null_mean.reindex(k)
    ax.plot(k, sm, color=INK2, lw=1.0, ls=(0, (3, 2)), zorder=2)
    ax.plot(k, t.observed_n_families, "o", ms=7, color=col, mec="white", mew=1.2, zorder=3)
    for kk, o, p in zip(k, t.observed_n_families, t.p_perm):
        if kk <= 4:
            star = " *" if p < 0.05 else ""
            ax.annotate(f"{int(o)}{star}", (kk, max(o, 0.003)), xytext=(7, 4), textcoords="offset points", fontsize=7,
                        color=INK)
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_ylim(-0.02, 6000)
    ax.set_yticks([0, 0.1, 1, 10, 100, 1000])
    ax.set_yticklabels(["0", "0.1", "1", "10", "100", "1000"])
    ax.set_xticks(range(1, 8))
    ax.set_xlim(0.6, 7.4)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_ylabel("Families significant in ≥ k origins")
    ax.set_title(title, loc="left", fontsize=9, fontweight="bold", color=INK, x=-0.17)
axC.set_xlabel("k (number of the seven origins)")
axE.plot([], [], "o", color=INK2, ms=6, label="observed")
axE.plot([], [], color=INK2, lw=1.2, label="chance (family shuffle)")
axE.add_patch(plt.Rectangle((0, 0), 0, 0, color="#d9d9d9", label="95% chance range"))
axE.plot([], [], color=INK2, lw=1.0, ls=(0, (3, 2)), label="chance, size-matched shuffle")
axE.legend(loc="upper right", fontsize=7, frameon=False)
fig.savefig(f"{FIG}/fig8_recurrence_ladder.png", dpi=200)
fig.savefig(f"{FIG}/fig8_recurrence_ladder.pdf")
print("ok")
