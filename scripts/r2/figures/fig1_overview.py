"""Figure 1. The eight dense panels: species by state, coded events, and cladograms with event branches.
Sources: data/dense/<R>/{species.tsv,transitions.tsv,guide_tree.nwk}; results/r2/panels/<R>/event_table.tsv (M node, seq-novel counts).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from common import (P, MM, STATE, UNI, SIMPLE, COMPLEX, REVERSAL, CONFIRM, INK, INK2, GRID, PANELS, PANEL_NAMES,
                    P1_EVENTS, read_newick, layout_cladogram, save, fig_label)

KIND_ORDER = ["origin", "origin_complex", "nested", "reversal", "other"]
KIND_LAB = {"origin": "unicellular → simple (origin)", "origin_complex": "unicellular → complex (animals)",
            "nested": "simple → complex", "reversal": "reversal", "other": "other (complex → complex)"}


def kind_of(row):
    if row.kind in ("origin", "origin_complex", "nested", "reversal"):
        return row.kind
    if row.kind == "T_event":
        return "other"
    return None  # sensitivity (alternative) codings are not drawn


def event_colour(k, eid):
    if eid in P1_EVENTS:
        return CONFIRM
    return {"origin": SIMPLE, "origin_complex": COMPLEX, "nested": COMPLEX, "reversal": REVERSAL, "other": INK2}[k]


def load():
    sp, ev = {}, {}
    for R in PANELS:
        s = pd.read_csv(f"{P}/data/dense/{R}/species.tsv", sep="\t")
        sp[R] = s
        t = pd.read_csv(f"{P}/data/dense/{R}/transitions.tsv", sep="\t")
        e = pd.read_csv(f"{P}/results/r2/panels/{R}/event_table.tsv", sep="\t")[["event_id", "M", "seq_novel"]]
        t = t.drop(columns=["M"]).merge(e, on="event_id", how="left")
        t["k"] = t.apply(kind_of, axis=1)
        ev[R] = t
    return sp, ev


def draw_tree(ax, R, sp, ev):
    root = read_newick(f"{P}/data/dense/{R}/guide_tree.nwk")
    tips = layout_cladogram(root)
    state = dict(zip(sp[R].code, sp[R].mc_state))
    e = ev[R][ev[R].k.notna()]
    emap = {m: (k, eid, n) for m, k, eid, n in zip(e.M, e.k, e.event_id, e.seq_novel)}
    for n in root.walk():
        if n.children:
            ys = [c.y for c in n.children]
            ax.plot([n.x, n.x], [min(ys), max(ys)], color="#9a9a9a", lw=0.5, solid_capstyle="butt", zorder=1)
        if n.parent is not None:
            x0 = n.parent.x
            if n.name in emap:
                k, eid, nov = emap[n.name]
                col = event_colour(k, eid)
                ax.plot([x0, n.x], [n.y, n.y], color=col, lw=2.6 if eid in P1_EVENTS else 1.8,
                        solid_capstyle="butt", zorder=3)
            else:
                ax.plot([x0, n.x], [n.y, n.y], color="#9a9a9a", lw=0.5, solid_capstyle="butt", zorder=1)
    xs = [t.x for t in tips]
    for t in tips:
        ax.scatter(t.x + 0.35, t.y, s=7, color=STATE.get(state.get(t.name, ""), "white"), edgecolor="none", zorder=4)
    ax.set_xlim(-0.3, max(xs) + 0.8)
    ax.set_ylim(-1, len(tips))
    ax.axis("off")
    n_ev = int(e.k.notna().sum())
    ax.set_title(f"{PANEL_NAMES[R]}\n{len(tips)} proteomes, {n_ev} events", fontsize=6.5, color=INK, pad=2, loc="left")


def main():
    sp, ev = load()
    fig = plt.figure(figsize=(180 * MM, 178 * MM))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.78, 1.0, 1.0], hspace=0.36, wspace=0.12,
                          left=0.03, right=0.99, top=0.965, bottom=0.065)
    gtop = gs[0, :].subgridspec(1, 2, wspace=0.08, width_ratios=[1, 1.08])
    # ---- (a) species by state ----
    axa = fig.add_subplot(gtop[0, 0])
    axa.set_position([0.10, axa.get_position().y0, 0.37, axa.get_position().height])
    y = np.arange(len(PANELS))[::-1]
    for yi, R in zip(y, PANELS):
        c = sp[R].mc_state.value_counts()
        left = 0
        for st in ["unicellular", "simple", "complex"]:
            v = int(c.get(st, 0))
            if v:
                axa.barh(yi, v, left=left, color=STATE[st], height=0.72, edgecolor="white", linewidth=0.8)
                axa.text(left + v / 2, yi, str(v), ha="center", va="center", fontsize=5.8,
                         color="white" if st != "unicellular" else INK)
            left += v
        axa.text(left + 0.6, yi, f"{left}", va="center", fontsize=6, color=INK2)
    axa.set_yticks(y)
    axa.set_yticklabels([PANEL_NAMES[R] for R in PANELS])
    axa.set_xlabel("species (proteomes) per panel")
    axa.set_xlim(0, 58)
    axa.set_xticks(range(0, 45, 10)); axa.spines["bottom"].set_bounds(0, 40)
    axa.tick_params(axis="y", length=0)
    axa.spines["left"].set_visible(False)
    axa.legend(handles=[Patch(color=STATE[s], label=s) for s in STATE], loc="center right", ncol=1,
               bbox_to_anchor=(1.02, 0.5), handlelength=1.0, title="species state", title_fontsize=6.5)
    axa.legend_._legend_box.align = "left"
    axa.set_title("Species by multicellular state", loc="left", color=INK)
    # ---- (b) events by kind ----
    axb = fig.add_subplot(gtop[0, 1])
    axb.set_position([0.53, axb.get_position().y0, 0.46, axb.get_position().height])
    for yi, R in zip(y, PANELS):
        e = ev[R][ev[R].k.notna()]
        left = 0
        for k in KIND_ORDER:
            sub = e[e.k == k]
            if k == "nested":
                parts = [("conf", sub[sub.event_id.isin(P1_EVENTS)]), ("nested", sub[~sub.event_id.isin(P1_EVENTS)])]
            else:
                parts = [(k, sub)]
            for kk, ss in parts:
                v = len(ss)
                if not v:
                    continue
                col = CONFIRM if kk == "conf" else event_colour(kk, "")
                axb.barh(yi, v, left=left, color=col, height=0.72, edgecolor="white", linewidth=0.8,
                         hatch="////" if kk == "origin_complex" else None)
                axb.text(left + v / 2, yi, str(v), ha="center", va="center", fontsize=5.8,
                         color="white" if kk in ("origin_complex", "nested", "other") else INK)
                left += v
        axb.text(left + 0.2, yi, f"{left}", va="center", fontsize=6, color=INK2)
    plt.rcParams["hatch.color"] = "white"
    axb.set_yticks(y)
    axb.set_yticklabels([])
    axb.tick_params(axis="y", length=0)
    axb.spines["left"].set_visible(False)
    axb.set_xlabel("coded events per panel (alternative codings excluded)")
    axb.set_xlim(0, 29)
    axb.set_xticks(range(0, 16, 2)); axb.spines["bottom"].set_bounds(0, 14)
    hand = [Patch(color=SIMPLE, label=KIND_LAB["origin"]),
            Patch(facecolor=COMPLEX, hatch="////", edgecolor="white", label=KIND_LAB["origin_complex"]),
            Patch(color=COMPLEX, label=KIND_LAB["nested"]),
            Patch(color=CONFIRM, label="simple → complex, confirmatory"),
            Patch(color=REVERSAL, label=KIND_LAB["reversal"]),
            Patch(color=INK2, label=KIND_LAB["other"])]
    axb.legend(handles=hand, loc="center right", bbox_to_anchor=(1.02, 0.5), handlelength=1.0, fontsize=6)
    axb.set_title("Coded events by kind", loc="left", color=INK)
    # ---- (c) cladograms ----
    tree_axes = []
    for i, R in enumerate(PANELS):
        ax = fig.add_subplot(gs[1 + i // 4, i % 4])
        draw_tree(ax, R, sp, ev)
        tree_axes.append(ax)
    # shared legend for trees
    lh = [Line2D([], [], color=SIMPLE, lw=1.8, label="origin"),
          Line2D([], [], color=COMPLEX, lw=1.8, label="to complex (animals: from unicellular)"),
          Line2D([], [], color=CONFIRM, lw=2.6, label="confirmatory simple → complex"),
          Line2D([], [], color=REVERSAL, lw=1.8, label="reversal"),
          Line2D([], [], color=INK2, lw=1.8, label="other (complex → complex)"),
          Line2D([], [], marker="o", ls="", color=UNI, ms=3, label="unicellular species"),
          Line2D([], [], marker="o", ls="", color=SIMPLE, ms=3, label="simple species"),
          Line2D([], [], marker="o", ls="", color=COMPLEX, ms=3, label="complex species")]
    fig.legend(handles=lh, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0), handlelength=1.6, columnspacing=1.2)
    fig_label(fig, "a", 0.005, 0.995)
    fig_label(fig, "b", 0.49, 0.995)
    fig_label(fig, "c", 0.005, 0.705)
    save(fig, "fig1_overview")

    # numbers for the caption
    tot_sp = sum(len(sp[R]) for R in PANELS)
    tot_ev = sum(int(ev[R].k.notna().sum()) for R in PANELS)
    tot_all = sum(len(ev[R]) for R in PANELS)
    print("species", tot_sp, "events drawn", tot_ev, "all coded incl. alternatives", tot_all)
    for R in PANELS:
        print(R, sp[R].mc_state.value_counts().to_dict(), ev[R].k.value_counts(dropna=False).to_dict())


if __name__ == "__main__":
    main()
