#!/usr/bin/env python
"""Figures (PNG + PDF) for a pooled tag: per-panel trees with event branches and sequence-novel counts,
forest plot of event rate ratios (P1 highlighted), simple->complex recurrent-domain heatmap.

    figures.py --tag <pooled tag>        # -> results/r2/panels/pooled/<tag>/figures/
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import OUT, log  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from mcell.phylo import read_tree  # noqa: E402

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#9a998f", "#e4e3dc"
KIND_COL = {"nested": "#2a78d6", "origin": "#1baf7a", "origin_complex": "#1baf7a", "reversal": "#eda100",
            "sensitivity": "#9a998f", "T_event": "#4a3aa7"}
P1_COL = "#eb6834"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "text.color": INK, "font.family": "DejaVu Sans"})


def save(fig, d, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(d, f"{name}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def tree_fig(R, br, ER, d):
    t = read_tree(os.path.join(OUT, R, "tree_bl.nwk"))
    b = br[br["panel"] == R].set_index("branch_child")
    ev = ER[ER["panel"] == R]
    evk = {r.M: r for r in ev.itertuples()}
    leaves = t.get_leaf_names()
    y = {n: i for i, n in enumerate(leaves)}
    x = {}
    for n in t.traverse("preorder"):
        x[n.name] = 0.0 if n.up is None else x[n.up.name] + n.dist
    for n in t.traverse("postorder"):
        if not n.is_leaf():
            y[n.name] = np.mean([y[c.name] for c in n.children])
    xmax = max(x.values())
    fig, ax = plt.subplots(figsize=(8, 0.22 * len(leaves) + 1.2))
    for n in t.traverse():
        if n.up is None:
            continue
        e = evk.get(n.name)
        col, lw = INK2, 0.9
        if e is not None:
            col = P1_COL if e.p1_role == "confirmatory" else KIND_COL.get(e.kind, INK2)
            lw = 2.4
        elif n.name in b.index and b.at[n.name, "inside_crown"]:
            col = MUTED
        ax.plot([x[n.up.name], x[n.name]], [y[n.name], y[n.name]], color=col, lw=lw, solid_capstyle="round")
        if not n.is_leaf():
            ys = [y[c.name] for c in n.children]
            ax.plot([x[n.name], x[n.name]], [min(ys), max(ys)], color=INK2, lw=0.8)
        if n.name in b.index and np.isfinite(b.at[n.name, "seq_novel"]):
            ax.text((x[n.up.name] + x[n.name]) / 2, y[n.name] - 0.18, f"{int(b.at[n.name, 'seq_novel'])}",
                    ha="center", va="bottom", fontsize=6, color=INK if e is not None else INK2)
    for lf in leaves:
        ax.text(x[lf] + xmax * 0.01, y[lf], lf, va="center", fontsize=7)
    ax.set_ylim(len(leaves), -1)
    ax.set_xlim(-xmax * 0.02, xmax * 1.25)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("substitutions per site (LG+G4, fixed topology)")
    handles = [plt.Line2D([], [], color=P1_COL, lw=2.4, label="P1 confirmatory event")]
    for k, c in (("nested", "simple→complex"), ("origin", "origin"), ("reversal", "reversal"),
                 ("sensitivity", "sensitivity coding"), ("T_event", "other (T_event)")):
        if (ev["kind"] == k).any() or (k == "origin" and ev["kind"].str.startswith("origin").any()):
            handles.append(plt.Line2D([], [], color=KIND_COL[k], lw=2.4, label=c))
    handles.append(plt.Line2D([], [], color=MUTED, lw=0.9, label="inside multicellular crown"))
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=7)
    ax.set_title(f"{R}: sequence-novel gains per branch (numbers); event branches coloured", loc="left", fontsize=9)
    save(fig, d, f"tree_{R}")


def forest(ER, meta, d):
    e = ER[ER["log_rr"].notna() & ER["kind"].isin(["nested", "origin", "origin_complex", "reversal"])].copy()
    order = {"nested": 0, "origin_complex": 1, "origin": 2, "reversal": 3}
    e["o"] = e["kind"].map(order)
    e = e.sort_values(["o", "panel", "event_id"])
    rows = list(e.itertuples())
    pooled = [(k, lab) for k, lab in (("p1_meta", "P1 confirmatory (mod. HK)"), ("secondary_all_nested", "all simple→complex"),
                                      ("secondary_nested_internal_only", "simple→complex, internal only"),
                                      ("secondary_one_per_panel", "one effect per panel"), ("secondary_T_only", "T-only coding"),
                                      ("descriptive_origins", "origins (descriptive)")) if meta.get(k)]
    n = len(rows) + len(pooled) + 1
    fig, ax = plt.subplots(figsize=(7, 0.2 * n + 1.2))
    for i, r in enumerate(rows):
        col = P1_COL if r.p1_role == "confirmatory" else KIND_COL.get(r.kind, INK2)
        ax.plot([r.rr_low, r.rr_high], [i, i], color=col, lw=1.4, solid_capstyle="round")
        if r.is_tip:   # hollow marker: tip branch (species-specific orphans; extrapolated from internal controls)
            ax.plot([r.rr], [i], "o", ms=4.5, mfc="white", mec=col, mew=1.2)
        else:
            ax.plot([r.rr], [i], "o", ms=4.5 if r.p1_role == "confirmatory" else 3.5, color=col, mec="white", mew=0.8)
    yt = [f"{r.event_id}{' (tip)' if r.is_tip else ''}" for r in rows]
    base = len(rows) + 1
    for j, (k, lab) in enumerate(pooled):
        m = meta[k]
        yy = base + j
        ax.fill([m["RR_low"], m["RR"], m["RR_high"], m["RR"]], [yy, yy - 0.3, yy, yy + 0.3], color=P1_COL if k == "p1_meta" else INK2)
        if m.get("PI_RR_low") not in (None, "None") and m.get("PI_RR_low") == m.get("PI_RR_low"):
            ax.plot([m["PI_RR_low"], m["PI_RR_high"]], [yy, yy], color=MUTED, lw=0.8, ls=":")
        yt.append(f"{lab}, k={m['k']}")
    ax.set_yticks(list(range(len(rows))) + [base + j for j in range(len(pooled))])
    ax.set_yticklabels(yt, fontsize=6.5)
    ax.axvline(1, color=MUTED, lw=0.8)
    ax.set_xscale("log")
    ax.set_ylim(n, -1)
    ax.grid(axis="x", color=GRID, lw=0.5)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("rate ratio: observed / expected sequence-novel gains (95% CI; pooled NB GLM on controls)")
    ax.set_title("Event rate ratios (orange = P1 confirmatory; blue = simple→complex; green = origins; yellow = reversals;\n"
                 "hollow = tip branch, extrapolated from internal controls)", loc="left", fontsize=9)
    save(fig, d, "forest_event_rate_ratios")


def heatmap(pdir, d):
    f = os.path.join(pdir, "w3", "w3_heatmap_matrix.tsv")
    if not os.path.exists(f):
        return
    M = pd.read_csv(f, sep="\t", index_col=0)
    if M.empty:
        return
    names = M.pop("name")
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#f3f2ee", "#7fb2ec", "#eda100", "#2a78d6"][: 4])
    fig, ax = plt.subplots(figsize=(0.45 * M.shape[1] + 3, 0.18 * M.shape[0] + 1.5))
    ax.imshow(M.to_numpy(), cmap=cmap, vmin=0, vmax=3, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(M.columns, rotation=60, ha="right", fontsize=6.5)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels([f"{n} ({i})" for i, n in zip(M.index, names)], fontsize=6)
    ax.set_xticks(np.arange(-.5, M.shape[1]), minor=True)
    ax.set_yticks(np.arange(-.5, M.shape[0]), minor=True)
    ax.grid(which="minor", color="white", lw=1.5)
    ax.tick_params(which="minor", length=0)
    hs = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ("#7fb2ec", "#eda100", "#2a78d6")]
    ax.legend(hs, ["sequence-novel family", "duplicated HOG", "both"], loc="upper left", bbox_to_anchor=(1.01, 1),
              frameon=False, fontsize=7)
    ax.set_title("Pfam domains recurrent (≥2 events) in simple→complex events", loc="left", fontsize=9)
    save(fig, d, "heatmap_recurrent_domains")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    a = ap.parse_args(argv)
    pdir = os.path.join(OUT, "pooled", a.tag)
    d = os.path.join(pdir, "figures")
    if os.path.exists(d):
        raise SystemExit(f"refusing to overwrite {d}")
    os.makedirs(d)
    br = pd.read_csv(os.path.join(pdir, "branch_residuals.tsv"), sep="\t", index_col=0)
    ER = pd.read_csv(os.path.join(pdir, "event_residuals.tsv"), sep="\t", keep_default_na=True)
    ER["p1_role"] = ER["p1_role"].fillna("")
    meta = json.load(open(os.path.join(pdir, "p1_and_meta.json")))
    for R in sorted(br["panel"].unique()):
        tree_fig(R, br, ER, d)
    forest(ER, meta, d)
    heatmap(pdir, d)
    log(f"figures in {d}")


if __name__ == "__main__":
    main()
