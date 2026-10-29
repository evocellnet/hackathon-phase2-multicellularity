#!/usr/bin/env python
"""Report figure 7 (EXPLORATORY): exploratory results and what survived the post-hoc skeptical checks.

Every number drawn is read from a committed table; nothing is typed in.
  (a) new multi-domain Pfam architectures per primary stem: percentile u among matched controls (NB-GLM residual)
      under four species/model variants; leave-one-transition-out and completeness checks as text
      results/robustness/exploratory/{arch_per_transition.tsv, arch_checks.tsv}
      (base values = results/exploratory/euk/architectures/per_transition.tsv, asserted equal)
  (b) co-option index C (Pfam-expanded) per transition vs the matched-control median; subsampled to 3-5 species
      results/exploratory/euk/cooption/per_transition.tsv, results/robustness/exploratory/coopt_subsample_summary.tsv
  (c) convergent Pfam-domain loss: domains lost in >= 2 transitions, observed vs null mean, full data and after removals
      results/robustness/exploratory/{loss_convergence_checks.tsv, loss_convergence_nonciliary.tsv}
  (d) structure of sequence-novel gains: no usable model / structurally old, stems vs matched controls, Wilson CIs
      results/structure/{per_stem_swissprot.tsv, cmh_swissprot.tsv}

Usage: python scripts/report/fig_exploratory.py  -> figures/fig7_exploratory.{png,svg}
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_results import (BLUE, CTRL, GREEN, GRID, INK, INK2, MUTED, NAME, ORANGE, TID, fmt_p, save,  # noqa: E402
                         style_rows)

P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.path.join(P, "results")
ROB = os.path.join(R, "robustness", "exploratory")


def tsv(*parts, **kw):
    return pd.read_csv(os.path.join(*parts), sep="\t", **kw)


def head(ax, letter, title, sub, pad=30):
    ax.set_title(f"{letter}   {title}", fontsize=10, fontweight="bold", color=INK, pad=pad, loc="left")
    ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=7.9, color=INK2, va="bottom", linespacing=1.3)


# ----------------------------------------------------------------------------- (a)
def panel_a(ax):
    apt = tsv(ROB, "arch_per_transition.tsv")
    chk = tsv(ROB, "arch_checks.tsv", keep_default_na=False, na_values=[""])
    chk = chk[(chk["statistic"] == "new_arch_multi") & (chk["version"] == "residual")]
    base = tsv(R, "exploratory", "euk", "architectures", "per_transition.tsv")
    base = base[base["statistic"] == "new_arch_multi"].set_index("transition_id")
    a0 = apt[(apt["species_set"] == "all") & (apt["variant"] == "base_GLM")].set_index("transition")
    assert np.allclose(a0.loc[TID, "u_resid"], base.loc[TID, "u"]), "arch u mismatch"

    def test(sp, var, dropped=""):
        r = chk[(chk["species_set"] == sp) & (chk["variant"] == var) & (chk["dropped"].fillna("") == dropped)]
        assert len(r) == 1, (sp, var, dropped)
        return r.iloc[0]

    series = [("all", "base_GLM", "all species", BLUE, "D", 52),
              ("all", "GLM+mean_busco_desc", "+ BUSCO completeness covariate", ORANGE, "o", 34),
              ("genome_only", "base_GLM", "genome-based species only", GREEN, "^", 40),
              ("busco80", "base_GLM", "species with BUSCO ≥ 80% only", INK2, "s", 32)]
    style_rows(ax, len(TID))
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    offs = np.linspace(-0.24, 0.24, len(series))
    hs = []
    for (sp, var, lab, col, mk, sz), o in zip(series, offs):
        d = apt[(apt["species_set"] == sp) & (apt["variant"] == var)].set_index("transition").reindex(TID)
        for i, t in enumerate(TID):
            u = d.at[t, "u_resid"]
            if pd.isna(u):
                ax.text(0.03, i + o, "not evaluable", fontsize=6.6, color=MUTED, va="center", style="italic")
                continue
            ax.scatter(u, i + o, s=sz, marker=mk, color=col, edgecolor="white", lw=0.6, zorder=4)
        r = test(sp, var)
        hs.append(Line2D([], [], ls="", marker=mk, ms=6.5, color=col, mec="white",
                         label=f"{lab}:  T = {r['T_obs']:.2f}, p = {fmt_p(r['p'])} (n = {int(r['n'])})"))
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("u: percentile of the stem among its 20 matched controls")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, -0.11), frameon=False, fontsize=7.6,
              handletextpad=0.3)
    loo = chk[(chk["species_set"] == "all") & (chk["variant"] == "base_GLM") & (chk["dropped"].fillna("") != "")]
    b80m = test("busco80", "base_GLM", "T_metazoa")
    head(ax, "a", "New domain architectures: survives some checks, not completeness",
         f"Leave one transition out (all species): every p ≤ {loo['p'].max():.3f}.  "
         f"BUSCO ≥ 80% without Animals: p = {b80m['p']:.2f}.\n"
         "Statistic: residual count of new multi-domain Pfam architectures (NB-GLM on control branches)")


# ----------------------------------------------------------------------------- (b)
def panel_b(ax):
    pt = tsv(R, "exploratory", "euk", "cooption", "per_transition.tsv")
    pf = pt[pt["version"] == "pfam_expanded"].set_index("transition").reindex(TID)
    hg = pt[pt["version"] == "hog_toolkit"].set_index("transition").reindex(TID)
    sub = tsv(ROB, "coopt_subsample_summary.tsv")
    sub = sub[sub["version"] == "pfam_expanded"]
    style_rows(ax, len(TID))
    trf = blended_transform_factory(ax.transAxes, ax.transData)
    ms = sorted(sub["m"].unique())
    for i, t in enumerate(TID):
        r = pf.loc[t]
        ax.scatter(r["median_C_controls"], i, s=120, marker="|", color=INK, lw=2.0, zorder=3)
        ax.scatter(r["C"], i, s=56, marker="D", color=BLUE, edgecolor="white", lw=0.8, zorder=4)
        ax.text(1.03, i, f"{r['u']:.2f}", transform=trf, fontsize=8, va="center", color=INK)
        ax.text(1.15, i, f"{hg.at[t, 'u']:.2f}", transform=trf, fontsize=8, va="center", color=INK2)
        s = sub[sub["transition"] == t].set_index("m")
        for j, m in enumerate(ms):
            if m not in s.index:
                continue
            y = i + 0.17 + 0.1 * j
            ax.plot([s.at[m, "C_q025"], s.at[m, "C_q975"]], [y, y], color=ORANGE, lw=1.4, zorder=3)
            ax.scatter(s.at[m, "C_median"], y, s=16, color=ORANGE, zorder=4)
            ax.scatter(s.at[m, "median_C_controls_size_matched"], y, s=40, marker="|", color=MUTED, lw=1.4,
                       zorder=3)
    ax.text(1.03, -0.95, "u", transform=trf, fontsize=8.3, fontweight="bold", color=INK2, va="center")
    ax.text(1.15, -0.95, "u HOG", transform=trf, fontsize=8.3, fontweight="bold", color=INK2, va="center")
    ax.set_xlabel("co-option index C (share of expanded Pfam families already present at P)")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    nsub = int(sub["n_sub"].iloc[0])
    hs = [Line2D([], [], ls="", marker="D", ms=6.5, color=BLUE, mec="white", label="transition, all species"),
          Line2D([], [], ls="", marker="|", ms=10, mew=2, color=INK, label="median of 20 matched controls"),
          Line2D([], [], color=ORANGE, lw=1.4, marker="o", ms=4,
                 label=f"subsampled to {ms[0]}–{ms[-1]} species (median, 95% of {nsub} draws)"),
          Line2D([], [], ls="", marker="|", ms=8, mew=1.4, color=MUTED, label="size-matched control median")]
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, -0.11), frameon=False, fontsize=7.6, ncol=2,
              columnspacing=1.0, handletextpad=0.4)
    head(ax, "b", "Co-option index: lower in Animals and Streptophytes, survives subsampling",
         "Pfam-based index; u HOG = the same index from HOG toolkit families (no effect there).\n"
         "Lower C = more expanded families new at the stem rather than co-opted from P")


# ----------------------------------------------------------------------------- (c)
def panel_c(ax):
    c = tsv(ROB, "loss_convergence_checks.tsv", keep_default_na=False, na_values=[""])
    nc = tsv(ROB, "loss_convergence_nonciliary.tsv", keep_default_na=False, na_values=[""])

    def row(analysis, dropped, sp):
        r = c[(c["analysis"] == analysis) & (c["dropped"].fillna("") == dropped) & (c["species_set"] == sp)]
        assert len(r) == 1, (analysis, dropped, sp)
        r = r.iloc[0]
        return r["G2_obs"], r["G2_permute_null_mean"], r["G2_permute_p"], int(r["n_transitions"])

    ncf = nc[nc["dropped"].fillna("") == ""].iloc[0]
    rows = [("All 7 transitions, all species", *row("baseline_reproduction", "", "all")),
            ("Without Dictyostelids", *row("LOTO", "T_dictyostelia", "all")),
            ("Without Dictyostelids and Red algae", *row("LOTO2", "T_dictyostelia;T_rhodophyta", "all")),
            ("Genome-based species only", *row("species_filter", "T_dictyostelia", "genome_only")),
            ("Genomes with BUSCO ≥ 80% only", *row("species_filter", "T_dictyostelia,T_rhodophyta",
                                                  "genome_busco80_hq")),
            (f"Cilium-linked domains removed ({int(ncf['n_domains_excluded'])})", ncf["G2_obs"],
             ncf["G2_permute_null_mean"], ncf["G2_permute_p"], int(ncf["n_transitions"])),
            ("  … and null keeps each transition's loss count", ncf["G2_obs"], ncf["fm_uniform_null_mean"],
             ncf["fm_uniform_p"], int(ncf["n_transitions"]))]
    n = len(rows)
    ax.set_ylim(n - 0.5, -0.6)
    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{r[0]}  (n = {r[4]})" for r in rows], fontsize=8)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    trf = blended_transform_factory(ax.transAxes, ax.transData)
    for i, (lab, obs, null, p, nt) in enumerate(rows):
        ax.axhline(i, color=GRID, lw=0.6, zorder=0)
        ax.plot([null, obs], [i, i], color=MUTED, lw=1.3, zorder=2)
        ax.scatter(null, i, s=48, facecolor="white", edgecolor=INK2, lw=1.2, zorder=3)
        ax.scatter(obs, i, s=56, marker="D", color=BLUE, edgecolor="white", lw=0.8, zorder=4)
        ax.text(1.03, i, fmt_p(p), transform=trf, fontsize=8, va="center", color=INK)
    ax.text(1.03, -0.95, "p", transform=trf, fontsize=8.3, fontweight="bold", color=INK2, va="center")
    ax.set_xlim(0, None)
    ax.set_xlabel("Pfam domains lost on ≥ 2 transition stems")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    hs = [Line2D([], [], ls="", marker="D", ms=6.5, color=BLUE, mec="white", label="observed"),
          Line2D([], [], ls="", marker="o", ms=6.5, mfc="white", mec=INK2, label="null mean (permutation)")]
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, -0.11), frameon=False, fontsize=7.6, ncol=2)
    head(ax, "c", "Convergent domain loss: the excess disappears",
         "The excess comes from Dictyostelids, transcriptome proteomes and cilium domains,\n"
         "i.e. the known flagellum loss already tested in H12")


# ----------------------------------------------------------------------------- (d)
def panel_d(ax1, ax2):
    ps = tsv(R, "structure", "per_stem_swissprot.tsv")
    cmh = tsv(R, "structure", "cmh_swissprot.tsv").set_index("measure")
    order = TID + ["POOLED"]
    for ax, meas, lab, letter in ((ax1, "frac_no_usable_model", "no usable model (%)", "d"),
                                  (ax2, "frac_old_of_searched", "structurally old (% of usable)", "")):
        d = ps[ps["measure"] == meas].set_index("transition_id").reindex(order)
        n = len(order)
        ax.set_ylim(n - 0.5, -0.6)
        ax.set_yticks(range(n))
        ax.set_yticklabels([NAME.get(t, "Pooled") for t in order] if ax is ax1 else [])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        for i, t in enumerate(order):
            r = d.loc[t]
            ax.axhline(i, color=GRID, lw=0.6, zorder=0)
            bold = t == "POOLED"
            for o, pre, col, mk in ((-0.15, "stem", BLUE, "D"), (0.15, "control", CTRL, "o")):
                ax.plot([100 * r[f"{pre}_ci_lo"], 100 * r[f"{pre}_ci_hi"]], [i + o, i + o], color=col,
                        lw=2.4 if bold else 1.5, solid_capstyle="round", zorder=2)
                ax.scatter(100 * r[f"{pre}_frac"], i + o, s=46 if bold else 30, marker=mk,
                           color=col if pre == "stem" else INK2, edgecolor="white", lw=0.6, zorder=3)
            if ax is ax2:
                ax.text(101, i, f"{int(r['stem_k'])}/{int(r['stem_n'])} vs "
                        f"{int(r['control_k'])}/{int(r['control_n'])}", fontsize=6.8, color=INK2, va="center")
        ax.set_xlim(0, 100)
        ax.set_xlabel(lab)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        m = cmh.loc[meas]
        ax.text(0, 1.02, f"Mantel–Haenszel OR {m['mh_OR']:.2f} ({m['mh_ci_lo']:.2f}–{m['mh_ci_hi']:.2f}), "
                f"p = {fmt_p(m['cmh_p'])}", transform=ax.transAxes, fontsize=7.6, color=INK2, va="bottom")
    ax1.set_title("d   Structure of sequence-novel gains: few testable, stems not older", fontsize=10,
                  fontweight="bold", color=INK, pad=30, loc="left")
    hs = [Line2D([], [], color=BLUE, lw=1.5, marker="D", ms=5.5, mec="white", label="stem gains"),
          Line2D([], [], color=CTRL, lw=1.5, marker="o", ms=5.5, mfc=INK2, mec="white",
                 label="matched-control gains"),
          ]
    ax1.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.02, -0.11), frameon=False, fontsize=7.6, ncol=2)
    nper = ps[(ps["measure"] == "frac_no_usable_model") & (ps["transition_id"] != "POOLED")]["stem_n"].unique()
    ax2.text(0, -0.15, f"lines: Wilson 95% CI; numbers: stem k/n vs control k/n\n"
             f"({'/'.join(str(int(x)) for x in nper)} sampled gains per stem and per control set)", transform=ax2.transAxes,
             fontsize=7.4, color=INK2, va="top")


def main():
    fig = plt.figure(figsize=(16.5, 12.6))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.42, left=0.13, right=0.95, top=0.86, bottom=0.1)
    panel_a(fig.add_subplot(gs[0, 0]))
    panel_b(fig.add_subplot(gs[0, 1]))
    panel_c(fig.add_subplot(gs[1, 0]))
    sg = gs[1, 1].subgridspec(1, 2, wspace=0.34)
    panel_d(fig.add_subplot(sg[0]), fig.add_subplot(sg[1]))
    fig.add_artist(FancyBboxPatch((0.012, 0.935), 0.105, 0.032, boxstyle="round,pad=0.004,rounding_size=0.006",
                                  transform=fig.transFigure, facecolor=INK2, edgecolor="none"))
    fig.text(0.0645, 0.951, "EXPLORATORY", ha="center", va="center", fontsize=11, fontweight="bold", color="white")
    fig.text(0.128, 0.951, "Exploratory results and what survived skeptical checks", fontsize=13.5,
             fontweight="bold", color=INK, va="center")
    fig.text(0.128, 0.922, "Not pre-registered and not in the Holm family; the checks were designed after seeing "
             "these results, so they can only weaken a claim.", fontsize=8.8, color=INK2, va="center")
    save(fig, "fig7_exploratory")


if __name__ == "__main__":
    main()
