#!/usr/bin/env python
"""Report figures 3-6 (results). Every number drawn is read from a committed result table; nothing is typed in.

  Fig 3  primary tests: per transition, the stem value against its matched control values, with u, T, p, Holm p
         results/primary/euk/{holm.tsv, primary_tests.tsv, per_transition_*.tsv, branch_values.tsv, matched_sets.tsv},
         results/report/h1h2_branch_residuals.tsv (H1/H2 residuals of all branches; scripts/report/export_h1h2_residuals.py),
         results/h9/{h9_transitions.tsv, h9_controls.tsv, h9_tests.tsv}
  Fig 4  gain uncertainty: (a) gain-evidence ladder by branch class, (b) per primary stem HOG gains, sequence-novel
         gains and the NB-GLM expectation, (c) inferred/true ratios on simulations
         results/gainladder/euk/branch_gain_ladder.tsv, results/primary/euk/per_transition_H1H2.tsv,
         results/validation_hogparser/metrics.tsv, results/calibration_tree/class_summary.tsv
  Fig 5  the 16 eukaryotic positive controls: expectation, outcome, failure category
         results/positive_controls/euk/{controls.tsv, failure_diagnosis.md}, data/positive_controls.tsv,
         data/toolkit_pfam.tsv + data/raw/pfam/Pfam-A.prereg_subset.hmm (Pfam names only)
  Fig 6  rates and flagellar modules: H9 r_i with 95% CI vs controls; H12 module presence at P and M
         results/h9/*, results/primary/euk/{h12_module_coverage.tsv, per_transition_H12.tsv},
         data/transitions.tsv (P and M nodes), data/h12_expected_direction_frozen.tsv

Each per-transition percentile u recomputed from the plotted control values is asserted equal to the committed u.

Usage: python scripts/report/fig_results.py [3 4 5 6]   -> figures/fig{3..6}_*.{png,svg}
"""
from __future__ import annotations

import os
import re
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

P = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = os.path.join(P, "results")
OUT = os.path.join(P, "figures")

# Style shared with fig_trees.py / fig_workflow.py (reference categorical slots 1-3, validated all-pairs).
INK, INK2, MUTED = "#1b1b1b", "#52514e", "#8f8d88"
GRID, FILL = "#e4e2dc", "#f4f3ef"
CTRL = "#b4b2ab"  # control values (recessive)
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "svg.fonttype": "none", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlelocation": "left",
})

# Transition order and plain names as in Figure 1 / Table 1.
TRANS = [("T_metazoa", "Animals"), ("T_fungi_hyphae", "Fungal hyphae"), ("T_dictyostelia", "Dictyostelids"),
         ("T_streptophyta", "Streptophytes"), ("T_rhodophyta", "Red algae"),
         ("T_px_filaments", "Brown-algal lineage (PX)"), ("T_volvocine_tgv", "Volvocine algae")]
TID = [t for t, _ in TRANS]
NAME = dict(TRANS)


def tsv(*parts, **kw):
    return pd.read_csv(os.path.join(R, *parts), sep="\t", **kw)


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote figures/{stem}.png/.svg")


def u_score(s, others):
    o = np.asarray(others, float)
    return float(((o < s).sum() + 0.5 * (o == s).sum() + 0.5) / (o.size + 1))


def fmt_p(p):
    return f"{p:.3f}" if p >= 0.001 else f"{p:.1e}"


def style_rows(ax, n, labels=True):
    ax.set_ylim(n - 0.5, -0.6)
    ax.set_yticks(range(n))
    ax.set_yticklabels([NAME[t] for t in TID] if labels else [])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    for i in range(n):
        ax.axhline(i, color=GRID, lw=0.6, zorder=0)


# =============================================================================
# Figure 3
# =============================================================================

def fig3_data():
    """{test: dict(rows=[(tid, stem_value, control_values, u, note)], xlabel, title, ref)}."""
    holm = tsv("primary", "euk", "holm.tsv").set_index("hypothesis")
    pt = tsv("primary", "euk", "primary_tests.tsv")
    pt = pt[(pt["set"] == "primary") & (pt["analysis"] == "primary")].set_index("hypothesis")
    bv = tsv("primary", "euk", "branch_values.tsv", index_col=0)
    ms = tsv("primary", "euk", "matched_sets.tsv")
    ms = ms[ms["set"] == "primary"]
    res = tsv("report", "h1h2_branch_residuals.tsv", index_col=0)
    out = {}
    spec = {
        "H1H2": ("per_transition_H1H2.tsv", res["resid_s"], "resid_s",
                 "NB residual  log[(y + 0.5) / (ŷ + 0.5)]", "H1/H2  Sequence-novel gains", 0.0, "two-sided"),
        "H4": ("per_transition_H4.tsv", bv["H4_s"], "y",
               "s = logit(class share of duplicated HOGs) − logit(share at parent)", "H4  Transcription-factor duplications", None,
               "one-sided, greater"),
        "H5": ("per_transition_H5.tsv", bv["H5_s"], "y",
               "s = logit(class share of duplicated HOGs) − logit(share at parent)", "H5  Adhesion/ECM + signalling duplications",
               None, "one-sided, greater"),
        "H12": ("per_transition_H12.tsv", bv["H12_f"], "y",
                "signed module loss  w · f  (w = +1 loss expected, −1 not)", "H12  Flagellar-module loss", 0.0,
                "one-sided, greater"),
    }
    for h, (fn, vals, stemcol, xl, title, ref, alt) in spec.items():
        per = tsv("primary", "euk", fn)
        per = per[per["set"] == "primary"].set_index("transition_id")
        rows = []
        for t in TID:
            if t not in per.index or pd.isna(per.at[t, "u"]):
                why = "not applicable (no flagellum at P)" if h == "H12" else "no value"
                rows.append((t, np.nan, [], np.nan, why))
                continue
            r = per.loc[t]
            w = float(r["w"]) if h == "H12" else 1.0
            ctrl = ms[(ms["hypothesis"] == h) & (ms["transition_id"] == t)]["control"]
            cv = w * vals.reindex(ctrl).dropna().to_numpy(float)
            sv = w * float(r[stemcol])
            assert np.isclose(sv, w * vals[r["stem"]], atol=1e-4), (h, t)
            u = u_score(sv, cv)
            assert np.isclose(u, r["u"], atol=1e-4), (h, t, u, r["u"])
            rows.append((t, sv, cv, float(r["u"]), ""))
        out[h] = dict(rows=rows, xlabel=xl, title=title, ref=ref, T=pt.at[h, "T_obs"], p=holm.at[h, "p"],
                      holm=holm.at[h, "holm_p"], n=int(holm.at[h, "n"]), k=int(holm.at[h, "k"]), alt=alt,
                      scheme=holm.at[h, "scheme"])
    # H9
    h9t = tsv("h9", "h9_transitions.tsv")
    h9t = h9t[h9t["set"] == "primary"].set_index("transition_id")
    h9c = tsv("h9", "h9_controls.tsv")["r_oriented"].to_numpy(float)
    h9s = tsv("h9", "h9_tests.tsv").set_index("analysis").loc["primary"]
    rows = []
    for t in TID:
        sv = float(h9t.at[t, "r"])
        u = u_score(sv, h9c)
        assert np.isclose(u, h9t.at[t, "u_percentile"], atol=1e-4), (t, u)
        rows.append((t, sv, h9c, u, ""))
    out["H9"] = dict(rows=rows, xlabel="log rate ratio r  (P→M vs P→sister; median over genes)",
                     title="H9  Rate of protein evolution", ref=0.0, T=h9s["T_obs"], p=holm.at["H9", "p"],
                     holm=holm.at["H9", "holm_p"], n=int(holm.at["H9", "n"]), k=int(holm.at["H9", "k"]),
                     alt="two-sided", scheme=holm.at["H9", "scheme"])
    assert np.isclose(h9s["p"], holm.at["H9", "p"], atol=1e-6)
    return out


def fig3():
    D = fig3_data()
    order = ["H1H2", "H4", "H5", "H12", "H9"]
    fig, axs = plt.subplots(2, 3, figsize=(14.2, 8.6), gridspec_kw=dict(hspace=0.55, wspace=0.42))
    jit = np.linspace(-0.22, 0.22, 20)
    for ax, h in zip(axs.flat, order):
        d = D[h]
        style_rows(ax, len(TID), labels=True)
        if d["ref"] is not None:
            ax.axvline(d["ref"], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
        tr = blended_transform_factory(ax.transAxes, ax.transData)
        for i, (t, sv, cv, u, note) in enumerate(d["rows"]):
            if np.isnan(sv):
                ax.text(0.5, i, note, transform=tr, ha="center", va="center", fontsize=7.6, color=MUTED,
                        style="italic")
                continue
            ax.scatter(cv, i + jit[:len(cv)], s=13, color=CTRL, edgecolor="none", zorder=2)
            ax.scatter([sv], [i], s=78, marker="D", color=BLUE, edgecolor="white", lw=1.0, zorder=4)
            ax.text(1.03, i, f"{u:.2f}", transform=tr, ha="left", va="center", fontsize=8, color=INK)
        ax.text(1.03, 1.0, "u", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, color=INK2,
                fontweight="bold")
        ax.set_xlabel(d["xlabel"], fontsize=7.8)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_title(d["title"], fontsize=10, fontweight="bold", color=INK, pad=30)
        ax.text(0, 1.06, f"T = {d['T']:.3f}   p = {fmt_p(d['p'])}   Holm p = {fmt_p(d['holm'])}   "
                f"(n = {d['n']}, k = {d['k']})", transform=ax.transAxes, fontsize=8.2, color=INK2, va="bottom")
        if ax not in axs[:, 0]:
            ax.set_yticklabels([])
    # key panel
    k = axs.flat[5]
    k.axis("off")
    hs = [Line2D([], [], ls="", marker="D", ms=8, color=BLUE, mec="white", label="transition stem"),
          Line2D([], [], ls="", marker="o", ms=4.5, color=CTRL, label="matched control branches"),
          Line2D([], [], color=MUTED, lw=0.8, ls=(0, (3, 2)), label="zero (no difference)")]
    k.legend(handles=hs, loc="upper left", frameon=False, fontsize=8.8, bbox_to_anchor=(0.0, 1.0))
    txt = ("u: mid-rank percentile of the stem among its\n"
           "own matched controls (0.5 = typical, 1 = above all).\n\n"
           "T: mean u over the n transitions; k controls\n"
           "per transition (k = 20 nearest of the 30-branch\n"
           "pool; H9: the 10 pool branches with two\n"
           "≥2-tip child clades, oriented at random).\n\n"
           "p: permutation of stem/control labels;\n"
           "H1/H2: parametric bootstrap (two-sided);\n"
           "H4, H5, H12 one-sided; H9 two-sided.\n"
           "Holm p: adjusted across the 5 primary tests.")
    k.text(0.0, 0.70, txt, transform=k.transAxes, va="top", fontsize=8.3, color=INK2, linespacing=1.35)
    fig.suptitle("Primary tests: stems against their matched controls; none is significant after Holm correction",
                 x=0.075, y=0.985, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    save(fig, "fig3_primary_tests")


# =============================================================================
# Figure 4
# =============================================================================

CLASS_LABEL = {"stem": "Transition stems", "unicellular_clade": "Within unicellular clades",
               "multicellular_clade": "Within multicellular clades", "mixed_clade": "Mixed clades",
               "root_child": "Root children"}


def fig4():
    lad = tsv("gainladder", "euk", "branch_gain_ladder.tsv")
    lad = lad[~lad["child_is_leaf"].astype(bool)]
    g = lad.groupby("branch_class").agg(n=("branch_child", "size"), hog=("hog_gains", "sum"),
                                        seq=("seq_novel", "sum"))
    g = g.reindex([c for c in CLASS_LABEL if c in g.index])
    tot_hog, tot_seq = int(lad["hog_gains"].sum()), int(lad["seq_novel"].sum())
    per = tsv("primary", "euk", "per_transition_H1H2.tsv")
    per = per[per["set"] == "primary"].set_index("transition_id").reindex(TID)
    met = tsv("validation_hogparser", "metrics.tsv")
    met = met[(met["scope"] == "internal_branches") & (met["preset"] == "fast")]
    cal = tsv("calibration_tree", "class_summary.tsv")
    cal = cal[cal["branch_class"].isin(["control_U", "primary_stem_edge"])]

    fig = plt.figure(figsize=(15.0, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.12, 0.9], wspace=0.55)

    # (a) ladder by branch class
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(g))
    h = 0.36
    ax.barh(y - h / 2, g["hog"], height=h, color=ORANGE, label="HOG gains (level 1)", zorder=3)
    ax.barh(y + h / 2, g["seq"], height=h, color=BLUE, label="sequence-novel gains (level 3)", zorder=3)
    xmax = g["hog"].max()
    for i, (c, r) in enumerate(g.iterrows()):
        ax.text(r["hog"] + xmax * 0.015, i - h / 2, f"{int(r['hog']):,}", va="center", fontsize=7.6, color=INK2)
        ax.text(r["seq"] + xmax * 0.015, i + h / 2, f"{int(r['seq']):,}  ({r['seq'] / r['hog']:.0%})",
                va="center", fontsize=7.6, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{CLASS_LABEL[c]}\n({int(n)} branches)" for c, n in g["n"].items()], fontsize=8.2)
    ax.set_ylim(len(g) - 0.4, -0.6)
    ax.set_xlim(0, xmax * 1.32)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("gains summed over internal branches")
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.legend(loc="center right", frameon=False, fontsize=7.8)
    ax.set_title("a   Gain-evidence ladder by branch class", fontsize=10, fontweight="bold", color=INK,
                 pad=32)
    ax.text(0, 1.02, f"{len(lad)} internal branches in the ladder table:\n{tot_hog:,} HOG gains, {tot_seq:,} sequence-novel "
            f"({tot_seq / tot_hog:.0%})", transform=ax.transAxes, fontsize=7.8, color=INK2, va="bottom")

    # (b) per primary stem
    ax = fig.add_subplot(gs[1])
    style_rows(ax, len(TID))
    for i, t in enumerate(TID):
        r = per.loc[t]
        ax.plot([r["seq_novel"], r["hog_gains"]], [i, i], color=MUTED, lw=1.2, zorder=2)
        ax.scatter(r["hog_gains"], i, s=46, color=ORANGE, edgecolor="white", lw=0.8, zorder=4)
        ax.scatter(r["seq_novel"], i, s=46, color=BLUE, edgecolor="white", lw=0.8, zorder=4)
        ax.scatter(r["yhat"], i, s=110, marker="|", color=INK, lw=2.0, zorder=5)
    ax.set_xscale("log")
    ax.set_xlim(8, 1500)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.grid(axis="x", which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("gains on the stem (log scale; path stems summed over edges)")
    hs = [Line2D([], [], ls="", marker="o", ms=6.5, color=ORANGE, mec="white", label="HOG gains"),
          Line2D([], [], ls="", marker="o", ms=6.5, color=BLUE, mec="white", label="sequence-novel gains (y)"),
          Line2D([], [], ls="", marker="|", ms=11, mew=2, color=INK, label="NB-GLM expectation (ŷ)")]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, fontsize=7.6,
              handletextpad=0.3, columnspacing=1.0)
    ax.set_title("b   Primary stems vs the control-branch model", fontsize=10, fontweight="bold", color=INK,
                 pad=32)
    ax.text(0, 1.02, "ŷ: negative-binomial model fitted on the 30 control branches", transform=ax.transAxes,
            fontsize=7.8, color=INK2, va="bottom")

    # (c) simulations
    ax = fig.add_subplot(gs[2])
    qs = [("gains", "Gains"), ("losses", "Losses"), ("dups_total", "Duplications")]
    sims = sorted(met["simulation"].unique())
    off = np.linspace(-0.12, 0.12, len(sims))
    for j, (q, _) in enumerate(qs):
        m = met[met["quantity"] == q].set_index("simulation")
        for o, s in zip(off, sims):
            ax.scatter(j - 0.14 + o, m.at[s, "ratio"], s=34, color=INK2, edgecolor="white", lw=0.6, zorder=4)
    cal_cols = {"losses": ("losses_true", "losses_inferred"), "dups_total": ("dups_true", "dups_inferred"),
                "gains": ("hog_gains_true", "hog_gains_inferred")}
    coff = np.linspace(-0.1, 0.1, len(cal))
    for j, (q, _) in enumerate(qs):
        tcol, icol = cal_cols[q]
        if q == "gains":
            tt, ii = cal[tcol].sum(), cal[icol].sum()
            ax.text(j + 0.2, 1.0 * 1.18, f"tree: {int(ii)} inferred\nvs {int(tt)} true\n(ratio undefined)",
                    ha="center", va="bottom", fontsize=6.9, color=INK2)
            continue
        for o, (_, r) in zip(coff, cal.iterrows()):
            ax.scatter(j + 0.18 + o, r[icol] / r[tcol], s=34, marker="s", facecolor="white", edgecolor=INK,
                       lw=1.0, zorder=4)
    ax.axhline(1, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_yscale("log")
    ax.set_ylim(0.8, 3.4)
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.8, 1, 1.5, 2, 2.5, 3]))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
    ax.set_xticks(range(len(qs)))
    ax.set_xticklabels([lab for _, lab in qs], fontsize=8.6)
    ax.set_xlim(-0.55, len(qs) - 0.45)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_ylabel("inferred / true events (log scale)")
    hs = [Line2D([], [], ls="", marker="o", ms=6, color=INK2, mec="white",
                 label=f"{len(sims)} simulations, internal branches"),
          Line2D([], [], ls="", marker="s", ms=6, mfc="white", mec=INK,
                 label="calibration on our tree:\nstem + control branches, 2 settings")]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=1, frameon=False, fontsize=7.6)
    ax.set_title("c   Gains are inflated most", fontsize=10, fontweight="bold", color=INK, pad=32)
    ax.text(0, 1.02, "Ancestra 'fast' preset vs simulated truth", transform=ax.transAxes, fontsize=7.8,
            color=INK2, va="bottom")
    fig.suptitle("Gain counts are the least reliable event type",
                 x=0.07, y=1.07, ha="left", fontsize=12.5, fontweight="bold", color=INK)
    save(fig, "fig4_gain_uncertainty")


# =============================================================================
# Figure 5
# =============================================================================

# Verdict codes used in results/positive_controls/euk/failure_diagnosis.md, grouped into the report's three causes.
VERDICT = {"iv": "Ancestra\ninference", "i": "Expectation\nor threshold", "ii": "Annotation\n(Pfam)"}
CATS = ["Ancestra\ninference", "Expectation\nor threshold", "Annotation\n(Pfam)"]
EVENT = {"gain": "gain", "co-option_present_at_P": "present at P", "co-option_present_at_P+loss":
         "present at P + loss", "co-option_present_at_P+loss+duplication": "present at P + loss + dup.",
         "expansion": "expansion", "duplication+gain": "duplication + gain", "null_no_expansion": "no expansion (null)"}


def pfam_names():
    names = {}
    tk = pd.read_csv(os.path.join(P, "data", "toolkit_pfam.tsv"), sep="\t")
    names.update(dict(zip(tk["pfam_acc"], tk["pfam_name"])))
    name = None
    with open(os.path.join(P, "data", "raw", "pfam", "Pfam-A.prereg_subset.hmm")) as f:
        for line in f:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("ACC "):
                acc = line.split()[1].split(".")[0]
                names[acc] = name
            elif line.startswith("DESC ") and name.startswith("DUF"):  # unnamed family: add its description
                names[acc] = f"{name} ({line[5:].strip().lower()})"
    return names


def failure_categories():
    md = open(os.path.join(R, "positive_controls", "euk", "failure_diagnosis.md")).read()
    out = {}
    for line in md.splitlines():
        m = re.match(r"\|\s*(C\d\d)\b[^|]*\|[^|]*\|[^|]*\|([^|]*)\|", line)
        if not m:
            continue
        codes = re.findall(r"\((i{1,3}|iv)\)", m.group(2))
        unknown = set(codes) - set(VERDICT)
        if unknown:
            sys.exit(f"unmapped verdict code(s) {unknown} for {m.group(1)}")
        out[m.group(1)] = {VERDICT[c] for c in codes}
    return out


def fig5():
    c = tsv("positive_controls", "euk", "controls.tsv").set_index("control")
    d = pd.read_csv(os.path.join(P, "data", "positive_controls.tsv"), sep="\t").set_index("control_id")
    d = d[d["run"] == "euk"].reindex(c.index)
    names = pfam_names()
    cats = failure_categories()
    n = len(c)
    fig = plt.figure(figsize=(14.0, 0.44 * n + 1.9))
    ax = fig.add_axes([0.01, 0.1, 0.98, 0.78])
    ax.set_xlim(0, 1)
    ax.set_ylim(n - 0.5, -1.4)
    ax.axis("off")
    X = dict(id=0.0, fam=0.035, ev=0.215, node=0.355, conf=0.665, out=0.735)
    cx = [0.84, 0.905, 0.97]
    hdr = dict(fontsize=8.6, fontweight="bold", color=INK, va="center")
    ax.text(X["id"], -0.95, "ID", **hdr)
    ax.text(X["fam"], -0.95, "Gene family (Pfam)", **hdr)
    ax.text(X["ev"], -0.95, "Expected event", **hdr)
    ax.text(X["node"], -0.95, "Expected node or branch", **hdr)
    ax.text(X["conf"], -0.95, "Confidence", **hdr)
    ax.text(X["out"], -0.95, "Outcome", **hdr)
    for x, lab in zip(cx, CATS):
        ax.text(x, -0.95, lab, ha="center", fontsize=8.2, fontweight="bold", color=INK, va="center",
                linespacing=1.1)
    ax.text((cx[0] + cx[-1]) / 2, -1.55, "Cause of failure", ha="center", fontsize=8.6, color=INK2, va="center")
    ax.plot([cx[0] - 0.03, cx[-1] + 0.03], [-1.38, -1.38], color=MUTED, lw=0.8)
    ax.plot([0, 1], [-0.5, -0.5], color=INK2, lw=0.9)
    for i, (cid, r) in enumerate(c.iterrows()):
        if i % 2 == 0:
            ax.add_patch(Rectangle((-0.005, i - 0.5), 1.01, 1, facecolor=FILL, edgecolor="none", zorder=0))
        acc = str(d.at[cid, "pfam_accessions"])
        if isinstance(d.at[cid, "anchor_protein"], str):
            fam = d.at[cid, "anchor_protein"].split(" (")[0].replace("Neurospora crassa ", "N. crassa ")
        elif acc.startswith("@toolkit"):
            fam = "toolkit classes: " + ", ".join(x.split(":")[1].replace("_", " ") for x in acc.split(","))
        else:
            fam = ", ".join("+".join(names.get(a, a) for a in grp.split("+")) for grp in acc.split(","))
        node = textwrap.fill(str(d.at[cid, "expected_node_description"]), 58)
        rec = str(r["recovered"]).strip().lower() == "yes"
        kw = dict(fontsize=7.9, color=INK, va="center")
        ax.text(X["id"], i, cid, fontweight="bold", **kw)
        ax.text(X["fam"], i, textwrap.fill(fam, 32), linespacing=1.05, **kw)
        ax.text(X["ev"], i, EVENT.get(d.at[cid, "expected_event"], d.at[cid, "expected_event"]), **kw)
        ax.text(X["node"], i, node, linespacing=1.05, **kw)
        ax.text(X["conf"], i, r["confidence"], **dict(kw, color=INK2))
        ax.scatter(X["out"] + 0.006, i, s=60, marker="o", facecolor=BLUE if rec else "white",
                   edgecolor=BLUE if rec else ORANGE, lw=1.4, zorder=3)
        ax.text(X["out"] + 0.018, i, "recovered" if rec else "not recovered", **kw)
        if not rec and cid not in cats:
            sys.exit(f"{cid} failed but has no row in failure_diagnosis.md")
        for x, lab in zip(cx, CATS):
            if lab in cats.get(cid, set()):
                ax.scatter(x, i, s=70, marker="s", color=ORANGE, edgecolor="white", zorder=3)
            elif not rec:
                ax.scatter(x, i, s=10, marker="o", color=GRID, zorder=2)
    n_rec = int((c["recovered"].str.strip().str.lower() == "yes").sum())
    hi_fail = [x for x, r in c.iterrows() if r["confidence"] == "high" and r["recovered"].strip().lower() != "yes"]
    fig.text(0.01, 0.955, f"Positive controls: {n_rec} of {n} expected events recovered; every failure traced "
             "to a cause", fontsize=12.5, fontweight="bold", color=INK, ha="left")
    fig.text(0.01, 0.915, f"Pre-registered criterion: ≥ 11 of {n} recovered and every high-confidence failure "
             f"({', '.join(hi_fail)}) explained.  P = parent node of a transition stem.", fontsize=8.4,
             color=INK2, ha="left")
    save(fig, "fig5_positive_controls")


# =============================================================================
# Figure 6
# =============================================================================

MOD_LABEL = {"IFT-B": "IFT-B", "BBSome": "BBSome", "axoneme": "Axoneme", "transition_zone": "Transition\nzone"}


def fig6():
    h9t = tsv("h9", "h9_transitions.tsv")
    h9t = h9t[h9t["set"] == "primary"].set_index("transition_id").reindex(TID)
    h9c = tsv("h9", "h9_controls.tsv")
    holm = tsv("primary", "euk", "holm.tsv").set_index("hypothesis")
    tr = pd.read_csv(os.path.join(P, "data", "transitions.tsv"), sep="\t").set_index("transition_id")
    cov = tsv("primary", "euk", "h12_module_coverage.tsv")
    cov = cov[cov["module"] != "module"]
    cov["present"] = cov["present"].astype(str) == "True"
    tot = cov.groupby("module")["n_components_total"].max()
    cov = cov.set_index(["level", "module"])
    exp = pd.read_csv(os.path.join(P, "data", "h12_expected_direction_frozen.tsv"), sep="\t").set_index(
        "transition_id")["expected"]
    p12 = tsv("primary", "euk", "per_transition_H12.tsv")
    p12 = p12[p12["set"] == "primary"].set_index("transition_id")
    mods = ["IFT-B", "BBSome", "axoneme", "transition_zone"]

    fig = plt.figure(figsize=(14.6, 5.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.2)
    # (a) H9
    ax = fig.add_subplot(gs[0])
    style_rows(ax, len(TID))
    ax.axvline(0, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    cv = h9c["r_oriented"].to_numpy(float)
    jit = np.linspace(-0.2, 0.2, len(cv))
    trf = blended_transform_factory(ax.transAxes, ax.transData)
    for i, t in enumerate(TID):
        r = h9t.loc[t]
        ax.scatter(cv, i + jit, s=13, color=CTRL, edgecolor="none", zorder=2)
        ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=BLUE, lw=2.2, solid_capstyle="round", zorder=3)
        ax.scatter(r["r"], i, s=70, marker="D", color=BLUE, edgecolor="white", lw=1.0, zorder=4)
        ax.text(1.03, i, f"{r['u_percentile']:.2f}", transform=trf, fontsize=8, va="center", color=INK)
        ax.text(1.12, i, f"{int(r['n_genes'])}", transform=trf, fontsize=8, va="center", color=INK2)
    ax.text(1.03, -0.95, "u", transform=trf, fontsize=8.5, fontweight="bold", va="center", color=INK2)
    ax.text(1.12, -0.95, "genes", transform=trf, fontsize=8.5, fontweight="bold", va="center", color=INK2)
    ax.set_xlabel("log rate ratio r  (P→M vs P→sister tips; median over BUSCO genes)")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    sec = ax.secondary_xaxis("top", functions=(np.exp, np.log))
    sec.set_xticks([0.7, 0.8, 0.9, 1.0, 1.2, 1.4])
    sec.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
    sec.tick_params(labelsize=7.6, colors=INK2)
    sec.spines["top"].set_color(MUTED)
    sec.set_xlabel("rate ratio", fontsize=7.8, color=INK2)
    hs = [Line2D([], [], color=BLUE, lw=2.2, marker="D", ms=7, mec="white", label="stem r with 95% CI"),
          Line2D([], [], ls="", marker="o", ms=4.5, color=CTRL, label=f"{len(cv)} control branches (oriented r)")]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=7.8)
    ax.set_title("a   H9: stem rate shifts fall within the control range", fontsize=10, fontweight="bold",
                 color=INK, pad=64)
    ax.text(0, 1.165, f"p = {fmt_p(holm.at['H9', 'p'])}, Holm p = {fmt_p(holm.at['H9', 'holm_p'])} "
            f"(n = {int(holm.at['H9', 'n'])}, k = {int(holm.at['H9', 'k'])})", transform=ax.transAxes,
            fontsize=8.2, color=INK2, va="bottom")

    # (b) H12 grid
    ax = fig.add_subplot(gs[1])
    ncol = 2 * len(mods)
    ax.set_xlim(-0.6, ncol + 2.6)
    ax.set_ylim(len(TID) - 0.5, -0.6)
    ax.axis("off")
    for j, m in enumerate(mods):
        x0 = 2 * j
        ax.text(x0 + 0.5, -1.3, MOD_LABEL[m], ha="center", va="center", fontsize=8.2, fontweight="bold",
                color=INK, linespacing=1.0)
        ax.text(x0, -0.75, "P", ha="center", va="center", fontsize=8, color=INK2)
        ax.text(x0 + 1, -0.75, "M", ha="center", va="center", fontsize=8, color=INK2)
    ax.text(ncol + 0.35, -1.05, "Expected\non stem", ha="center", va="center", fontsize=8.2, fontweight="bold",
            color=INK, linespacing=1.0)
    ax.text(ncol + 1.75, -1.05, "Loss f\n(H12 u)", ha="center", va="center", fontsize=8.2, fontweight="bold",
            color=INK, linespacing=1.0)
    for i, t in enumerate(TID):
        # rows share panel a's y-axis (same order and limits), so the transition names are not repeated
        Pn, Mn = tr.at[t, "P"], tr.at[t, "M"]
        for j, m in enumerate(mods):
            for k, node in enumerate((Pn, Mn)):
                x = 2 * j + k
                if (node, m) in cov.index:
                    nc, pres = int(cov.at[(node, m), "n_components"]), bool(cov.at[(node, m), "present"])
                else:
                    nc, pres = 0, False
                ax.add_patch(Rectangle((x - 0.44, i - 0.4), 0.88, 0.8, facecolor=BLUE if pres else "white",
                                       edgecolor=BLUE if pres else MUTED, lw=0.9))
                ax.text(x, i, f"{nc}/{int(tot[m])}", ha="center", va="center", fontsize=7.4,
                        color="white" if pres else INK2)
        e = str(exp.get(t, "")).replace("_", " ").replace("not applicable", "n/a")
        ax.text(ncol + 0.35, i, e, ha="center", va="center", fontsize=8, color=INK)
        if t in p12.index and pd.notna(p12.at[t, "u"]):
            ax.text(ncol + 1.75, i, f"{p12.at[t, 'f_b']:.2f} ({p12.at[t, 'u']:.2f})", ha="center", va="center",
                    fontsize=8, color=INK)
        else:
            ax.text(ncol + 1.75, i, "excluded", ha="center", va="center", fontsize=7.8, color=MUTED,
                    style="italic")
    hs = [Rectangle((0, 0), 1, 1, facecolor=BLUE, edgecolor=BLUE, label="module present (≥ 50% of components)"),
          Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=MUTED, label="absent")]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False, fontsize=7.8)
    ax.set_title("b   H12: flagellar modules at P and M", fontsize=10, fontweight="bold", color=INK, pad=64,
                 loc="left", x=0.0)
    ax.text(0.0, 1.165, f"p = {fmt_p(holm.at['H12', 'p'])}, Holm p = {fmt_p(holm.at['H12', 'holm_p'])} "
            f"(n = {int(holm.at['H12', 'n'])}, k = {int(holm.at['H12', 'k'])}); cell: components found / total",
            transform=ax.transAxes, fontsize=8.2, color=INK2, va="bottom")
    fig.suptitle("Rates and flagellar modules at the seven transitions", x=0.07, y=1.12, ha="left",
                 fontsize=12.5, fontweight="bold", color=INK)
    save(fig, "fig6_rates_flagella")


if __name__ == "__main__":
    which = sys.argv[1:] or ["3", "4", "5", "6"]
    for w in which:
        {"3": fig3, "4": fig4, "5": fig5, "6": fig6}[w]()
