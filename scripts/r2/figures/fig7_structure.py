"""Figure 7. Protein structure predicted from sequence: accuracy, use inside Ancestra, structural novelty, divergence.
Sources:
 (a) results/r2/w5_3di_validation/per_protein_accuracy.tsv, accuracy_meta.txt
 (b) results/r2/w5/r1/xval_agreement_summary.tsv, controls_{default,3di,3di_guard}.tsv
 (c) results/r2/verification/w5/work/glmm_results.txt, glmm_results_exh.txt (independent verifier's mixed models)
 (d) results/r2/w6/panel_tests.tsv
Colours: Okabe-Ito (validated set in common.py).
"""
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import P, MM, INK, INK2, GRID, save, panel_label

V = f"{P}/results/r2/w5_3di_validation"
W5 = f"{P}/results/r2/w5/r1"
VER = f"{P}/results/r2/verification/w5/work"
W6 = f"{P}/results/r2/w6"

C_DEF = "#555555"     # default / unadjusted / real-data neutral
C_STR = "#56B4E9"     # structure mode
C_GRD = "#0072B2"     # structure mode + guard
C_ADJ = "#0072B2"     # fully adjusted
C_EXH = "#E69F00"     # stricter (exhaustive) search
C_BR = "#CC79A7"      # branch-level
C_REAL = "#0072B2"
C_SHUF = "#E69F00"
C_SPLIT = "#BDBDBD"


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def panel_a(ax):
    d = pd.read_csv(f"{V}/per_protein_accuracy.tsv", sep="\t")
    meta = dict(l.strip().split("\t") for l in open(f"{V}/accuracy_meta.txt") if "\t" in l)
    chance = float(meta["independent_baseline_identity"])
    major = float(meta["most_common_state_frequency"])
    edges = [50, 150, 300, 500, 1000, 3000]
    labs = ["50–150", "150–300", "300–500", "500–1,000", "1,000–3,000\n(windowed)"]
    x = np.arange(len(labs))
    for src, col, name, off in [("sp", C_GRD, "Swiss-Prot proteins (mostly animals)", -0.09),
                                ("pdb", C_EXH, "other protist proteins", 0.09)]:
        ys, ns = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = d[(d.src == src) & (d.L > lo) & (d.L <= hi)] if lo > 50 else d[(d.src == src) & (d.L >= lo) & (d.L <= hi)]
            ys.append(s.nagree.sum() / s.nvalid.sum()); ns.append(len(s))
        ax.plot(x + off, ys, "-o", color=col, ms=4.5, lw=1.4, label=name, mec="white", mew=0.6)
    allv = d.nagree.sum() / d.nvalid.sum()
    ax.axhline(major, color=INK2, lw=0.8, ls="--")
    ax.axhline(chance, color=INK2, lw=0.8, ls=":")
    ax.text(len(labs) - 0.6, major + 0.012, f"always the commonest letter ({major:.2f})", ha="right", fontsize=6, color=INK2)
    ax.text(len(labs) - 0.6, chance + 0.012, f"chance from letter frequencies ({chance:.2f})", ha="right", fontsize=6, color=INK2)
    ax.axhspan(0.41, 0.65, color="#e9e9e9", zorder=0, lw=0)
    ax.text(-0.35, 0.425, "published range for\nsequence-predicted letters", fontsize=5.8, color=INK2, va="bottom")
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=6)
    ax.set_xlim(-0.45, len(labs) - 0.55)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel("protein length (amino acids)")
    ax.set_ylabel("share of residues with the same\nstructural letter as AlphaFold")
    ax.set_title(f"Structure letters predicted from sequence\nvs AlphaFold-derived (overall {allv:.3f}; 4,981 proteins)", loc="left", color=INK)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), handlelength=1.6)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)


def panel_b(ax):
    s = pd.read_csv(f"{W5}/xval_agreement_summary.tsv", sep="\t")
    s = s[(s.direction == "ancestra->eukaryoma") & (s.event == "gain")]
    groups = [("Metazoa_stem", "animal\nstem"), ("Choanozoa_stem", "stem of\nanimals +\nchoanoflag."),
              ("both_stems", "both\nstems"), ("other_internal", "all other\ninternal\nbranches")]
    variants = [("default", "sequence only (default)", C_DEF), ("3di", "with structure", C_STR),
                ("3di_guard", "structure + guard", C_GRD)]
    x = np.arange(len(groups)); w = 0.25
    for j, (v, lab, col) in enumerate(variants):
        for i, (g, _) in enumerate(groups):
            r = s[(s.variant == v) & (s.group == g)].iloc[0]
            p = r.same / r.n_mapped; lo, hi = wilson(r.same, r.n_mapped)
            ax.bar(x[i] + (j - 1) * w, p * 100, width=w * 0.92, color=col, label=lab if i == 0 else None, edgecolor="white", lw=0.5)
            ax.errorbar(x[i] + (j - 1) * w, p * 100, yerr=[[100 * (p - lo)], [100 * (hi - p)]], color=INK, lw=0.7, capsize=0)
    ax.set_xticks(x); ax.set_xticklabels([g[1] for g in groups], fontsize=6)
    ax.set_ylabel("new families placed on the same\nbranch by EukaryOMA (%)")
    ax.set_ylim(0, 52)
    ax.set_title("Structure mode inside Ancestra\n(animals-and-relatives reconstruction)", loc="left", color=INK)
    # deep-node positive controls
    rec = {}
    for v, _, _ in variants:
        c = pd.read_csv(f"{W5}/controls_{v}.tsv", sep="\t")
        dn = c[c.deep_node == True]  # noqa: E712
        rec[v] = (int((dn.recovered == "yes").sum()), len(dn))
    txt = "Deep-node positive controls recovered:\n" + "   ".join(f"{lab.split(' (')[0]} {rec[v][0]}/{rec[v][1]}" for v, lab, _ in variants)
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, fontsize=5.8, color=INK, va="top")
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.87), handlelength=1.2, ncol=2, columnspacing=0.8, fontsize=5.6)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)


def parse_glmm(path):
    t = open(path).read()
    num = r"([\d.]+) \[([\d.]+), ([\d.]+)\]"
    out = {}
    out[("all", "unadj")] = re.search(r"GLMM unadjusted \(random.*?event OR " + num, t).groups()
    out[("all", "full")] = re.search(r"GLMM full: event OR " + num, t).groups()
    m = re.search(r"GLMM full by kind: nested OR " + num + r".*?origin OR " + num, t).groups()
    out[("nested", "full")], out[("origin", "full")] = m[:3], m[3:]
    m = re.search(r"GLMM unadjusted by kind: nested OR " + num + r".*?origin OR " + num, t).groups()
    out[("nested", "unadj")], out[("origin", "unadj")] = m[:3], m[3:]
    out[("all", "branch")] = re.search(r"; adjusted event OR " + num, t.split("Branch-level")[1]).groups()
    out[("nested", "branch")] = re.search(r"nested OR " + num, t.split("branch-level adjusted by kind")[1]).groups()
    return {k: tuple(float(x) for x in v) for k, v in out.items()}


def panel_c(ax):
    g = parse_glmm(f"{VER}/glmm_results.txt")
    ge = parse_glmm(f"{VER}/glmm_results_exh.txt")
    rows = [("all", "all events\n(19)"), ("nested", "simple-to-\ncomplex (7)"),
            ("origin", "origins of\nsimple (12)")]
    models = [("unadj", g, "unadjusted", C_DEF, "o"),
              ("full", g, "adjusted (length, disorder, family size, depth, search space)", C_ADJ, "s"),
              ("full", ge, "adjusted, stricter (exhaustive) search", C_EXH, "D"),
              ("branch", g, "adjusted, branches as units", C_BR, "^")]
    yc = np.arange(len(rows))[::-1] * 1.0
    offs = np.linspace(0.27, -0.27, len(models))
    for (key, lab), y0 in zip(rows, yc):
        for (m, src, mlab, col, mk), o in zip(models, offs):
            if (key, m) not in src:
                continue
            est, lo, hi = src[(key, m)]
            ax.plot([lo, hi], [y0 + o] * 2, color=col, lw=1.2)
            ax.plot(est, y0 + o, mk, color=col, ms=4.5, mec="white", mew=0.5,
                    label=mlab if key == "all" else None)
    ax.axvline(1, color=INK2, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xticks([0.8, 1, 1.5, 2, 3]); ax.set_xticklabels(["0.8", "1", "1.5", "2", "3"])
    ax.set_xlim(0.75, 3.6)
    ax.set_yticks(yc); ax.set_yticklabels([r[1] for r in rows], fontsize=6.3)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(rows) - 0.4 + 1.25)
    ax.set_xlabel("odds ratio: 'no structural match outside the clade',\nevent branches vs control branches (log scale)")
    ax.set_title("Are sequence-novel families more often\nstructurally unmatched at events?", loc="left", color=INK)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=1, handlelength=1.0, fontsize=5.6)
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)


def panel_d(ax):
    d = pd.read_csv(f"{W6}/panel_tests.tsv", sep="\t")
    d = d[d.analysis == "primary"]
    panels = [("R1", "animals and relatives"), ("R2", "fungi"), ("R4_curated", "volvocine algae")]
    pipes = [("real", "predicted structure (real data)", C_REAL, "o"),
             ("shuffled3di", "control: structure letters shuffled", C_SHUF, "D"),
             ("aa_splithalf", "control: two halves of the sequence", C_SPLIT, "s")]
    x = np.arange(len(panels)); offs = [-0.22, 0, 0.22]
    for (pp, lab, col, mk), o in zip(pipes, offs):
        for i, (pn, _) in enumerate(panels):
            r = d[(d.panel == pn) & (d.pipeline == pp)].iloc[0]
            ax.plot([x[i] + o] * 2, [r.b_effect_lo, r.b_effect_hi], color=col, lw=1.2)
            ax.plot(x[i] + o, r.b_effect, mk, color=col, ms=5, mec=INK2 if pp == "aa_splithalf" else "white", mew=0.5,
                    label=lab if i == 0 else None)
    ax.axhline(0, color=INK2, lw=0.8, ls="--")
    ax.set_xticks(x); ax.set_xticklabels([p[1] for p in panels], fontsize=6.3)
    ax.set_xlim(-0.5, len(panels) - 0.5)
    ax.set_ylim(-0.4, 0.5)
    ax.set_ylabel("crown minus unicellular branches\n(log structure-per-sequence change)")
    ax.set_title("Does structure change faster per sequence\nchange inside multicellular clades?", loc="left", color=INK)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), handlelength=1.0, fontsize=5.8)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)


def main():
    fig = plt.figure(figsize=(180 * MM, 150 * MM))
    gs = fig.add_gridspec(2, 2, wspace=0.38, hspace=0.55, left=0.11, right=0.985, top=0.95, bottom=0.1,
                          width_ratios=[1, 1])
    axa = fig.add_subplot(gs[0, 0]); panel_a(axa); panel_label(axa, "a", x=-0.2)
    axb = fig.add_subplot(gs[0, 1]); panel_b(axb); panel_label(axb, "b", x=-0.2)
    gs2 = gs[1, 0].subgridspec(1, 1)
    axc = fig.add_subplot(gs2[0, 0]); panel_c(axc); panel_label(axc, "c", x=-0.2)
    axd = fig.add_subplot(gs[1, 1]); panel_d(axd); panel_label(axd, "d", x=-0.2)
    save(fig, "fig7_structure")


if __name__ == "__main__":
    main()
