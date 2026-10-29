"""Figure 3. Sequence-novel gene families at simple-to-complex events: the confirmatory test, the tip artefact,
and the internal-only estimate.

Sources
 (a) results/r2/verification/panels/verify_stats.json + verify_event_residuals.tsv  -- SPEC-LITERAL values (no
     branch-length floor; deviation D-PI-1 in results/r2/panels/DEVIATIONS.txt), produced by
     scripts/r2/verification/panels/verify_stats.py.
 (b) results/r2/panels/pooled/final_8panels/{event_residuals.tsv, branch_residuals.tsv, event_class_comparators.tsv}
     (implemented model; the floor changes one amoebozoan control only).
 (c) event_residuals.tsv tip-control columns (tip-control model, R4/R6 only; unaffected by the floor),
     event_class_comparators.tsv (tip controls), spec-literal internal components for the two internal P1 events.
 (d) event_class_comparators.tsv (internal branches, implemented model) and event_residuals.tsv per event.
Spec-literal versions of the class pools in (b)-(d) are recomputed here from verify_event_residuals.tsv and printed.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

from common import (P, MM, UNI, SIMPLE, COMPLEX, REVERSAL, CONFIRM, INK, INK2, GRID, P1_EVENTS, save, fig_label, reml)

POOL = f"{P}/results/r2/panels/pooled/final_8panels"
VER = f"{P}/results/r2/verification/panels"
NAMES = {"R4_astrephomene": "Astrephomene (tip)", "R4_volvox_carteri": "Volvox carteri (tip)",
         "R4_volvox_africanus": "V. africanus/reticuliferus", "R6_fucophycidae": "Fucophycidae brown algae",
         "R2_pezizomycotina": "Pezizomycotina fungi", "R2_agaricomycetes": "Agaricomycetes fungi",
         "R3_embryophyta": "land plants", "R5_florideophyceae": "Florideophyceae red algae",
         "R8_heterocyst": "heterocystous Nostocales"}
XLIM = (0.02, 800)
XT = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300]


def logx(ax):
    ax.set_xscale("log")
    ax.set_xlim(*XLIM)
    ax.set_xticks(XT)
    ax.set_xticklabels([f"{t:g}" for t in XT])
    ax.minorticks_off()
    ax.axvline(1, color=INK2, lw=0.7, ls=(0, (2, 2)), zorder=0)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def ci_line(ax, y, rr, lo, hi, col, filled=True, ms=4.5, lw=1.1, marker="o"):
    ax.plot([max(lo, XLIM[0]), min(hi, XLIM[1])], [y, y], color=col, lw=lw, solid_capstyle="butt", zorder=2)
    if hi > XLIM[1]:
        ax.annotate("", xy=(XLIM[1], y), xytext=(XLIM[1] / 1.6, y), arrowprops=dict(arrowstyle="->", color=col, lw=lw))
    if lo < XLIM[0]:
        ax.annotate("", xy=(XLIM[0], y), xytext=(XLIM[0] * 1.6, y), arrowprops=dict(arrowstyle="->", color=col, lw=lw))
    ax.plot(rr, y, marker, color=col, mfc=col if filled else "white", mec=col, mew=1.0, ms=ms, zorder=3)


def diamond(ax, y, rr, lo, hi, col, h=0.32):
    ax.fill([lo, rr, hi, rr], [y, y + h, y, y - h], color=col, zorder=3, lw=0)


def main():
    vs = json.load(open(f"{VER}/verify_stats.json"))
    vr = pd.read_csv(f"{VER}/verify_event_residuals.tsv", sep="\t").set_index("event_id")
    er = pd.read_csv(f"{POOL}/event_residuals.tsv", sep="\t").set_index("event_id")
    br = pd.read_csv(f"{POOL}/branch_residuals.tsv", sep="\t")
    cc = pd.read_csv(f"{POOL}/event_class_comparators.tsv", sep="\t")

    def comp(cls, branch, model_prefix):
        r = cc[(cc.event_class == cls) & (cc.branch == branch) & (cc.model.str.startswith(model_prefix))].iloc[0]
        return r.RR, r.RR_low, r.RR_high, int(r.k)

    fig = plt.figure(figsize=(180 * MM, 215 * MM))
    gs = fig.add_gridspec(4, 1, height_ratios=[5.2, 4.6, 5.0, 10.2], hspace=0.62,
                          left=0.25, right=0.74, top=0.965, bottom=0.088)

    # ---------- (a) confirmatory events, spec-literal ----------
    ax = fig.add_subplot(gs[0])
    ys = [5, 4, 3, 2]
    for y, e in zip(ys, P1_EVENTS):
        d = vs["p1"][e]
        tip = bool(vr.loc[e, "tip"])
        ci_line(ax, y, d["rr"], d["lo"], d["hi"], CONFIRM, filled=not tip)
        ax.text(1.02, y, f"{d['rr']:.2f}   p = {d['p_perm']:.3f}", transform=ax.get_yaxis_transform(), va="center",
                fontsize=6, color=INK)
    m = vs["p1_meta"]
    diamond(ax, 0.8, m["RR"], m["lo"], m["hi"], CONFIRM)
    ax.text(1.02, 0.8, f"{m['RR']:.2f}   p = {vs['stouffer_p']:.4f}", transform=ax.get_yaxis_transform(), va="center",
            fontsize=6, color=INK, fontweight="bold")
    ax.set_yticks(ys + [0.8])
    ax.set_yticklabels([NAMES[e] for e in P1_EVENTS] + ["pooled (4 events)"])
    ax.set_ylim(0, 5.8)
    logx(ax)
    ax.set_xlabel("rate ratio, observed / expected sequence-novel families")
    ax.set_title("Confirmatory test: four simple → complex events against internal control branches", loc="left", color=INK)
    ax.text(1.02, 5.65, "RR, one-sided permutation p", transform=ax.get_yaxis_transform(), fontsize=5.8, color=INK2)

    # ---------- (b) tip artefact ----------
    ax = fig.add_subplot(gs[1])
    tipc = br[br.tip_control.astype(str) == "True"]
    rr_tc = np.exp(tipc.s.values)
    groups = [("unicellular tips (no event)", rr_tc, UNI, None),
              ("origin tips", "origin", SIMPLE, comp("origin", "tip", "primary")),
              ("simple → complex tips", "simple->complex", COMPLEX, comp("simple->complex", "tip", "primary")),
              ("reversal tips", "reversal", REVERSAL, comp("reversal", "tip", "primary"))]
    rng = np.random.default_rng(11)
    for gi, (lab, src, col, pooled) in enumerate(groups):
        y = len(groups) - gi
        if isinstance(src, str):
            ev = cc[(cc.event_class == src) & (cc.branch == "tip") & (cc.model.str.startswith("primary"))].events.iloc[0].split(";")
            vals = np.exp(er.loc[ev, "log_rr"].values)
            names = ev
        else:
            vals = src; names = [None] * len(src)
        jit = rng.uniform(-0.17, 0.17, len(vals))
        for v, j, nm in zip(vals, jit, names):
            if nm in P1_EVENTS:
                ax.plot(v, y + 0.28 + j * 0.3, "o", mfc="white", mec=CONFIRM, mew=1.2, ms=4.5, zorder=4)
            else:
                ax.plot(v, y + 0.28 + j * 0.3, "o", color=col if col != UNI else "#9a9a9a", ms=2.6, alpha=0.9, zorder=3)
        if pooled is None:
            gm = float(np.exp(np.mean(np.log(vals))))
            ax.plot([gm, gm], [y - 0.32, y + 0.05], color=INK, lw=1.4)
            ax.text(1.02, y, f"geometric mean {gm:.1f}, n = {len(vals)}", transform=ax.get_yaxis_transform(), fontsize=5.8, color=INK, va="center")
        else:
            rr, lo, hi, k = pooled
            diamond(ax, y - 0.18, rr, lo, hi, col, h=0.16)
            ax.text(1.02, y, f"pooled {rr:.1f} [{lo:.1f}, {hi:.1f}], n = {k}", transform=ax.get_yaxis_transform(), fontsize=5.8, color=INK, va="center")
    ax.set_yticks([len(groups) - i for i in range(len(groups))])
    ax.set_yticklabels([g[0] for g in groups])
    ax.set_ylim(0.4, len(groups) + 0.75)
    logx(ax)
    ax.set_xlabel("rate ratio against internal control branches")
    ax.set_title("Every single-species (tip) branch scores high against internal controls, whatever its biology", loc="left", color=INK)
    gm_tc = float(np.exp(np.mean(np.log(rr_tc))))

    # ---------- (c) like with like ----------
    ax = fig.add_subplot(gs[2])
    rows = []
    for e in ["R4_astrephomene", "R4_volvox_carteri"]:
        lr, v = er.loc[e, "tip_log_rr"], er.loc[e, "tip_var"]
        rows.append((NAMES[e], np.exp(lr), np.exp(lr - 1.96 * np.sqrt(v)), np.exp(lr + 1.96 * np.sqrt(v)), CONFIRM, False,
                     f"p = {er.loc[e, 'tip_p_perm']:.2f}"))
    # confirmatory test recomputed with tip controls for the two tips; internal events keep spec-literal components
    y_ = [er.loc[e, "tip_log_rr"] for e in P1_EVENTS[:2]] + [vr.loc[e, "log_rr"] for e in P1_EVENTS[2:]]
    v_ = [er.loc[e, "tip_var"] for e in P1_EVENTS[:2]] + [vr.loc[e, "var"] for e in P1_EVENTS[2:]]
    mt = reml(y_, v_)
    ps = [er.loc[e, "tip_p_perm"] for e in P1_EVENTS[:2]] + [vr.loc[e, "p_perm"] for e in P1_EVENTS[2:]]
    pst = 1 - stats.norm.cdf(stats.norm.ppf(1 - np.array(ps)).sum() / 2)
    rows.append(("confirmatory test, tip controls", mt["RR"], mt["lo"], mt["hi"], CONFIRM, True, f"p = {pst:.2f}"))
    a = comp("simple->complex", "tip", "tip controls")
    rows.append((f"simple → complex tips (n = {a[3]})", a[0], a[1], a[2], COMPLEX, True, ""))
    a = comp("reversal", "tip", "tip controls")
    rows.append((f"reversal tips (n = {a[3]})", a[0], a[1], a[2], REVERSAL, True, ""))
    yy = list(range(len(rows), 0, -1))
    for y, (lab, rr, lo, hi, col, pooled, ptxt) in zip(yy, rows):
        if pooled:
            diamond(ax, y, rr, lo, hi, col, h=0.25)
        else:
            ci_line(ax, y, rr, lo, hi, col, filled=False)
        ax.text(1.02, y, f"{rr:.2f} [{lo:.2f}, {hi:.2f}] {ptxt}", transform=ax.get_yaxis_transform(), va="center",
                fontsize=5.8, color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(0.4, len(rows) + 0.6)
    logx(ax)
    ax.set_xlim(0.1, 10)
    ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"])
    ax.set_xlabel("rate ratio against unicellular tip controls")
    ax.set_title("Compared tip with tip, the excess disappears", loc="left", color=INK)

    # ---------- (d) internal branches only ----------
    ax = fig.add_subplot(gs[3])
    ev_sc = cc[(cc.event_class == "simple->complex") & (cc.branch == "internal")].events.iloc[0].split(";")
    ev_sc = sorted(ev_sc, key=lambda e: -er.loc[e, "rr"])
    n = len(ev_sc)
    yy = list(range(n + 3, 3, -1))
    for y, e in zip(yy, ev_sc):
        col = CONFIRM if e in P1_EVENTS else COMPLEX
        ci_line(ax, y, er.loc[e, "rr"], er.loc[e, "rr_low"], er.loc[e, "rr_high"], col, ms=3.8, lw=0.9)
        ax.text(1.02, y, f"{er.loc[e, 'rr']:.2f} [{er.loc[e, 'rr_low']:.2f}, {er.loc[e, 'rr_high']:.1f}]",
                transform=ax.get_yaxis_transform(), va="center", fontsize=5.8, color=INK2)
    pools = [("simple → complex, pooled", comp("simple->complex", "internal", "primary"), COMPLEX, 2.6),
             ("origins, pooled", comp("origin", "internal", "primary"), SIMPLE, 1.4),
             ("reversals, pooled", comp("reversal", "internal", "primary"), REVERSAL, 0.4)]
    for lab, (rr, lo, hi, k), col, y in pools:
        if hi > 100 or lo < 0.1:
            ax.plot([0.1, 100], [y, y], color=col, lw=1.1)
            ax.plot(rr, y, "D", color=col, ms=5, zorder=3)
            ax.text(0.105, y + 0.25, "95% CI beyond axis", fontsize=5.3, color=INK2, va="bottom")
        else:
            diamond(ax, y, rr, lo, hi, col, h=0.3)
        ax.text(1.02, y, f"{rr:.2f} [{lo:.2f}, {hi:.3g}], n = {k}", transform=ax.get_yaxis_transform(), va="center",
                fontsize=5.8, color=INK)
    ax.set_yticks(yy + [p[3] for p in pools])
    ax.set_yticklabels([NAMES[e] for e in ev_sc] + [p[0] for p in pools])
    ax.set_ylim(-0.3, n + 3.7)
    logx(ax)
    ax.set_xlim(0.1, 100)
    ax.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.axvline(2.2, color=INK2, lw=0.7, ls=(0, (4, 2)), zorder=0)
    ax.text(2.2, n + 3.55, "2.2-fold bias bar", fontsize=5.5, color=INK2, ha="center", va="bottom")
    ax.set_xlabel("rate ratio against internal control branches")
    ax.set_title("Internal branches only (exploratory): simple → complex events vs origins", loc="left", color=INK)

    leg = [Line2D([], [], marker="o", ls="", color=INK2, mfc=INK2, ms=4, label="internal branch"),
           Line2D([], [], marker="o", ls="", color=INK2, mfc="white", ms=4, label="tip (single-species) branch"),
           Line2D([], [], marker="D", ls="", color=INK2, ms=4, label="pooled estimate (95% CI)"),
           Line2D([], [], marker="s", ls="", color=CONFIRM, ms=5, label="confirmatory event"),
           Line2D([], [], marker="s", ls="", color=COMPLEX, ms=5, label="other simple → complex"),
           Line2D([], [], marker="s", ls="", color=SIMPLE, ms=5, label="origin"),
           Line2D([], [], marker="s", ls="", color=REVERSAL, ms=5, label="reversal")]
    fig.legend(handles=leg, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.005), fontsize=6.2, handletextpad=0.3,
               columnspacing=1.0)
    for lab, a_ in zip("abcd", fig.axes[:4]):
        fig_label(fig, lab, 0.01, a_.get_position().y1 + 0.03)
    save(fig, "fig3_novelty_events")

    # ---- numbers for caption / cross-check ----
    print("P1 spec-literal:", {e: (round(vs['p1'][e]['rr'], 2), round(vs['p1'][e]['p_perm'], 3)) for e in P1_EVENTS},
          "pooled", round(m["RR"], 2), round(m["lo"], 2), round(m["hi"], 1), "Stouffer p", round(vs["stouffer_p"], 4))
    print("non-event unicellular tip controls n", len(rr_tc), "geomean", round(gm_tc, 2))
    print("P1 with tip controls (spec-literal internal parts): RR", round(mt["RR"], 2), round(mt["lo"], 2), round(mt["hi"], 2),
          "Stouffer p", round(pst, 3), "ps", ps)
    for cls, br_ in [("simple->complex", "internal"), ("origin", "internal"), ("reversal", "internal"),
                     ("reversal", "tip"), ("origin", "tip"), ("simple->complex", "tip")]:
        evs = cc[(cc.event_class == cls) & (cc.branch == br_) & cc.model.str.startswith("primary")].events.iloc[0].split(";")
        r = reml(vr.loc[evs, "log_rr"], vr.loc[evs, "var"])
        print(f"spec-literal {cls:16s} {br_:8s} k={r['k']:2d} RR {r['RR']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}]   implemented",
              [round(x, 2) for x in comp(cls, br_, 'primary')[:3]])


if __name__ == "__main__":
    main()
