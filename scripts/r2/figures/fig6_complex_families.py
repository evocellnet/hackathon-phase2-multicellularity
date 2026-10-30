"""Figure 6. What is gained and duplicated at simple-to-complex events: example domain families.
Sources: results/r2/panels/<R>/w3_event_families.tsv (primary panels; R4, not the curated sensitivity panel),
results/r2/panels/pooled/final_8panels/w3/{w3_event_summary.tsv, w3_recurrent_domains.tsv}.
A family is counted for a domain if at least one of its representative proteins carries that Pfam domain
(transposon-domain families excluded). 'Gained' = sequence-novel family originating on the event branch;
'duplicated' = family with at least one duplication on the event branch.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

from common import P, MM, COMPLEX, CONFIRM, INK, INK2, GRID, P1_EVENTS, save, fig_label

R = f"{P}/results/r2/panels"
EVENTS = [  # (event_id, panel, label, tip)
    ("R4_volvox_carteri", "R4", "Volvox carteri", True),
    ("R4_volvox_africanus", "R4", "V. africanus/reticuliferus", False),
    ("R4_astrephomene", "R4", "Astrephomene", True),
    ("R4_pleodorina", "R4", "Pleodorina", True),
    ("R8_heterocyst", "R8", "heterocystous Nostocales", False),
    ("R3_embryophyta", "R3", "land plants", False),
    ("R5_florideophyceae", "R5", "Florideophyceae", False),
    ("R6_fucophycidae", "R6", "Fucophycidae", False),
    ("R2_pezizomycotina", "R2", "Pezizomycotina", False),
    ("R2_agaricomycetes", "R2", "Agaricomycetes", False),
]
DOMAINS = [  # (group, pfam, label)
    ("volvocine\nextracellular matrix", "PF12499", "pherophorin-like (DUF3707)"),
    ("volvocine\nextracellular matrix", "PF05548", "matrix metalloprotease (Peptidase_M11)"),
    ("volvocine\nextracellular matrix", "PF00211", "guanylate cyclase"),
    ("two-component\nsignalling", "PF00072", "response regulator"),
    ("two-component\nsignalling", "PF00512", "histidine kinase (HisKA)"),
    ("two-component\nsignalling", "PF02518", "histidine kinase ATPase (HATPase_c)"),
    ("two-component\nsignalling", "PF01590", "GAF domain"),
    ("two-component\nsignalling", "PF00805", "pentapeptide repeat"),
    ("protein kinases", "PF00069", "protein kinase (Pkinase)"),
    ("protein kinases", "PF07714", "Tyr/Ser-Thr kinase (PK_Tyr_Ser-Thr)"),
    ("plant receptors", "PF13855", "leucine-rich repeat (LRR_8)"),
    ("plant receptors", "PF00646", "F-box"),
    ("transcription\nfactors", "PF00096", "C2H2 zinc finger"),
    ("transcription\nfactors", "PF00170", "bZIP"),
    ("transcription\nfactors", "PF00172", "Zn2Cys6 (Zn_clus)"),
    ("fungal\nmetabolism", "PF07690", "MFS transporter"),
    ("fungal\nmetabolism", "PF00106", "short-chain dehydrogenase"),
    ("fungal\nmetabolism", "PF00734", "cellulose-binding (CBM_1)"),
]
KINASES = {"PF00069", "PF07714"}


def load():
    fr = []
    for pan in sorted({e[1] for e in EVENTS}):
        d = pd.read_csv(f"{R}/{pan}/w3_event_families.tsv", sep="\t", dtype=str)
        fr.append(d)
    d = pd.concat(fr, ignore_index=True)
    d = d[d.event_id.isin([e[0] for e in EVENTS]) & (d.transposon_domain != "True")]
    d["pf"] = d.pfam.fillna("").apply(lambda s: set(x for x in s.split(";") if x))
    return d


def main():
    d = load()
    ev_ids = [e[0] for e in EVENTS]
    M = np.zeros((len(DOMAINS), len(EVENTS)), int)
    for j, eid in enumerate(ev_ids):
        sub = d[d.event_id == eid]
        for i, (_, pf, _) in enumerate(DOMAINS):
            M[i, j] = int(sub.pf.apply(lambda s: pf in s).sum())
    # kinase in sequence-novel families at the 7 internal events (the recurrence test unit)
    internal = [e for e in EVENTS if not e[3]]
    kin = []
    for eid, _, lab, _ in internal:
        sub = d[(d.event_id == eid) & (d.kind == "seq_novel")]
        kin.append(int(sub.pf.apply(lambda s: bool(s & {"PF00069"})).sum()))
    summ = pd.read_csv(f"{R}/pooled/final_8panels/w3/w3_event_summary.tsv", sep="\t").set_index("event_id")
    rec = pd.read_csv(f"{R}/pooled/final_8panels/w3/w3_recurrent_domains.tsv", sep="\t")
    kr = rec[(rec.units == "seq_novel") & (rec.pfam == "PF00069")].iloc[0]

    fig = plt.figure(figsize=(180 * MM, 175 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.3, 1.0], width_ratios=[1.0, 1.0], hspace=0.42, wspace=0.75,
                          left=0.2, right=0.9, top=0.86, bottom=0.07)

    # ---------- (a) heatmap ----------
    ax = fig.add_subplot(gs[0, :])
    pos = ax.get_position(); ax.set_position([0.44, pos.y0 + 0.02, 0.46, pos.height - 0.05])
    cmap = LinearSegmentedColormap.from_list("seq", ["#f4f8fb", "#9cc3e0", COMPLEX, "#00365a"])
    vmax = max(M.max(), 1)
    im = ax.imshow(M, aspect="auto", cmap=cmap, norm=PowerNorm(0.5, vmin=0, vmax=vmax), interpolation="nearest")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=5.8,
                        color="white" if v > vmax * 0.25 else INK)
            else:
                ax.text(j, i, "·", ha="center", va="center", fontsize=6, color="#b5b5b5")
    ax.set_xticks(range(len(EVENTS)))
    ax.set_xticklabels([e[2] + (" (tip)" if e[3] else "") for e in EVENTS], fontsize=6, rotation=40, ha="left",
                       rotation_mode="anchor")
    for t, e in zip(ax.get_xticklabels(), EVENTS):
        if e[0] in P1_EVENTS:
            t.set_color(CONFIRM); t.set_fontweight("bold")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=2)
    for t, e in zip(ax.get_xticklabels(), EVENTS):
        if e[0] in P1_EVENTS:
            t.set_color(CONFIRM); t.set_fontweight("bold")
    ax.set_yticks(range(len(DOMAINS)))
    ax.set_yticklabels([x[2] for x in DOMAINS], fontsize=6.2)
    ax.tick_params(axis="y", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    # white separators
    ax.set_xticks(np.arange(-0.5, len(EVENTS)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(DOMAINS)), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    # group brackets on the left
    groups = []
    for i, (g, _, _) in enumerate(DOMAINS):
        if not groups or groups[-1][0] != g:
            groups.append([g, i, i])
        else:
            groups[-1][2] = i
    for g, a, b in groups:
        ax.annotate("", xy=(-0.64, a - 0.35), xytext=(-0.64, b + 0.35), xycoords=("axes fraction", "data"),
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
        ax.text(-0.655, (a + b) / 2, g, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.2,
                color=INK)
        if b < len(DOMAINS) - 1:
            ax.axhline(b + 0.5, color="white", lw=3)
    # event-type brackets above (volvocine germ-soma / others)
    for x0, x1, t in [(-0.4, 3.4, "volvocine germ–soma"), (3.6, 9.4, "other simple → complex (internal branches)")]:
        yb = len(DOMAINS) - 0.5 + 0.35
        ax.plot([x0, x1], [yb, yb], color=INK2, lw=0.8, clip_on=False)
        ax.text((x0 + x1) / 2, yb + 0.25, t, ha="center", va="top", fontsize=6.2, color=INK)
    p_ = ax.get_position()
    cax = fig.add_axes([p_.x1 + 0.015, p_.y0 + p_.height * 0.25, 0.012, p_.height * 0.75])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("families gained or duplicated on the event branch\nthat carry the domain (count)", fontsize=6)
    cb.ax.tick_params(labelsize=5.8)
    cb.outline.set_visible(False)

    # ---------- (b) kinases: gained vs duplicated families ----------
    ax = fig.add_subplot(gs[1, 0])
    dup = []
    for eid, _, _, _ in internal:
        sub = d[(d.event_id == eid) & (d.kind == "dup")]
        dup.append(int(sub.pf.apply(lambda s_: "PF00069" in s_).sum()))
    y = np.arange(len(internal))
    h = 0.38
    cols = [CONFIRM if e[0] in P1_EVENTS else COMPLEX for e in internal]
    ax.barh(y - h / 2, kin, height=h, color=cols, edgecolor="white", linewidth=0.5)
    ax.barh(y + h / 2, dup, height=h, color=cols, alpha=0.4, edgecolor="white", linewidth=0.5)
    for yi, a, b in zip(y, kin, dup):
        ax.text(a + 0.6, yi - h / 2, str(a), va="center", fontsize=5.6, color=INK)
        ax.text(b + 0.6, yi + h / 2, str(b), va="center", fontsize=5.6, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([e[2] for e in internal], fontsize=6)
    for t, e in zip(ax.get_yticklabels(), internal):
        if e[0] in P1_EVENTS:
            t.set_color(CONFIRM)
    ax.invert_yaxis()
    ax.set_xlabel("families with a protein-kinase domain (count)")
    ax.set_title(f"Kinases: one new family at each of {int(kr.k_events)} of 7 events\n(recurrence p = {kr.p_perm:.3f}, q = {kr.q_bh:.2f})",
                 loc="left", color=INK)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=INK2, label="new (sequence-novel)"), Patch(color=INK2, alpha=0.4, label="duplicated")],
              loc="lower right", fontsize=5.8, handlelength=1.0)
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xlim(0, 45)

    # ---------- (c) most novel families carry no domain ----------
    ax = fig.add_subplot(gs[1, 1])
    evs = [e for e in EVENTS]
    tot = np.array([summ.loc[e[0], "seq_novel_families"] for e in evs])
    wp = np.array([summ.loc[e[0], "seq_novel_with_pfam"] for e in evs])
    frac = 100 * wp / tot
    y = np.arange(len(evs))
    ax.barh(y, frac, color=[CONFIRM if e[0] in P1_EVENTS else COMPLEX for e in evs], height=0.7,
            edgecolor="white", linewidth=0.6)
    for yi, f, a, b in zip(y, frac, wp, tot):
        ax.text(f + 1, yi, f"{a:,} of {b:,}", va="center", fontsize=5.6, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([e[2] + (" (tip)" if e[3] else "") for e in evs], fontsize=6)
    for t, e in zip(ax.get_yticklabels(), evs):
        if e[0] in P1_EVENTS:
            t.set_color(CONFIRM)
    ax.invert_yaxis()
    ax.set_xlim(0, 75)
    ax.set_xlabel("new families with any Pfam domain (%)")
    ax.set_title("Most new families have no known domain", loc="left", color=INK)
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)

    fig_label(fig, "a", 0.01, 0.985)
    fig_label(fig, "b", 0.01, 0.335)
    fig_label(fig, "c", 0.54, 0.335)
    save(fig, "fig6_complex_families")
    print(pd.DataFrame(M, index=[x[2] for x in DOMAINS], columns=ev_ids).to_string())
    print("kinase seq_novel per internal event", dict(zip([e[0] for e in internal], kin)), "dup", dup, "recurrence row", kr.to_dict())
    print("frac with pfam", dict(zip(ev_ids, np.round(frac, 1))))


if __name__ == "__main__":
    main()
