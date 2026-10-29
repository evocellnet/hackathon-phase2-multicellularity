"""Shared style, palette and helpers for the report v2 figures (report/figures_v2/).

Run any figure script with the project's 'ancestra' env, e.g.
    conda_envs/ancestra/bin/python scripts/r2/figures/fig1_overview.py
or all of them with scripts/r2/figures/make_all.sh.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import optimize, stats  # noqa: E402

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
FIGDIR = os.path.join(P, "report", "figures_v2")
os.makedirs(FIGDIR, exist_ok=True)

MM = 1 / 25.4  # inches per mm

# Okabe-Ito based, validated (dataviz validate_palette.js, all pairs, light mode).
# Colour encodes the multicellular state; an event branch takes the colour of the state it leads to.
UNI = "#BDBDBD"       # unicellular
SIMPLE = "#009E73"    # simple multicellular   (origins)
COMPLEX = "#0072B2"   # complex multicellular  (simple-to-complex events)
REVERSAL = "#E69F00"  # reversals (back towards unicellularity / simpler state)
CONFIRM = "#CC79A7"   # the four pre-registered (confirmatory) simple-to-complex events
STATE = {"unicellular": UNI, "simple": SIMPLE, "complex": COMPLEX}
INK = "#1a1a1a"
INK2 = "#555555"
GRID = "#e5e5e5"
# method identity (Fig 2)
ANC = "#3a3a3a"
EOMA = "#56B4E9"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

PANEL_NAMES = {"R1": "Holozoa", "R2": "Fungi", "R3": "Streptophyta", "R4": "Volvocales",
               "R5": "Red algae", "R6": "Stramenopiles", "R7": "Amoebozoa", "R8": "Cyanobacteria"}
PANELS = list(PANEL_NAMES)
P1_EVENTS = ["R4_astrephomene", "R4_volvox_carteri", "R4_volvox_africanus", "R6_fucophycidae"]


def panel_label(ax, s, x=-0.02, y=1.02, **kw):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", color=INK, **kw)


def fig_label(fig, s, x, y):
    fig.text(x, y, s, fontsize=9, fontweight="bold", va="top", ha="left", color=INK)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=300)
    plt.close(fig)
    print("wrote", os.path.join(FIGDIR, name) + ".{png,pdf}")


# ---------------- Newick ----------------
class Node:
    __slots__ = ("name", "children", "parent", "bl", "x", "y")

    def __init__(self):
        self.name, self.children, self.parent, self.bl, self.x, self.y = "", [], None, 0.0, 0.0, 0.0

    def tips(self):
        if not self.children:
            return [self]
        out = []
        for c in self.children:
            out += c.tips()
        return out

    def walk(self):
        yield self
        for c in self.children:
            yield from c.walk()


def read_newick(path):
    s = open(path).read().strip().rstrip(";")
    root = Node(); cur = root; i = 0; tok = ""
    stack = []

    def flush(node, t):
        if ":" in t:
            n, b = t.split(":", 1)
            node.name = n.strip(); node.bl = float(b)
        else:
            node.name = t.strip()

    while i < len(s):
        ch = s[i]
        if ch == "(":
            n = Node(); n.parent = cur; cur.children.append(n); stack.append(cur); cur = n
        elif ch == ",":
            flush(cur, tok); tok = ""
            par = stack[-1]; n = Node(); n.parent = par; par.children.append(n); cur = n
        elif ch == ")":
            flush(cur, tok); tok = ""
            cur = stack.pop()
        else:
            tok += ch
        i += 1
    flush(cur, tok)
    return root


def layout_cladogram(root):
    """Tips aligned on the right at x = max depth; internal x = min(child x) - 1; y = tip order."""
    tips = root.tips()
    for k, t in enumerate(tips):
        t.y = len(tips) - 1 - k

    def depth(n):
        return 0 if n.parent is None else depth(n.parent) + 1
    maxd = max(depth(t) for t in tips)

    def setx(n):
        if not n.children:
            n.x = maxd
        else:
            for c in n.children:
                setx(c)
            n.x = min(c.x for c in n.children) - 1
            n.y = float(np.mean([c.y for c in n.children]))
    setx(root)
    return tips


# ---------------- random-effects meta-analysis (copied from scripts/r2/verification/panels/verify_stats.py) ----------------
def reml(y, v):
    y = np.asarray(y, float); v = np.asarray(v, float); k = len(y)

    def nll(t2):
        w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
        return 0.5 * (np.log(v + t2).sum() + np.log(w.sum()) + (w * (y - mu) ** 2).sum())
    r = optimize.minimize_scalar(nll, bounds=(0, 100), method="bounded", options={"xatol": 1e-10})
    t2 = r.x if nll(r.x) < nll(0) else 0.0
    w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
    q = (w * (y - mu) ** 2).sum() / (k - 1); se = np.sqrt(max(1, q) / w.sum()); t = stats.t.ppf(0.975, k - 1)
    return dict(k=k, tau2=t2, mu=mu, RR=np.exp(mu), lo=np.exp(mu - t * se), hi=np.exp(mu + t * se),
                p_one=1 - stats.t.cdf(mu / se, k - 1))
