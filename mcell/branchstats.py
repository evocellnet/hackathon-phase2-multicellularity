"""Branch-level statistics for the transition tests (DESIGN.md §4.6, binding).

Everything that turns per-branch gene-content numbers (from :mod:`mcell.hogs`)
into the primary test results lives here:

1. **Tree and branch table** -- :func:`build_branch_table`: one row per branch
   ``parent -> child`` (indexed by the child name) with depth, descendant and
   sister tip counts, terminal / root-child flags, branch length, transcriptome
   fractions and mean BUSCO of descendants and sisters, the node states, the
   branch class (stem, nested stem, within-multicellular, within-unicellular,
   grey-zone, other), the transition it carries and membership of the control
   pool U with the reason for every exclusion.
2. **Matching** -- :func:`match_controls`: z-scored covariates, the
   k = min(20, |U|) Mahalanobis-nearest controls of every stem, SMD table.
3. **Relabelling test** -- :func:`relabel_test` (u = mid-rank percentile,
   T = mean over transitions, 10,000 draws, seed 11). Default scheme
   'permute' (DESIGN v3.1 D6); the v3 draw-without-replacement scheme is
   ``scheme='spec'`` (sensitivity; anti-conservative when the N_i overlap).
4. **Effect sizes** -- :func:`effect_table` (log O/E per transition with the
   5-95 % control range), :func:`random_effects` (REML, DL fallback),
   :func:`holm`, :func:`leave_one_out`.
5. **Statistics s_b** -- :func:`count_statistic` and
   :func:`count_residual_statistic` (H1/H2, NB-GLM residual of D6;
   :func:`residual_bootstrap_test`), path stems of D7 (:func:`path_sum`),
   :func:`compositional_statistic` (H4/H5), :func:`module_loss_statistic` and
   :func:`h12_signs` (H12), :func:`h9_values` (H9 plug-in).
6. **Power** -- :func:`simulate_power`, :func:`mde`, :func:`min_attainable_p`.
7. **Registry** -- :func:`append_registry`.

Branch identifiers are the names of the child nodes (all nodes of the tree must
be named and unique). A transition is identified by ``transition_id``; its stem
is the branch ending at M. See docs/branchstats.md for the mapping to §4.6.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import math
import os
import subprocess
import warnings
from dataclasses import dataclass, field
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import optimize, stats

__all__ = [
    "SEED", "N_DRAWS", "COVARIATES", "BRANCH_CLASSES",
    "load_tree", "branch_table", "node_states", "build_branch_table", "stem_branches",
    "MatchResult", "match_controls", "regression_adjusted_expectation",
    "u_score", "RelabelResult", "relabel_test", "min_attainable_p", "count_relabellings",
    "effect_table", "random_effects", "holm", "leave_one_out",
    "count_statistic", "path_sum", "path_edges", "compositional_statistic",
    "RESID_COVARIATES", "fit_count_glm", "count_residual_statistic", "residual_bootstrap_test", "branch_lengths_from_tree", "class_hogs_from_members",
    "module_loss_statistic", "h12_signs", "h9_control_branches", "h9_values",
    "simulate_power", "mde",
    "REGISTRY_COLUMNS", "append_registry", "md5sums", "git_commit",
]

SEED = 11
N_DRAWS = 10_000
K_MAX = 20
SMD_FLAG = 0.5
#: matching covariates of §4.6 (log d is added for H4/H5 via ``extra``)
COVARIATES = ("depth", "log_bl", "log_n_desc", "log_n_sister", "frac_tx_desc", "frac_tx_sister")
BRANCH_CLASSES = ("stem", "nested_stem", "within_multicellular", "within_unicellular",
                  "grey_zone", "other")
_MC = {"simple", "complex"}


# =============================================================================
# 1. tree and branch table
# =============================================================================

def load_tree(newick):
    """Rooted ete3 tree from a Newick string or file (format 1: named internal nodes).

    All nodes must be named and names must be unique."""
    import ete3

    if isinstance(newick, ete3.TreeNode):
        t = newick
    else:
        txt = newick
        if isinstance(newick, (str, os.PathLike)) and os.path.exists(os.fspath(newick)):
            with open(os.fspath(newick)) as fh:
                txt = fh.read()
        t = ete3.Tree(txt.strip(), format=1)
    names = [n.name for n in t.traverse()]
    if any(not nm for nm in names):
        raise ValueError("all tree nodes must be named")
    dup = pd.Series(names)[pd.Series(names).duplicated()].tolist()
    if dup:
        raise ValueError(f"duplicate node names: {dup[:5]}")
    return t


def _traits_index(traits: pd.DataFrame | None) -> pd.DataFrame:
    if traits is None:
        return pd.DataFrame(columns=["data_type", "busco_C", "mc_state"])
    tr = traits.copy()
    if "code" in tr.columns:
        tr = tr.set_index("code")
    return tr


def branch_table(tree, traits: pd.DataFrame | None = None, branch_lengths=None) -> pd.DataFrame:
    """Per-branch structure and covariates (one row per non-root node = branch into it).

    Parameters
    ----------
    tree : ete3 tree (see :func:`load_tree`).
    traits : species table with ``code``, ``data_type`` ('genome'/'transcriptome'),
        optional ``busco_C``. Tips without a row are left out of the fractions.
    branch_lengths : None (column NaN; ``log_bl`` is then dropped from matching),
        ``"tree"`` (the Newick lengths of ``tree``), a mapping child name ->
        length, or a tree with lengths (ete3 tree, path or Newick) mapped onto
        ``tree`` by :func:`branch_lengths_from_tree`.

    Columns: parent, child, depth (edges from the root to the child), n_desc
    (tips below the child), n_sister (tips below the other children of the
    parent, summed for polytomies), is_terminal, root_child, branch_length,
    log_bl, log_n_desc, log_n_sister, frac_tx_desc, frac_tx_sister (fraction of
    transcriptome tips among tips with a known data type), mean_busco_desc.
    """
    tr = _traits_index(traits)
    if branch_lengths is not None and not isinstance(branch_lengths, (str, Mapping)):
        branch_lengths = branch_lengths_from_tree(tree, branch_lengths)      # ete3 tree with lengths
    elif isinstance(branch_lengths, str) and branch_lengths != "tree":
        branch_lengths = branch_lengths_from_tree(tree, branch_lengths)      # path / Newick
    is_tx = {}
    if "data_type" in tr.columns:
        is_tx = {k: (str(v).lower() == "transcriptome") for k, v in tr["data_type"].items()
                 if isinstance(v, str)}
    busco = tr["busco_C"].to_dict() if "busco_C" in tr.columns else {}
    tips = {}
    for n in tree.traverse("postorder"):
        tips[n.name] = [n.name] if n.is_leaf() else [x for c in n.children for x in tips[c.name]]

    def frac(names):
        known = [is_tx[x] for x in names if x in is_tx]
        return float(np.mean(known)) if known else np.nan

    rows = []
    depth = {tree.name: 0}
    for n in tree.traverse("preorder"):
        if n.up is None:
            continue
        depth[n.name] = depth[n.up.name] + 1
        sis = [x for s in n.up.children if s is not n for x in tips[s.name]]
        if branch_lengths is None:
            bl = np.nan
        elif isinstance(branch_lengths, str) and branch_lengths == "tree":
            bl = float(n.dist)
        else:
            bl = float(branch_lengths.get(n.name, np.nan))
        b = [busco.get(x, np.nan) for x in tips[n.name]]
        b = [float(x) for x in b if x is not None and not pd.isna(x)]
        rows.append(dict(
            branch=n.name, parent=n.up.name, child=n.name, depth=depth[n.name],
            n_desc=len(tips[n.name]), n_sister=len(sis), is_terminal=n.is_leaf(),
            root_child=n.up.up is None, branch_length=bl,
            frac_tx_desc=frac(tips[n.name]), frac_tx_sister=frac(sis),
            mean_busco_desc=float(np.mean(b)) if b else np.nan))
    bt = pd.DataFrame(rows).set_index("branch")
    with np.errstate(divide="ignore", invalid="ignore"):
        bt["log_bl"] = np.log(bt["branch_length"].where(bt["branch_length"] > 0))
    bt["log_n_desc"] = np.log(bt["n_desc"])
    bt["log_n_sister"] = np.log(bt["n_sister"].where(bt["n_sister"] > 0))
    missing = [x for x in tips[tree.name] if x not in is_tx]
    bt.attrs["tips_without_data_type"] = missing
    return bt


def _parse_M(m) -> list[str]:
    if m is None or (isinstance(m, float) and math.isnan(m)):
        return []
    return [x.strip() for x in str(m).split(",") if x.strip()]


def _role(from_state: str, to_state: str) -> str:
    f, t = str(from_state).strip(), str(to_state).strip()
    if t in _MC:
        return "stem" if f == "unicellular" else "nested_stem"
    if t == "grey":
        return "grey_stem"
    if t == "reversal":
        return "reversal_stem"
    return "other_stem"


def node_states(tree, traits: pd.DataFrame | None, transitions: pd.DataFrame) -> pd.DataFrame:
    """State of every node: 'unicellular', 'multicellular', 'grey' or 'reversal'.

    Rule: every node at or below the crown M of a transition takes the state the
    transition leads to (simple/complex -> 'multicellular', grey -> 'grey',
    reversal -> 'reversal'); deeper crowns override shallower ones (so a reversal
    clade nested in a multicellular clade is 'reversal'). Nodes below no crown are
    'unicellular', except tips whose trait ``mc_state`` is simple/complex
    ('multicellular', flagged ``uncoded``) or grey ('grey'). Transitions whose M
    is missing (dropped) are ignored; M given as several tips makes each tip a
    crown.
    """
    tr = _traits_index(traits)
    nodes = {n.name: n for n in tree.traverse()}
    depth = {}
    for n in tree.traverse("preorder"):
        depth[n.name] = 0 if n.up is None else depth[n.up.name] + 1
    state = {nm: "unicellular" for nm in nodes}
    crown_tid = {nm: "" for nm in nodes}
    crowns = []
    for _, r in transitions.iterrows():
        for m in _parse_M(r.get("M")):
            if m not in nodes:
                warnings.warn(f"{r['transition_id']}: M node {m!r} not in tree; ignored")
                continue
            crowns.append((depth[m], m, r["transition_id"], _role(r["from_state"], r["to_state"])))
    crowns.sort(key=lambda x: (x[0], x[1]))
    st_of_role = {"stem": "multicellular", "nested_stem": "multicellular", "grey_stem": "grey",
                  "reversal_stem": "reversal"}
    for _, m, tid, role in crowns:
        s = st_of_role.get(role)
        if s is None:
            continue
        for x in nodes[m].traverse():
            state[x.name] = s
            crown_tid[x.name] = tid
    note = {nm: "" for nm in nodes}
    trait_state = tr["mc_state"].to_dict() if "mc_state" in tr.columns else {}
    for nm, n in nodes.items():
        if not n.is_leaf():
            continue
        ts = trait_state.get(nm)
        if not isinstance(ts, str):
            continue
        if crown_tid[nm] == "":
            if ts in _MC:
                state[nm] = "multicellular"
                note[nm] = "uncoded_multicellular_tip"
            elif ts == "grey":
                state[nm] = "grey"
                note[nm] = "grey_tip_without_transition"
        elif state[nm] == "multicellular" and ts == "unicellular":
            note[nm] = "unicellular_tip_inside_multicellular_crown"
    return pd.DataFrame({"state": state, "crown_transition": crown_tid, "state_note": note})


def build_branch_table(tree, traits: pd.DataFrame | None, transitions: pd.DataFrame,
                       branch_lengths=None) -> pd.DataFrame:
    """:func:`branch_table` + node states + branch class + control pool U.

    Added columns
    -------------
    state_parent, state_child : node states (:func:`node_states`).
    branch_class : one of :data:`BRANCH_CLASSES`. Transition stems by role:
        unicellular -> simple/complex = 'stem'; simple/complex -> simple/complex =
        'nested_stem'; -> grey = 'grey_zone'; reversal stems = 'other'.
        Otherwise: both ends multicellular = 'within_multicellular'; child grey =
        'grey_zone'; both ends unicellular = 'within_unicellular'; anything else
        (inside reversal clades, uncoded multicellular tips) = 'other'.
    transition_id, role : transition(s) whose M is the child ('' if none).
    primary, primary_transitions : stem/nested stem of a primary transition, i.e.
        ``analysis_set == 'primary'`` when that column exists, else
        ``primary_eligible == 1`` (several transitions sharing one M are listed;
        a warning is raised).
    adjacent_to_stem : edge into P or edge P -> sister, for every coded
        transition (P = first node of ``stem_path`` if given, else parent of M;
        also the P named in the table when it differs). For a path stem
        (``stem_path`` "P>X>M", DESIGN v3.1 D7) the edge into P and every child
        edge of P and of each intermediate node are adjacent.
    stem_path : on the stem row (child = M), the path edges (child names, ';'),
        e.g. "X;M"; intermediate path edges get the stem's class, its
        transition_id and role 'path_edge'. ``branch_length`` of a path stem is
        the sum over its edges (``edge_branch_length`` keeps the single edge);
        all other covariates are those of the last edge.
    in_U, U_exclusion : control-pool membership and the first failed rule
        (order: class, terminal, sister<2, root_child, adjacent).
    """
    bt = branch_table(tree, traits, branch_lengths)
    ns = node_states(tree, traits, transitions)
    bt["state_parent"] = ns.loc[bt["parent"], "state"].to_numpy()
    bt["state_child"] = ns.loc[bt["child"], "state"].to_numpy()
    bt["state_note"] = ns.loc[bt["child"], "state_note"].to_numpy()
    nodes = {n.name: n for n in tree.traverse()}

    tid_of: dict[str, list[str]] = {}
    role_of: dict[str, list[str]] = {}
    primary_of: dict[str, list[str]] = {}
    path_of: dict[str, list[str]] = {}
    path_edge_of: dict[str, tuple] = {}
    adjacent: set[str] = set()
    # When several transitions share an M node (e.g. the primary volvocine path stem and an exploratory
    # single-edge coding of the same clade), the row with the highest analysis priority defines the stem path:
    # process rows primary -> secondary -> exploratory -> reversal -> excluded (stable within a set).
    if "analysis_set" in transitions.columns:
        _prio = {"primary": 0, "secondary": 1, "exploratory": 2, "reversal": 3, "excluded": 4}
        transitions = (transitions.assign(_prio=transitions["analysis_set"].map(_prio).fillna(9))
                       .sort_values("_prio", kind="stable").drop(columns="_prio"))
    for _, r in transitions.iterrows():
        Ms = [m for m in _parse_M(r.get("M")) if m in nodes]
        role = _role(r["from_state"], r["to_state"])
        if "analysis_set" in transitions.columns:
            pe = str(r.get("analysis_set", "")).strip() == "primary"
        else:
            pe = str(r.get("primary_eligible", 0)).strip() in ("1", "1.0", "True", "true")
        for m in Ms:
            tid_of.setdefault(m, []).append(r["transition_id"])
            role_of.setdefault(m, []).append(role)
            if pe and role in ("stem", "nested_stem"):
                primary_of.setdefault(m, []).append(r["transition_id"])
                if nodes[m].is_leaf():
                    warnings.warn(f"{r['transition_id']}: primary_eligible but M={m} is a leaf "
                                  "(DESIGN §2 requires an internal M)")
            pm = nodes[m].up
            if pm is None:
                continue
            path = _parse_path(r.get("stem_path"), m, nodes, r["transition_id"]) if len(Ms) == 1 else None
            if path is None:
                path = [pm.name, m]
            path_of.setdefault(m, path[1:])
            for x in path[1:-1]:
                path_edge_of.setdefault(x, (r["transition_id"], role))
            ps = [nodes[x] for x in path[:-1]]
            Pn = r.get("P")
            if isinstance(Pn, str) and Pn in nodes and nodes[Pn] not in ps:
                ps.append(nodes[Pn])
            for p in ps:
                if p.up is not None and p.name not in path[1:]:
                    adjacent.add(p.name)                 # edge into P
                for s in p.children:
                    if s.name != m and s.name not in path[1:]:
                        adjacent.add(s.name)             # P (or path node) -> sister
    cls = []
    for b, r in bt.iterrows():
        roles = role_of.get(b)
        sp, sc = r["state_parent"], r["state_child"]
        if roles or b in path_edge_of:
            ro = roles[0] if roles else path_edge_of[b][1]
            c = {"stem": "stem", "nested_stem": "nested_stem", "grey_stem": "grey_zone"}.get(ro, "other")
        elif sp == "multicellular" and sc == "multicellular":
            c = "within_multicellular"
        elif sc == "grey":
            c = "grey_zone"
        elif sp == "unicellular" and sc == "unicellular":
            c = "within_unicellular"
        else:
            c = "other"
        cls.append(c)
    bt["branch_class"] = cls
    for b, c in zip(bt.index, cls):
        if c in ("stem", "nested_stem") and b not in path_edge_of:
            want = "unicellular" if c == "stem" else "multicellular"
            if bt.at[b, "state_parent"] != want:
                warnings.warn(f"{b}: {c} but parent state is {bt.at[b, 'state_parent']!r} "
                              f"(from_state implies {want!r})")
    bt["transition_id"] = [";".join(tid_of[b]) if b in tid_of else
                           (path_edge_of[b][0] if b in path_edge_of else "") for b in bt.index]
    bt["role"] = [";".join(role_of[b]) if b in role_of else
                  ("path_edge" if b in path_edge_of else "") for b in bt.index]
    bt["stem_path"] = [";".join(path_of[b]) if b in path_of else "" for b in bt.index]
    bt["edge_branch_length"] = bt["branch_length"]
    for b, pe in path_of.items():
        if len(pe) > 1 and b in bt.index:
            bt.at[b, "branch_length"] = float(bt.loc[pe, "edge_branch_length"].sum(min_count=len(pe)))
    with np.errstate(divide="ignore", invalid="ignore"):
        bt["log_bl"] = np.log(bt["branch_length"].where(bt["branch_length"] > 0))
    bt["primary"] = [b in primary_of for b in bt.index]
    bt["primary_transitions"] = [";".join(primary_of.get(b, [])) for b in bt.index]
    for b, t in primary_of.items():
        if len(t) > 1:
            warnings.warn(f"primary transitions {t} share the stem {b}")
    for b, t in tid_of.items():
        if len(t) > 1:
            warnings.warn(f"transitions {t} share M = {b}")
    bt["adjacent_to_stem"] = [b in adjacent for b in bt.index]
    reason = []
    for b, r in bt.iterrows():
        if r["branch_class"] != "within_unicellular":
            reason.append("class:" + r["branch_class"])
        elif r["is_terminal"]:
            reason.append("terminal")
        elif r["n_sister"] < 2:
            reason.append("sister_tips<2")
        elif r["root_child"]:
            reason.append("root_child")
        elif r["adjacent_to_stem"]:
            reason.append("adjacent_to_stem")
        else:
            reason.append("")
    bt["U_exclusion"] = reason
    bt["in_U"] = bt["U_exclusion"] == ""
    return bt


def _parse_path(sp, m, nodes, tid):
    """[P, ..., M] from a 'P>X>M' string, or None (absent/invalid -> warning)."""
    if not isinstance(sp, str) or ">" not in sp:
        return None
    path = [x.strip() for x in sp.split(">")]
    ok = (path[-1] == m and all(x in nodes for x in path)
          and all(nodes[path[i + 1]].up is nodes[path[i]] for i in range(len(path) - 1)))
    if not ok:
        warnings.warn(f"{tid}: stem_path {sp!r} is not a parent->child chain ending at M={m}; "
                      "using parent(M) -> M")
        return None
    return path


def path_edges(bt: pd.DataFrame, branch: str) -> list:
    """Edges (child names) of the stem ending at ``branch``; [branch] if not a path."""
    if "stem_path" in bt.columns and isinstance(bt.at[branch, "stem_path"], str) and bt.at[branch, "stem_path"]:
        return bt.at[branch, "stem_path"].split(";")
    return [branch]


def _paths(bt) -> dict:
    """{stem branch: [edges]} for path stems (>= 2 edges) in ``bt``."""
    if bt is None or "stem_path" not in getattr(bt, "columns", []):
        return {}
    return {b: p.split(";") for b, p in bt["stem_path"].items() if isinstance(p, str) and ";" in p}


def stem_branches(bt: pd.DataFrame, primary_only: bool = True,
                  roles: Iterable[str] = ("stem", "nested_stem")) -> pd.Series:
    """transition_id -> stem branch id, for stems (and nested stems) in ``bt``."""
    out = {}
    roles = set(roles)
    for b, r in bt.iterrows():
        if not r["transition_id"]:
            continue
        prim = set(r["primary_transitions"].split(";")) if "primary_transitions" in r else set()
        for tid, ro in zip(r["transition_id"].split(";"), r["role"].split(";")):
            if primary_only and tid not in prim:
                continue
            if ro in roles:
                out[tid] = b
    return pd.Series(out, dtype=object).sort_index()


# =============================================================================
# 2. matching
# =============================================================================

@dataclass
class MatchResult:
    """Matched control sets. ``sets[tid]`` lists control branch ids, nearest first."""
    stems: dict
    sets: dict
    distances: dict
    covariates: list
    dropped_covariates: dict
    smd: pd.DataFrame
    k: int
    pool: list

    def flagged(self) -> pd.DataFrame:
        return self.smd[self.smd["flag"]]


def _covariate_frame(bt: pd.DataFrame, covariates, extra) -> pd.DataFrame:
    X = pd.DataFrame(index=bt.index)
    for c in covariates:
        X[c] = bt[c] if c in bt.columns else np.nan
    if extra is not None:
        extra = pd.DataFrame(extra) if isinstance(extra, pd.Series) else extra
        for c in extra.columns:
            X[c] = extra[c].reindex(bt.index)
    return X.astype(float)


def match_controls(bt: pd.DataFrame, stems: Mapping[str, str] | pd.Series | None = None,
                   pool: Sequence[str] | None = None, covariates: Sequence[str] = COVARIATES,
                   extra: pd.DataFrame | pd.Series | None = None, k: int = K_MAX,
                   smd_threshold: float = SMD_FLAG) -> MatchResult:
    """k = min(k, |U|) Mahalanobis-nearest controls in the pool for every stem.

    Covariates are z-scored with the pool mean and SD; the Mahalanobis metric
    uses the pool covariance of the z-scores (pseudo-inverse). Covariates that
    are entirely missing or constant in the pool are dropped (``dropped_covariates``).
    Remaining missing values raise. Ties in distance are broken by branch id.
    ``extra`` adds columns (e.g. ``log_d`` for H4/H5, indexed by branch).

    SMD (per stem and covariate, raw units) = (x_stem - mean x_N) / SD_pool;
    the rows with ``transition_id == 'ALL'`` compare the mean over stems with the
    mean over all matched controls (with multiplicity). ``flag`` = |SMD| > 0.5.
    """
    if stems is None:
        stems = stem_branches(bt)
    stems = dict(stems)
    if not stems:
        raise ValueError("no stems to match")
    pool = list(bt.index[bt["in_U"]]) if pool is None else list(pool)
    if not pool:
        raise ValueError("empty control pool")
    X = _covariate_frame(bt, covariates, extra)
    sb = list(stems.values())
    dropped = {}
    use = []
    for c in X.columns:
        xp = X.loc[pool, c]
        if xp.isna().all():
            dropped[c] = "missing"
        elif xp.std(ddof=1) == 0 or len(pool) < 2:
            dropped[c] = "constant_in_pool"
        else:
            use.append(c)
    if not use:
        raise ValueError("no usable matching covariates")
    sub = X.loc[pool + sb, use]
    if sub.isna().any().any():
        bad = sub[sub.isna().any(axis=1)]
        raise ValueError(f"missing covariate values for {list(bad.index)[:10]} ({use})")
    mu = X.loc[pool, use].mean()
    sd = X.loc[pool, use].std(ddof=1)
    Z = (X[use] - mu) / sd
    Zp = Z.loc[pool].to_numpy()
    S = np.atleast_2d(np.cov(Zp, rowvar=False))
    Si = np.linalg.pinv(S)
    kk = min(int(k), len(pool))
    sets, dists = {}, {}
    pool_arr = np.asarray(pool, dtype=object)
    for tid, b in stems.items():
        d = Zp - Z.loc[b].to_numpy()
        d2 = np.einsum("ij,jk,ik->i", d, Si, d)
        order = np.lexsort((pool_arr.astype(str), np.round(d2, 12)))[:kk]
        sets[tid] = [str(x) for x in pool_arr[order]]
        dists[tid] = np.sqrt(np.maximum(d2[order], 0))
    rows = []
    for c in use:
        for tid, b in stems.items():
            xs, xc = X.loc[b, c], X.loc[sets[tid], c].mean()
            rows.append(dict(transition_id=tid, stem=b, covariate=c, stem_value=xs,
                             control_mean=xc, smd=(xs - xc) / sd[c]))
        xs = X.loc[sb, c].mean()
        xc = np.mean([X.loc[x, c] for t in stems for x in sets[t]])
        rows.append(dict(transition_id="ALL", stem="", covariate=c, stem_value=xs,
                         control_mean=xc, smd=(xs - xc) / sd[c]))
    smd = pd.DataFrame(rows)
    smd["flag"] = smd["smd"].abs() > smd_threshold
    return MatchResult(stems=stems, sets=sets, distances=dists, covariates=use,
                       dropped_covariates=dropped, smd=smd, k=kk, pool=pool)


def regression_adjusted_expectation(y: pd.Series, bt: pd.DataFrame, match: MatchResult,
                                    extra=None) -> pd.Series:
    """Optional §4.6 / review C1: E_i = mean_c y_c * exp(beta'(x_i - x_c)), beta from a
    Poisson log-link GLM of y on the matching covariates over the pool."""
    import statsmodels.api as sm

    X = _covariate_frame(bt, match.covariates, extra)[match.covariates]
    Xp = sm.add_constant(X.loc[match.pool].to_numpy(), has_constant="add")
    fit = sm.GLM(y.loc[match.pool].to_numpy(float), Xp, family=sm.families.Poisson()).fit()
    beta = np.asarray(fit.params)[1:]
    out = {}
    for tid, b in match.stems.items():
        C = match.sets[tid]
        adj = np.exp((X.loc[b].to_numpy() - X.loc[C].to_numpy()) @ beta)
        out[tid] = float(np.mean(y.loc[C].to_numpy(float) * adj))
    return pd.Series(out, name="E_adjusted")


# =============================================================================
# 3. relabelling test
# =============================================================================

def u_score(s_b: float, s_others) -> float:
    """u(b, N) = (#{c: s_c < s_b} + 0.5 #{c: s_c = s_b} + 0.5) / (|N \\ b| + 1).

    ``s_others`` are the values of N \\ b (b itself already removed)."""
    o = np.asarray(s_others, dtype=float)
    return float(((o < s_b).sum() + 0.5 * (o == s_b).sum() + 0.5) / (o.size + 1))


@dataclass
class RelabelResult:
    T_obs: float
    u: pd.Series
    T_null: np.ndarray
    p: float
    p_mcse: float
    alternative: str
    n: int
    k: dict
    n_draws: int
    seed: int
    scheme: str
    acceptance_rate: float
    n_relabellings: float
    p_floor_mc: float
    p_floor_exact: float
    dropped: dict = field(default_factory=dict)

    @property
    def min_attainable_p(self) -> float:
        return max(self.p_floor_mc, self.p_floor_exact)

    def summary(self) -> dict:
        ks = list(self.k.values())
        return dict(T_obs=self.T_obs, p=self.p, p_mcse=self.p_mcse, n=self.n,
                    k=int(np.median(ks)) if ks else 0, k_min=min(ks) if ks else 0,
                    alternative=self.alternative, n_draws=self.n_draws, seed=self.seed,
                    scheme=self.scheme, min_attainable_p=self.min_attainable_p,
                    p_floor_mc=self.p_floor_mc, p_floor_exact=self.p_floor_exact,
                    acceptance_rate=self.acceptance_rate)


def count_relabellings(sets: Sequence[Sequence], max_n: int = 22) -> float:
    """Number of tuples (c_1..c_n), c_i in N_i, all distinct (rectangular permanent,
    DP over subsets of strata). NaN if n > max_n."""
    n = len(sets)
    if n == 0:
        return 1.0
    if n > max_n:
        return float("nan")
    member: dict = {}
    for i, s in enumerate(sets):
        for c in set(s):
            member.setdefault(c, []).append(i)
    full = 1 << n
    dp = np.zeros(full, dtype=float)
    dp[0] = 1.0
    allm = np.arange(full)
    for c, idx in member.items():
        new = dp.copy()
        for i in idx:
            bit = 1 << i
            wo = allm[(allm & bit) == 0]
            new[wo | bit] += dp[wo]
        dp = new
    return float(dp[full - 1])


def min_attainable_p(k: Sequence[int] | Sequence[Sequence], n_draws: int = N_DRAWS) -> dict:
    """Floors of the relabelling p-value.

    ``k`` is either the matched sets themselves (exact count of relabellings with
    distinct controls) or the set sizes |N_i| (then sets are assumed disjoint,
    count = prod |N_i|, an upper bound). Returns p_floor_mc = 1/(n_draws+1),
    p_floor_exact = 1/(M+1) with M the number of relabellings (the smallest
    p a full enumeration could give), and min_attainable_p = their maximum.
    """
    if len(k) and not np.isscalar(k[0]):
        M = count_relabellings(k)
    else:
        M = float(np.prod([float(x) for x in k])) if len(k) else 1.0
    fm = 1.0 / (n_draws + 1)
    fe = 1.0 / (M + 1) if np.isfinite(M) else 0.0
    return dict(n_relabellings=M, p_floor_mc=fm, p_floor_exact=fe, min_attainable_p=max(fm, fe))


def _sets_and_values(values, match, stems):
    if isinstance(match, MatchResult):
        stems = dict(match.stems)
        sets = {t: list(s) for t, s in match.sets.items()}
    else:
        if stems is None:
            raise ValueError("stems mapping required when match is a plain mapping")
        stems = dict(stems)
        sets = {t: list(s) for t, s in dict(match).items()}
    v = pd.Series(values, dtype=float) if not isinstance(values, pd.Series) else values.astype(float)
    return v, stems, sets


def relabel_test(values, match, stems: Mapping[str, str] | None = None,
                 alternative: str = "two-sided", n_draws: int = N_DRAWS, seed: int = SEED,
                 scheme: str = "permute", signs: Mapping[str, float] | None = None,
                 transitions: Sequence[str] | None = None, max_tries: int = 10**8) -> RelabelResult:
    """Stratified relabelling test of §4.6.

    Parameters
    ----------
    values : per-branch statistic s_b (Series/mapping branch id -> value) covering
        stems and controls. NaN stems are dropped from the test, NaN controls are
        removed from their sets (both reported in ``dropped``).
    match : :class:`MatchResult` or mapping transition_id -> list of control ids
        (then ``stems`` maps transition_id -> stem branch id).
    alternative : 'greater' (stems higher), 'less', or 'two-sided' (|T - 1/2|).
    scheme :
        ``'permute'`` (default, primary since DESIGN v3.1 D6) -- see below.
        ``'spec'`` (sensitivity) -- literal v3 §4.6: c_i drawn in N_i, distinct across
        i (duplicate draws rejected), T_d = mean_i u(c_i, N_i \\ c_i).
        Details of ``'permute'``: permutation of the values over all branches involved
        (stems and every control of any N_i), T recomputed exactly as T_obs
        (stems vs their own N_i). Exact under exchangeability of those branches;
        unlike 'spec' it stays calibrated when the N_i overlap (see
        docs/branchstats.md, "Calibration").
    signs : optional per-transition multiplier w_i (H12): stratum i compares
        w_i * s.
    transitions : optional subset/order of transition ids to use.

    'spec' draws: c_i uniform in N_i independently for every i; draws in which
    two transitions get the same control are rejected (uniform over relabellings
    with distinct controls). Both schemes: p = (1 + #{T_d >= T_obs}) / (n_draws + 1)
    (``>=`` replaced by ``<=`` for 'less' and by |T-1/2| comparisons for
    two-sided); Monte-Carlo SE = sqrt(p (1-p) / n_draws).
    """
    if alternative not in ("greater", "less", "two-sided"):
        raise ValueError("alternative must be 'greater', 'less' or 'two-sided'")
    if scheme not in ("spec", "permute"):
        raise ValueError("scheme must be 'spec' or 'permute'")
    v, stems, sets = _sets_and_values(values, match, stems)
    tids = list(transitions) if transitions is not None else list(sets.keys())
    signs = {} if signs is None else dict(signs)
    dropped = {}
    use = []
    for t in tids:
        b = stems[t]
        if b not in v.index or pd.isna(v.get(b)):
            dropped[t] = "stem_value_missing"
            continue
        s = [c for c in sets[t] if c in v.index and not pd.isna(v[c])]
        if len(s) < len(sets[t]):
            dropped[f"{t}:controls"] = len(sets[t]) - len(s)
        if len(s) < 2:
            dropped[t] = "fewer_than_2_controls"
            continue
        if b in s:
            raise ValueError(f"stem {b} of {t} is in its own control set")
        sets[t] = s
        use.append(t)
    n = len(use)
    if n == 0:
        raise ValueError("no transition with a stem value and >= 2 controls")
    stem_ids = [stems[t] for t in use]
    if len(set(stem_ids)) < n:
        raise ValueError("two transitions share a stem branch")
    all_ctrl = sorted({c for t in use for c in sets[t]})
    if set(all_ctrl) & set(stem_ids):
        raise ValueError("a stem branch is a control of another transition")
    ids = all_ctrl + stem_ids
    gid = {x: i for i, x in enumerate(ids)}
    vals = v.reindex(ids).to_numpy(float)
    w = np.array([float(signs.get(t, 1.0)) for t in use])
    S = [np.array([gid[c] for c in sets[t]]) for t in use]
    stem_g = np.array([gid[b] for b in stem_ids])
    kdict = {t: len(sets[t]) for t in use}

    u_obs = np.array([u_score(w[i] * vals[stem_g[i]], w[i] * vals[S[i]]) for i in range(n)])
    T_obs = float(u_obs.mean())

    # leave-self-out u of every control in every stratum ('spec')
    if scheme == "spec":
        Uloo = []
        for i in range(n):
            x = w[i] * vals[S[i]]
            less = (x[None, :] < x[:, None]).sum(1)
            eq = (x[None, :] == x[:, None]).sum(1) - 1
            Uloo.append((less + 0.5 * eq + 0.5) / len(x))
    rng = np.random.default_rng(seed)
    if scheme == "spec":
        T_null, acc_rate = _draw_spec(rng, S, Uloo, n_draws, max_tries)
    else:
        T_null, acc_rate = _draw_permute(rng, vals, stem_g, S, w, n_draws), 1.0
    eps = 1e-12
    if alternative == "greater":
        hits = np.sum(T_null >= T_obs - eps)
    elif alternative == "less":
        hits = np.sum(T_null <= T_obs + eps)
    else:
        hits = np.sum(np.abs(T_null - 0.5) >= abs(T_obs - 0.5) - eps)
    p = float((1 + hits) / (n_draws + 1))
    if scheme == "spec":
        floors = min_attainable_p([sets[t] for t in use], n_draws)
    else:
        M = float(math.perm(len(ids), n))
        floors = dict(n_relabellings=M, p_floor_mc=1 / (n_draws + 1), p_floor_exact=1 / (M + 1))
    return RelabelResult(
        T_obs=T_obs, u=pd.Series(u_obs, index=use, name="u"), T_null=T_null, p=p,
        p_mcse=float(math.sqrt(p * (1 - p) / n_draws)), alternative=alternative, n=n, k=kdict,
        n_draws=n_draws, seed=seed, scheme=scheme, acceptance_rate=float(acc_rate),
        n_relabellings=floors["n_relabellings"], p_floor_mc=floors["p_floor_mc"],
        p_floor_exact=floors["p_floor_exact"], dropped=dropped)


def _draw_spec(rng, S, Uloo, n_draws, max_tries):
    n = len(S)
    kk = np.array([len(x) for x in S])
    out, got, tries = [], 0, 0
    batch = max(1000, n_draws)
    acc_rate = 1.0
    while got < n_draws:
        J = np.stack([rng.integers(0, kk[i], size=batch) for i in range(n)], axis=1)
        tries += batch
        if n > 1:
            G = np.stack([S[i][J[:, i]] for i in range(n)], axis=1)
            J = J[~(np.diff(np.sort(G, axis=1), axis=1) == 0).any(axis=1)]
        if len(J):
            out.append(np.mean([Uloo[i][J[:, i]] for i in range(n)], axis=0))
            got += len(J)
        acc_rate = got / tries
        if tries >= max_tries and got < n_draws:
            raise RuntimeError(f"relabelling: only {got} accepted draws after {tries} tries "
                               "(control sets overlap too much for distinct draws)")
        if got < n_draws:
            batch = int(min(max_tries, max(1000, 1.2 * (n_draws - got) / max(acc_rate, 1e-6))))
    return np.concatenate(out)[:n_draws], acc_rate


def _T_matrix(V, stem_g, S, w):
    """T for every row of the value matrix V (columns = global branch index)."""
    T = np.zeros(V.shape[0])
    for i in range(len(stem_g)):
        st = w[i] * V[:, [stem_g[i]]]
        c = w[i] * V[:, S[i]]
        T += ((c < st).sum(1) + 0.5 * (c == st).sum(1) + 0.5) / (c.shape[1] + 1)
    return T / len(stem_g)


def _draw_permute(rng, vals, stem_g, S, w, n_draws, chunk: int = 5000):
    out = []
    for a in range(0, n_draws, chunk):
        m = min(chunk, n_draws - a)
        V = rng.permuted(np.tile(vals, (m, 1)), axis=1)
        out.append(_T_matrix(V, stem_g, S, w))
    return np.concatenate(out)


# =============================================================================
# 4. effect sizes, pooling, multiplicity, leave-one-out
# =============================================================================

def effect_table(values, match, stems: Mapping[str, str] | None = None,
                 scale: str = "log_ratio", signs: Mapping[str, float] | None = None,
                 q: tuple = (0.05, 0.95)) -> pd.DataFrame:
    """Per-transition effect sizes.

    ``scale='log_ratio'`` (counts, H1/H2): est = log((y_i+0.5)/(mean_N y + 0.5)).
    ``scale='difference'`` (statistics that can be negative: H4/H5 logits, H9 r,
    H12): est = w_i * (s_i - mean_N s).
    The same transformation applied to every control c in N_i gives the control
    effects e_c; ``pi_low``/``pi_high`` = their 5 % / 95 % quantiles (the §4.6
    "prediction interval") and the within-transition variance is
    var = var(e_c, ddof=1) * (1 + 1/k) (spread of a single branch around the
    matched mean, plus the uncertainty of the mean). Zero variances are floored at
    the smallest positive variance in the table (or 1e-6).
    """
    v, stems, sets = _sets_and_values(values, match, stems)
    signs = {} if signs is None else dict(signs)
    rows = []
    for t, C in sets.items():
        b = stems[t]
        yc = v.reindex(C).dropna().to_numpy(float)
        y = v.get(b, np.nan)
        k = len(yc)
        if pd.isna(y) or k < 2:
            rows.append(dict(transition_id=t, stem=b, k=k, y=y))
            continue
        m = yc.mean()
        wi = float(signs.get(t, 1.0))
        if scale == "log_ratio":
            est = math.log((y + 0.5) / (m + 0.5))
            ec = np.log((yc + 0.5) / (m + 0.5))
        elif scale == "difference":
            est = wi * (y - m)
            ec = wi * (yc - m)
        else:
            raise ValueError("scale must be 'log_ratio' or 'difference'")
        var = float(np.var(ec, ddof=1) * (1 + 1 / k))
        rows.append(dict(transition_id=t, stem=b, k=k, y=float(y), expected=float(m), est=est,
                         var=var, pi_low=float(np.quantile(ec, q[0])),
                         pi_high=float(np.quantile(ec, q[1])),
                         percentile=u_score(wi * y, wi * yc)))
    df = pd.DataFrame(rows).set_index("transition_id")
    if "var" in df.columns:
        pos = df["var"][df["var"] > 0]
        floor = float(pos.min()) if len(pos) else 1e-6
        df["var_floored"] = df["var"].notna() & (df["var"] <= 0)
        df.loc[df["var_floored"], "var"] = floor
    return df


def _reml_negll(tau2, y, v):
    w = 1.0 / (v + tau2)
    mu = np.sum(w * y) / np.sum(w)
    return 0.5 * (np.sum(np.log(v + tau2)) + math.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))


def random_effects(est, var, method: str = "REML", ci_method: str = "z", alpha: float = 0.05) -> dict:
    """Random-effects pooling of per-transition estimates.

    method 'REML' (restricted likelihood maximised over tau^2 >= 0 with a bounded
    scalar optimiser, then checked against the boundary tau^2 = 0) or 'DL'
    (DerSimonian-Laird). If REML fails (optimiser error or non-finite result) the
    DL estimate is returned with ``method='DL(fallback)'``.
    ci_method 'z' (Wald, as metafor's default) or 'hksj' (Hartung-Knapp-Sidik-
    Jonkman, t with n-1 df; recommended for few transitions).
    Returns mu, se, ci_low, ci_high, z, p, tau2, Q, df, I2, n, method, ci_method.
    """
    y = np.asarray(est, dtype=float)
    v = np.asarray(var, dtype=float)
    ok = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[ok], v[ok]
    n = len(y)
    if n == 0:
        raise ValueError("no finite estimates")
    w0 = 1 / v
    mu0 = np.sum(w0 * y) / np.sum(w0)
    Q = float(np.sum(w0 * (y - mu0) ** 2))
    cdl = np.sum(w0) - np.sum(w0 ** 2) / np.sum(w0)
    tau2_dl = max(0.0, (Q - (n - 1)) / cdl) if n > 1 and cdl > 0 else 0.0
    used = method
    if method.upper() == "DL" or n < 2:
        tau2 = tau2_dl
        used = "DL" if method.upper() == "DL" else method
    elif method.upper() == "REML":
        try:
            ub = max(10.0 * (np.var(y, ddof=1) + v.max()), 1e-8)
            r = optimize.minimize_scalar(_reml_negll, bounds=(0.0, ub), args=(y, v),
                                         method="bounded", options={"xatol": 1e-12, "maxiter": 2000})
            tau2 = float(r.x)
            if _reml_negll(0.0, y, v) <= _reml_negll(tau2, y, v):
                tau2 = 0.0
            if not np.isfinite(tau2):
                raise FloatingPointError
            used = "REML"
        except Exception:  # pragma: no cover - defensive
            tau2, used = tau2_dl, "DL(fallback)"
    else:
        raise ValueError("method must be 'REML' or 'DL'")
    w = 1 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    se = float(math.sqrt(1 / np.sum(w)))
    if ci_method == "z":
        crit = stats.norm.ppf(1 - alpha / 2)
        z = mu / se
        p = 2 * stats.norm.sf(abs(z))
    elif ci_method == "hksj":
        qh = float(np.sum(w * (y - mu) ** 2) / (n - 1)) if n > 1 else 1.0
        se = float(math.sqrt(qh / np.sum(w)))
        crit = stats.t.ppf(1 - alpha / 2, max(n - 1, 1))
        z = mu / se
        p = 2 * stats.t.sf(abs(z), max(n - 1, 1))
    else:
        raise ValueError("ci_method must be 'z' or 'hksj'")
    I2 = max(0.0, (Q - (n - 1)) / Q) if Q > 0 else 0.0
    return dict(mu=mu, se=se, ci_low=mu - crit * se, ci_high=mu + crit * se, z=float(z),
                p=float(p), tau2=float(tau2), Q=Q, df=n - 1, I2=I2, n=n, method=used,
                ci_method=ci_method)


def holm(pvals) -> np.ndarray:
    """Holm step-down adjusted p-values (NaN kept in place)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    m = len(q)
    if m == 0:
        return out
    order = np.argsort(q, kind="mergesort")
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, (m - r) * q[i])
        adj[i] = min(1.0, run)
    out[ok] = adj
    return out


def leave_one_out(transitions: Sequence[str], fn: Callable[[list], Mapping | float],
                  groups: Mapping[str, str] | None = None) -> pd.DataFrame:
    """Re-run ``fn(kept_transition_ids)`` leaving out one transition at a time, or
    (with ``groups``: transition_id -> supergroup) one supergroup at a time.

    ``fn`` returns a dict (or a number, stored as ``value``); the first row
    (``left_out='none'``) is the full analysis. Use e.g.
    ``fn = lambda keep: random_effects(eff.loc[keep, 'est'], eff.loc[keep, 'var'])``
    or ``lambda keep: relabel_test(s, match, transitions=keep).summary()``.
    """
    tids = list(transitions)
    if groups is None:
        units = [(t, [t]) for t in tids]
    else:
        g = pd.Series({t: groups.get(t, "NA") for t in tids})
        units = [(gname, list(g.index[g == gname])) for gname in sorted(g.unique())]
    rows = []
    for name, drop in [("none", [])] + units:
        keep = [t for t in tids if t not in drop]
        if not keep:
            continue
        try:
            r = fn(keep)
            r = dict(r) if isinstance(r, Mapping) else {"value": r}
        except Exception as e:  # report, do not abort the sweep
            r = {"error": f"{type(e).__name__}: {e}"}
        r.update(left_out=name, n_kept=len(keep))
        rows.append(r)
    df = pd.DataFrame(rows)
    return df[["left_out", "n_kept"] + [c for c in df.columns if c not in ("left_out", "n_kept")]]


# =============================================================================
# 5. per-branch statistics s_b
# =============================================================================

def _branches_index(events: pd.DataFrame, branches) -> pd.Index:
    if branches is None:
        return pd.Index(sorted(events["branch_child"].unique()))
    return pd.Index(branches.index if isinstance(branches, pd.DataFrame) else list(branches))


def path_sum(values: pd.Series, bt: pd.DataFrame) -> pd.Series:
    """Replace the value of every path stem (``bt.stem_path`` with >= 2 edges) by
    the sum over its edges (DESIGN v3.1 D7). Other branches are unchanged."""
    out = values.copy()
    for b, edges in _paths(bt).items():
        out[b] = float(values.reindex(edges).sum(min_count=len(edges)))
    return out


def count_statistic(events: pd.DataFrame, event: str = "gain", hogs: Iterable[str] | None = None,
                    branches=None, weight: str | None = None) -> pd.Series:
    """H1/H2: number of HOG events of one kind per branch (``branch_child``).

    ``hogs`` restricts to a set of hog_ids (e.g. the sequence-novel gains of the
    gain-evidence ladder). By default one HOG = one event (``weight=None``); pass
    ``weight='n_events'`` to count duplication copies. At-origin duplications are
    excluded for ``event='dup'``. Branches without events get 0. If ``branches``
    is a branch table with path stems, their counts are summed over the path edges.
    """
    e = events[events["event"] == event]
    if event == "dup" and "at_origin" in e.columns:
        e = e[~e["at_origin"].astype(bool)]
    if hogs is not None:
        e = e[e["hog_id"].isin(set(hogs))]
    if weight is None:
        s = e.groupby("branch_child")["hog_id"].nunique()
    else:
        s = e.groupby("branch_child")[weight].sum()
    idx = _branches_index(events, branches)
    s = s.reindex(idx, fill_value=0).astype(float).rename(f"n_{event}")
    if isinstance(branches, pd.DataFrame):
        s = path_sum(s, branches)
    return s


def class_hogs_from_members(members: pd.DataFrame, genes_in_class: Iterable[str],
                            min_frac: float = 0.5) -> set:
    """HOG ids (any level) whose members carry a class domain in >= ``min_frac`` of
    genes. ``members``: long table with columns hog_id, gene (e.g. concatenated
    :meth:`mcell.hogs.HogResult.members` over levels)."""
    g = set(genes_in_class)
    f = members.assign(_in=members["gene"].isin(g)).groupby("hog_id")["_in"].mean()
    return set(f.index[f >= min_frac])


def _logit(p):
    return math.log(p / (1 - p))


def _lifter(hogs: pd.DataFrame | None):
    """Function (hog_ids, level) -> ids of their ancestors at ``level`` (a HOG
    without an ancestor there keeps its own id). Needs hog_id, level, parent_hog_id."""
    if hogs is None or "parent_hog_id" not in hogs.columns:
        def lift(ids, level):
            return set(ids)
        return lift
    lvl = dict(zip(hogs["hog_id"], hogs["level"].astype(str)))
    par = dict(zip(hogs["hog_id"], hogs["parent_hog_id"]))

    def lift(ids, level):
        out = set()
        for h in ids:
            x = h
            while x is not None and isinstance(x, str) and lvl.get(x) != level:
                x = par.get(x)
            out.add(x if isinstance(x, str) else h)
        return out
    return lift


def _parents_and_edges(events, bt):
    """(branch -> u, branch -> list of edges); u of a path stem = its first node."""
    if bt is not None:
        par = bt["parent"].to_dict()
    else:
        par = events.drop_duplicates("branch_child").set_index("branch_child")["branch_parent"].to_dict()
    edges = {b: [b] for b in par}
    for b, pe in _paths(bt).items():
        edges[b] = pe
        par[b] = bt.at[pe[0], "parent"]
    return par, edges


def compositional_statistic(events: pd.DataFrame, hogs: pd.DataFrame, class_hogs: Iterable[str],
                            bt: pd.DataFrame | None = None, min_d: int = 10,
                            exclude: Iterable[str] = ()) -> pd.DataFrame:
    """H4/H5: s_b = logit((a + 0.5)/(d + 1)) - logit(pi_u) for every branch b = u -> v.

    D_b = HOGs at u with >= 1 (non-origin) duplication on b, d = |D_b|,
    a = |D_b ∩ C|; pi_u = share of the HOGs present at u (``hogs.level == u``)
    that are in C. If pi_u is 0 or 1 it is replaced by (n_C + 0.5)/(n_u + 1)
    (``pi_smoothed``). s is NaN when d < ``min_d`` (``reason='d<10'``) or the
    branch is in ``exclude`` (e.g. WGD-flagged; ``reason='excluded'``).
    Path stems (``bt.stem_path``): D = union over the path edges of the
    duplicated HOGs, each mapped to its ancestor at the path's first node u
    (HOGs born on the path keep their own id); pi_u from that first node.
    Returns a table indexed by branch (child) with parent (= u), d, a, n_u, c_u,
    pi_u, s, log_d, reason. ``bt`` (branch table) supplies parents and paths;
    without it parents are taken from ``events``.
    """
    C = set(class_hogs)
    dup = events[(events["event"] == "dup")]
    if "at_origin" in dup.columns:
        dup = dup[~dup["at_origin"].astype(bool)]
    par, edges = _parents_and_edges(events, bt)
    lift = _lifter(hogs)
    lv = hogs["level"].astype(str)
    n_u = lv.value_counts()
    c_u = lv[hogs["hog_id"].isin(C)].value_counts()
    D = dup.groupby("branch_child")["hog_id"].agg(lambda x: set(x))
    excl = set(exclude)
    rows = []
    for b, u in par.items():
        if not u:
            continue
        if len(edges[b]) == 1:
            Db = D.get(b, set())
        else:
            Db = set()
            for e in edges[b]:
                Db |= lift(D.get(e, set()), u)
        d, a = len(Db), len(Db & C)
        nu, cu = int(n_u.get(u, 0)), int(c_u.get(u, 0))
        smooth = False
        if nu == 0:
            pi = np.nan
        else:
            pi = cu / nu
            if pi <= 0 or pi >= 1:
                pi, smooth = (cu + 0.5) / (nu + 1), True
        reason = ""
        if b in excl or any(e in excl for e in edges[b]):
            reason = "excluded"
        elif d < min_d:
            reason = f"d<{min_d}"
        elif not np.isfinite(pi):
            reason = "no_hogs_at_u"
        s = np.nan if reason else _logit((a + 0.5) / (d + 1)) - _logit(pi)
        rows.append(dict(branch=b, parent=u, n_edges=len(edges[b]), d=d, a=a, n_u=nu, c_u=cu,
                         pi_u=pi, pi_smoothed=smooth, s=s,
                         log_d=math.log(d) if d > 0 else np.nan, reason=reason))
    return pd.DataFrame(rows).set_index("branch")


def module_loss_statistic(events: pd.DataFrame, module_hogs: pd.DataFrame,
                          module_ref_size: Mapping[str, float], bt: pd.DataFrame | None = None,
                          min_present: float = 0.5, hogs: pd.DataFrame | None = None) -> pd.DataFrame:
    """H12: f_b = mean over modules of (module HOGs present at u lost on b) /
    (module HOGs present at u).

    ``module_hogs``: columns level, hog_id, module -- the HOGs present at each node
    that belong to each flagellar module. A module is skipped at u when its number
    of HOGs at u is < ``min_present`` x ``module_ref_size[module]`` (the module's
    reference HOG/family count, e.g. its size in data/h12_sets.tsv). f_b is NaN if
    every module is skipped. Path stems: module HOGs present at P (first node)
    lost on any path edge (losses on later edges are mapped back to their
    ancestor at P through ``hogs``, which is then required). The signed
    statistic s_b = w_i f_b is applied in :func:`relabel_test` /
    :func:`effect_table` through ``signs`` (see :func:`h12_signs`), because w
    depends on the transition whose stratum b is in.
    """
    loss = events[events["event"] == "loss"]
    lost = loss.groupby("branch_child")["hog_id"].agg(lambda x: set(x))
    mh = module_hogs.assign(level=module_hogs["level"].astype(str))
    by_level = {lvl: g.groupby("module")["hog_id"].agg(lambda x: set(x)) for lvl, g in mh.groupby("level")}
    par, edges = _parents_and_edges(events, bt)
    if any(len(e) > 1 for e in edges.values()) and hogs is None:
        raise ValueError("path stems need the hogs table to map losses back to P")
    lift = _lifter(hogs)
    rows = []
    for b, u in par.items():
        if not u:
            continue
        mods = by_level.get(u, pd.Series(dtype=object))
        if len(edges[b]) == 1:
            L = lost.get(b, set())
        else:
            L = set()
            for e in edges[b]:
                L |= lift(lost.get(e, set()), u)
        fr = {}
        for m, ref in module_ref_size.items():
            H = mods.get(m, set()) if len(mods) else set()
            if len(H) == 0 or len(H) < min_present * ref:
                continue
            fr[m] = len(H & L) / len(H)
        rows.append(dict(branch=b, parent=u, n_modules=len(fr),
                         f=float(np.mean(list(fr.values()))) if fr else np.nan,
                         **{f"f_{m}": fr.get(m, np.nan) for m in module_ref_size}))
    return pd.DataFrame(rows).set_index("branch")


# ---- H1/H2 residual count statistic (DESIGN v3.1, D6) ------------------------

#: covariates of the H1/H2 NB-GLM (fitted on U only)
RESID_COVARIATES = ("log_bl", "log_n_desc", "log_n_sister", "frac_tx_desc")


def _irls(X, y, alpha, beta=None, max_iter=100, tol=1e-9):
    if beta is None:
        beta = np.zeros(X.shape[1])
        beta[0] = math.log(y.mean() + 0.1)
    for it in range(max_iter):
        eta = np.clip(X @ beta, -30, 30)
        mu = np.exp(eta)
        W = mu / (1 + alpha * mu)
        z = eta + (y - mu) / mu
        XtW = X.T * W
        try:
            new = np.linalg.solve(XtW @ X, XtW @ z)
        except np.linalg.LinAlgError:
            new = np.linalg.lstsq(XtW @ X, XtW @ z, rcond=None)[0]
        if not np.all(np.isfinite(new)):
            return beta, False
        if np.max(np.abs(new - beta)) < tol * (1 + np.max(np.abs(beta))):
            return new, True
        beta = new
    return beta, False


def _nb_ll(alpha, y, mu):
    from scipy.special import gammaln
    r = 1.0 / alpha
    return float(np.sum(gammaln(y + r) - gammaln(r) - gammaln(y + 1)
                        + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu) + 1e-300)))


def fit_count_glm(X, y, family: str = "nb", max_outer: int = 50, start: dict | None = None) -> dict:
    """NB2 GLM (log link; alpha by maximum likelihood, alternating with IRLS) or
    Poisson. If the NB fit fails (non-convergence, non-finite values, alpha at
    the upper bound) it falls back to Poisson with quasi-dispersion
    phi = Pearson X^2/(n-p) (``method='poisson_quasi'``, ``fallback_reason``).
    X must contain the intercept column; ``start`` (a previous result) warm-starts
    beta and alpha. Returns beta, alpha, phi, method, converged, fallback_reason."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, p = X.shape
    reason = ""
    if family == "nb":
        try:
            if start is not None and start.get("method") == "nb":
                beta, la = np.asarray(start["beta"], float).copy(), math.log(max(start["alpha"], 1e-5))
            else:
                beta, ok = _irls(X, y, 0.0)
                la = math.log(max(_nb_dispersion(y), 1e-3))
            conv = False
            ll_old = -np.inf
            for _ in range(max_outer):
                beta, ok = _irls(X, y, math.exp(la), beta)
                mu = np.exp(np.clip(X @ beta, -30, 30))
                lo, hi = (-12, 4) if start is None else (max(-12, la - 2), min(4, la + 2))
                r = optimize.minimize_scalar(lambda t: -_nb_ll(math.exp(t), y, mu), bounds=(lo, hi),
                                             method="bounded", options={"xatol": 1e-8})
                d = abs(r.x - la)
                la = r.x
                ll = -float(r.fun)
                # converged when alpha (log scale) or the likelihood stops moving
                # (a flat profile near alpha -> 0 makes log alpha drift slowly)
                if ok and (d < 1e-6 or abs(ll - ll_old) < 1e-9 * (1 + abs(ll))):
                    conv = True
                    break
                ll_old = ll
            beta, ok = _irls(X, y, math.exp(la), beta)
            if not (conv and ok):
                reason = "nb_not_converged"
            elif la > 3.9:
                reason = "nb_alpha_at_bound"
            elif not np.all(np.isfinite(beta)):
                reason = "nb_nonfinite"
            else:
                return dict(beta=beta, alpha=math.exp(la), phi=1.0, method="nb", converged=True,
                            fallback_reason="")
        except Exception as e:  # pragma: no cover - defensive
            reason = f"nb_error:{type(e).__name__}"
    beta, ok = _irls(X, y, 0.0)
    mu = np.exp(np.clip(X @ beta, -30, 30))
    phi = float(np.sum((y - mu) ** 2 / mu) / max(n - p, 1))
    if reason:
        warnings.warn(f"NB-GLM failed ({reason}); Poisson with quasi-dispersion used")
    return dict(beta=beta, alpha=0.0, phi=phi, method="poisson_quasi" if reason else "poisson",
                converged=bool(ok), fallback_reason=reason)


def count_residual_statistic(y: pd.Series, bt: pd.DataFrame, pool: Sequence[str] | None = None,
                             covariates: Sequence[str] = RESID_COVARIATES, loo: bool = True,
                             family: str = "nb") -> pd.DataFrame:
    """H1/H2 statistic of DESIGN v3.1 (D6): s_b = log((y_b + 0.5)/(yhat_b + 0.5)).

    yhat from an NB-GLM ``y ~ log_bl + log_n_desc + log_n_sister + frac_tx_desc``
    fitted on the control pool U only (``bt.in_U`` or ``pool``); every pool
    branch gets a leave-one-out prediction (``loo=True``), all other branches
    the full-pool fit. ``y`` must already hold path-summed counts for path
    stems (:func:`count_statistic` with the branch table does this) and ``bt``
    the path branch lengths (:func:`build_branch_table`). Covariates constant
    or missing in the pool (e.g. the 1.0 branch-length placeholder) are dropped;
    branches with a missing covariate get NaN. Fit details (method, alpha,
    fallbacks, covariates used) are in ``.attrs``.
    """
    pool = list(bt.index[bt["in_U"]]) if pool is None else list(pool)
    X = bt.reindex(columns=list(covariates)).astype(float)
    use, dropped = [], {}
    for c in covariates:
        xp = X.loc[pool, c]
        if xp.isna().all():
            dropped[c] = "missing"
        elif xp.std(ddof=1) == 0 or xp.isna().any():
            dropped[c] = "constant_in_pool" if xp.std(ddof=1) == 0 else "missing_in_pool"
        else:
            use.append(c)
    mu_, sd_ = X.loc[pool, use].mean(), X.loc[pool, use].std(ddof=1)
    Z = (X[use] - mu_) / sd_
    D = np.column_stack([np.ones(len(Z)), Z.to_numpy()])
    D = pd.DataFrame(D, index=bt.index)
    yy = y.reindex(bt.index).astype(float)
    pool = [b for b in pool if np.isfinite(yy[b])]
    full = fit_count_glm(D.loc[pool].to_numpy(), yy[pool].to_numpy(), family)
    fits = [dict(left_out="", **{k: v for k, v in full.items() if k != "beta"})]
    yhat = pd.Series(np.exp(np.clip(D.to_numpy() @ full["beta"], -30, 30)), index=bt.index)
    how = pd.Series("full", index=bt.index)
    if loo:
        for b in pool:
            keep = [x for x in pool if x != b]
            f = fit_count_glm(D.loc[keep].to_numpy(), yy[keep].to_numpy(), family)
            fits.append(dict(left_out=b, **{k: v for k, v in f.items() if k != "beta"}))
            yhat[b] = float(np.exp(np.clip(D.loc[b].to_numpy() @ f["beta"], -30, 30)))
            how[b] = "loo"
    bad = Z.isna().any(axis=1)
    yhat[bad] = np.nan
    out = pd.DataFrame({"y": yy, "yhat": yhat, "s": np.log((yy + 0.5) / (yhat + 0.5)), "fit": how})
    fits = pd.DataFrame(fits)
    out.attrs.update(covariates=use, dropped_covariates=dropped, method=full["method"],
                     alpha=full["alpha"], phi=full["phi"], n_pool=len(pool),
                     n_fallback=int((fits["method"] == "poisson_quasi").sum()), fits=fits)
    return out


def residual_bootstrap_test(y: pd.Series, bt: pd.DataFrame, match, stems=None,
                            alternative: str = "two-sided", n_boot: int = 999, seed: int = SEED,
                            covariates: Sequence[str] = RESID_COVARIATES, family: str = "nb",
                            signs=None) -> dict:
    """Parametric-bootstrap p for T on the H1/H2 residual statistic.

    The permutation null (``relabel_test`` on s_b) ignores that all stems share
    the error of the fitted GLM; when stems lie outside the covariate range of
    U (extrapolation) their residuals move together and the permutation test is
    anti-conservative (docs/branchstats.md). Here H0 is "stems follow the model
    fitted on U": counts of every stem and pool branch are simulated from the
    full-pool fit (NB with alpha-hat; Poisson-quasi fits use an NB with
    alpha = (phi - 1)/mean(mu)), the residual statistic (with LOO) and T are
    recomputed, p = (1 + #{T* at least as extreme}) / (n_boot + 1).
    Returns T_obs, p, p_mcse, n_boot, T_boot, and the observed residual table.
    """
    obs = count_residual_statistic(y, bt, covariates=covariates, family=family)
    r0 = relabel_test(obs["s"], match, stems=stems, n_draws=1, signs=signs)
    T_obs = r0.T_obs
    pool = list(bt.index[bt["in_U"]])
    cols = obs.attrs["covariates"]
    X = bt[cols].astype(float)
    Z = (X - X.loc[pool].mean()) / X.loc[pool].std(ddof=1)
    D = np.column_stack([np.ones(len(Z)), Z.to_numpy()])
    yy = y.reindex(bt.index).astype(float)
    ok = np.isfinite(yy.to_numpy()) & ~Z.isna().any(axis=1).to_numpy()
    pidx = [bt.index.get_loc(b) for b in pool if ok[bt.index.get_loc(b)]]
    full = fit_count_glm(D[pidx], yy.to_numpy()[pidx], family)
    mu = np.exp(np.clip(D @ full["beta"], -30, 30))
    alpha = full["alpha"] if full["method"] == "nb" else max((full["phi"] - 1) / mu[pidx].mean(), 0.0)
    rng = np.random.default_rng(seed)
    Tb = np.empty(n_boot)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i in range(n_boot):
            ys = pd.Series(_rnb(rng, mu, alpha), index=bt.index, dtype=float).where(ok)
            sb = count_residual_statistic(ys, bt, covariates=covariates, family=family)["s"]
            Tb[i] = relabel_test(sb, match, stems=stems, n_draws=1, signs=signs).T_obs
    eps = 1e-12
    if alternative == "greater":
        hits = np.sum(Tb >= T_obs - eps)
    elif alternative == "less":
        hits = np.sum(Tb <= T_obs + eps)
    else:
        hits = np.sum(np.abs(Tb - 0.5) >= abs(T_obs - 0.5) - eps)
    p = float((1 + hits) / (n_boot + 1))
    return dict(T_obs=T_obs, p=p, p_mcse=math.sqrt(p * (1 - p) / n_boot), n_boot=n_boot,
                T_boot=Tb, residuals=obs, alpha=alpha, method=full["method"])


def branch_lengths_from_tree(guide_tree, bl_tree=None, placeholder: float = 1.0) -> dict:
    """Branch lengths for the guide-tree branches, read from a tree with lengths
    (e.g. data/species_tree_bl_euk.nwk; ete3 tree, Newick string or path).

    Branches are matched by their descendant tip set (or its complement, for an
    unrooted / differently rooted tree); unmatched branches get NaN (warning).
    If ``bl_tree`` is None or the path does not exist yet, every branch gets
    ``placeholder`` (warning); the constant is then dropped from the GLM.
    """
    import ete3
    names = [n.name for n in guide_tree.traverse() if n.up is not None]
    if bl_tree is None or (isinstance(bl_tree, (str, os.PathLike)) and not isinstance(bl_tree, ete3.TreeNode)
                           and not str(bl_tree).lstrip().startswith("(") and not os.path.exists(bl_tree)):
        warnings.warn(f"no branch-length tree ({bl_tree}); placeholder {placeholder} used")
        return {nm: placeholder for nm in names}
    if not isinstance(bl_tree, ete3.TreeNode):
        txt = open(bl_tree).read() if os.path.exists(str(bl_tree)) else str(bl_tree)
        try:
            bl_tree = ete3.Tree(txt.strip(), format=1)
        except Exception:
            bl_tree = ete3.Tree(txt.strip(), format=0)
    alltips = frozenset(guide_tree.get_leaf_names())
    split = {}
    for n in bl_tree.traverse():
        if n.up is not None:
            split.setdefault(frozenset(n.get_leaf_names()) & alltips, float(n.dist))
    out, miss = {}, []
    for n in guide_tree.traverse():
        if n.up is None:
            continue
        t = frozenset(n.get_leaf_names())
        v = split.get(t, split.get(alltips - t))
        if v is None:
            miss.append(n.name)
            v = np.nan
        out[n.name] = v
    if miss:
        warnings.warn(f"{len(miss)} guide-tree branches without a length, e.g. {miss[:3]}")
    return out


def h12_signs(expected: Mapping[str, str]) -> dict:
    """transition_id -> w (+1 where loss is expected, -1 otherwise)."""
    return {t: (1.0 if str(e).strip().lower() in ("loss", "loss_expected", "1", "+1") else -1.0)
            for t, e in expected.items()}


def h9_control_branches(tree, bt: pd.DataFrame, min_tips: int = 2) -> list:
    """Branches in U whose child node has >= 2 child clades with >= ``min_tips`` tips."""
    nodes = {n.name: n for n in tree.traverse()}
    out = []
    for b in bt.index[bt["in_U"]]:
        big = [c for c in nodes[b].children if len(c) >= min_tips]
        if len(big) >= 2:
            out.append(b)
    return out


def h9_values(r_stems: Mapping[str, float], stems: Mapping[str, str],
              r_controls: Mapping[str, float], orient: bool = True, seed: int = SEED) -> pd.Series:
    """H9 plug-in: per-branch values for :func:`relabel_test`.

    ``r_stems``: transition_id -> r_i (precomputed median log rate ratio, P->M vs
    P->S). ``r_controls``: control branch -> r_c computed in an arbitrary fixed
    orientation of its two child clades; with ``orient=True`` each r_c is
    multiplied by a random sign (seeded; branches processed in sorted order) to
    give the random orientation of §4.6.
    """
    out = {stems[t]: float(r) for t, r in r_stems.items()}
    keys = sorted(r_controls)
    sg = np.random.default_rng(seed).choice([-1.0, 1.0], size=len(keys)) if orient else np.ones(len(keys))
    for kx, s in zip(keys, sg):
        out[kx] = s * float(r_controls[kx])
    return pd.Series(out, dtype=float)


# =============================================================================
# 6. power / minimum detectable effect
# =============================================================================

def _nb_dispersion(y) -> float:
    y = np.asarray(y, dtype=float)
    m, s2 = y.mean(), y.var(ddof=1)
    return max((s2 - m) / (m ** 2), 1e-8) if m > 0 else 1e-8


def _rnb(rng, mean, alpha):
    mean = np.asarray(mean, dtype=float)
    if alpha <= 1e-8:
        return rng.poisson(mean)
    r = 1.0 / alpha
    return rng.negative_binomial(r, r / (r + np.maximum(mean, 1e-12)))


def simulate_power(y: pd.Series, match: MatchResult, effects: Sequence[float] = (1.0, 1.5, 2.0, 3.0),
                   n_sim: int = 200, alphas: Sequence[float] = (0.01, 0.05),
                   alternative: str = "greater", n_draws: int = 2000, seed: int = SEED,
                   scheme: str = "permute", mean_mode: str = "observed",
                   dispersion: float | None = None) -> pd.DataFrame:
    """Power of the relabelling test for a planted multiplicative effect on stems.

    Counts are negative binomial with dispersion alpha (var = m + alpha m^2)
    estimated by moments from the observed control counts ``y`` over the pool
    (or given). Means: ``mean_mode='observed'`` -> each control keeps its observed
    count as its mean and stem i has mean effect * mean_{N_i} y; ``'pooled'`` ->
    every control has the pool mean and stems effect * pool mean. Each replicate
    re-draws all counts and runs :func:`relabel_test` (``n_draws`` draws).
    Returns one row per effect with power at each alpha.
    """
    rng = np.random.default_rng(seed)
    ctrl = sorted({c for s in match.sets.values() for c in s})
    yc = y.reindex(match.pool).dropna().to_numpy(float)
    a = _nb_dispersion(yc) if dispersion is None else float(dispersion)
    if mean_mode == "observed":
        base_c = y.reindex(ctrl).to_numpy(float)
        base_s = np.array([y.reindex(match.sets[t]).mean() for t in match.stems])
    elif mean_mode == "pooled":
        base_c = np.full(len(ctrl), yc.mean())
        base_s = np.full(len(match.stems), yc.mean())
    else:
        raise ValueError("mean_mode must be 'observed' or 'pooled'")
    stem_ids = [match.stems[t] for t in match.stems]
    rows = []
    for eff in effects:
        ps = np.empty(n_sim)
        for r in range(n_sim):
            vals = pd.Series(np.concatenate([_rnb(rng, base_c, a), _rnb(rng, eff * base_s, a)]),
                             index=ctrl + stem_ids, dtype=float)
            res = relabel_test(vals, match, alternative=alternative, n_draws=n_draws,
                               seed=int(rng.integers(2**31)), scheme=scheme)
            ps[r] = res.p
        row = dict(effect=eff, n=len(stem_ids), k=match.k, dispersion=a, n_sim=n_sim,
                   n_draws=n_draws, scheme=scheme, alternative=alternative)
        for al in alphas:
            row[f"power_{al:g}"] = float(np.mean(ps <= al))
        rows.append(row)
    return pd.DataFrame(rows)


def mde(power: pd.DataFrame, target: float = 0.8) -> dict:
    """Smallest planted effect reaching ``target`` power, per alpha (log-linear
    interpolation between the simulated effects; NaN if never reached)."""
    out = {}
    df = power.sort_values("effect")
    for c in [c for c in df.columns if c.startswith("power_")]:
        e, pw = df["effect"].to_numpy(float), df[c].to_numpy(float)
        val = np.nan
        for i in range(len(e)):
            if pw[i] >= target:
                if i == 0 or pw[i - 1] >= target:
                    val = e[i]
                else:
                    f = (target - pw[i - 1]) / (pw[i] - pw[i - 1])
                    val = float(np.exp(np.log(e[i - 1]) + f * (np.log(e[i]) - np.log(e[i - 1]))))
                break
        out[c.replace("power_", "alpha_")] = val
    return out


# =============================================================================
# 7. analysis registry
# =============================================================================

REGISTRY_COLUMNS = ["timestamp", "hypothesis", "statistic", "n", "k", "estimate", "ci_low", "ci_high",
                    "p", "p_mcse", "holm_p", "tier", "input_md5s", "git_commit",
                    "seed", "n_draws", "scheme", "notes"]


def md5sums(paths: Iterable[str]) -> str:
    """'basename:md5;...' for the given files (in the given order)."""
    out = []
    for p in paths:
        h = hashlib.md5()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        out.append(f"{os.path.basename(p)}:{h.hexdigest()}")
    return ";".join(out)


def git_commit(path: str | None = None) -> str:
    """HEAD commit of the repository containing ``path`` (or the cwd); 'NA' if none."""
    d = path or os.getcwd()
    if os.path.isfile(d):
        d = os.path.dirname(d)
    try:
        r = subprocess.run(["git", "-C", d or ".", "rev-parse", "HEAD"], capture_output=True,
                           text=True, timeout=10)
        return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else "NA"
    except Exception:
        return "NA"


def append_registry(path: str, row: Mapping) -> dict:
    """Append one row to the analysis registry TSV (results/analysis_registry.tsv).

    Keys must be among :data:`REGISTRY_COLUMNS` (unknown keys raise). ``timestamp``
    (UTC ISO) and ``git_commit`` are filled when absent; ``input_md5s`` may be a
    list of file paths (hashed with :func:`md5sums`) or a pre-formatted string;
    ``tier`` must be 'primary' or 'exploratory' (or 'secondary'). The header is
    written when the file is new and checked otherwise. Returns the written row.
    """
    bad = set(row) - set(REGISTRY_COLUMNS)
    if bad:
        raise ValueError(f"unknown registry fields: {sorted(bad)}")
    r = {c: row.get(c, "") for c in REGISTRY_COLUMNS}
    if not r["timestamp"]:
        r["timestamp"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    if not r["git_commit"]:
        r["git_commit"] = git_commit(os.path.dirname(os.path.abspath(path)))
    if isinstance(r["input_md5s"], (list, tuple)):
        r["input_md5s"] = md5sums(r["input_md5s"])
    if r["tier"] not in ("primary", "exploratory", "secondary"):
        raise ValueError("tier must be 'primary', 'exploratory' or 'secondary'")
    if not r["hypothesis"]:
        raise ValueError("hypothesis is required")

    def fmt(x):
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return "NA"
        if isinstance(x, float):
            return f"{x:.6g}"
        return str(x).replace("\t", " ").replace("\n", " ")

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    new = not os.path.exists(path) or os.path.getsize(path) == 0
    if not new:
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
        if head != REGISTRY_COLUMNS:
            raise ValueError(f"{path}: header differs from REGISTRY_COLUMNS")
    with open(path, "a") as fh:
        if new:
            fh.write("\t".join(REGISTRY_COLUMNS) + "\n")
        fh.write("\t".join(fmt(r[c]) for c in REGISTRY_COLUMNS) + "\n")
    return r
