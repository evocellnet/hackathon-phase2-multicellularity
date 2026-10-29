#!/usr/bin/env python
"""W5 task A: Ancestra structure-rescue variants vs default on panel R1 (Holozoa).

Steps (one per call; each refuses to overwrite its outputs):
  cache   <v>   parse the variant's OrthoXML (mcell.hogs) + pre-registered Pfam subset -> analysis cache
  controls <v>  round-1 positive controls relevant to R1 (scripts/r2/panels/controls_panel: same rules,
                aliases Fungi->Holomycota, Apusomonadida_Opisthokonta->Eukaryota; no gain ladder)
  xvalprep <v>  restrict the variant's OrthoXML to R1 species present in EukaryOMA, strip level
                annotations, parse on the pruned R1 tree (or v = eukaryoma: the EukaryOMA OrthoXML
                already filtered to our 125 species in runs/eukaryoma_xval, restricted further)
  xval          family-level agreement (classify_gain / classify_dup of scripts/eukaryoma/xval_analysis.py,
                unchanged) on every internal branch of the pruned tree, both directions, each variant
  counts        per-branch gains/losses/dups, singletons, runtime per variant
Variants: default = runs/r2/R1; 3di = runs/r2/R1_3di; 3di_guard = runs/r2/R1_3di_guard.
Outputs: results/r2/w5/r1/  (work files in results/r2/w5/work/r1/)
"""
import json, os, sys, time, warnings

ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts/eukaryoma"))
sys.path.insert(0, os.path.join(ROOT, "scripts/r2/panels"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

VARIANTS = {"default": "runs/r2/R1", "3di": "runs/r2/R1_3di", "3di_guard": "runs/r2/R1_3di_guard"}
OUT = os.path.join(ROOT, "results/r2/w5/r1")
WK = os.path.join(ROOT, "results/r2/w5/work/r1")
TREE_FULL = os.path.join(ROOT, "runs/r2/R1/stages/input/guide_tree.nwk")
TREE_PR = os.path.join(WK, "xval", "tree_pruned.nwk")
STEMS = {"Metazoa": "Choanozoa>Metazoa", "Choanozoa": "Filozoa>Choanozoa"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def guard(p):
    if os.path.exists(p):
        sys.exit(f"refusing to overwrite {p}")


def cpaths(v):
    rd = os.path.join(ROOT, VARIANTS[v])
    return dict(orthoxml=os.path.join(rd, "hogs.orthoxml"), tree=os.path.join(rd, "stages/input/guide_tree.nwk"),
                cache=os.path.join(WK, f"cache_{v}"), toolkit=os.path.join(ROOT, "data/toolkit_pfam.tsv"),
                h12_sets=os.path.join(ROOT, "data/h12_sets.tsv"),
                pfam=os.path.join(ROOT, "results/r2/panels/R1/work/pfam_prereg"),
                transitions=os.path.join(ROOT, "data/transitions.tsv"))


def step_cache(v):
    from mcell import analysis as an
    cp = cpaths(v)
    guard(os.path.join(cp["cache"], "DONE.json"))
    an.build_cache("R1", paths=cp, levels="all", log=log)


def step_controls(v):
    from mcell import analysis as an
    from mcell import controls as pc
    import controls_panel as cpn
    os.makedirs(OUT, exist_ok=True)
    dst = os.path.join(OUT, f"controls_{v}.tsv")
    guard(dst)
    cp = cpaths(v)
    cache = an.load_cache("R1", paths=cp)
    st = pd.read_csv(os.path.join(ROOT, "data/dense/R1/species.tsv"), sep="\t")
    sp = st[["code", "data_type", "busco_C", "mc_state"]].copy()
    sp["low_quality"] = (sp["busco_C"] < 50).astype(int)
    cpn.PanelCtx.alias = cpn.ALIASES["R1"]
    ctx = cpn.PanelCtx(cache=cache, species=sp, transitions=pd.read_csv(
        os.path.join(ROOT, "data/dense/R1/transitions.tsv"), sep="\t"), paths=cp)
    ref = pd.read_csv(os.path.join(ROOT, "results/positive_controls/euk/controls.tsv"), sep="\t").set_index("control")
    rows = []
    for cid in cpn.RELEVANT["R1"]:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = pc.EUK_RULES[cid](ctx)
        except KeyError as e:
            s = pc.Score(cid, "not_applicable", f"missing node/species in panel: {e}")
        rows.append(dict(variant=v, control=cid, deep_node=cid in ("C03", "C06", "C07"),
                         round1_recovered=ref.loc[cid, "recovered"] if cid in ref.index else "",
                         recovered=s.recovered, components=";".join(f"{k}={x}" for k, x in s.components.items()),
                         observed=s.observed, notes=s.notes))
    # diagnostics for the deep-node controls (origin / duplication branches of the PF lineages)
    from positive_controls import diag
    dg = pd.DataFrame([diag(v, cache, pf) for pf in ("PF01049", "PF00362", "PF00907")])
    pd.DataFrame(rows).to_csv(dst, sep="\t", index=False)
    dg.to_csv(os.path.join(OUT, f"controls_diag_{v}.tsv"), sep="\t", index=False)
    log(f"{v}: " + ", ".join(f"{r['control']}:{r['recovered']}" for r in rows))


def r1_eoma_species():
    sm = pd.read_csv(os.path.join(ROOT, "results/eukaryoma_xval/species_map.tsv"), sep="\t")
    r1 = pd.read_csv(os.path.join(ROOT, "data/dense/R1/species.tsv"), sep="\t")["code"]
    return sorted(set(sm["code"]) & set(r1))


def step_xvalprep(v):
    import xval_lib as xl
    from mcell.hogs import parse_hogs
    os.makedirs(os.path.join(WK, "xval"), exist_ok=True)
    keep = r1_eoma_species()
    if not os.path.exists(TREE_PR):
        t = xl.prune_tree(xl.load_tree(TREE_FULL), set(keep))
        xl.write_newick(t, TREE_PR + ".tmp")
        os.replace(TREE_PR + ".tmp", TREE_PR)
        log(f"pruned tree {len(keep)} species: {open(TREE_PR).read().strip()}")
    out = os.path.join(WK, "xval", "parsed", v)
    guard(os.path.join(out, "run_info.json"))
    xml = os.path.join(WK, "xval", f"{v}.orthoxml.gz")
    src = (os.path.join(ROOT, "runs/eukaryoma_xval/eukaryoma.orthoxml.gz") if v == "eukaryoma"
           else os.path.join(ROOT, VARIANTS[v], "hogs.orthoxml"))
    if not os.path.exists(xml):
        st = xl.filter_orthoxml(src, xml + ".tmp.gz", {c: c for c in keep}, strip_levels=True, log=log)
        os.replace(xml + ".tmp.gz", xml)
        log(f"{v} filter: {st}")
    t0 = time.time()
    res = parse_hogs(xml, newick=TREE_PR)
    bad = res.check_conservation()
    res.write(out, fmt="parquet")
    json.dump(dict(stats=res.stats, conservation_violations=len(bad), wall_s=time.time() - t0, xml=xml),
              open(os.path.join(out, "run_info.json"), "w"), indent=1, default=str)
    log(f"{v}: parsed {res.stats['n_genes']:,} genes; violations {len(bad)}")


def step_xval():
    from xval_analysis import M, classify_gain, classify_dup, load_rename
    from mcell.hogs import load_species_tree
    dst = os.path.join(OUT, "xval_agreement_per_branch.tsv")
    guard(dst)
    os.makedirs(OUT, exist_ok=True)
    keep = r1_eoma_species()
    tree = load_species_tree(TREE_PR, keep)
    ren = load_rename()
    P = os.path.join(WK, "xval", "parsed")
    eo = M("eukaryoma", tree, ren, d=os.path.join(P, "eukaryoma"))
    variants = [v for v in VARIANTS if os.path.exists(os.path.join(P, v, "run_info.json"))]
    ev_rows = []
    internal = [i for i in range(len(tree.names)) if not tree.is_leaf[i] and tree.parent[i] >= 0]
    for v in variants:
        an_ = M(v, tree, d=os.path.join(P, v))
        for ci in internal:
            c, top = tree.names[ci], int(tree.parent[ci])
            role = "stem" if c in STEMS else "other"
            for src, dst_, direc in ((an_, eo, "ancestra->eukaryoma"), (eo, an_, "eukaryoma->ancestra")):
                e = src.events[src.events["branch_child"] == c]
                for r in e[(e["event"] == "gain") & (e["n_species"] >= 2)].itertuples():
                    res = classify_gain(src, dst_, r.family_id, ci, {ci}, top)
                    ev_rows.append(dict(variant=v, branch=c, role=role, direction=direc, event="gain", **res))
                for r in e[(e["event"] == "dup") & ~e["at_origin"]].itertuples():
                    ids = [x for x in r.child_hog_ids.split(";") if x]
                    res = classify_dup(src, dst_, ids, ci, top)
                    ev_rows.append(dict(variant=v, branch=c, role=role, direction=direc, event="dup", **res))
        log(f"{v}: classified ({len(ev_rows):,} rows so far)")
    E = pd.DataFrame(ev_rows)
    E.to_csv(os.path.join(OUT, "xval_events.tsv.gz"), sep="\t", index=False)
    classes = ["same", "deeper", "shallower", "other", "unmapped"]
    A = E.groupby(["variant", "branch", "role", "direction", "event"])["cls"].value_counts().unstack(fill_value=0)
    A = A.reindex(columns=classes, fill_value=0)
    A["n"] = A.sum(axis=1)
    A["n_mapped"] = A["n"] - A["unmapped"]
    A["frac_same"] = (A["same"] / A["n_mapped"]).round(4)
    A = A.reset_index()
    A.to_csv(dst, sep="\t", index=False)
    rows = []
    for (v, direc, evt), s in A.groupby(["variant", "direction", "event"]):
        for g, ss in (("Metazoa_stem", s[s.branch == "Metazoa"]), ("Choanozoa_stem", s[s.branch == "Choanozoa"]),
                      ("both_stems", s[s.role == "stem"]), ("other_internal", s[s.role == "other"])):
            same, nm = int(ss["same"].sum()), int(ss["n_mapped"].sum())
            rows.append(dict(variant=v, direction=direc, event=evt, group=g, n_events=int(ss["n"].sum()),
                             n_mapped=nm, same=same, frac_same=round(same / nm, 4) if nm else np.nan,
                             deeper=int(ss["deeper"].sum()), shallower=int(ss["shallower"].sum())))
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(OUT, "xval_agreement_summary.tsv"), sep="\t", index=False)
    print(S.to_string())


def step_counts():
    from mcell.hogs import load_species_tree  # noqa: F401
    dst = os.path.join(OUT, "counts_per_branch.tsv")
    guard(dst)
    rows, run = [], []
    for v in VARIANTS:
        cdir = os.path.join(WK, f"cache_{v}")
        if not os.path.exists(os.path.join(cdir, "DONE.json")):
            continue
        b = pd.read_csv(os.path.join(cdir, "branch_summary.tsv"), sep="\t", keep_default_na=False)
        b.insert(0, "variant", v)
        rows.append(b)
        genes = pd.read_parquet(os.path.join(cdir, "genes.parquet"))
        fam_n = genes.groupby("family_id").size()
        rd = os.path.join(ROOT, VARIANTS[v])
        man = json.load(open(os.path.join(rd, "manifest.json")))
        stg = {k: x.get("wall_s") for k, x in man["stages"].items()}
        resc = man["stages"].get("rescue", {}).get("stats", {})
        merge_stats = man["stages"]["merge"].get("stats", {})
        run.append(dict(variant=v, n_genes=len(genes), n_families=len(fam_n),
                        n_singleton_genes=int((fam_n == 1).sum()),
                        frac_singleton=round(float((fam_n == 1).sum() / len(genes)), 5),
                        n_multi_gene_families=int((fam_n > 1).sum()),
                        stage_wall_s=json.dumps(stg), stage_skipped=json.dumps(
                            {k: x.get("skipped") for k, x in man["stages"].items()}),
                        rescue_stats=json.dumps({k: x for k, x in resc.items() if not isinstance(x, (list, dict))}),
                        rescue_guard=json.dumps(merge_stats.get("rescue_guard", {}), default=str)[:2000]))
    B = pd.concat(rows, ignore_index=True)
    os.makedirs(OUT, exist_ok=True)
    B.to_csv(dst, sep="\t", index=False)
    pd.DataFrame(run).to_csv(os.path.join(OUT, "run_summary.tsv"), sep="\t", index=False)
    piv = B.pivot_table(index="branch_child", columns="variant", values=["gains", "losses", "dups"], aggfunc="sum")
    piv.to_csv(os.path.join(OUT, "counts_wide.tsv"), sep="\t")
    print(pd.DataFrame(run).drop(columns=["rescue_guard"]).to_string())
    print(piv.loc[["Metazoa", "Choanozoa", "Filozoa", "Holozoa", "Choanoflagellata"]].to_string())


if __name__ == "__main__":
    step, *args = sys.argv[1:]
    {"cache": step_cache, "controls": step_controls, "xvalprep": step_xvalprep, "xval": step_xval,
     "counts": step_counts}[step](*args)
