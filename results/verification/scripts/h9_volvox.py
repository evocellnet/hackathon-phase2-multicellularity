"""Independent H9 r_i (own code; ete3 only for Newick parsing)."""
import sys, os, glob, numpy as np, ete3, pandas as pd
P = sys.argv[1]
guide = ete3.Tree(f"{P}/data/guide_tree_euk.nwk", format=1)
gtips = {n.name: set(n.get_leaf_names()) for n in guide.traverse()}
def esm(n):   # equal-split mean path length to tips
    if n.is_leaf(): return 0.0
    return sum(c.dist + esm(c) for c in n.children) / len(n.children)
def r_gene(t, Pn, Mn, Sn):
    present = set(t.get_leaf_names())
    tipsP, tipsM = gtips[Pn], gtips[Mn]
    tipsS = set().union(*(gtips[s] for s in Sn))
    if len(present - tipsP) < 1: return None, "no_outgroup"
    mt = present & tipsM
    if len(mt) < 2: return None, "lt2_M"
    mrca = t.get_common_ancestor(list(mt))
    if set(mrca.get_leaf_names()) != mt: return None, "M_not_monophyletic"
    # binding rule (DESIGN 4.6): MRCA of the present M tips must be M itself, i.e. present tips span both
    # children of guide node M (checked on the guide tree, independent of gene-tree node labels)
    gM = guide.search_nodes(name=Mn)[0]
    if sum(1 for ch in gM.children if set(ch.get_leaf_names()) & mt) < 2: return None, "MRCA_not_M"
    if len(present & tipsS) < 1: return None, "no_S"
    # P node in gene tree: MRCA of present (M ∪ S ∪ intermediate) tips under guide P
    pt = present & tipsP
    Pnode = t.get_common_ancestor(list(pt))
    # stem: sum of edges from mrca(M) up to Pnode
    stem, n = 0.0, mrca
    while n is not Pnode:
        stem += n.dist; n = n.up
    s_kids = [c for c in Pnode.children if set(c.get_leaf_names()) <= tipsS]
    if not s_kids: return None, "no_S_child"
    LM = stem + esm(mrca)
    LS = np.mean([c.dist + esm(c) for c in s_kids])
    if LM <= 0 or LS <= 0: return None, "zero"
    return float(np.log(LM / LS)), ""
def unit(Pn, Mn, Sn, files):
    rs, why = {}, {}
    for f in files:
        g = os.path.basename(f)[:-7]
        t = ete3.Tree(open(f).read().strip(), format=1)
        r, w = r_gene(t, Pn, Mn, Sn)
        if r is None: why[w] = why.get(w, 0) + 1
        else: rs[g] = r
    return rs, why
files = sorted(glob.glob(f"{P}/runs/phylo/euk/genetrees/*.bl.nwk"))
print("gene trees:", len(files))
rep = pd.read_csv(f"{P}/results/h9/h9_gene_table.tsv", sep="\t")
out = {}
for name, (Pn, Mn, Sn) in {"primary T_volvocine_tgv (P=Chlamydomonas_TGV, M=Gonium_Volvocaceae, S=Chlamydomonas)": ("Chlamydomonas_TGV", "Gonium_Volvocaceae", ["Chlamydomonas"]),
                           "buggy single-edge alternative (P=TGV_clade, M=Gonium_Volvocaceae, S=Tetrsoc)": ("TGV_clade", "Gonium_Volvocaceae", ["Tetrsoc"]),
                           "check T_metazoa": ("Choanozoa", "Metazoa", ["Choanoflagellata"])}.items():
    rs, why = unit(Pn, Mn, Sn, files)
    v = np.array(list(rs.values()))
    rng = np.random.default_rng(7)
    boots = np.median(v[rng.integers(0, len(v), (2000, len(v)))], axis=1)
    print(f"{name}: n_genes={len(v)} r_i={np.median(v):.5f} CI~[{np.quantile(boots,.025):.4f},{np.quantile(boots,.975):.4f}] excluded={why}")
    tid = {"Chlamydomonas_TGV": "T_volvocine_tgv", "Choanozoa": "T_metazoa"}.get(Pn)
    if tid:
        rr = rep[(rep.unit == tid) & rep.included].set_index("gene").r
        common = sorted(set(rr.index) & set(rs))
        print(f"   reported: n_genes={len(rr)} r_i={rr.median():.5f}; same gene set: {set(rr.index)==set(rs)}; max |diff r_g| = {max(abs(rr[g]-rs[g]) for g in common):.2e}")
