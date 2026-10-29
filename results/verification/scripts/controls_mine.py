"""Independent check of positive controls C01, C05, C10 from Pfam tables + OrthoXML (lxml) + raw DIAMOND hits."""
import sys, gzip, glob, os, json, collections
sys.setrecursionlimit(200000)
from lxml import etree
import ete3, pandas as pd
P, S = sys.argv[1], sys.argv[2]
NS = "{http://orthoXML.org/2011/}"
ACC = {"PF00110", "PF03514", "PF00102", "PF00017", "PF07714"}
hits = collections.defaultdict(set)          # acc -> proteins
dom = collections.defaultdict(list)          # protein -> [(env_from, acc)] for SH2/PK
sp_with = collections.defaultdict(set)
for f in sorted(glob.glob(f"{P}/data/annotations/pfam/euk/*.pfam.tsv.gz")):
    sp = os.path.basename(f).split(".")[0]
    with gzip.open(f, "rt") as fh:
        next(fh)
        for line in fh:
            x = line.split("\t")
            if x[1] in ACC:
                hits[x[1]].add(x[0]); sp_with[x[1]].add(sp)
                if x[1] in ("PF00017", "PF07714"):
                    dom[x[0]].append((int(x[5]), x[1]))
print("pfam species files:", len(glob.glob(f'{P}/data/annotations/pfam/euk/*.pfam.tsv.gz')), {a: len(v) for a, v in hits.items()})
guide = ete3.Tree(f"{P}/runs/euk_fast/stages/input/guide_tree.nwk", format=1)
tips = {n.name: set(n.get_leaf_names()) for n in guide.traverse()}
spc = lambda g: g.split("_", 1)[0]
res = {}
# ---- C05 (a): SH2 before PK_Tyr_Ser-Thr by envelope start
arch = collections.Counter()
for prot, d in dom.items():
    sh2 = [e for e, a in d if a == "PF00017"]; pk = [e for e, a in d if a == "PF07714"]
    if sh2 and pk and min(sh2) < max(pk):
        s = spc(prot)
        for cl in ("Metazoa", "Choanoflagellata", "Filasterea"):
            if s in tips[cl]: arch[cl] += 1
res["C05a_SH2_before_PK_proteins"] = dict(arch)
# ---- C01(b), C10(b): outside species with hits
res["C01b_outside_Metazoa_species_PF00110"] = sorted(sp_with["PF00110"] - tips["Metazoa"])
res["C10b_outside_Anydrophyta_species_PF03514"] = sorted(sp_with["PF03514"] - tips["Anydrophyta"])
# ---- OrthoXML pass
tin, tout, name2id, id2name, parent_of, depth = {}, {}, {}, {}, {}, {}
gene_sp, gene_prot = {}, {}
cnt = [0]
def walk(el, par, d):
    t = int(el.get("id")); nm = el.get("name"); name2id[nm] = t; id2name[t] = nm; parent_of[t] = par; depth[t] = d
    tin[t] = cnt[0]; cnt[0] += 1
    for ch in el:
        if ch.tag == NS + "taxon": walk(ch, t, d + 1)
    tout[t] = cnt[0]
anc_eq = lambda a, b: tin[a] <= tin[b] and tout[b] <= tout[a]
internal = set()
def build(el):
    tag = el.tag[len(NS):]
    if tag == "geneRef":
        gid = int(el.get("id")); return ("g", name2id[gene_sp[gid]], [], gene_prot[gid], el.get("id"))
    kids = [build(c) for c in el if c.tag in (NS + "orthologGroup", NS + "paralogGroup", NS + "geneRef")]
    return ("o" if tag == "orthologGroup" else "p", int(el.get("taxonId")), kids, None, el.get("id"))
def genes(v, acc=None):
    if v[0] == "g": return [v[3]]
    out = []
    for w in v[2]: out += genes(w)
    return out
def crossed(u, v):
    """internal species nodes X crossed by gene-tree edge u->v (u may be ('virt', t))."""
    ku, tu = u[0], u[1]; kv, tv = v[0], v[1]
    # X on path from tv up to tu
    xs = []
    x = tv
    while x is not None:
        below = (anc_eq(x, tv) and x != tv) if kv == "p" else anc_eq(x, tv)
        above = anc_eq(tu, x) if ku in ("p", "virt") else (anc_eq(tu, x) and tu != x)
        if below and above and x in internal: xs.append(x)
        if not above: break
        x = parent_of[x]
    return xs
def edges(u, v):
    yield u, v
    for w in v[2]:
        yield from edges(v, w)
pfhog = lambda mem, acc: len(mem) >= 2 and sum(m in hits[acc] for m in mem) >= 0.5 * len(mem)
c01, c10, c05b = [], [], []
fam_c10 = set(); fam_c01 = set()
ctx = etree.iterparse(f"{P}/runs/euk_fast/hogs.orthoxml", events=("end",), huge_tree=True)
for ev, el in ctx:
    tag = el.tag
    if tag == NS + "species":
        nm = el.get("name")
        for g in el.iter(NS + "gene"):
            gene_sp[int(g.get("id"))] = nm; gene_prot[int(g.get("id"))] = g.get("protId")
        el.clear()
    elif tag == NS + "taxonomy":
        for ch in el:
            if ch.tag == NS + "taxon": walk(ch, None, 0)
        internal = {t for t in id2name if any(parent_of.get(c) == t for c in id2name)}
        MET, ANY, CHO = name2id["Metazoa"], name2id["Anydrophyta"], name2id["Choanozoa"]
        el.clear()
    elif tag in (NS + "orthologGroup", NS + "paralogGroup") and el.getparent() is not None and el.getparent().tag == NS + "groups":
        fid = el.get("id")
        root = build(el)
        allg = genes(root)
        rel = [a for a in ("PF00110", "PF03514", "PF00102") if any(g in hits[a] for g in allg)]
        if rel:
            k, t = root[0], root[1]
            top = [(("virt", t, [], None, None), root)] if k == "o" else [(root, w) for w in root[2]]
            gained_at = t
            for (u0, v0) in top:
                for (u, v) in edges(u0, v0):
                    xs = crossed(u, v)
                    if not xs: continue
                    mem = genes(v)
                    for x in xs:
                        if gained_at == MET and "PF00110" in rel and pfhog(mem, "PF00110"):
                            c01.append((fid, id2name[x], v[4])); fam_c01.add(fid)
                        if gained_at == ANY and "PF03514" in rel and pfhog(mem, "PF03514"):
                            c10.append((fid, id2name[x], v[4])); fam_c10.add(fid)
                        if x == CHO and "PF00102" in rel and pfhog(mem, "PF00102"):
                            sps = {spc(m) for m in mem}
                            if sps & tips["Metazoa"] and sps & (tips["Filozoa"] - tips["Metazoa"]):
                                c05b.append((fid, v[4]))
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
res["C01a_PF-HOG_PF00110_origin_Metazoa_stem"] = dict(n_hog_rows=len(c01), n_families=len(fam_c01), families=sorted(fam_c01))
res["C10a_PF-HOG_PF03514_origin_Anydrophyta_stem"] = dict(n_hog_rows=len(c10), n_families=len(fam_c10), families=sorted(fam_c10))
res["C05b_PF-HOG_PF00102_at_Choanozoa_with_Metazoa_and_nonmetazoan_Filozoa"] = dict(n=len(c05b), n_families=len({f for f, _ in c05b}))
# ---- C10 (c): best outside-clade DIAMOND hit domain for reps of gained PF-HOG families
reps = pd.read_parquet(f"{P}/runs/gainladder/euk/reps.parquet")
r = reps[reps.family_id.isin(fam_c10 | fam_c01)]
q2f = dict(zip(r.gene, r.family_id))
eukc = set(tips[guide.get_tree_root().name])
best = {}
for fn in ("hits_pass1.tsv.gz", "hits_pass2.tsv.gz"):
    with gzip.open(f"{P}/runs/gainladder/euk/{fn}", "rt") as fh:
        for line in fh:
            q = line.split("\t", 1)[0]
            if q in q2f:
                x = line.rstrip("\n").split("\t"); s = spc(x[1]); f = q2f[q]
                clade = tips["Anydrophyta"] if f in fam_c10 else tips["Metazoa"]
                if s in clade: continue
                b = float(x[5])
                if q not in best or b > best[q][0]:
                    best[q] = (b, s, "euk" if s in eukc else "prok", float(x[4]))
per_fam = collections.defaultdict(list)
for q, v in best.items(): per_fam[q2f[q]].append(v)
res["C10c_best_outside_hits"] = {f: sorted(per_fam.get(f, []), reverse=True)[:3] for f in sorted(fam_c10)}
res["C10c_families_with_rep_best_hit_prok"] = sum(any(v[2] == "prok" for v in per_fam.get(f, [])) for f in fam_c10)
res["C01_families_with_any_outside_hit(seq-novel if 0)"] = {f: len(per_fam.get(f, [])) for f in sorted(fam_c01)}
print(json.dumps(res, indent=1))
json.dump(res, open(f"{S}/controls_mine.json", "w"), indent=1)
