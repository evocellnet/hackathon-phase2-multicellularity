"""Independent per-branch HOG gain/loss/dup counts from Ancestra OrthoXML (lxml only).
Semantics (OrthoXML/pyHAM style): orthologGroup = speciation at taxonId; paralogGroup = duplication on the
branch directly above its taxonId; gene = at its species leaf. Ancestral gene at X = gene-tree edge u->v with
v at-or-below X and u strictly above X. Gain on p->c = top-level group whose top lies on p->c.
"""
import sys, json, collections
sys.setrecursionlimit(200000)
from lxml import etree
NS = "{http://orthoXML.org/2011/}"
fn = sys.argv[1]
targets = [tuple(x.split(">")) for x in sys.argv[2].split(",")]   # parent>child
out = sys.argv[3]

# taxonomy -> intervals
tin, tout, name2id, id2name, parent_of = {}, {}, {}, {}, {}
gene_sp = {}
counter = [0]
def walk_tax(el, par):
    tid = int(el.get("id")); nm = el.get("name")
    name2id[nm] = tid; id2name[tid] = nm; parent_of[tid] = par
    tin[tid] = counter[0]; counter[0] += 1
    for ch in el:
        if ch.tag == NS + "taxon":
            walk_tax(ch, tid)
    tout[tid] = counter[0]
def anc_eq(a, b):  # a ancestor of or equal to b
    return tin[a] <= tin[b] and tout[b] <= tout[a]

sp_leaf = {}   # species name -> taxon id
res = {t: collections.Counter() for t in targets}
tgt_ids = None

gained = collections.defaultdict(list)
dupfams = collections.defaultdict(set)
cur = [None]
def build(el):
    """return node tuple (kind, taxon, children)"""
    tag = el.tag[len(NS):]
    if tag == "geneRef":
        return ("g", sp_leaf[gene_sp[int(el.get("id"))]], [])
    kids = [build(c) for c in el if c.tag in (NS + "orthologGroup", NS + "paralogGroup", NS + "geneRef")]
    t = int(el.get("taxonId"))
    if tag == "orthologGroup":
        tr = [p.get("value") for p in el if p.tag == NS + "property" and p.get("name") == "TaxRange"]
        if tr and name2id.get(tr[0]) != t:
            res_meta["taxrange_mismatch"] += 1
        return ("o", t, kids)
    return ("p", t, kids)

def at_or_below(v, X):
    k, t, _ = v
    if k == "p":
        return anc_eq(X, t) and X != t
    return anc_eq(X, t)
def strictly_above(u, X):   # u is a node tuple or ('virt', t) for virtual edge above a top ortho
    k, t = u[0], u[1]
    if k in ("p", "virt"):
        return anc_eq(t, X)
    return anc_eq(t, X) and t != X

def crossings(u, v, X):
    """ancestral genes at X among edge u->v and its descendants (first crossings on each path)."""
    if at_or_below(v, X) and strictly_above(u, X):
        return [(u, v)]
    out = []
    # descend only if v is above X (strictly above, or para on branch above X... general: X in subtree region of v)
    if strictly_above(v, X) or (v[0] == "p" and v[1] == X):
        for w in v[2]:
            out.extend(crossings(v, w, X))
    return out

res_meta = collections.Counter()
def process(root):
    k, t, kids = root
    for (p, c) in tgt_ids:
        r = res[(id2name[p], id2name[c])]
        # gains
        if k == "o" and t == c:
            r["gains"] += 1; gained[f"{id2name[p]}>{id2name[c]}"].append(cur[0])
        if k == "p" and t == c:
            r["gains"] += 1; r["dups_at_origin"] += len(kids) - 1; gained[f"{id2name[p]}>{id2name[c]}"].append(cur[0])
        virt = ("virt", t) if k == "o" else ("p_above", None)
        if k == "o":
            top_edges = [(("virt", t), root)]
        else:
            top_edges = [(root, w) for w in kids]   # para top: children edges start on branch above T
        anc_p = []
        for (u, v) in top_edges:
            anc_p.extend(crossings(u, v, p))
        r["n_parent"] += len(anc_p)
        nchild = []
        for (u, v) in top_edges:
            nchild.extend(crossings(u, v, c))
        r["n_child"] += len(nchild)
        hog_dup = 0
        for (u, v) in anc_p:
            kk = len(crossings(u, v, c))
            if kk == 0:
                r["losses"] += 1
            elif kk > 1:
                r["dups"] += kk - 1; hog_dup = 1; dupfams[f"{id2name[p]}>{id2name[c]}"].add(cur[0])
        if k == "p" and t == c and len(kids) > 1:
            hog_dup = hog_dup  # dups at origin counted separately
        r["fam_with_dup"] += hog_dup
        r["anc_with_dup"] += sum(1 for (u, v) in anc_p if len(crossings(u, v, c)) > 1)

ctx = etree.iterparse(fn, events=("end",), huge_tree=True)
nfam = 0
for ev, el in ctx:
    tag = el.tag
    if tag == NS + "species":
        nm = el.get("name")
        for g in el.iter(NS + "gene"):
            gene_sp[int(g.get("id"))] = nm
        el.clear()
    elif tag == NS + "taxonomy":
        for ch in el:
            if ch.tag == NS + "taxon":
                walk_tax(ch, None)
        for nm, tid in name2id.items():
            sp_leaf[nm] = tid
        tgt_ids = [(name2id[p], name2id[c]) for p, c in targets]
        for p, c in tgt_ids:
            assert parent_of[c] == p, (id2name[p], id2name[c])
        el.clear()
    elif tag in (NS + "orthologGroup", NS + "paralogGroup") and el.getparent() is not None and el.getparent().tag == NS + "groups":
        cur[0] = el.get("id"); process(build(el)); nfam += 1
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
out_d = {f"{p}>{c}": dict(v) for (p, c), v in res.items()}
out_d["_meta"] = dict(n_top_groups=nfam, n_genes=len(gene_sp), **res_meta)
json.dump(out_d, open(out, "w"), indent=1)
json.dump({"gained": gained, "dup_families": {k: sorted(v) for k, v in dupfams.items()}}, open(out.replace(".json", "_lists.json"), "w"))
print(json.dumps(out_d, indent=1))
