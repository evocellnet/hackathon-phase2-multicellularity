"""Task 2a: independently restrict both OrthoXMLs to the 125 shared species,
strip all level annotations (TaxRange properties, taxonId on groups, <taxonomy>),
drop empty groups and collapse single-child groups. EukaryOMA species are
renamed to our species codes (via the OMA X-code prefix of protId).

usage: v2_filter_orthoxml.py ancestra|eukaryoma
"""
import gzip, sys
from xml.sax.saxutils import quoteattr
import pandas as pd
from lxml import etree

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
W = f"{P}/results/r2/verification/eukaryoma/work"
which = sys.argv[1]
smap = pd.read_csv(f"{P}/results/eukaryoma_xval/species_map.tsv", sep="\t")
smap = smap[smap.eukaryoma_code.notna()]
keep_codes = set(smap.code)
x2code = dict(zip(smap.eukaryoma_code, smap.code))
src = f"{P}/runs/euk_fast/hogs.orthoxml" if which == "ancestra" else f"{P}/data/raw/eukaryoma/oma-hogs.orthoXML.gz"
out = gzip.open(f"{W}/{which}.filtered.orthoxml.gz", "wt", compresslevel=3)
out.write('<?xml version="1.0" encoding="UTF-8"?>\n<orthoXML xmlns="http://orthoXML.org/2011/" origin="verif" originVersion="1" version="0.5">\n')

def L(tag):
    return tag.rsplit("}", 1)[-1]

kept_gene_ids = set()
n_sp = n_genes = n_fam = 0

def prune(el):
    """return serialised pruned subtree or None; collapses unary groups."""
    t = L(el.tag)
    if t == "geneRef":
        gid = el.get("id")
        return ("g", f'<geneRef id="{gid}"/>') if gid in kept_gene_ids else None
    if t not in ("orthologGroup", "paralogGroup"):
        return None
    kids = [k for k in (prune(c) for c in el) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        return kids[0]
    return ("grp", f"<{t}>" + "".join(k[1] for k in kids) + f"</{t}>")

fh = gzip.open(src, "rb") if src.endswith(".gz") else open(src, "rb")
groups_started = False
depth = 0
for ev, el in etree.iterparse(fh, events=("start", "end"), huge_tree=True):
    t = L(el.tag)
    if ev == "start":
        if t in ("orthologGroup", "paralogGroup"):
            if not groups_started:
                out.write("<groups>\n")
                groups_started = True
            depth += 1
        continue
    if t == "species":
        genes = [g for g in el.iter() if L(g.tag) == "gene"]
        if which == "ancestra":
            code = el.get("name")
        else:
            xs = {g.get("protId")[:5] for g in genes}
            assert len(xs) == 1, xs
            code = x2code.get(xs.pop())
        if code in keep_codes:
            n_sp += 1
            out.write(f'<species name="{code}"><database name="d" version="1"><genes>\n')
            for g in genes:
                kept_gene_ids.add(g.get("id"))
                out.write(f'<gene id="{g.get("id")}" protId={quoteattr(g.get("protId"))}/>\n')
                n_genes += 1
            out.write("</genes></database></species>\n")
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
    elif t in ("orthologGroup", "paralogGroup"):
        depth -= 1
        if depth == 0:
            r = prune(el)
            if r is not None and r[0] == "grp":
                out.write(r[1] + "\n")
                n_fam += 1
            el.clear()
            while el.getprevious() is not None:
                del el.getparent()[0]
out.write("</groups>\n</orthoXML>\n")
out.close()
print(which, "species", n_sp, "genes", n_genes, "multi-gene families", n_fam)
