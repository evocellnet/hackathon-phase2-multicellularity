"""Diagnosis: count, in the producer's filtered OrthoXMLs, unary groups and in particular
orthologGroups whose only child is a paralogGroup (the parser merges such a paralogGroup into
its parent, erasing the duplication)."""
import gzip, sys
from lxml import etree
P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
for name in ("eukaryoma", "ancestra"):
    n_unary = n_o_p = 0
    fh = gzip.open(f"{P}/runs/eukaryoma_xval/{name}.orthoxml.gz", "rb")
    for ev, el in etree.iterparse(fh, events=("end",), huge_tree=True):
        t = el.tag.rsplit("}", 1)[-1]
        if t in ("orthologGroup", "paralogGroup"):
            kids = [k for k in el if k.tag.rsplit("}", 1)[-1] in ("orthologGroup", "paralogGroup", "geneRef")]
            if len(kids) == 1:
                n_unary += 1
                if t == "orthologGroup" and kids[0].tag.endswith("paralogGroup"):
                    n_o_p += 1
            if el.getparent() is not None and el.getparent().tag.endswith("groups"):
                el.clear()
                while el.getprevious() is not None:
                    del el.getparent()[0]
        elif t == "gene":
            el.clear()
    print(name, "unary groups", n_unary, "orthologGroup->single paralogGroup", n_o_p, flush=True)
