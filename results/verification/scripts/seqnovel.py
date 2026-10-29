import gzip, json, sys, collections
import pandas as pd, ete3
P, S = sys.argv[1], sys.argv[2]
lists = json.load(open(f"{S}/events_mine_lists.json"))["gained"]
tree = ete3.Tree(f"{P}/runs/euk_fast/stages/input/guide_tree.nwk", format=1)
targets = sys.argv[3].split(",")
reps = pd.read_parquet(f"{P}/runs/gainladder/euk/reps.parquet")
searched = set(open(f"{P}/runs/gainladder/euk/searched_genes.txt").read().split())
q2fam, fam2br, clade = {}, {}, {}
for k in targets:
    p, c = k.split(">")
    node = tree.search_nodes(name=c)[0]
    clade[k] = set(node.get_leaf_names())
    fams = set(lists[k])
    for f in fams: fam2br[f] = k
    r = reps[reps.family_id.isin(fams)]
    bad = r[~r.species.isin(clade[k])]
    print(k, "gained families", len(fams), "reps", len(r), "reps outside clade (should be 0):", len(bad),
          "reps searched", r.gene.isin(searched).sum(), "families w/o searched rep", len(fams - set(r[r.gene.isin(searched)].family_id)))
    for g, f in zip(r.gene, r.family_id):
        q2fam.setdefault(g, set()).add(f)
outside = collections.Counter(); outside_tx = collections.Counter()
nlines = 0
for fn in ("hits_pass1.tsv.gz", "hits_pass2.tsv.gz"):
    with gzip.open(f"{P}/runs/gainladder/euk/{fn}", "rt") as fh:
        for line in fh:
            q, s = line.split("\t", 2)[:2]
            if q in q2fam:
                nlines += 1
                sp = s.split("_", 1)[0]
                for f in q2fam[q]:
                    if sp not in clade[fam2br[f]]:
                        outside[f] += 1
res = {}
for k in targets:
    fams = set(lists[k])
    nov = [f for f in fams if outside[f] == 0]
    res[k] = dict(hog_gains=len(fams), seq_novel=len(nov), non_novel=len(fams) - len(nov))
    json.dump(sorted(nov), open(f"{S}/novel_{k.replace('>','__')}.json", "w"))
print("hit lines for target queries:", nlines)
print(json.dumps(res, indent=1))
json.dump(res, open(f"{S}/seqnovel_mine.json", "w"), indent=1)
