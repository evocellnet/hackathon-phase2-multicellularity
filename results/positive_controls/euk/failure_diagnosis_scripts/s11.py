exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
ev=cache.events; H=cache.hogs
g=ev[ev.event=="gain"]
print("gain events:",len(g)," at_origin True:",int(g.at_origin.sum()))
# what HOG ids do gain events refer to? are they root HOGs?
ids=set(H.hog_id[H.parent_idx<0])
print("gain hog_id is a root HOG:", g.hog_id.isin(ids).mean(), " child ids root:", g.child_hog_ids.astype(str).map(lambda s: all(x in ids for x in s.split(';'))).mean())
print(g.head(3)[["branch_parent","branch_child","hog_id","child_hog_ids","at_origin"]].to_string())
# origin column of hogs vs parent_idx<0
print("H.origin == (parent_idx<0):", ((H.origin)==(H.parent_idx<0)).mean())
# C12 raw counts
import gzip
for s in ["Volvcar","Chlarei"]:
    d=pd.read_csv(f"data/annotations/pfam_prereg/euk/{s}.pfam.tsv.gz",sep="\t")
    print(s,"raw PF12499 proteins:",d[d.pfam_acc.str.startswith("PF12499")].protein.nunique(), "proteome size", sum(1 for l in open(f"data/proteomes/euk/{s}.fa") if l.startswith(">")))
