exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
ev=cache.events; H=cache.hogs
met=ctx.clade("Metazoa"); cho=ctx.clade("Choanoflagellata")
for br in ["Choanozoa","Metazoa"]:
    d=ev[(ev.event=="dup")&(ev.branch_child==br)]
    kids=[k for s in d.child_hog_ids.astype(str) for k in s.split(";") if k in L.row_of.index]
    rows=[int(L.row_of[k]) for k in kids]
    ms=L.member_species(rows)
    nmet=sum(1 for r in rows if set(ms[r])<=met); ncho=sum(1 for r in rows if set(ms[r])<=cho)
    nboth=sum(1 for r in rows if set(ms[r])&met and set(ms[r])&cho)
    print(f"dup events on stem({br}): {len(d)} events, n_events sum {int(d.n_events.sum())}; child HOGs {len(rows)}: metazoan-only {nmet}, choano-only {ncho}, both {nboth}")
# losses on Choanoflagellata stem whose HOG has only metazoan members and whose parent dup is on Choanozoa stem
l=ev[(ev.event=="loss")&(ev.branch_child=="Choanoflagellata")]
print("loss events on stem(Choanoflagellata):",int(l.n_events.sum()))
lm=ev[(ev.event=="loss")&(ev.branch_child=="Metazoa")]
print("loss events on stem(Metazoa):",int(lm.n_events.sum()))
