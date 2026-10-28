exec(open("/tmp/claude-209372/-scratch-spascare-agentic-ai-workspaces-phase2-multicellularity/c136af28-0d96-5ed4-9d5c-35c3f102e2f8/scratchpad/load.py").read())
import collections
hol=ctx.clade("Holozoa"); op=ctx.clade("Opisthokonta")
for hid in ["Eukaryota|HOG0058630","Eukaryota|HOG0060409"]:
    r=L.row_of[hid]; ms=L.member_species([r])[r]
    c=collections.Counter("Holozoa" if s in hol else ("Holomycota" if s in op else "non-Opisthokonta") for s in ms)
    no=collections.Counter(s for s in ms if s not in op)
    print(hid, len(ms), dict(c), "non-opisthokont species:", dict(no))
