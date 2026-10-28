import pyhmmer, sys, gzip, collections
P="/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
want={"PF01049","PF00362","PF00907","PF01414","PF06816","PF07684"}
hmms=[]
with pyhmmer.plan7.HMMFile(f"{P}/data/raw/pfam/Pfam-A.prereg_subset.hmm") as f:
    for h in f:
        acc=h.accession.decode() if isinstance(h.accession,bytes) else h.accession
        if acc.split(".")[0] in want: hmms.append(h)
print("hmms",[ (h.name, h.accession, h.cutoffs.gathering) for h in hmms], file=sys.stderr)
sps=open("obazoa.txt").read().split()
alpha=pyhmmer.easel.Alphabet.amino()
out=open("relaxed_hits.tsv","w")
out.write("species\tprotein\tpfam\tseq_evalue\tseq_score\tbest_dom_score\tbest_dom_ievalue\tGA_seq\tGA_dom\tpasses_GA\n")
for s in sps:
    with pyhmmer.easel.SequenceFile(f"{P}/data/proteomes/euk/{s}.fa", digital=True, alphabet=alpha) as sf:
        seqs=sf.read_block()
    for top in pyhmmer.hmmsearch(hmms, seqs, cpus=4, E=0.01, domE=0.01):
        q=top.query; acc=(q.accession.decode() if isinstance(q.accession,bytes) else q.accession).split(".")[0]
        ga=q.cutoffs.gathering
        for hit in top:
            if hit.evalue>0.01: continue
            bd=max(hit.domains, key=lambda d:d.score)
            name=hit.name.decode() if isinstance(hit.name,bytes) else hit.name
            passes = hit.score>=ga[0] and any(d.score>=ga[1] for d in hit.domains)
            out.write(f"{s}\t{name}\t{acc}\t{hit.evalue:.2e}\t{hit.score:.1f}\t{bd.score:.1f}\t{bd.i_evalue:.2e}\t{ga[0]}\t{ga[1]}\t{int(passes)}\n")
    out.flush()
print("done",file=sys.stderr)
